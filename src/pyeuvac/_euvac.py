import numpy as np
import xarray as xr
import pyeuvac._misc as _m


class Euvac:
    '''
    EUVAC model class.
    '''
    def __init__(self):

        self._bands_dataset, self._lines_dataset, self._full_dataset = _m.get_euvac()

        self._full_f74113 = np.array(self._full_dataset['F74113'], dtype=np.float64).reshape(37, 1)
        self._full_ai = np.array(self._full_dataset['Ai'], dtype=np.float64).reshape(37, 1)

        self._bands_f74113 = np.array(self._bands_dataset['F74113'], dtype=np.float64).reshape(20,1)
        self._bands_ai = np.array(self._bands_dataset['Ai'], dtype=np.float64).reshape(20,1)

        self._lines_f74113 = np.array(self._lines_dataset['F74113'], dtype=np.float64).reshape(17,1)
        self._lines_ai = np.array(self._lines_dataset['Ai'], dtype=np.float64).reshape(17,1)

    def _get_p(self, i, f107, f107avg):
        '''
        Method for getting the P parameter. For each pair of values f107 and f107a, the value
        P = (f107 + f107 avg) / 2. - 80. is calculated, and then a matrix is constructed,
        each column of which is the value P.
        :param i: number of rows.
        :param f107: single value of the daily index F10.7 (in s.f.u.) or an array of such values.
        :param f107avg: a single value of the index F10.7 (in s.f.u.) averaged over 81 days or an array of such values.
        :return: numpy array.
        '''

        if f107.size != f107avg.size:
            raise Exception(f'The number of F10.7 and F10.7_avg values does not match. f107 contained {f107.size} '
                            f'elements, f107avg contained {f107avg.size} elements.')

        p = (f107 + f107avg) / 2. - 80.
        p_0 = p[:]
        for j in range(i-1):
            p_0 = np.vstack((p_0, p))

        return p_0

    def _predict(self, matrix_ai, matrix_f74113, vector_x, correction):
        pai = matrix_ai * vector_x + 1.0

        if correction:
            pai[pai < 0.8] = 0.8

        return matrix_f74113 * pai * 1.e9

    def _check_types(self, *proxies):
        if not all([isinstance(x, (float, int, list, np.ndarray)) for x in proxies]):
            raise TypeError(f'Only float, int, list and np.ndarray types are allowed. f107 was {type(proxies[0])}, '
                            f'f107avg was {type(proxies[1])}')
        return True

    def get_spectral_bands(self, *, f107, f107avg, correction=False):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 20 intervals
        of the spectrum of the interval 10-105 nm.
        :param f107: single value of the daily index F10.7 (in s.f.u.) or an array of such values.
        :param f107avg: a single value of the index F10.7 (in s.f.u.) averaged over 81 days or an array of such values.
        :param correction: parameter for applying flux factor correction. False by default.
        :return: xarray Dataset [euv_flux_spectra, lband, uband].
        '''

        bands = 20

        if self._check_types(f107, f107avg):

            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(bands,f107, f107avg)
            spectra = self._predict(self._bands_ai, self._bands_f74113, p, correction)

            res = np.zeros((spectra.shape[1], spectra.shape[1], spectra.shape[0]))
            for i in range(spectra.shape[1]):
                res[i, i, :] = spectra[:, i]

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), res),
                                         'lband': ('band_number', self._bands_dataset['lband'].values),
                                         'uband': ('band_number', self._bands_dataset['uband'].values)},
                              coords={'F107': f107,
                                      'F107AVG':  f107avg,
                                      'band_center': self._bands_dataset['center'].values,
                                      'band_number': np.arange(bands)})

    def get_spectral_lines(self, *, f107, f107avg, correction=False):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 17 lines
        of the spectrum of the interval 10-105 nm.
        :param f107: single value of the daily index F10.7 (in s.f.u.) or an array of such values.
        :param f107avg: a single value of the index F10.7 (in s.f.u.) averaged over 81 days or an array of such values.
        :param correction: parameter for applying flux factor correction. False by default.
        :return: xarray Dataset [euv_flux_spectra, wavelength].
        '''

        lines = 17

        if self._check_types(f107, f107avg):

            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(lines, f107, f107avg)
            spectra = self._predict(self._lines_ai, self._lines_f74113, p, correction)

            res = np.zeros((spectra.shape[1], spectra.shape[1], spectra.shape[0]))
            for i in range(spectra.shape[1]):
                res[i, i, :] = spectra[:, i]

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'line_wavelength'), res),
                                         'wavelength': ('line_number', self._lines_dataset['lambda'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'line_wavelength': self._lines_dataset['lambda'].values,
                                      'line_number': np.arange(lines)})

    def get_spectra(self, *, f107, f107avg, correction=False):
        '''
        Model calculation method. Combines the get_spectra_bands() and get_spectral_lines() methods.
        :param f107: single value of the daily index F10.7 (in s.f.u.) or an array of such values.
        :param f107avg: a single value of the index F10.7 (in s.f.u.) averaged over 81 days or an array of such values.
        :param correction: parameter for applying flux factor correction. False by default.
        :return: xarray Dataset [euv_flux_spectra, lband, uband], xarray Dataset [euv_flux_spectra, wavelength].
        '''

        return (self.get_spectral_bands(f107=f107, f107avg=f107avg, correction=correction),
                self.get_spectral_lines(f107=f107, f107avg=f107avg, correction=correction))

    def predict(self, f107, f107avg, correction=False):
        data = 37

        if self._check_types(f107, f107avg):

            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(data, f107, f107avg)
            spectra = self._predict(self._full_ai, self._full_f74113, p, correction)

            res = np.zeros((spectra.shape[1], spectra.shape[1], spectra.shape[0]))
            for i in range(spectra.shape[1]):
                res[i, i, :] = spectra[:, i]

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), res),
                                         'lband': ('band_number', self._full_dataset['lband'].values),
                                         'uband': ('band_number', self._full_dataset['uband'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'band_center': self._full_dataset['center'].values,
                                      'band_number': np.arange(data)})