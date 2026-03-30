import numpy as np
import xarray as xr
import sparse
import pyeuvac._misc as _m


class Euvac:
    '''
    EUVAC model class.
    '''
    def __init__(self):
        self._bands_coeffs, self._lines_coeffs, self._full_coeffs = _m.get_euvac_coeffs()

    @classmethod
    def scale_SI(cls, input):
        return input * 1.e4

    @classmethod
    def unscale_SI(cls, input):
        return input * 1.e-4

    def _get_p(self, i, f107, f107avg):
        if f107.size != f107avg.size:
            raise Exception(f'The number of F10.7 and F10.7_avg values does not match. f107 contained {f107.size} '
                            f'elements, f107avg contained {f107avg.size} elements.')

        p = (f107 + f107avg) / 2. - 80.
        p_0 = p[:]
        for j in range(i-1):
            p_0 = np.vstack((p_0, p))

        return p_0

    def _predict(self, matrix_ai, matrix_f74113, vector_x, autoscale_input=False):
        pai = matrix_ai * vector_x + 1.0
        pai[pai < 0.8] = 0.8

        return matrix_f74113 * pai * 1.e13 if autoscale_input else matrix_f74113 * pai * 1.e9


    def _check_types(self, *proxies):
        if not all([isinstance(x, (float, int, list, np.ndarray)) for x in proxies]):
            raise TypeError(f'Only float, int, list and np.ndarray types are allowed. f107 was {type(proxies[0])}, '
                            f'f107avg was {type(proxies[1])}')
        return True

    def get_spectral_bands(self, *, f107, f107avg, autoscale_input=False):
        if self._check_types(f107, f107avg):
            f107 = np.array(f107, dtype=np.float64).reshape(-1, )
            f107avg = np.array(f107avg, dtype=np.float64).reshape(-1, )

            bands = 20
            cols = len(f107)

            p = self._get_p(bands,f107, f107avg)
            f74113 = np.array(self._bands_coeffs['F74113'], dtype=np.float64).reshape(20, 1)
            ai = np.array(self._bands_coeffs['Ai'], dtype=np.float64).reshape(20, 1)
            spectra = self._predict(ai, f74113, p, autoscale_input).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, bands) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, bands) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(bands) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, bands))

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), s),
                                         'lband': ('band_number', self._bands_coeffs['lband'].data),
                                         'uband': ('band_number', self._bands_coeffs['uband'].data)},
                              coords={'F107': f107,
                                      'F107AVG':  f107avg,
                                      'band_center': self._bands_coeffs['center'].values,
                                      'band_number': np.arange(bands)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'lband': 'lower boundary of wavelength interval',
                                     'uband': 'upper boundary of wavelength interval'})

    def get_spectral_lines(self, *, f107, f107avg, autoscale_input=False):
        if self._check_types(f107, f107avg):
            f107 = np.array(f107, dtype=np.float64).reshape(-1, )
            f107avg = np.array(f107avg, dtype=np.float64).reshape(-1, )

            lines = 17
            cols = len(f107)

            p = self._get_p(lines, f107, f107avg)
            f74113 = np.array(self._lines_coeffs['F74113'], dtype=np.float64).reshape(17, 1)
            ai = np.array(self._lines_coeffs['Ai'], dtype=np.float64).reshape(17, 1)
            spectra = self._predict(ai, f74113, p, autoscale_input).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, lines) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, lines) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(lines) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, lines))

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'line_wavelength'), s),
                                         'wavelength': ('line_number', self._lines_coeffs['lambda'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'line_wavelength': self._lines_coeffs['lambda'].values,
                                      'line_number': np.arange(lines)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'wavelength': 'the wavelength of a discrete line'})

    def get_spectra(self, *, f107, f107avg, correction=False):
        return (self.get_spectral_bands(f107=f107, f107avg=f107avg, autoscale_input=correction),
                self.get_spectral_lines(f107=f107, f107avg=f107avg, autoscale_input=correction))

    def predict(self, *, f107, f107avg, correction=False):
        if self._check_types(f107, f107avg):
            f107 = np.array(f107, dtype=np.float64).reshape(-1, )
            f107avg = np.array(f107avg, dtype=np.float64).reshape(-1, )

            data = 37
            cols = len(f107)

            p = self._get_p(data, f107, f107avg)
            f74113 = np.array(self._full_coeffs['F74113'], dtype=np.float64).reshape(37, 1)
            ai = np.array(self._full_coeffs['Ai'], dtype=np.float64).reshape(37, 1)

            spectra = self._predict(ai, f74113, p, correction).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, data) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, data) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(data) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, data))

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), s),
                                         'lband': ('band_number', self._full_coeffs['lband'].values),
                                         'uband': ('band_number', self._full_coeffs['uband'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'band_center': self._full_coeffs['center'].values,
                                      'band_number': np.arange(data)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'lband': 'lower boundary of wavelength interval',
                                     'uband': 'upper boundary of wavelength interval'})
