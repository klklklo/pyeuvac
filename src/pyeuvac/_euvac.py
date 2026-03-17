import numpy as np
import xarray as xr
import sparse
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
        bands = 20
        cols = len(f107)

        if self._check_types(f107, f107avg):
            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(bands,f107, f107avg)
            spectra = self._predict(self._bands_ai, self._bands_f74113, p, correction).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, bands) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, bands) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(bands) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, bands))


            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), s)},
                              coords={'F107': f107,
                                      'F107AVG':  f107avg,
                                      'band_center': self._bands_dataset['center'].values,
                                      'band_number': np.arange(bands)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'lband': 'lower boundary of wavelength interval',
                                     'uband': 'upper boundary of wavelength interval'})

    def get_spectral_lines(self, *, f107, f107avg, correction=False):
        lines = 17
        cols = len(f107)

        if self._check_types(f107, f107avg):
            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(lines, f107, f107avg)
            spectra = self._predict(self._lines_ai, self._lines_f74113, p, correction).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, lines) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, lines) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(lines) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, lines))

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'line_wavelength'), s),
                                         'wavelength': ('line_number', self._lines_dataset['lambda'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'line_wavelength': self._lines_dataset['lambda'].values,
                                      'line_number': np.arange(lines)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'wavelength': 'the wavelength of a discrete line'})

    def get_spectra(self, *, f107, f107avg, correction=False):
        return (self.get_spectral_bands(f107=f107, f107avg=f107avg, correction=correction),
                self.get_spectral_lines(f107=f107, f107avg=f107avg, correction=correction))

    def predict(self, *, f107, f107avg, correction=False):
        data = 37
        cols = len(f107)

        if self._check_types(f107, f107avg):
            f107 = np.array([f107], dtype=np.float64) if isinstance(f107, (type(None), int, float)) \
                else np.array(f107, dtype=np.float64)
            f107avg = np.array([f107avg], dtype=np.float64) if isinstance(f107avg, (type(None), int, float)) \
                else np.array(f107avg, dtype=np.float64)

            p = self._get_p(data, f107, f107avg)
            spectra = self._predict(self._full_ai, self._full_f74113, p, correction).flatten(order='F')

            coords = np.array([np.hstack([x for x in [np.repeat(i, data) for i in range(cols)]]),
                               np.hstack([x for x in [np.repeat(i, data) for i in range(cols)]]),
                               np.hstack([x for x in [np.arange(data) for _ in range(cols)]])])

            s = sparse.COO(coords, spectra, shape=(cols, cols, data))

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F107', 'F107AVG', 'band_center'), s),
                                         'lband': ('band_number', self._full_dataset['lband'].values),
                                         'uband': ('band_number', self._full_dataset['uband'].values)},
                              coords={'F107': f107,
                                      'F107AVG': f107avg,
                                      'band_center': self._full_dataset['center'].values,
                                      'band_number': np.arange(data)},
                              attrs={'F10.7 units': '10^-22 · W · m^-2 · Hz^-1',
                                     'F10.7 81-day average units': '10^-22 · W · m^-2 · Hz^-1',
                                     'spectrum units': 'photons · m^-2 · s^-1',
                                     'wavelength units': 'nm',
                                     'euv_flux_spectra': 'modeled EUV solar irradiance',
                                     'lband': 'lower boundary of wavelength interval',
                                     'uband': 'upper boundary of wavelength interval'})
