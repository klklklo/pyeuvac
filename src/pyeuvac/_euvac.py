import numpy as np
import xarray as xr
import pyeuvac._misc as _m


class Euvac:
    '''
    EUVAC model class.
    '''
    def __init__(self):
        # The equation of the model Fi = F74113 * (1 + Ai * (P - 80)) <=> Fi = F74113 + F74113 * Ai * (P - 80)
        # In the form of a matrix product: F = (F74113 F74113*Ai) x (1 X)^T, where X = (P - 80)
        # Therefore _bands_coeffs and _lines_coeffs are represented by matrices (F74113 F74113*Ai)

        self._bands_dataset, self._lines_dataset = _m.get_euvac()
        self._bands_coeffs = np.vstack((np.array(self._bands_dataset['F74113'], dtype=np.float64),
                                        np.array(self._bands_dataset['F74113']) *
                                        np.array(self._bands_dataset['Ai'], dtype=np.float64))).transpose()

        self._lines_coeffs = np.vstack((np.array(self._lines_dataset['F74113'], dtype=np.float64),
                                        np.array(self._lines_dataset['F74113']) *
                                        np.array(self._lines_dataset['Ai'], dtype=np.float64))).transpose()

    def _get_p(self, f107, f107avg):
        '''
        Method for preparing data. It creates a two-dimensional array, the first column of which is filled with ones,
        the second with the values of P = (F10.7 + F10.7avg) / 2.
        :param f107: single value of the daily index F10.7 or an array of such values.
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values.
        :return: numpy-array for model calculation.
        '''

        f107 = np.array(f107, dtype=float)
        f107avg = np.array(f107avg, dtype=float)

        if f107.size != f107avg.size:
            raise Exception(f'The number of F10.7 and F10.7_avg values does not match. f107 contained {f107.size} '
                            f'elements, f107avg contained {f107avg.size} elements.')
        return np.vstack([[1., x] for x in (f107 + f107avg) / 2. - 80.])

    def _check_types(self, f107, f107avg):
        if not isinstance(f107, (float, int, list, np.ndarray) or not isinstance(f107avg, (float, int, list, np.ndarray))):
            raise TypeError(f'Only float, int, list and np.ndarray. f107 was {type(f107)}, f107avg was {type(f107avg)}')

        if type(f107) != type(f107avg):
            raise TypeError(f'f107 and f107avg types must be equal. f107 was {type(f107).__name__}, '
                            f'f107avg was {type(f107avg).__name__}')
        return True

    def get_spectral_bands(self, *, f107, f107avg):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 20 intervals
        of the spectrum of the interval 10-105 nm.
        :param f107: single value of the daily index F10.7 or an array of such values.
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values.
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center].
        '''

        if self._check_types(f107, f107avg):

            if isinstance(f107, (int, float)):
                f107 = np.array(f107, dtype=float).reshape(1,)
                f107avg = np.array(f107avg, dtype=float).reshape(1,)

            p = self._get_p(f107, f107avg)
            spectra = np.dot(self._bands_coeffs, p.T)

            res = np.zeros((spectra.shape[1], spectra.shape[1], spectra.shape[0]))
            for i in range(spectra.shape[1]):
                res[i, i, :] = spectra[:, i]

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F10.7', 'F10.7AVG', 'band_center'), res),
                                         'lband': ('band_number', self._bands_dataset['lband'].values),
                                         'uband': ('band_number', self._bands_dataset['uband'].values),
                                         'center': ('band_number', self._bands_dataset['center'].values)
                                         },
                              coords={'band_center': self._bands_dataset['center'].values,
                                      'F10.7': f107,
                                      'F10.7AVG':  f107avg,
                                      'band_number': np.arange(20)})

    def get_spectral_lines(self, *, f107, f107avg):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 17 lines
        of the spectrum of the interval 10-105 nm.
        :param f107: single value of the daily index F10.7 or an array of such values.
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values.
        :return: xarray Dataset [euv_flux_spectra, line_lambda].
        '''

        if self._check_types(f107, f107avg):

            if isinstance(f107, (int, float)):
                f107 = np.array(f107, dtype=float).reshape(1,)
                f107avg = np.array(f107avg, dtype=float).reshape(1,)

            p = self._get_p(f107, f107avg)
            spectra = np.dot(self._lines_coeffs, p.T)

            res = np.zeros((spectra.shape[1], spectra.shape[1], spectra.shape[0]))
            for i in range(spectra.shape[1]):
                res[i, i, :] = spectra[:, i]

            return xr.Dataset(data_vars={'euv_flux_spectra': (('F10.7', 'F10.7AVG', 'lambda'), res),
                                         'line_lambda': ('line_number', self._lines_dataset['lambda'].values)},
                              coords={'lambda': self._lines_dataset['lambda'].values,
                                      'line_number': np.arange(17),
                                      'F10.7': f107,
                                      'F10.7AVG': f107avg})
    def get_spectra(self, *, f107, f107avg):
        '''
        Model calculation method. Combines the get_spectra_bands() and get_spectral_lines() methods.
        :param f107: single value of the daily index F10.7 or an array of such values.
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values.
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center], xarray Dataset [euv_flux_spectra, line_lambda].
        '''

        return self.get_spectral_bands(f107=f107, f107avg=f107avg), self.get_spectral_lines(f107=f107, f107avg=f107avg)