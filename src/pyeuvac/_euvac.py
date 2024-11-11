import numpy as np
import xarray as xr
import pyeuvac._misc as _m


class Euvac:
    '''
    Euvac model class. Wavelength range 5-105 nm
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
        the second with the values of P = (F10.7 + F10.7avg) / 2
        :param f107: single value of the daily index F10.7 or an array of such values
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values
        :return: numpy-array for model calculation
        '''

        if isinstance(f107, (float, int)) and isinstance(f107avg, (float, int)):
            return np.array([1., (f107 + f107avg) / 2. - 80.])

        if isinstance(f107, (float, int)) and not isinstance(f107avg, (float, int)) \
                or not isinstance(f107, (float, int)) and isinstance(f107avg, (float, int)):
            raise TypeError(f'F10.7 and F10.7 must be same type. f107 was {type(f107)}, f107avg was {type(f107avg)}')
        f107 = np.array(f107)
        f107avg = np.array(f107avg)
        if len(f107) != len(f107avg):
            raise Exception(f'The number of F10.7 and F10.7_avg values does not match. f107 contained {len(f107)} '
                            f'elements, f107avg contained {len(f107avg)} elements ')
        return np.vstack([[1., x] for x in (f107 + f107avg) / 2. - 80.])

    def get_spectral_bands(self, f107, f107avg):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 20 intervals
        of the spectrum of the interval 10-105 nm
        :param f107: single value of the daily index F10.7 or an array of such values
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center]
        '''

        p = self._get_p(f107, f107avg)
        res = np.dot(self._bands_coeffs, p.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('band_center', ('F10.7', 'F10.7_avg')), res),
                                     'lband': ('band_number', self._bands_dataset['lband'].values),
                                     'uband': ('band_number', self._bands_dataset['uband'].values),
                                     'center': ('band_number', self._bands_dataset['center'].values)},
                          coords={'band_center': self._bands_dataset['center'].values,
                                  'F10.7': f107,
                                  'F10.7_avg': f107avg,
                                  'band_number': np.arange(20)})

    def get_spectral_lines(self, f107, f107avg):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 17 lines
        of the spectrum of the interval 10-105 nm
        :param f107: single value of the daily index F10.7 or an array of such values
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values
        :return: xarray Dataset [euv_flux_spectra, line_lambda]
        '''

        x = self._get_p(f107, f107avg)
        res = np.dot(self._lines_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('lambda', ('F10.7', 'F10.7_avg')), res),
                                     'line_lambda': ('line_number', self._lines_dataset['lambda'].values)},
                          coords={'lambda': self._lines_dataset['lambda'].values,
                                  'line_number': np.arange(17),
                                  'F10.7': f107,
                                  'F10.7_avg': f107avg,})

    def get_spectra(self, f107, f107avg):
        '''
        Model calculation method. Combines the get_spectra_bands() and get_spectral_lines() methods
        :param f107: single value of the daily index F10.7 or an array of such values
        :param f107avg: a single value of the F10.7 index averaged over 81 days or an array of such values
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center], xarray Dataset [euv_flux_spectra, line_lambda]
        '''
        return self.get_spectral_bands(f107, f107avg), self.get_spectral_lines(f107, f107avg)
