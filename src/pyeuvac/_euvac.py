import numpy as np
import xarray as xr
import pyeuvac._misc as _m

class EUVAC:
    '''
    EUVAC model class. Wavelength range 5-105 nm
    '''
    def __init__(self):
        # The equation of the model Fi = F74113 * (1 + Ai * (P - 80)) <=> Fi = F74113 + F74113 * Ai * (P - 80)
        # In the form of a matrix product: F = (F74113 F74113*Ai) x (1 X)^T, where X = (P - 80)
        # Therefore _bands_coeffs and _lines_coeffs are represented by matrices (F74113 F74113*Ai)

        self._bands_dataset, self._lines_dataset = _m.get_euvac()
        self._bands_coeffs = np.vstack((np.array(self._bands_dataset['F74113'], dtype=np.float64),
                                        np.array(self._bands_dataset['F74113']) * np.array(self._bands_dataset['Ai'], dtype=np.float64))).transpose()

        self._lines_coeffs = np.vstack((np.array(self._lines_dataset['F74113'], dtype=np.float64),
                                        np.array(self._lines_dataset['F74113']) * np.array(self._lines_dataset['Ai'], dtype=np.float64))).transpose()

    def _get_P(self, list_of_F):
        '''
        Method for preparing data. It creates a two-dimensional array, the first column of which is filled with ones,
        the second with the values of P = (F10.7 + F10.7A) / 2
        :param list_of_F: tuple (F10.7, F10.7A) or list of these tuples
        :return: numpy-array for model calculation
        '''

        if isinstance(list_of_F, tuple):
            P = np.array(sum(list_of_F)/2.)
            return np.array([1., P-80])[None, :]

        tmp = np.array([sum(i) / 2. for i in list_of_F], dtype=np.float64)
        tmp = tmp.reshape((tmp.size, 1))
        array = np.ones((tmp.size, 1), dtype=np.float64)
        return np.hstack([array, tmp-80])

    def get_spectra_bands(self, P):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 20 intervals
        of the spectrum of the interval 10-105 nm
        :param P: tuple (F10.7, F10.7A) or list of these tuples
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center]
        '''
        x = self._get_P(P)
        res = np.dot(self._bands_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('band_center', 'P'), res),
                                     'lband': ('band_number', self._bands_dataset['lband'].values),
                                     'uband': ('band_number', self._bands_dataset['uband'].values),
                                     'center': ('band_number', self._bands_dataset['center'].values)},
                          coords={'band_center': self._bands_dataset['center'].values,
                                  'P': x[:, 1] + 80,
                                  'band_number': np.arange(20)})

    def get_spectra_lines(self, P):
        '''
        Model calculation method. Returns the values of radiation fluxes in all 17 lines
        of the spectrum of the interval 10-105 nm
        :param P: tuple (F10.7, F10.7A) or list of these tuples
        :return: xarray Dataset [euv_flux_spectra]
        '''
        x = self._get_P(P)
        res = np.dot(self._lines_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('lambda', 'P'), res)},
                          coords={'lambda': self._lines_dataset['lambda'].values,
                                  'P': x[:, 1]+80})

    def get_spectra(self, P):
        '''
        Model calculation method. Combines the get_spectra_bands() and get_spectra_lines() methods
        :param P: tuple (F10.7, F10.7A) or list of these tuples
        :return: xarray Dataset [euv_flux_spectra, lband, uband, center], xarray Dataset [euv_flux_spectra]
        '''
        return self.get_spectra_bands(P), self.get_spectra_lines(P)
