import numpy as np
import pandas as pd
import xarray as xr
import pyeuvac._misc as _m

class EUVAC:
    def __init__(self):
        self._bands_dataset, self._lines_dataset = _m.get_euvac()
        self._bands_coeffs = np.vstack((np.array(self._bands_dataset['F74113'], dtype=np.float64),
                                        np.array(self._bands_dataset['F74113']) * np.array(self._bands_dataset['Ai'], dtype=np.float64))).transpose()

        self._lines_coeffs = np.vstack((np.array(self._lines_dataset['F74113'], dtype=np.float64),
                                        np.array(self._lines_dataset['F74113']) * np.array(self._lines_dataset['Ai'], dtype=np.float64))).transpose()

    def _get_P(self, list_of_F):
        if isinstance(list_of_F, tuple):
            array = np.array([1., sum(list_of_F)/2. - 80], dtype=np.float64)
            return array.reshape((1, 2))


        tmp = np.array([sum(i) / 2. for i in list_of_F], dtype=np.float64)
        tmp = tmp.reshape((tmp.size, 1))
        array = np.ones((tmp.size, 1), dtype=np.float64)
        return tmp, np.hstack([array, tmp-80])

    def get_spectra_bands(self, P):
        p, x = self._get_P(P)
        res = np.dot(self._bands_coeffs, x.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('band_center', 'P'), res),
                                     'lband': ('band_number', self._bands_dataset['lband'].values),
                                     'uband': ('band_number', self._bands_dataset['uband'].values),
                                     'center': ('band_number', self._bands_dataset['center'].values)},
                          coords={'band_center': self._bands_dataset['center'].values,
                                  'P': p[:, 0],
                                  'band_number': np.arange(20)})

    def get_spectra_lines(self, P):
        p, x = self._get_P(P)
        res = np.dot(self._lines_coeffs, p.T)
        return xr.Dataset(data_vars={'euv_flux_spectra': (('lambda', 'P'), res)},
                          coords={'lambda': self._lines_dataset['lambda'].values,
                                  'P': p[:, 0],})

    def get_spectra(self, P):
        return self.get_spectra_bands(P), self.get_spectra_lines(P)
