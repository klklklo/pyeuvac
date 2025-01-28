import functools
import xarray as xr
from importlib_resources import files


@functools.cache
def _read_coeffs(file):
    return xr.open_dataset(files('pyeuvac._coeffs').joinpath(file))

def get_euvac():
    return (_read_coeffs('euvac_bands_coeffs.nc').copy(), _read_coeffs('euvac_lines_coeffs.nc').copy(),
            _read_coeffs('euvac_coeffs_full.nc').copy())
