import functools
import xarray as xr
from importlib_resources import files


@functools.cache
def read_coeffs(file):
    return xr.open_dataset(files('pyeuvac._coeffs').joinpath(file))

def get_euvac():
    return read_coeffs('euvac_model_bands.nc').copy(), read_coeffs('euvac_model_lines.nc').copy(),
