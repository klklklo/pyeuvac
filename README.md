# pyeuvac
<!--Basic information-->
pyeuvac is a Python3 implementation of the extra ultraviolet (EUV) flux model described by P. G. Richards, 
J. A. Fennelly, D. G. Torr. This EUV model provides fluxes in the range 5-105 nm, divided into 20 intervals 
of 5 nm width and into 17 separate lines.

If you use pyeuvac or Euvac model directly or indirectly, please, cite in your research the following paper:

1. Richards, P. G., J. A. Fennelly, and D. G. Torr (1994), Euvac: A solar EUV Flux Model for aeronomic calculations, 
J. Geophys. Res., 99(A5), 8981-8992. https://doi.org/10.1029/94JA00518

# User's guide

<!--Users guide-->

## Installation

The following command is used to install the package:

```
python -m pip install pyeuvac
```

pyeuvac is the name of the package.

## Usage example

The pyeuvac package contains one class Euvac which has 3 methods: get_spectral_bands() for calculating the spectrum 
over intervals, get_spectral_lines() for calculating the spectrum along individual lines, and get_spectra() combining both methods.

1. get_spectral_bands()
```
# importing a package with the alias p
import pyeuvac as pe
# creating an instance of the Euvt2021 class
ex = pe.Euvac()
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200 (P = 200 as an example of the Richards et al.) using get_spectral_bands()
spectra = ex.get_spectral_bands((200,200))
# output the resulting EUV-spectra
print(spectra['euv_flux_spectra'])


<xarray.DataArray 'euv_flux_spectra' (band_center: 20, P: 1)> Size: 160B
array([[ 2.642448  ],
       [ 0.83475   ],
       ...
       [ 2.2567559 ],
       [ 3.762175  ]])
Coordinates:
  * band_center  (band_center) float64 160B 7.5 12.5 17.5 ... 92.5 97.5 102.5
  * P            (P) float64 8B 200.0
```

If you need to calculate the spectrum for several P values, pass them using a list:
```
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200, F10.7 = 190 and F10.7A = 210, F10.7 = 200 and F10.7A = 220
spectra = ex.get_spectral_bands([(200,200), (190,210), (200,220)])
# output the resulting EUV-spectra
print(spectra['euv_flux_spectra'])


<xarray.DataArray 'euv_flux_spectra' (band_center: 20, P: 3)> Size: 480B
array([[ 2.642448  ,  2.642448  ,  2.762652  ],
       [ 0.83475   ,  0.83475   ,  0.8668125 ],
       ...            ...          ...
       [ 2.2567559 ,  2.2567559 ,  2.32190223],
       [ 3.762175  ,  3.762175  ,  3.87010625]])
Coordinates:
  * band_center  (band_center) float64 160B 7.5 12.5 17.5 ... 92.5 97.5 102.5
  * P            (P) float64 24B 200.0 200.0 210.0
```

2. get_spectral_lines()

This method calculates the EUV flux in the range of 5-105 nm, divided into 20 intervals with a width of 5 nm each. 

Input parameters:
- f107 - single value of the daily index F10.7 (in s.f.u.);
- f107avg - 81-day average F10.7 value (in s.f.u.).

f107 and f107avg can be represented by lists for calculating spectra for several values of f107 and f107avg.
In this case, the lengths of these lists should be the same.

Output parameters:
- xarray dataset
```
<xarray.Dataset> Size: 492B
Dimensions:           (lambda: 17, ('F10.7', 'F10.7_avg'): 1, line_number: 17,
                       F10.7: 1, F10.7_avg: 1)
Coordinates:
  * lambda            (lambda) float64 136B 25.63 28.41 30.33 ... 102.6 103.2
  * line_number       (line_number) int32 68B 0 1 2 3 4 5 ... 11 12 13 14 15 16
  * F10.7             (F10.7) float64 8B <F10.7 values>
  * F10.7_avg         (F10.7_avg) float64 8B <F10.7 81-day average values>
Dimensions without coordinates: ('F10.7', 'F10.7_avg')
Data variables:
    euv_flux_spectra  (lambda, ('F10.7', 'F10.7_avg')) float64 136B <EUV flux spectra values>
    line_lambda       (line_number) float64 136B 25.63 28.41 ... 102.6 103.2
```

Below is an example of spectrum calculation using get_spectra_lines() method

```

```


If you need to calculate the spectrum for several P values, pass them using a list:
```
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200, F10.7 = 190 and F10.7A = 210, F10.7 = 200 and F10.7A = 220
spectra = ex.get_spectral_lines([(200,200), (190,210), (200,220)])
# output the resulting EUV-spectra
print(spectra['euv_flux_spectra'])


<xarray.DataArray 'euv_flux_spectra' (lambda: 17, P: 3)> Size: 408B
array([[0.61318   , 0.61318   , 0.625945  ],
       [3.679536  , 3.679536  , 3.968664  ],
       ...          ...         ...
       [5.676986  , 5.676986  , 5.8584015 ],
       [3.4313916 , 3.4313916 , 3.5423409 ]])
Coordinates:
  * lambda   (lambda) float64 136B 25.63 28.41 30.33 30.38 ... 97.7 102.6 103.2
  * P        (P) float64 24B 200.0 200.0 210.0

```

3. get_spectra()

This method combines the get_spectral_bands() and get_spectral_lines() methods. The method returns a tuple of 
xarray Dataset (lines, bands), the first element is the flux in intervals, the second is the flux in individual lines.

```
# importing a package with the alias p
import pyeuvac as pe
# creating an instance of the Euvt2021 class
ex = pe.Euvac()
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200 (P = 200 as an example of the Richards et al.) using get_spectra()
spectra = ex.get_spectra((200,200))
# output the resulting EUV-spectra
print(spectra)


(<xarray.Dataset> Size: 968B
Dimensions:           (band_center: 20, P: 1, band_number: 20)
Coordinates:
  * band_center       (band_center) float64 160B 7.5 12.5 17.5 ... 97.5 102.5
  * P                 (P) float64 8B 200.0
  * band_number       (band_number) int64 160B 0 1 2 3 4 5 ... 14 15 16 17 18 19
Data variables:
    euv_flux_spectra  (band_center, P) float64 160B 2.642 0.8348 ... 2.257 3.762
    lband             (band_number) float64 160B 5.0 10.0 15.0 ... 95.0 100.0
    uband             (band_number) float64 160B 10.0 15.0 20.0 ... 100.0 105.0
    center            (band_number) float64 160B 7.5 12.5 17.5 ... 97.5 102.5, <xarray.Dataset> Size: 280B
Dimensions:           (lambda: 17, P: 1)
Coordinates:
  * lambda            (lambda) float64 136B 25.63 28.41 30.33 ... 102.6 103.2
  * P                 (P) float64 8B 200.0
Data variables:
    euv_flux_spectra  (lambda, P) float64 136B 0.6132 3.68 3.2 ... 5.677 3.431)
```


