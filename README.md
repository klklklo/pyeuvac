# pyeuvac
<!--Basic information-->
pyeuvac is a Python3 implementation of the extra ultraviolet (EUV) flux model EUVAC described by P. G. Richards, 
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

## Euvac

The pyeuvac package contains one class Euvac which has 3 methods:
- get_spectral_bands() for calculating the spectrum over intervals;
- get_spectral_lines() for calculating the spectrum along individual lines;
- get_spectra() for calculating the spectrum in a wavelength interval and in an individual wavelength.

All methods of the class have two input parameters:
- F<sub>10.7</sub> - daily value of the F<sub>10.7</sub> solar activity index (in s.f.u.);
- F<sub>10.7A</sub> - the average value of the F<sub>10.7</sub> solar activity index for 81 days.

Parameters can only be of types float, int, list and numpy.array. If there are several parameters, they are passed as 
lists with an equal number of elements.
Calculated EUV spectrum in units of cm<sup>-2</sup> s<sup>-1</sup>.

## Usage example

1. get_spectral_bands()

Method for calculating spectrum in 5 nm wide intervals from 5-105 nm range. Method returns xarray Dataset class object.

Input parameters:
- f107 - single value of the daily index F<sub>10.7</sub> (in s.f.u.);
- f107avg - 81-day average F<sub>10.7</sub> value (in s.f.u.).

```
<xarray.Dataset> Size: 736B
Dimensions:           (F10.7: 1, F10.7AVG: 1, band_center: 20, band_number: 20)
Coordinates:
  * F10.7             (F10.7) float64 8B <input F10.7 values>
  * F10.7AVG          (F10.7AVG) float64 8B <input F10.7A values>
  * band_center       (band_center) float64 160B 7.5 12.5 17.5 ... 97.5 102.5
  * band_number       (band_number) int32 80B 0 1 2 3 4 5 ... 14 15 16 17 18 19
Data variables:
    euv_flux_spectra  (F10.7, F10.7AVG, band_center) float64 160B <output spectrum>
    lband             (band_number) float64 160B 5.0 10.0 15.0 ... 95.0 100.0
    uband             (band_number) float64 160B 10.0 15.0 20.0 ... 100.0 105.0
```
The resulting spectrum is contained in a three-dimensional array with dimensions (F<sub>10.7</sub>, F<sub>10.7A</sub>, euv_spectrum)

Below is an example of calculating the spectrum with input parameters F<sub>10.7</sub> = F<sub>10.7A</sub> = 200 s.f.u.
```
# importing a package with the alias pe
import pyeuvac as pe
# creating an instance of the Euvac class
example = pe.Euvac()
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200 (P = 200 as an example of the Richards et al.) using get_spectral_bands()
spectrum = example.get_spectral_bands(f107=200., f107avg=200.)
# output the resulting EUV-spectra
print(spectrum['euv_flux_spectra'])


<xarray.DataArray 'euv_flux_spectra' (F10.7: 1, F10.7AVG: 1, band_center: 20)> Size: 160B
array([[[ 2.642448  ,  0.83475   , 12.504     , 10.3367336 ,
          7.01157116,  3.56471   ,  1.69090256,  0.72348547,
          0.976695  ,  0.92705019,  0.51372157,  0.826272  ,
          0.34776092,  0.22137   ,  1.19157903,  2.5642565 ,
          5.945697  ,  4.793988  ,  2.2567559 ,  3.762175  ]]])
Coordinates:
  * F10.7        (F10.7) float64 8B 200.0
  * F10.7AVG     (F10.7AVG) float64 8B 200.0
  * band_center  (band_center) float64 160B 7.5 12.5 17.5 ... 92.5 97.5 102.5
```

If you need to calculate the spectrum for several F<sub>10.7</sub> and F<sub>10.7A</sub> values, pass them using a list.
The number of values in the lists must be equal.
```
# calculate the spectrum values at F10.7 = [200., 210., 220.] and F10.7A = [200., 210., 220.]
spectrum = example.get_spectral_bands(f107=[200., 210., 220.], f107avg=[200., 210., 220.])
# output the resulting EUV-spectra for the first pair of F10.7 and F10.7A
print(spectrum['euv_flux_spectra'][0,0,:])


<xarray.DataArray 'euv_flux_spectra' (band_center: 20)> Size: 160B
array([ 2.642448  ,  0.83475   , 12.504     , 10.3367336 ,  7.01157116,
        3.56471   ,  1.69090256,  0.72348547,  0.976695  ,  0.92705019,
        0.51372157,  0.826272  ,  0.34776092,  0.22137   ,  1.19157903,
        2.5642565 ,  5.945697  ,  4.793988  ,  2.2567559 ,  3.762175  ])
Coordinates:
    F10.7        float64 8B 200.0
    F10.7AVG     float64 8B 200.0
  * band_center  (band_center) float64 160B 7.5 12.5 17.5 ... 92.5 97.5 102.5
```

2. get_spectral_lines()

Method for calculating spectrum in 17 separate lines from the range 5-105 nm. Method returns xarray Dataset class object.

Input parameters:
- f107 - single value of the daily index F10.7 (in s.f.u.);
- f107avg - 81-day average F10.7 value (in s.f.u.).

f107 and f107avg can be represented by lists for calculating spectra for several values of F<sub>10.7</sub> and F<sub>10.7A</sub>.
In this case, the lengths of these lists should be the same.

Output parameters:
- xarray dataset
``` 
<xarray.Dataset> Size: 492B
Dimensions:           (F10.7: 1, F10.7AVG: 1, line_wavelength: 17,
                       line_number: 17)
Coordinates:
  * F10.7             (F10.7) float64 8B <input F10.7 values>
  * F10.7AVG          (F10.7AVG) float64 8B <input F10.7A values>
  * line_wavelength   (line_wavelength) float64 136B 25.63 28.41 ... 102.6 103.2
  * line_number       (line_number) int32 68B 0 1 2 3 4 5 ... 11 12 13 14 15 16
Data variables:
    euv_flux_spectra  (F10.7, F10.7AVG, line_wavelength) float64 136B <output spectrum>
    wavelength        (line_number) float64 136B 25.63 28.41 ... 102.6 103.2
```

Below is an example of spectrum calculation using get_spectra_lines() method
```
# importing a package with the alias pe
import pyeuvac as pe
# creating an instance of the Euvac class
example = pe.Euvac()
# calculate the spectrum values at F10.7 = 200 and F10.7A = 200 (P = 200 as an example of the Richards et al.) using get_spectral_bands()
spectrum = example.get_spectral_lines(f107=200., f107avg=200.)
# output the resulting EUV-spectra
print(spectrum['euv_flux_spectra'])


<xarray.DataArray 'euv_flux_spectra' (F10.7: 1, F10.7AVG: 1, line_wavelength: 17)> Size: 136B
array([[[0.61318   , 3.679536  , 3.2       , 9.6599724 , 1.1641526 ,
         0.55071116, 1.00224288, 2.05612492, 1.55873   , 2.22441   ,
         0.49140144, 0.24854   , 0.6596096 , 0.977886  , 6.4812176 ,
         5.676986  , 3.4313916 ]]])
Coordinates:
  * F10.7            (F10.7) float64 8B 200.0
  * F10.7AVG         (F10.7AVG) float64 8B 200.0
  * line_wavelength  (line_wavelength) float64 136B 25.63 28.41 ... 102.6 103.2
```

If you need to calculate the spectrum for several P values, pass them using a list:
```
# calculate the spectrum values at F10.7 = [200., 210., 220.] and F10.7A = [200., 210., 220.]
spectrum = example.get_spectral_lines(f107=[200., 210., 220.], f107avg=[200., 210., 220.])
# output the resulting EUV-spectra for the first pair of F10.7 and F10.7A
print(spectrum['euv_flux_spectra'][0,0,:])


<xarray.DataArray 'euv_flux_spectra' (line_wavelength: 17)> Size: 136B
array([0.61318   , 3.679536  , 3.2       , 9.6599724 , 1.1641526 ,
       0.55071116, 1.00224288, 2.05612492, 1.55873   , 2.22441   ,
       0.49140144, 0.24854   , 0.6596096 , 0.977886  , 6.4812176 ,
       5.676986  , 3.4313916 ])
Coordinates:
    F10.7            float64 8B 200.0
    F10.7AVG         float64 8B 200.0
  * line_wavelength  (line_wavelength) float64 136B 25.63 28.41 ... 102.6 103.2
```

3. get_spectra()

This method combines the get_spectral_bands() and get_spectral_lines() methods. The method returns a tuple of 
xarray Dataset (lines, bands), the first element is the flux in intervals, the second is the flux in individual lines.
