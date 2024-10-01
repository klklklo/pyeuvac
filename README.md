# pynusinov
<!--Basic information-->
pynusinov is a Python3 implementation of models of the ultraviolet radiation spectra of the Sun described by A.A. Nusinov. 
The EUV model describes variations in the 5–105 nm spectral region, which are responsible for the ionization of the main components of the earth’s atmosphere.
The FUV model describes the flux changes in the 115–242 nm region, which determines heating of the upper atmosphere and the dissociation of molecular oxygen.
The input parameter for both models is the intensity of the photon flux in the Lyman-alpha line, which has been measured for decades. 
Using this parameter allows you to calculate solar radiation fluxes for any period of time.

If you use pynusinov or Nusinov's EUV/FUV models directly or indirectly, please, cite in your research the following paper:

1. Nusinov, A.A., Kazachevskaya, T.V., Katyushina, V.V. - Solar Extreme and Far Ultraviolet Radiation Modeling for Aeronomic
Calculations. Remote Sens. 2021, 13, 1454. https://doi.org/10.3390/rs13081454

# User's guide

<!--Users guide-->

## Installation

The following command is used to install the package:

```
python -m pip install pynusinov
```

pynusinov is the name of the package.

The package contains two classes: Euvt2021 and Fuvt2021.

## Fuvt2021

Implementation of the Nusinov model for calculating the spectrum of far ultraviolet radiation from the Sun (FUV)
in the wavelength range 115-242 nm. The model is based on the idea of a linear dependence of radiation fluxes in
1 nm wide intervals on the intensity in the Lyman-alpha hydrogen line (l = 121.6 nm).

Input parameters:
- flow in the Lyman-alpha line Nla (in units of 10^15 m^-2 * s^-1). You can set one or more Nla values.
Use a list to pass multiple values.

Output parameters:
- xarray dataset

```
<xarray.Dataset> Size: 6kB
Dimensions:                (band_center: 127, lyman_alpha_composite: 1, band_number: 127)
Coordinates:
  * band_center            (band_center) float64 1kB 115.5 116.5 ... 240.5 241.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B <Nl input values>
  * band_number            (band_number) int64 1kB 0 1 2 3 4 ... 123 124 125 126
Data variables:
    fuv_flux_spectra       (band_center, lyman_alpha_composite) float64 1kB 1...
    lband                  (band_number) int64 1kB 115 116 117 ... 239 240 241
    uband                  (band_number) int64 1kB 116 117 118 ... 240 241 242
    fuv_band_width         (band_number) float64 1kB 1.0 1.0 1.0 ... 1.0 1.0 1.0
```

## Usage example

- import the pynusinov package;
- create an instance of the Fuvt2021 class;
- perform calculations with the created instance.

The following is an example of performing the described steps:

```
>>> # importing a package with the alias p
>>> import pynusinov as p
>>> # creating an instance of the Fuvt2021 class
>>> ex = p.Fuvt2021()
>>> # calculate the spectra values at Nla = 3.31 (10^15) using get_spectra()
>>> spectra = ex.get_spectra(3.31)
>>> # output the resulting FUV-spectra
>>> print(spectra['fuv_flux_spectra'])

<xarray.DataArray 'fuv_flux_spectra' (band_center: 127, lyman_alpha_composite: 1)> Size: 1kB
array([[1.0226240e+13],
       [1.3365010e+13],
...
       [4.5222314e+16],
       [5.3300029e+16]])
Coordinates:
  * band_center            (band_center) float64 1kB 115.5 116.5 ... 240.5 241.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B 3.31
```

If you need to calculate the spectrum for several Na values, pass them using a list:

```
>>> # calculate the spectrum values at Nl_1 = 3.31 (10^15) and Nl_2 = 7.12 (10^15) using get_spectra()
>>> spectra = ex.get_spectra([3.31, 7.12])
>>> # output the resulting FUV-spectrum
>>> print(spectra['fuv_flux_spectra'])

<xarray.DataArray 'fuv_flux_spectra' (band_center: 127, lyman_alpha_composite: 2)> Size: 2kB
array([[1.0226240e+13, 1.7099480e+13],
       [1.3365010e+13, 1.7826520e+13],
...
       [4.5222314e+16, 4.7239328e+16],
       [5.3300029e+16, 5.5418008e+16]])
Coordinates:
  * band_center            (band_center) float64 1kB 115.5 116.5 ... 240.5 241.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 16B 3.31 7.12
```

## Euvt2021

Implementation of the Nusinov model for calculating the spectra of the extreme ultraviolet radiation of the Sun (EUV)
in the wavelength range of 10-105 nm. This model calculates the ultraviolet spectra for an individual wavelength or 
a wavelength interval. The model is based on the idea of a linear dependence of radiation fluxes in intervals
of unequal width on the intensity in the HeI helium line (l = 58.4 nm). 

Input parameters:
- the flow in the HeI line Nl (in units of 10^15 m^-2 * s^-1)

Output parameters:
- xarray dataset

For calculations of the model by interval wavelength and by wavelength interval xarray is different:

```
# wavelength interval
<xarray.Dataset> Size: 968B
Dimensions:                (band_center: 20, lyman_alpha_composite: 1,
                            band_number: 20)
Coordinates:
  * band_center            (band_center) float64 160B 7.5 12.5 ... 97.5 102.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B <Nl input values>
  * band_number            (band_number) int64 160B 0 1 2 3 4 ... 15 16 17 18 19
Data variables:
    euv_flux_spectra       (band_center, lyman_alpha_composite) float64 160B ...
    lband                  (band_number) int64 160B 5 10 15 20 ... 85 90 95 100
    uband                  (band_number) int64 160B 10 15 20 25 ... 95 100 105
    center                 (band_number) float64 160B 7.5 12.5 ... 97.5 102.5


# wavelength line
<xarray.Dataset> Size: 264B
Dimensions:                (line: 16, lyman_alpha_composite: 1)
Coordinates:
  * line                   (line) float64 128B 25.6 28.4 30.4 ... 102.6 103.2
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B <Nl input values>
Data variables:
    euv_flux_spectra       (line, lyman_alpha_composite) float64 128B ...
```

## Usage example

This class contains two methods for calculating the spectrum:
- get_spectral_bands() for calculating the spectrum in a wavelength interval;
- get_spectral_lines() for calculating the spectrum for an individual wavelength.

The steps of work are similar to the steps described for the Fuvt2021 class. 

Below is an example of working with the Fuvt2021 class:

1. get_spectral_lines()
```
>>> # importing a package with the alias p
>>> import pynusinov as p
>>> # creating an instance of the Euvt2021 class
>>> ex = p.Euvt2021()
>>> # calculate the spectrum values at Nl = 3.31 (10^15) using get_spectral_lines()
>>> spectra = ex.get_spectral_lines(3.31)
>>> # output the resulting EUV-spectra
>>> print(spectra['euv_flux_spectra'])

<xarray.DataArray 'euv_flux_spectra' (line: 16, lyman_alpha_composite: 1)> Size: 128B
array([[ 1.07475700e+13],
       [-3.48013400e+11],
...   
       [ 3.01426805e+13],
       [ 5.22986620e+12]])
Coordinates:
  * line                   (line) float64 128B 25.6 28.4 30.4 ... 102.6 103.2
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B 3.31
```

If you need to calculate the spectrum for several Na values, pass them using a list:

```
>>> # calculate the spectrum values at Nl_1 = 3.31 (10^15) and Nl_2 = 7.12 (10^15) using get_spectral_lines()
>>> spectra = ex.get_spectral_lines([3.31, 7.12])
>>> # output the resulting EUV-spectrum
>>> print(spectra['euv_flux_spectra'])

<xarray.DataArray 'euv_flux_spectra' (line: 16, lyman_alpha_composite: 2)> Size: 256B
array([[ 1.07475700e+13,  6.92348800e+13],
       [-3.48013400e+11,  1.29777664e+13],
...
       [ 3.01426805e+13,  9.21014720e+13],
       [ 5.22986620e+12,  1.51018048e+13]])
Coordinates:
  * line                   (line) float64 128B 25.6 28.4 30.4 ... 102.6 103.2
  * lyman_alpha_composite  (lyman_alpha_composite) float64 16B 3.31 7.12
```

2. get_spectral_bands()
```
>>> # importing a package with the alias p
>>> import pynusinov as p
>>> # creating an instance of the Euvt2021 class
>>> ex = p.Euvt2021()
>>> # calculate the spectrum values at Nl = 3.31 (10^15) using get_spectral_bands()
>>> spectra = ex.get_spectral_bands(3.31)
>>> # output the resulting EUV-spectra
>>> print(spectra['euv_flux_spectra'])

<xarray.DataArray 'euv_flux_spectra' (band_center: 20, lyman_alpha_composite: 1)> Size: 160B
array([[2.52122700e+12],
       [2.59186240e+12],
...
       [5.73289352e+13],
       [9.57620734e+13]])
Coordinates:
  * band_center            (band_center) float64 160B 7.5 12.5 ... 97.5 102.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 8B 3.31
```

If you need to calculate the spectrum for several Na values, pass them using a list:

```
>>> # calculate the spectrum values at Nl_1 = 3.31 (10^15) and Nl_2 = 7.12 (10^15) using get_spectral_bands()
>>> spectra = ex.get_spectral_bands([3.31, 7.12])
>>> # output the resulting EUV-spectrum
>>> print(spectra['euv_flux_spectra'])

<xarray.DataArray 'euv_flux_spectra' (band_center: 20, lyman_alpha_composite: 2)> Size: 320B
array([[2.52122700e+12, 3.44494080e+13],
       [2.59186240e+12, 2.14175296e+13],
...
       [5.73289352e+13, 1.07909581e+14],
       [9.57620734e+13, 2.62794074e+14]])
Coordinates:
  * band_center            (band_center) float64 160B 7.5 12.5 ... 97.5 102.5
  * lyman_alpha_composite  (lyman_alpha_composite) float64 16B 3.31 7.12

```



