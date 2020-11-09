# ACROBAT


[![CodeFactor](https://www.codefactor.io/repository/github/tylerjackoliver/acrobat/badge?s=24e1157e0259c7ca0a3a435459b55a10881e7481)](https://www.codefactor.io/repository/github/tylerjackoliver/acrobat)
[![Contributors][contributors-shield]][contributors-url]
[![Issues](https://img.shields.io/github/issues/tylerjackoliver/ACROBAT)](https://github.com/tylerjackoliver/ACROBAT/issues)
![License](https://img.shields.io/github/license/tylerjackoliver/ACROBAT)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/tylerjackoliver/ACROBAT">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>
  <h3 align="center">ACROBAT: bAllistic CaptuRe OrBit Analysis Tool</h3>

  <p align="center">
    This program identifies sets of points in the Elliptic-Restricted Three-body Problem that lead to temporary ballistic capture, as per <a href="https://doi.org/10.1007/s10569-014-9580-5">Z.F. Luo et. al., 2014</a>.
    <br />
    <a href="https://github.com/tylerjackoliver/ACROBAT/issues">Report Bug</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#building)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

As part of ongoing research into ballistic capture, this program is designed to recreate results on identifying ballistic capture sets for generic systems in the Elliptic-Restricted Three-body Problem (ER3BP).
The program is designed to mirror the latest research on ballistic capture by Luo & Topputo, where the equations of motion are studied in a planet-centered inertial frame under the general equations of n-body motion for three bodies;
 this formulation makes the identification of stopping criteria far easier than traditional approaches in Polar coordinates under Levi-Civita regularisation found in e.g. Topputo, 2012.

This program leverages CPU parallelism to reduce program run-time.

### Changes from the previous literature

The previous literature used the [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009](https://doi.org/10.1007/s10569-010-9320-4) to obtain Right Ascension and Declination angles for the definition of the Body Mean Equator frame at a given epoch. This is labour-intensive, as often these are expansions of multiple nutation and precession angles and are unique for each target body. Instead, this code interfaces directly with the SPICE library to extract this information from binary PCK files that are published regularly by NAIF.

Specifically, the function `bods2c_c` ([here](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bods2c_c.html)) is used to convert from a body name or identification string to the integer ID representation contained in the given binary PCK file. The `bodeul_` function (F2C'd as part of the SPICE library and undocumented) can then be used to determine RA and DEC to the accuracy contained in the SPK files. The likelihood of mistakes is reduced, at the cost of program execution time being slightly -- _slightly_ -- higher. The FORTRAN documentation for the equivalent `BODEUL` function is available [here](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/bodeul.html).

### Built With

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [BOOST](https://www.boost.org/)
* [OpenMP](https://www.openmp.org/)
* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)
* [Intel MPI](https://software.intel.com/content/www/us/en/develop/tools/mpi-library.html), although any MPI implementation that supports the Intel C++ compiler will do.
* [CMake](https://cmake.org/)
* [SPICE](https://naif.jpl.nasa.gov/naif/)

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The dependencies given above are required prior to build. The program is currently built and tested with the Intel C++ compiler for better compiler-directed vectorisation and optimization support. GNU compilers have been tested to work, but are not officially supported.

Any compiler used must support the C++14 standard. At the moment, the program itself is not fully compliant with the Standard but will make use of C++14 intrinsics once a preliminary version is implemented.

### Building

1. Clone the repo
```sh
git clone https://github.com/tylerjackoliver/ACROBAT.git
```
2. Run the CMake wrapper
```sh
sh build.sh
```

### Running
The user should edit the parameters for their desired configuration in `Params.hpp`. Six unique parameters are required to instantiate all of the derived parameters for the problem:

* `HOST`, a `std::string` object containing either the NAIF designation (e.g. `"399"`) or common identifier (`"Earth"`) for the HOST body, i.e. the major body in the ER3BP.
* `TARGET`, a `std::string` object containing either the NAIF designation (e.g. `"499"`) or common identifier (`"Mars"`) for the TARGET body, i.e. the minor body in the ER3BP.
* `ECC`, a double-precision number containing the desired eccentricity of the capture orbits, [0.0, 1.0).
* `INC`, a double-precision number containing the inclination of the desired capture orbits, [0, 2 * pi).
* `LONGTD`, a double-precision number containing the longitude of the desired capture orbits, [0, 2 * pi).
* `EPOCH`, a double-precision number containing the epoch of the investigation, in ephemeris seconds past J2000. For abstract dates, the `str2et_c()` function may be used in the Parameters file.

The following derived parameters are then calculated from the values given above using the SPICE system. Note that this relies on the use of planetary constant kernels (PCKs), lightsecond kernels (LSKs) and planetary ephemerides (SPKs). An instruction on including these kernels may be found in the `spice/` subdirectory.

* `R`, the mean planetary radius of the `TARGET`.
* `RS`, the sphere of influence of the `TARGET` about the `HOST`, assuming a perfectly spherical SOI.
* `targetGM`, the gravitational parameter of the  `TARGET`.
* `hostGM`, the gravitational parameter of the `HOST`.
* `M`, the mean anomaly of the `TARGET` about the `HOST` at the given `EPOCH`.

<!-- CONTACT -->
## Contact

Jack Tyler - [@tylerjackoliver](https://twitter.com/tylerjackoliver) - jack.tyler@soton.ac.uk

Project Link: [https://github.com/tylerjackoliver/ACROBAT](https://github.com/tylerjackoliver/ACROBAT)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/tylerjackoliver/ACROBAT
[contributors-url]: https://img.shields.io/github/contributors/tylerjackoliver/ACROBAT
[issues-shield]: https://github.com/tylerjackoliver/ACROBAT/issues
[issues-url]: https://img.shields.io/github/issues/tylerjackoliver/ACROBAT
[license-shield]: ""
[license-url]: https://img.shields.io/github/license/tylerjackoliver/ACROBAT
[product-screenshot]: images/screenshot.png
