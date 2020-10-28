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
* [CMake](https://cmake.org/)
* [SPICE](https://naif.jpl.nasa.gov/naif/)

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

The dependencies given above are required prior to build. The program is currently built and tested with the Intel C++ compiler for better compiler-directed vectorisation and optimization support. GNU compilers have been tested to work, but are not officially supported.

Any compiler used must support the C++11 standard. At the moment, the program itself is not fully compliant with the Standard but will make use of C++11 intrinsics once a preliminary version is implemented.

### Building

1. Clone the repo
```sh
git clone https://github.com/tylerjackoliver/ACROBAT.git
```
2. Run the CMake wrapper
```sh
sh build.sh
```

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
