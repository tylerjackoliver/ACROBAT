# ACROBAT

[![CodeFactor](https://www.codefactor.io/repository/github/tylerjackoliver/acrobat/badge?s=24e1157e0259c7ca0a3a435459b55a10881e7481)](https://www.codefactor.io/repository/github/tylerjackoliver/acrobat)

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>
  <h3 align="center">ACROBAT: bAllistic CaptuRe OrBit Analysis Tool</h3>

  <p align="center">
    This program identifies sets of points in the Elliptic-Restricted Three-body Problem that lead to temporary ballistic capture, as per Z.F. Luo; F. Topputo; 2016.
    <br />
    <a href="https://github.com/tylerjackoliver/ACROBAT/issues">Report Bug</a>
    ·
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

As part of ongoing research into ballistic capture, this program is designed to recreate results on identifying ballistic capture sets for generic systems in the Elliptic-Restricted Three-body Problem (ER3BP).
The program is designed to mirror the latest research on ballistic capture by Luo & Topputo, where the equations of motion are studied in a planet-centered inertial frame under the general equations of n-body motion for three bodies;
 this formulation makes the identification of stopping criteria far easier than traditional approaches in Polar coordinates under Levi-Civita regularisation found in e.g. Topputo, 2012.

This program leverages CPU parallelism to reduce program run-time.

### Built With

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [BOOST](https://www.boost.org/)
* [OpenMP](https://www.openmp.org/)
* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)
* [CMake](https://cmake.org/)

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

Your Name - [@tylerjackoliver](https://twitter.com/tylerjackoliver) - jack.tyler@soton.ac.uk

Project Link: [https://github.com/tylerjackoliver/ACROBAT](https://github.com/tylerjackoliver/ACROBAT)

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/tylerjackoliver/repo.svg?style=flat-square
[contributors-url]: https://github.com/tylerjackoliver/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/tylerjackoliver/repo.svg?style=flat-square
[forks-url]: https://github.com/tylerjackoliver/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/tylerjackoliver/repo.svg?style=flat-square
[stars-url]: https://github.com/tylerjackoliver/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/tylerjackoliver/repo.svg?style=flat-square
[issues-url]: https://github.com/tylerjackoliver/repo/issues
[license-shield]: https://img.shields.io/github/license/tylerjackoliver/repo.svg?style=flat-square
[license-url]: https://github.com/tylerjackoliver/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/tylerjackoliver
[product-screenshot]: images/screenshot.png
