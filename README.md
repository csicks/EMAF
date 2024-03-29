# EMAF: Fast Cryo-EM Image Alignment Algorithm using Power Spectrum Features
[![CMake](https://github.com/csicks/EMAF/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/csicks/EMAF/actions/workflows/cmake.yml)
[![DOI](https://img.shields.io/badge/doi-10.1021/acs.jcim.1c00745-blue.svg)](https://doi.org/10.1021/acs.jcim.1c00745)

## Introduction
- EMAF is an alignment algorithm program which aims at cryo-EM two-dimensional particle images. The program is implemented in C++ under C++14 standard. Some other alignment methods like Fourier-Mellin transform and XMIPP alignment are also re-implemented and included in this program. Both source code and complied program are available.
- This program has been tested on Ubuntu 14.04, Ubuntu 16.04 and Ubuntu 18.04. Other Linux system should also be able to run this program. For Windows, you may be able to use this program after installation of FFTW.

## Dependencies
- Since most utilities are implemented in our own source code, the only dependent package is FFTW. You may simply install FFTW3 using `sudo apt-get install libfftw3-dev` in ubuntu command line or compiling from source code following instructions on the [official website](http://fftw.org/).

## Usage
- You may directly use the complied program `EMAF` in Linux commandline. Please input `./EMAF -h` to see commandline instructions. If some error occurs like GLIBC error, you may install the required libraries or consider compile from source code following the below instructions.
- You may compile from source code using the following command in Linux system. CMake version 3.10 and above have been tested. You may change the [CMake File](./CMakeLists.txt) if other CMake version are used.

  ```bash
  cd /path/to/README
  mkdir build
  cd build
  cmake ..
  make
  ```

- Another possible way to use this algorithm is directly running code from main function in source code. This way may be suitable for people who would like to understand the code. Proper comments are provided in source code.
- By now, only MRC format and certain XMD format are supported as input file. XMD format is used for XMIPP relied programs and only tested under certain cases. Therefore, data is strongly recommended to be converted to MRC format when it is feed into this program.

## Suggestions
- When automatic estimate of filter parameters is not accurate, manual fine tuning may help.
- All particle images used in experiments are square images which means their heights equal their widths. Define the height of image as `s`.
  1. For synthetic GroEL datasets, EMAF-MD with low pass filter width `s/10` is suggested.
  2. For synthetic EMPIAR 10065 datasets, EMAF-MD with low pass filter width `s/20` is suggested.
  3. For synthetic EMPIAR 10398 datasets, EMAF-MD with low pass filter width `s/20` is suggested.
  4. For synthetic Apoferritin datasets, EMAF-MD with low pass filter width `s/15` is suggested.
  5. For real GroEL datasets, EMAF-MD with low pass filter width `s/10` is suggested.
  6. For real EMPIAR 10398 datasets, EMAF-MD with low pass filter width `s/20` is suggested.

## License
- EMAF is licensed under the GNU General Public License v3.0; you may not use this file except in compliance with the License. See file LICENSE for details.
- All materials are made available under the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License (CC BY-NC 4.0) license. You can find details at: https://creativecommons.org/licenses/by-nc/4.0/legalcode 

## Others
- Please refer to http://www.csbio.sjtu.edu.cn/bioinf/EMAF/ for used datasets and more information.
- All code is provided for research purpose and without any warranty. If you use code or ideas from this project for your research, please cite our paper:
  ```bibtex
  @Article{Chen2021a,
  author    = {Chen, Yu-Xuan and Xie, Rui and Yang, Yang and He, Lin and Feng, Dagan and Shen, Hong-Bin},
  journal   = {J. Chem. Inf. Model.},
  title     = {Fast Cryo-EM Image Alignment Algorithm Using Power Spectrum Features},
  year      = {2021},
  issn      = {1549-9596},
  month     = sep,
  number    = {9},
  pages     = {4795--4806},
  volume    = {61},
  comment   = {doi: 10.1021/acs.jcim.1c00745},
  doi       = {10.1021/acs.jcim.1c00745},
  publisher = {American Chemical Society},
  url       = {https://doi.org/10.1021/acs.jcim.1c00745},
  }
  ```
