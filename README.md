# EMAF: A Fast Cryo-EM Image Alignment Algorithm using Power Spectrum Features

## Introduction
- EMAF is an alignment algorithm program which aims at cryo-EM two-dimensional particle images. The program is implemented in C++ under C++11 standard. Some other alignment methods like Fourier-Mellin transform and XMIPP alignment are also re-implemented and included in this program. Both source code and complied program are available.
- This program has been tested on  Ubuntu 14.04, Ubuntu 16.04 and Ubuntu 18.04. Other Linux system should also be able to run this program. For Windows, you may be able to use this program after installation of FFTW.

## Dependencies
- Since most utilities are implemented in our own source code, the only dependent package is FFTW. You may simply install FFTW3 using "sudo apt-get install libfftw3-dev" in ubuntu command line or compiling from source code following instructions on the official website http://fftw.org/.

## Usage
- You may directly use the complied program "EMAF" in Linux commandline. Please input "./EMAF -h" to see commandline instructions. If some error occurs like GLIBC error, you may install the required libraries or consider compile from source code following the below instructions.
- You may compile from source code using the following command in Linux system. CMake version 3.10 and above have been tested. You may change the CMakeLists.txt if other CMake version are used.
  1. cd /path/to/this README.md
  2. mkdir build
  3. cd build
  4. cmake ..
  5. make
- Another possible way to use this algorithm is directly running code from main function in source code. This way may be suitable for people who would like to understand the code. Proper comments are provided in source code.
- This program temporally uses an INT variable to store the whole size of data (number_of_images * width * height). Then it might crash if the data size is too big. Consider using the BIGDATA version if data size is very large.
- By now, only MRC format and certain XMD format are supported as input file. XMD format is used for XMIPP relied programs and only tested under certain cases. Therefore, data is strongly recommended to be converted to MRC format when it is feed into this program.

## Suggestions
- When automatic estimate of filter parameters is not accurate, manual fine tuning may help.
- All particle images used in experiments are square images which means their heights equal their widths. Define the height of image as 's'.
  1. For synthetic GroEL datasets, EMAF-MD with low pass filter width s/10 is suggested.
  2. For synthetic EMPIAR 10065 datasets, EMAF-MD with low pass filter width s/20 is suggested.
  3. For synthetic EMPIAR 10398 datasets, EMAF-MD with low pass filter width s/20 is suggested.
  4. For synthetic Apoferritin datasets, EMAF-MD with low pass filter width s/15 is suggested.
  5. For real GroEL datasets, EMAF-MD with low pass filter width s/10 is suggested.
  6. For real EMPIAR 10398 datasets, EMAF-MD with low pass filter width s/20 is suggested.
