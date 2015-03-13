# C++ 1D/2D DWT Implementation for Win32 and Linux #

### Wavelet2d Libraries ###
  1. 1D DWT and IDWT Implementation (Two Modes)
  1. 2D DWT and IDWT Implementation (Two Modes)
  1. 1D SWT and ISWT Implementation ( Stationary Wavelet Transform)
  1. 2D SWT Implementation
  1. Implemented using FFTW3 Library
  1. Shared(.so) and static(.a) libraries for Linux
  1. Shared(.dll) and static(.a) libraries for Win32 GCC (MinGW).
  1. Shared(.dll) libraries for Microsoft VC++

> This implementation uses C++ vector objects to accept input, store and output data. This is done to overcome certain programming limitations associated with using arrays. Dynamic arrays were also an option but they require a lot of house-cleaning without having any advantages over vector objects. To learn more about vectors, please refer to excellent tutorial at [CPLUSPLUS](http://www.cplusplus.com/reference/stl/vector/). Vectors and Arrays can be converted back and forth so even if you are more comfortable with arrays, using vectors should not really be a problem.

### Changes from Previous Version ###

Focus is on speed.

1. FFTW3 Library is used to improve computation speed.

2. No GNUPLOT outputs are generated.

3. Periodic and Symmetric functions take same arguments. They will still be accessed by the same names, though.(dwt ,dwt\_sym etc.)

4. Focus on shared libraries. No static libraries for MSVC++ although MinGW and LINUX version will come with static as well as shared libraries.


wavelet2d - dynamic libraries( .so and .dll)

wavelet2s - static libraries (.a)

|**[List of Functions](http://code.google.com/p/wavelet1d/wiki/newfunc)**| Lists all Functions available in the Library|
|:-----------------------------------------------------------------------|:--------------------------------------------|
|**[Example Code 1](http://code.google.com/p/wavelet1d/wiki/new1DDWTdemo)**| 1D DWT/IDWT Demo Code|
|**[Example Code 2](http://code.google.com/p/wavelet1d/wiki/new1DAppx)**| 1D Signal Approximation Demo|
|**[Example Code 3](http://code.google.com/p/wavelet1d/wiki/new2DDWTdemo)**| 2D DWT/IDWT Demo Code|
|**[Example Code 4](http://code.google.com/p/wavelet1d/wiki/new2DAppx)**| Image Approximation Demo|
|**[Example Code 5](http://code.google.com/p/wavelet1d/wiki/new2DSWTdemo)**| 2D SWT Demo|
|**[Example Code 6](http://code.google.com/p/wavelet1d/wiki/new1DSWTDemo)**| 1D SWT Demo Code|

![https://lh4.googleusercontent.com/-4w-fWQyoieo/Thn_s7qb8-I/AAAAAAAAAGs/B0HFvLru6CU/s512/empire.jpg](https://lh4.googleusercontent.com/-4w-fWQyoieo/Thn_s7qb8-I/AAAAAAAAAGs/B0HFvLru6CU/s512/empire.jpg)

_J=2 Level Discrete Wavelet Transform of a 569X800 image using dwt\_2d\_sym function (Resized here for Display)_

Image Processing Note : I have not implemented any image class in C++ so I'm using OPENCV to handle images. It happens to be a rather bulky package and if you are not already working in image processing area, there are other more convenient options that may be used to handle simple image operations( loading, displaying, saving and converting them to c++ STL formats).
