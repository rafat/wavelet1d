I have added a DWT GUI for windows user :: [dyadwaves ](http://code.google.com/p/dyadwaves/). Any feedback will be appreciated.

# 1D/2D Discrete Wavelet Transform Implementation in C++ #

  1. 1D DWT and IDWT Implementation (Two Modes)
  1. 2D DWT and IDWT Implementation (Two Modes)
  1. 1D SWT and ISWT Implementation ( Stationary Wavelet Transform)
  1. 2D SWT Implementation
  1. Implemented using Danielson Lanczos FFT Algorithm

Libraries (Compiled using MinGW and MSVC++ 2010 compilers in Windows) are available in Downloads section while source code is available in the source "src" folder.

The libraries and source codes are available at http://code.google.com/p/wavelet1d/. I have tried to make this a more "user-friendly" implementation by keeping the functions as simple and descriptive as possible. Additionally, the 1D DWT/IDWT routines will output data that can be used easily with GNUPLOT. It causes redundancy and possibly some clutter if signal files are large but ,in my opinion, the trade-off is worth it. The accompanying script also works best only when DWT and IDWT routines are called in the same program.

I have used C++ vector objects to accept input, store and output data. This is done to overcome certain programming limitations associated with using arrays. Dynamic arrays were also an option but they require a lot of house-cleaning without having any advantages over vector objects. To learn more about vectors, please refer to excellent tutorial at http://www.cplusplus.com/reference/stl/vector/ Vectors and Arrays can be converted back and forth so even if you are more comfortable with arrays, using vectors should not really be a problem.

|**[List of Functions](http://code.google.com/p/wavelet1d/wiki/Functions)**| Lists all Functions available in the Library|
|:-------------------------------------------------------------------------|:--------------------------------------------|
|**[Example Code 1](http://code.google.com/p/wavelet1d/wiki/1DDWTdemo)**| 1D DWT/IDWT Demo Code|
|**[Example Code 2](http://code.google.com/p/wavelet1d/wiki/1DAppx)**| 1D Signal Approximation Demo|
|**[Example Code 3](http://code.google.com/p/wavelet1d/wiki/2DDWTdemo)**| 2D DWT/IDWT Demo Code|
|**[Example Code 4](http://code.google.com/p/wavelet1d/wiki/2DAppx)**| Image Approximation Demo|
|**[Example Code 5](http://code.google.com/p/wavelet1d/wiki/2DSWTdemo)**| 2D SWT Demo Code|
|**[Example Code 6](http://code.google.com/p/wavelet1d/wiki/1DSWTDemo)**| 1D SWT Demo Code|

![https://lh4.googleusercontent.com/-4w-fWQyoieo/Thn_s7qb8-I/AAAAAAAAAGs/B0HFvLru6CU/s512/empire.jpg](https://lh4.googleusercontent.com/-4w-fWQyoieo/Thn_s7qb8-I/AAAAAAAAAGs/B0HFvLru6CU/s512/empire.jpg)

_J=2 Level Discrete Wavelet Transform of a 569X800 image using dwt\_2d\_sym function (Resized here for Display)_

Image Processing Note : I have not implemented any image class in C++ as of yet so I'm using OPENCV to handle images. It happens to be a rather bulky package and if you are not already working in image processing area, there are other more convenient options that may be used to handle simple image operations( loading, displaying, saving and converting them to c++ STL formats).

