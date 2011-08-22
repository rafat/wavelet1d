Notes

1.dwt and dwt_sym functions are interchangeable in all the codes as they take the same arguments.

2.Please make sure that right header file is included. Header file for dynamic libraries is wavelet2d.h and for static libraries, wavelet2s.h is used.

3.I’m using OPENCV to handle images as i already use OPENCV for other image processing work. You may want to use some simpler image libraries as OPENCV is a full image processing suite and is very bulky or you can just use 2D matrices/build your own image classes. Regardless, DWT/IDWT operations are more important than the choice of libraries.

4.All updated example codes are available at
http://code.google.com/p/wavelet1d/source/browse/trunk/wavelet2d/demo/