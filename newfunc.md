# List of Functions for wavelet2d Library #

## 1D DWT/IDWT Functions ##

Both periodic and symmetric extension methods for Decimated DWT take exactly same input arguments.The difference is that Output vector in periodic extension has usually the same size (~N) as the input vector(N) while output vector in symmetric extension case has redundancies with length depending on size of filter used.

### Periodic Extension ###

**1. DWT:**  `void* dwt(vector<double> &sig, int J, string nm, vector<double> &dwt_output,vector<double> &flag, vector<int> &length )`

where,

_**sig**_ :: _Input Signal vector_

_**J**_ :: _Decomposition levels_

_**nm**_ :: _Wavelet Name_ _(See filtcoef for available wavelet families)_

_**length**_ :: _Lengths of respective approximation and detail vectors are stored in this integer vector._

_**dwt\_output**_ :: _Output of Discrete Wavelet Transform. It stores coefficients in following format:_

**`[A(J) D(J) D(J-1) ..... D(1)]`**

_where **A(J)** is the approximation coefficient vector at the Jth level while **D(n)** are the detail coefficient vectors at the nth level._**length**_contains the lengths of corresponding vectors. Last entry of the_**length**_vector is the length of the original signal._

_**flag**_ :: _Housekeeping vector. In this implementation it contains two values-_

_flag[0](0.md) is 0 if the signal is even and it is 1 if signal is odd and if it is made even by repeating the last value one more time_

_flag[1](1.md) - contains the decomposition levels._

_Housekeeping vector is a double vector as it was originally meant to store more values than it currently does_

![https://lh3.googleusercontent.com/-axX5BI71K7k/Tk7CvelHwrI/AAAAAAAAALI/Tc9h_KznCgM/s912/dwtperscreen.png](https://lh3.googleusercontent.com/-axX5BI71K7k/Tk7CvelHwrI/AAAAAAAAALI/Tc9h_KznCgM/s912/dwtperscreen.png)
_DWT stats (periodic extension) for an input signal of length 256_



**2. IDWT:** `void* idwt(vector<double> &dwtop,vector<double> &flag, string nm,vector<double> &idwt_output,vector<int> &length)`

where,

_**dwtop**_ :: _is the DWT vector_

_**flag**_ :: _Same Housekeeping function as obtained from the DWT function_

_**nm**_ :: _Wavelet Used_

_**idwt\_output**_ :: _Output of IDWT_

_**length**_ :: _Length vector obtained from the DWT computations_

### Symmetric Extension ###

**3. DWT:**  `void* dwt_sym(vector<double> &sig, int J, string nm, vector<double> &dwt_output,vector<double> &flag, vector<int> &length )`

where,

_**sig**_ :: _Input Signal vector_

_**J**_ :: _Decomposition levels_

_**nm**_ :: _Wavelet Name_ _(See filtcoef for available wavelet families)_

_**length**_ :: _Lengths of respective approximation and detail vectors are stored in this integer vector._

_**dwt\_output**_ :: _Output of Discrete Wavelet Transform. It stores coefficients in following format:_

**`[A(J) D(J) D(J-1) ..... D(1)]`**

_where **A(J)** is the approximation coefficient vector at the Jth level while **D(n)** are the detail coefficient vectors at the nth level._**length**_contains the lengths of corresponding vectors. Last entry of the_**length**_vector is the length of the original signal._

_**flag**_ :: _Housekeeping vector. In this implementation it contains two values-_

_flag[0](0.md) is 0 if the signal is even and it is 1 if signal is odd and if it is made even by repeating the last value one more time_

_flag[1](1.md) - contains the decomposition levels._

_Housekeeping vector is a double vector as it was originally meant to store more values than it currently does_

_**idwt\_output**_ :: _Output of the Inverse Discrete Wavelet Transform_

_**length**_ :: _Length Vector Obtained from the DWT computations_

![https://lh5.googleusercontent.com/-rooeVHG4pdc/Tk7FP_sPOTI/AAAAAAAAALQ/xWTasj0blNE/s912/dwtsymscreen.png](https://lh5.googleusercontent.com/-rooeVHG4pdc/Tk7FP_sPOTI/AAAAAAAAALQ/xWTasj0blNE/s912/dwtsymscreen.png)
_DWT stats (symmetric extension) for an input signal of length 256_

**4. IDWT:** `void* idwt_sym(vector<double> &dwtop,vector<double> &flag, string nm,vector<double> &idwt_output,vector<int> &length)`

where,

_**dwtop**_ :: _is the DWT vector_

_**flag**_ :: _Same Housekeeping function as obtained from the DWT function_

_**nm**_ :: _Wavelet Used_

_**idwt\_output**_ :: _Output of IDWT_

_**length**_ :: _Length vector obtained from the DWT computations_

## 1D SWT/ISWT Functions ##

**5. SWT:** `void* swt(vector<double> &sig, int J, string nm, vector<double> &swt_output, int &length)`

> All the coefficients are of equal lengths and that value is stored in _**length**_. _**swt\_output**_ stores value in the same format as _**dwt**_ and _**dwt\_sym**_ functions - Approximation coefficient vector at level J is stored at the beginning of the **swt\_output** vector followed by detail coefficients vectors at levels J, J-1,...., 1. The signal length has to be divisible by 2^J for reliable results. You can use signal extension (see below) to make the lengths compatible with the SWT.

_Two Level SWT Decomposition of a 247 length signal vector_

**6. ISWT:** `void* iswt(vector<double> &swtop,int J, string nm, vector<double> &iswt_output)`


_**swtop**_ is the output of SWT stage , _**J**_ - number of levels and _**nm**_ is the wavelet as before. Output of ISWT is stored in _**iswt\_output**_ vector.

## 2D DWT/IDWT Functions ##

As in 1D case, both periodic and symmetric extension methods for Decimated DWT take exactly same input arguments.The difference is that Output vector in periodic extension has usually the same size (~NXN) as the input vector(NXN) while output vector in symmetric extension case has redundancies with length/breadth depending on size of filter used.

### Periodic Extension ###

**7. 2D DWT:** `void* dwt_2d(vector<vector<double> > &origsig, int J, string nm, vector<double> &dwt_output,vector<double> &flag, vector<int> &length) `

_**origsig**_ :: _Input Image/Matrix_

_**J**_ :: _Number of Decomposition Levels_

_**nm**_ :: _Wavelet Name_

_**flag**_ :: _Stores values for IDWT function.Only flag[0](0.md) value is important as it contains decomposition levels._

_**dwt\_output**_ :: _1D vector that stores the output in the following format
**A(J) D<sub>h</sub>(J) D<sub>v</sub>(J) D<sub>d</sub>(J)  ..... D<sub>h</sub>(1) D<sub>v</sub>(1) D<sub>d</sub>(1)**_

where _**A(J)**_ is the approximation coefficient vector at the Jth level while _**D(n)**_ are the three detail coefficient vectors(horizontal,vertical and detail) at the nth level. It is important to remember that approximation and detail coefficients are actually two dimensional so we need a length vector that stores rows and columns values of each coefficient element. The length vector is given by _**length**_.

For example, the first element of _**output**_ vector is the approximation matrix stored as a vector and the first two elements of _**length**_ vectors are row and column values of the approximation matrix. In other words, a 300 element approximation matrix ( 15 rows X 20 columns) can be extracted from the 300 element approximation vector.

![https://lh3.googleusercontent.com/-mQm6JsFSDSE/Tk7eqaNPhLI/AAAAAAAAALY/CF70J_UuJOk/s912/dwt2dperscreen.png](https://lh3.googleusercontent.com/-mQm6JsFSDSE/Tk7eqaNPhLI/AAAAAAAAALY/CF70J_UuJOk/s912/dwt2dperscreen.png)
_2D DWT computation using periodic extension_

**8. 2D IDWT:** `void* idwt_2d(vector<double>  &dwtop,vector<double> &flag, string nm,vector<vector<double> > &idwt_output, vector<int> &length)`

where,

_**dwtop**_ :: _is the DWT vector_

_**flag**_ :: _Same Housekeeping function as obtained from the DWT function_

_**nm**_ :: _Wavelet Used_

_**idwt\_output**_ :: _Output of IDWT which should be defined to have the same number of rows and columns as the input image/matrix_

_**length**_ :: _Length vector obtained from the DWT computations_

### Symmetric Extension ###

**9. 2D DWT:** `void* dwt_2d_sym(vector<vector<double> > &origsig, int J, string nm, vector<double> &dwt_output,vector<double> &flag, vector<int> &length) `

_**origsig**_ :: _Input Image/Matrix_

_**J**_ :: _Number of Decomposition Levels_

_**nm**_ :: _Wavelet Name_

_**flag**_ :: _Stores values for IDWT function.Only flag[0](0.md) value is important as it contains decomposition levels._

_**dwt\_output**_ :: _1D vector that stores the output in the following format
**A(J) D<sub>h</sub>(J) D<sub>v</sub>(J) D<sub>d</sub>(J)  ..... D<sub>h</sub>(1) D<sub>v</sub>(1) D<sub>d</sub>(1)**_

where _**A(J)**_ is the approximation coefficient vector at the Jth level while _**D(n)**_ are the three detail coefficient vectors(horizontal,vertical and detail) at the nth level. It is important to remember that approximation and detail coefficients are actually two dimensional so we need a length vector that stores rows and columns values of each coefficient element. The length vector is given by _**length**_.

For example, the first element of _**output**_ vector is the approximation matrix stored as a vector and the first two elements of _**length**_ vectors are row and column values of the approximation matrix. In other words, a 300 element approximation matrix ( 15 rows X 20 columns) can be extracted from the 300 element approximation vector.

![https://lh3.googleusercontent.com/-H_5vbaS10FY/Tk7i3pyq5jI/AAAAAAAAALc/8ajg-QiCjpk/s912/dwt2dsymscreen.png](https://lh3.googleusercontent.com/-H_5vbaS10FY/Tk7i3pyq5jI/AAAAAAAAALc/8ajg-QiCjpk/s912/dwt2dsymscreen.png)
_2D DWT computation using symmetric extension_

**10. 2D IDWT:** `void* idwt_2d_sym(vector<double>  &dwtop,vector<double> &flag, string nm,vector<vector<double> > &idwt_output, vector<int> &length)`

where,

_**dwtop**_ :: _is the DWT vector_

_**flag**_ :: _Same Housekeeping function as obtained from the DWT function_

_**nm**_ :: _Wavelet Used_

_**idwt\_output**_ :: _Output of IDWT which should be defined to have the same number of rows and columns as the input image/matrix_

_**length**_ :: _Length vector obtained from the DWT computations_

## 2D SWT Function ##

**11. 2D SWT:** `void* swt_2d(vector<vector<double> > &sig,int J, string nm, vector<double> &swt_output)`

_**swt\_output**_ is a 1D vector which is arranged the same way as DWT output vector in the Decimated 2D cases above except that in this case all coefficients are of same size. This is a highly redundant transform as a three level decomposition of a 512X512 image results in 10 512X512 images - one approximation image and 9 detail images (three at each level).

## Convolution ##

**12. Convolution FFT\_ESTIMATE (Recommended):** `double convfft(vector<double> &a, vector<double> &b, vector<double> &c)`

Convolution function is pretty straightforward. _**a**_ and _**b**_ are input vectors and _**c**_ is the convolution output. _**convfft**_ uses FFT so it gives better results in most cases than the regular convolution which is implemented by _**convol**_ function.


**13. Convolution Direct** (Use it for only smaller vectors)**:** `double convol(vector<double> &a, vector<double> &b, vector<double> &c)`


**14. Convolution FFT\_MEASURE** (Recommended if you are going to perform convolutions of same length hundreds or thousands of time in one program)**:** `double convfftm(vector<double> &a, vector<double> &b, vector<double> &c)`

_**convfftm**_ is performed using MEASUREing capabilities of FFTW3 library so it is not recommended if you are going to convolve two vectors only once. This has some overhead but gives good results if multiple instances of same convolution are performed repeatedly.

## Wavelet Filters ##

**15. Filters:** `int filtcoef(string nm, vector<double> &lpd, vector<double> &hpd, vector<double> &lpr, vector<double> &hpr)`


_**nm:**_ Wavelet name.

_**lpd:**_ Low Pass Decomposition Filter Coefficients.

_**hpd:**_ High Pass Decomposition Filter Coefficients.

_**lpr:**_ Low Pass Reconstruction Filter Coefficients.

_**hpr:**_ High Pass Reconstruction Filter Coefficients.

All filters are `vector<double>` objects and can be obtained by specifying the wavelet name. Currently, following Wavelets are available:

_**Daubechies :**_ db1,db2,.., ,db15

_**Biorthogonal:**_ bior1.1 ,bior1.3 ,bior1.5 ,bior2.2 ,bior2.4 ,bior2.6 ,bior2.8 ,bior3.1 ,bior3.3 ,bior3.5 ,bior3.7 ,bior3.9 ,bior4.4 ,bior5.5 ,bior6.8

_**Coiflets:**_ coif1,coif2,coif3,coif4,coif5

_**Symmlets:**_ sym2,........, sym10

## 1D Vector Manipulation ##

**16. Downsampling:** `void downsamp(vector<double> &sig, int M, vector<double> &sig_d)`

_**sig**_ :: Signal to be Downsampled

_**M**_ :: Downsampling factor

_**sig\_d**_ :: Downsampled signal

**17. Upsampling:** `void upsamp(vector<double> &sig, int M, vector<double> &sig_u)`

_**sig**_ :: Signal to be Upsampled

_**M**_ :: Upsampling factor

_**sig\_u**_ :: Upsampled signal

**18. Periodic Extension:** `void* per_ext(vector<double> &sig, int a) `

_**per\_ext**_ periodically extends the signal _**sig**_ by value _**a**_ in either direction.

**19. Symmetric Extension:** `void* symm_ext(vector<double> &sig, int a) `

_**symm\_ext**_ symmetrically extensd the signal _**sig**_ by value _**a**_ in either direction. This function needs refinement as it doesn't gives good results for smaller vectors.

## 2D Vector Manipulation ##

**20. 2D Downsampling:** `void* downsamp2(vector<vector<double> > & vec1,vector<vector<double> > & vec2, int rows_dn, int cols_dn) `

_**vec1**_ is the input, _**rows\_dn**_ and _**cols\_dn**_ are row and column downsampling factors. _**vec2**_ is the downsampled matrix.

**21. 2D Upsampling:** `void* upsamp2(vector<vector<double> > & vec1,vector<vector<double> > & vec2, int rows_up, int cols_up)`

_**vec1**_ is the input, _**rows\_up**_ and _**cols\_up**_ are row and column upsampling factors. _**vec2**_ is the upsampled matrix.

**22. 2D Periodic Extension:** `void* per_ext2d(vector<vector<double> > &signal,vector<vector<double> > &temp2, int a)`

_**signal**_ is extended by **_a_** in all directions and the result is returned in 2D vector _**temp2**_

**23. 2D Symmetric Extension:** `void* symm_ext2d(vector<vector<double> > &signal,vector<vector<double> > &temp2, int a)`

_**signal**_ is extended by **_a_** in all directions and the result is returned in 2D vector _**temp2**_