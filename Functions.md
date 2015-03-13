_Note : Only first ten functions are available in the wave1d library and convfft is also not included. It is my recommendation that you use the latest library(wavelet.lib and libwavelet.a in wavelet-03.rar) even for 1D operations as 1D implementation is faster and more accurate._

# 1D Functions #

### 1. J-Level Discrete Wavelet Transform ###
**`void* dwt(vector<double> &, int ,string , vector<double> &, vector<double> &)`**

Example Usage : **dwt(sig, J, nm, dwt\_output,flag )**

where,

**sig :** is the input signal . It is a `vector<double>` object.

**J :** Number of DWT Levels. Integer.

**nm :** Name of Wavelet. See function filtcoef for more information. String.

**dwt\_ouput:** is a `vector<double>` object that will return the DWT output.

**flag:** is a vector that contains two `vector<double>` values.

`flag[0]` - contains the number of zeros if the signal is zeropadded.

`flag[1]`- contains the number of decomposition levels. (J)

flag is used by IDWT to help process the DWT output.

Additional Outputs : In addition to above, the function also generates .txt and .dat files.

GNUPLOT ready .dat files

**gnusig.dat :** contains signal data

**gnufilt.dat:** Contains filters data.[Four filters are arranged columnwise in the file](All.md)

**gnuout.dat :** Contains DWT coefficients. It stores approximation and detail coefficients columnwise in alternating fashion. For example, three level of approximation and detail coefficients will result in gnuout.dat having 7 columns. First column is index followed by alternating approximation and detail columns.

Text Files : These two files are used by the code to generate gnuout.dat file.

**appx.txt :** Values of Approximation Coefficients from each stage. Values are appended after each stage.

**det.txt :** Values of Detail Coefficients after each stage. Values are likewise appended after each stage.

**dwtout.txt:** Stores the DWT ouput in text format. It contains Approximation coefficients at Jth level, followed by detail coefficients at Jthe, J-1,.... level 1.

**flag.txt:** Contains values of flags in text format.

### 2. 1-Level Discrete Wavelet Transform ###

**`void* dwt1(string, vector<double> &, vector<double> &, vector<double> &)`**

Example Usage: **dwt1(nm,sig, appx\_sig, det\_sig)**

**nm :** Wavelet Name. String

**sig :** Input Signal. `Vector<double>` Object.

**appx\_sig:** Approximation Coefficients. `Vector<double>` object.

**det\_sig:** Detail Coefficients. `Vector<double>` object.

Same Functionality can be achieved by setting J=1 in the dwt function but dwt1 directly outputs detail and approximation coefficients for one level of decomposition. This function will not generate any of GNUPLOT files but may be preferable in any number of scenarios.

### 3. Dyadic Zero-Padding ###

**`void* dyadic_zpad_1d(vector<double> &)`**

Example Usage: **dyadic\_zpad\_1d(sig)**

It zeropads signal sig such that its length is now dyadic. For example a signal with length 195 will be zeropadded to 256.

### 4. Convolution ###

Recommended -- **`double convfft(vector<double> &, vector<double> &, vector<double> &)`**

**-or-**

**`double convol(vector<double> &, vector<double> &, vector<double> &)`**

Example Usage: **convfft(signal,lpd,cA\_undec)**

where, signal and lpd are two inputs that are convolved and the convolution output is stored in cA\_undec. All three are `vector<double>` objects.

Important- **convfft** will not modify the value of signal while convol will. I have not corrected this as **convfft** is the faster, better implementation and should be used instead of **convol**.

### 5. Filter Coefficients ###

**`int filtcoef(string , vector<double> &, vector<double> &, vector<double> &, vector<double> &)`**

Example Usage: **filtcoef(nm,lpd,hpd,lpr,hpr)**

**nm:** Wavelet name.

**lpd:** Low Pass Decomposition Filter Coefficients.

**hpd:** High Pass Decomposition Filter Coefficients.

**lpr:** Low Pass Reconstruction Filter Coefficients.

**hpr:** High Pass Reconstruction Filter Coefficients.

All filters are `vector<double>` objects and can be obtained by specifying the wavelet name. Currently, following Wavelets are available:

Daubechies : db1,db2,.., ,db15

Biorthogonal: bior1.1 ,bior1.3 ,bior1.5 ,bior2.2 ,bior2.4 ,bior2.6 ,bior2.8 ,bior3.1 ,bior3.3 ,bior3.5 ,bior3.7 ,bior3.9 ,bior4.4 ,bior5.5 ,bior6.8

Coiflets: coif1,coif2,coif3,coif4,coif5

### 6. Downsampling ###

**`void downsamp(vector<double> &, int , vector<double> &)`**

Example Usage : **downsamp(cA\_undec, D, cA)**

**cA\_undec:** Signal to be downsampled. `Vector<double>`

**D :** Downsampling Factor. Integer

**cA:** Downsampled Signal. `Vector<double>`

### 7. Upsampling ###

**`void upsamp(vector<double> &, int, vector<double> &)`**

Example Usage : **upsamp(cA, U, cA\_up)**

**cA:** Signal to be upsampled. `Vector<double>`

**U :** Upsampling Factor. Integer

**cA\_up:** Upsampled Signal. `Vector<double>`

### 8. J-Level Inverse Discrete Wavelet Transform ###

**`void* idwt(vector<double> &,vector<double> &, string , vector<double> &)`**

Example Usage: **idwt(dwt\_output, flag,nm,output)**

**dwt\_output:** Output of DWT which serves as Input to the IDWT. `vector<double>`

**flag:** Contains values of J and zeropadding. `vector<double>`

**nm :** Wavelet name. String

**output:** Reconstructed IDWT Output. `vector<double>`

Additional .dat output

**gnurecon.dat :** GNUPLOT ready reconstructed signal.

### 9. 1-Level Inverse Discrete Wavelet Transform ###

**`void* idwt1(string wname, vector<double> &, vector<double> &, vector<double> &)`**

Example Usage : **idwt1(nm,idwt\_output, app,detail)**

**nm:** Wavelet Name.

**idwt\_output:** Outputs IDWT.

**app:** Approximation Coefficients.

**detail:** Detail Coefficients.

### 10. GNUDWTPLOT ###

**`void* gnudwtplot(int)`**

Example Usage: **gnudwtplot(J)**

Generates GNUPLOT script gnudwt.gnu.txt that contains plots of signals, filters and wavelet coefficients. it works if it is called following DWT and IDWT operations. To draw the plots, just load the script in gnuplot. See wavedemo1.cpp for more details.



# 1D Symmetric Extension DWT Functions #

### 11. J-Level Symmetric Extension DWT ###

**`void* dwt_sym(vector<double> &, int ,string , vector<double> &,vector<double> &,vector<int> &, int );`**

Example Usage: **dwt\_sym (signal, J,nm, dwt\_output, flag,length, e)**

Same as **dwt** function except it contains two more arguments and the output is a 1D vector. **length** is an integer vector which contains lengths of approximation and detail coefficients. **e** is the length of symmetric extension in each direction. **dwt\_output** stores coefficients in following format:
**`[A(J) D(J) D(J-1) ..... D(1)]`**

> where **A(J)** is the approximation coefficient vector at the Jth level while **D(n)** are the detail coefficient vectors at the nth level. **length** contains the lengths of corresponding vectors. Last entry of the **length** vector is the length of the original signal.

### 12. J-Level Symmetric Extension IDWT ###

**`void* idwt_sym(vector<double> &,vector<double> &, string,vector<double> &, vector<int> &);`**

Example Usage: **idwt\_sym(dwtop,flag,nm,idwt\_output,length)**

**idwt\_sym** is identical to **idwt** above except that it contains one more argument in lengths of the coefficients vector **length**. It is the same vector described in **11**.

### 13. Symmetric Extension ###

**`void* symm_ext(vector<double> &, int )`**

Example Usage: **symm\_ext(signal,e)**

**symm\_ext** symmetrically extends the signal by length e in either direction.

# 1D Stationary Wavelet Transform #

### 14. J-Level Stationary Wavelet Transform ###

**`void* swt(vector<double> &, int , string , vector<double> &, int &) `**

Example Usage: **swt(sig,J,nm,swt\_output,length)**

This function works on dyadic length signal. All the coefficients are of equal lengths and that value is stored in **length**. **swt\_output** stores value in the same format as **11** - Approximation coefficient vector at level J is stored at the beginning of the **swt\_output** vector followed by detail coefficients vectors at levels J, J-1,...., 1. Additionally, ALL coefficient vectors(including redundant approximation vectors) are stored in **gnuout.dat** file. See **gnuswtplot** for more details.

### 15. J-Level Inverse Stationary Wavelet Transform ###

**`void* iswt(vector<double> &,int , string, vector<double> &)`**

Example Usage: **iswt(swtop,J,nm,iswt\_output)**

**swtop** is the output of SWT stage , J - number of levels and nm is the wavelet as before. Output of ISWT is stored in **iswt\_output** vector.

### 16. GNUSWTPLOT ###

**`void* gnuswtplot(int)`**

Example Usage: **gnuswtplot(J)**

Generates GNUPLOT script gnudwt.gnu.txt that contains plots of signals, filters and wavelet coefficients. it works if it is called following SWT and ISWT operations. To draw the plots, just load the script in gnuplot.

### 17. 1D Periodic Extension ###

**`void* per_ext(vector<double> &, int )`**

Example Usage: **per\_ext(sig,a)**

**per\_ext** periodically extends the signal **sig** by value **a** in either direction.

# 2D Functions #

There are three sets of 2D DWT options.

1. Regular decimated 2D DWT and IDWT. It works with dyadic length 2D signals and outputs 2D signals.(The DWT output is in the format of a dyadically decomposed image with cLL, the low pass output, in the top left corner) If the input signal isn't of dyadic length then it will be zeropadded and you can calculate output size by using **dwt\_output\_dim** function and remove zeros from the output by using the **zero\_remove** function. More on this will follow. For a N X N dydadic length signal, the output is N X N which is useful in many applications.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/Tbk_OaFazJI/AAAAAAAAABI/N9CQwURPJfs/ezw.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/Tbk_OaFazJI/AAAAAAAAABI/N9CQwURPJfs/ezw.png)

2. Symmetric Extension 2D DWT and IDWT signal. It works with signal of any length and DWT output is a 1D vector (just as in Matlab) which may make it easier to handle in certain situations. On the flip side, the output is not of the same size as input and so there are some redundancies. To display the image, you'll have to convert 1D vector into image format and also account for downsampling and convolution which result in non-linear variation of length across deocmposition levels. The functions **dispDWT** and **dwt\_output\_dim\_sym** will help in displaying the image. On the flip side, performing operations on 1D
DWT output vector should be straightforward.

3. 2D Stationary Wavelet Transform. Inverse Stationary Transform is not implemented as of yet but 2d SWT should be good for image analysis, pattern recognition and edge detection operations that don't require ISWT for the most part.

### 18. J-Level 2D Discrete Wavelet Transform ###

**`void* dwt_2d(vector<vector<double> > &, int , string , vector<vector<double> > &, vector<double> &)`**

Example Usage: **dwt\_2d(vec1, J, nm, dwt\_output,flag )**

**vec1** is the 2D input signal.

**J** is the Decomposition level.

**nm** gives the name of wavelet family. See **filtcoef** for more details.

**dwt\_output** is the 2D output coefficients. They are arranged as shown in the figure above. The dimensions can be obtained by using **dwt\_output\_dim** function. A code fragment may look like this-
```
  int rr1,cc1;
  string nm = "db4";
  // Finding DWT output dimensions as the DWT output is zeropadded
  dwt_output_dim(vec1, rr1, cc1 );// vec1 is the input image/matrix
  int J = 2;
  vector<double> flag;// Flag stores values including size of input and level of decomposition and this vector needs to be passed to idwt_2d function.
  vector<vector<double> >  dwt_output(rr1, vector<double>(cc1));
  // Computing 2d DWT ( vec1is the signal while dwt_output is the DWT   output)
  dwt_2d(vec1, J, nm, dwt_output,flag );
```

### 19. J-Level 2D Inverse Discrete Wavelet Transform ###

**`void* idwt_2d(vector<vector<double> > &,vector<double> &, string ,vector<vector<double> > &)`**

Example Usage: **idwt\_2d(dwt\_output,flag, nm ,idwt\_output)**

**idwt\_output** is a 2D vector object and may require zeropadding removal depending on the size of input image/matrix. Continuing from above, code fragment for IDWT may look like following

```
 // Finding IDWT
 vector<vector<double>  > final(rr1, vector<double>(cc1));
 idwt_2d(dwt_hold,flag, nm ,final);
 // Removing Zeropadding
 zero_remove(vec1,final);// vec1 is the input image/matrix
```

### 20. 2D Dyadic Zero Padding ###

**`void* dyadic_zpad_2d(vector<vector<double> > &,vector<vector<double> > &)`**

Example Usage: **dyadic\_zpad\_2d(sig,sig2)**

This function makes the signal of dyadic length by adding zeros. In this particular DWT implementation **dwt\_2d**, the output **sig2** will be a NXN matrix which may require a lot of zeros being padded to the signal. It may be desirable to use **dwt\_2d\_sym** function in this case which operates on images/matrices of any size without any zero-padding.

### 21. Get Output Dimension ###

**`void* dwt_output_dim(vector<vector<double> >&, int &, int & )`**

Example Usage: **dwt\_output\_dim(vec1, rr1, cc1 )**

Because of zero-padding, the output image may not be equal in size to the input image. This function returns the dimensions of output matrix. If you are using dyadic length image then output will have the same dimensions as the input. This function works only with **dwt\_2d** functions and not with other DWT implementations.

### 22. 2D Zeropad Removal ###

**`void* zero_remove(vector<vector<double> > &,vector<vector<double> > &) `**

Example Usage: **zero\_remove(vec1,final)**

**vec1** is the input signal while **final** is the zeropadded IDWT output. **zero\_remove** will remove zeros based on the size of the input image/matrix.

### 23. Get 2D DWT coefficients at level N ###

**`void* getcoeff2d(vector<vector<double> > &, vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<double> &, int &)`**

Example Usage: **getcoeff2d(dwtoutput,cH,cV,cD,flag,N)**

Given 2D DWT output **dwtoutput** and the Nth level of decomposition, **getcoeff2d** will return the three detail coefficients at the Nth level of decomposition. Works only with **dwt\_2d** function.

### 24. 1-Level 2D DWT ###

**`void* dwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,vector<vector<double> >  &, vector<vector<double> > &, vector<vector<double> > &)`**

Example Usage: **dwt2(nm,sig,cA,cH,cV,cD)**

Given a 2D signal **sig** and wavelet name **nm**, **dwt2** will return the approximation(**cA**) and detail coefficients(**cH,cV,cD**) after 1 level of decomposition.

### 25. 1-Level 2D IDWT ###

**`void* idwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,vector<vector<double> >  &, vector<vector<double> >  &, vector<vector<double> > &)`**

Example Usage: **idwt2(nm,idwt\_output, cA, cH,cV,cD)**

This is the exact inverse of **dwt2** function. Inputs are the four approximation and detail coefficients while the output is the singke stage **idwt\_output**.

### 26. 2D Downsampling ###

**`void* downsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int)`**

Example Usage: **downsamp2(vec1,vec2,rows\_dn,cols\_dn)**

**vec1** is the input, **rows\_dn** and **cols\_dn** are row and column downsampling factors. **vec2** is the downsampled matrix.

### 27. 2D Upsampling ###

**`void* upsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int)`**

Example Usage: **upsamp2(vec1,vec2,rows\_up,cols\_up)**

**vec1** is the input, **rows\_up** and **cols\_up** are row and column upsampling factors. **vec2** is the upsampled matrix.

# 2D Symmetric Extension DWT Functions #

### 28. J-Level Symmetric Extension 2D DWT ###

**`void* dwt_2d_sym(vector<vector<double> > &, int , string , vector<double> &, vector<double> & ,vector<int> &, int )`**

Example Usage: **dwt\_2d\_sym(vec1,J,nm,output,flag,length,e)**

**vec1** Input Image/Matrix

**J** Number of Decomposition Levels

**nm** Wavelet Name

**flag** Stores values for IDWT function

**output** 1D vector that stores the output in the following format
**`[A(J) D(J) D(J-1) ..... D(1)]`**

where **A(J)** is the approximation coefficient vector at the Jth level while **D(n)** are the detail coefficient vectors at the nth level. It is important to remember that approximation and detail coefficients are actually two dimensional so we need a length vector that stores rows and columns values of each coefficient element. The length vector is given by **length**.

For example, the first element of **output** vector is the approximation matrix stored as a vector and the first two elements of **length** vectors are row and column values of the approximation matrix. In other words, a 300 element approximation matrix ( 15 rows X 20 columns) can be extracted from the 300 element approximation vector.

**e** is the length of symmetric extension in all directions.

A code fragment is as following

```
 	string nm = "db4";

	// Finding 2D DWT Transform of the image using symmetric      extension algorithm
	// Extension is set to 3 (eg., int e = 3)

	vector<int> length;
	vector<double> output,flag;
	int J =3;
	int e=3;
	dwt_2d_sym(vec1,J,nm,output,flag,length,e);
```

This DWT implementation works on image/matrices of any size and no dyadic zero padding is required.

### 29. J-Level Symmetric Extension 2D IDWT ###

**`void* idwt_2d_sym(vector<double>  &,vector<double> &, string ,vector<vector<double> > &,vector<int> &)`**

Example Usage: **idwt\_2d\_sym(dwt\_output,flag, nm, idwt\_output,length)**

**idwt\_output** will have the same dimensions as the input image/matrix.

```
 	// Finding IDWT

	vector<vector<double> > idwt_output(rows, vector<double>(cols));// rows and cols are the rows and columns of the input image/matrix

	idwt_2d_sym( output,flag, nm, idwt_output,length);
```

> For more on other arguments, see 28.

### 30. 2D Symmetric Extension ###

**`void symm_ext2d(vector<vector<double> > &,vector<vector<double> > &, int )`**

Example Usage: **symm\_ext2d(origsig,sig, ext)**

**origsig** is the original signal, **sig** is the extended signal. It is extended by integer value **ext** in all directions.

### 31. Additional dwt\_2d\_sym Functions ###

**`void* dispDWT(vector<double> &,vector<vector<double> > &, vector<int> &, vector<int> &, int )`**

and

**`void* dwt_output_dim_sym(vector<int> &,vector<int> &, int )`**

> As we mentioned previously,DWT output of symmetric extension 2D DWT is a one dimensional vector. Additionally, the output isn't of dyadic length across different decomposition levels which makes the rectangular display of DWT output image a bit challenging and imprecise. **dispDWT** and **dwt\_output\_dim\_sym** are two functions that are used to return a rectangular set of coefficients arranged to display the DWT in traditional format. A code fragment that computes symmetric DWT and then rearranges coefficients is shown next

```
        string nm = "db4";

	// Finding 2D DWT Transform of the image using symmetric      extension algorithm
	// Extension is set to 3 (eg., int e = 3)

	vector<int> length;
	vector<double> output,flag;
	int J =3;
	int e=3;
	dwt_2d_sym(vec1,J,nm,output,flag,length,e);
        
        vector<int> length2;
	// This algorithm computes DWT of image of any given size. Together with convolution and
	// subsampling operations it is clear that subsampled images are of different length than
	// dyadic length images. In order to compute the "effective" size of DWT we do additional
	// calculations.
	dwt_output_dim_sym(length,length2,J);
	// length2 gives the integer vector that contains the size of subimages that will
	// combine to form the displayed output image. The last two entries of length2 gives the
	// size of DWT ( rows_n by cols_n)

	int siz = length2.size();
        // The last two values of length2 vector give the dimensions of the DWT output.
	int rows_n=length2[siz-2];
	int cols_n = length2[siz-1];

	vector<vector< double> > dwtdisp(rows_n, vector<double>(cols_n));
	dispDWT(output,dwtdisp, length ,length2, J);
```

> There are ,obviously, other ways to display the DWT output. A more precise and more efficient way could be to display output separately at each level instead of trying to bruteforce them into one image.

# 2D Stationary Wavelet Transform Functions #

### 32. J-Level 2D Stationary Wavelet Transform ###

**`void* swt_2d(vector<vector<double> > &,int , string , vector<double> &)`**

Example Usage: **swt\_2d(vec1,J,nm,output)**

**output** is a 1D vector which is exactly equal to 1D DWT output vector of **dwt\_2d\_sym** except that in this case all coefficients are of same size. If **vec1** isn't a (2<sup>M)</sup> X (2<sup>M)</sup> image/matrix , it will be zeropadded to dyadic dimensions with rows = columns.

For example, a 3 level decomposition of 256X256 image yields 10 256X256 sets of coefficients - one approximation coefficient set at level 3 and three detail coefficients at each level.

![https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbnv_x4dUbI/AAAAAAAAACk/9tvLG_-iZb4/apprx.jpg](https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbnv_x4dUbI/AAAAAAAAACk/9tvLG_-iZb4/apprx.jpg)

_Approximation Coefficients at J=3 (256X256 Image)_

![https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbnvFEiRhpI/AAAAAAAAACc/doxXP9gpcl0/s640/detail.jpg](https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbnvFEiRhpI/AAAAAAAAACc/doxXP9gpcl0/s640/detail.jpg)

_Detail Coefficients at levels 3,2,1 respectively (768X768 image shown as 640X640)_

### 33. 2D Periodic Extension ###

**`void per_ext2d(vector<vector<double> > &,vector<vector<double> > &, int )`**

Example Usage: **per\_ext2d(origsig,sig, ext)**

**origsig** is the original signal, **sig** is the extended signal. It is extended by integer value **ext** in all directions.

# FFT Functions #

### 34. 1D FFT and IFFT ###

**`void* fft(vector<complex<double> > &,int ,unsigned int)`**

Example Usage: **fft(inp,1,N) -or- fft(inp,-1,N)**

**fft** performs in-place FFT operation on complex vector **inp**.

Second argument is 1 for forward FFT and -1 for Inverse FFT.

Third argument gives the length of the output. If N is larger than the length of **inp** then it will be zeropadded, FFT will be performed and N-length complex output will be returned.

### 35. Frequency Response ###

**`void* freq(vector<double> &, vector<double> &)`**

Example Usage: **freq(sig,fre\_oup)**

This function returns frequency response **fre\_oup** of signal **sig**. The function will also output a freq.dat file that can be directly plotted in GNUPLOT.