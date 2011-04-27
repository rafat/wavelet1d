//============================================================================
// Name        : imagedemo_sym.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description :Image Approximation using symmetric extension DWT
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <algorithm>
#include "wavelet.h"
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

using namespace std;
using namespace cv;

void findthresh(vector<double> &vector1, int N, double& t){
	sort(vector1.begin(), vector1.end(), greater<double>());
	t = vector1.at(N-1);
}

void* maxval(vector<vector<double> > &arr, double &max){
	max = 0;
	for (unsigned int i =0; i < arr.size(); i++) {
		for (unsigned int j =0; j < arr[0].size(); j++) {
			if (max <= arr[i][j]){
				max = arr[i][j];
			}
		}
	}
	return 0;
}

void* maxval1(vector<double> &arr, double &max){
	max = 0;
	for (unsigned int i =0; i < arr.size(); i++) {
		if (max <= arr[i]){
			max = arr[i];
		}

	}
	return 0;
}


int main() {
	IplImage* img = cvLoadImage("lena512.bmp");
	if (!img){
		cout << " Can't read Image. Try Different Format." << endl;
		exit(1);
	}
	int height, width;
	height = img->height;
	width = img->width;
	int nc = img->nChannels;
	//   uchar* ptr2 =(uchar*) img->imageData;
	int pix_depth = img->depth;
	CvSize size;
	size.width =width;
	size.height=height;
	cout << "depth" << pix_depth <<  "Channels" << nc << endl;


	cvNamedWindow("Original Image", CV_WINDOW_AUTOSIZE);
	cvShowImage("Original Image", img);
	cvWaitKey();
	cvDestroyWindow("Original Image");
	cvSaveImage("orig.bmp",img);


	int rows =(int) height;
	int cols =(int) width;
	Mat matimg(img);

	vector<vector<double> > vec1(rows, vector<double>(cols));


	int k =1;
	for (int i=0; i < rows; i++) {
		for (int j =0; j < cols; j++){
			unsigned char temp;
			temp = ((uchar*) matimg.data + i * matimg.step)[j  * matimg.elemSize() + k ];
			vec1[i][j] = (double) temp;
		}

	}

	string nm = "db2";
	vector<double> l1,h1,l2,h2;
	filtcoef(nm,l1,h1,l2,h2);
	// unsigned int lf=l1.size();
	//  int rows_n =(int) (rows+ J*(lf-1));
	//  int cols_n =(int)  (cols + J * ( lf -1));

	// Finding 2D DWT Transform of the image using symetric extension algorithm
	// Extension is set to 0 (eg., int e = 0)

	vector<int> length;
	vector<double> output,flag;
	int J =6;
	int e=0;
	dwt_2d_sym(vec1,J,nm,output,flag,length,e);

	double max;
	vector<int> length2;
	// This algorithm computes DWT of image of any given size. Together with convolution and
	// subsampling operations it is clear that subsampled images are of different length than
	// dyadic length images. In order to compute the "effective" size of DWT we do additional
	// calculations.
	dwt_output_dim_sym(length,length2,J);
	// length2 is gives the integer vector that contains the size of subimages that will
	// combine to form the displayed output image. The last two entries of length2 gives the
	// size of DWT ( rows_n by cols_n)

	int siz = length2.size();
	int rows_n=length2[siz-2];
	int cols_n = length2[siz-1];

	vector<vector< double> > dwtdisp(rows_n, vector<double>(cols_n));
	dispDWT(output,dwtdisp, length ,length2, J);

	// dispDWT returns the 2D object dwtdisp which will be displayed using OPENCV's image
	// handling functions

	vector<vector<double> >  dwt_output= dwtdisp;

	// Stroing the DWT coefficients in two different vectors that will be used to approximate
	// Image with two different sets of chosen coefficients.

	vector<double> dwt_coef1;
	vector<double> dwt_coef2;

	dwt_coef1 = output;
	dwt_coef2 = output;

	maxval(dwt_output,max);// max value is needed to take care of overflow which happens because
	// of convolution operations performed on unsigned 8 bit images

	//Displaying Scaled Image
	// Creating Image in OPENCV
	IplImage *cvImg; // image used for output
	CvSize imgSize; // size of output image

	imgSize.width = cols_n;
	imgSize.height = rows_n;

	cvImg = cvCreateImage( imgSize, 8, 1 );
	// dwt_hold is created to hold the dwt output as further operations need to be
	// carried out on dwt_output in order to display scaled images.
	vector<vector<double> > dwt_hold(rows_n, vector<double>( cols_n));
	dwt_hold = dwt_output;
	// Setting coefficients of created image to the scaled DWT output values
	for (int i = 0; i < imgSize.height; i++ ) {
		for (int j = 0; j < imgSize.width; j++ ){
			if ( dwt_output[i][j] <= 0.0){
				dwt_output[i][j] = 0.0;
			}
			if ( i <= (length2[0]) && j <= (length2[1]) ) {
				((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
						(char) ( (dwt_output[i][j] / max) * 255.0);
			} else {
				((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
						(char) (dwt_output[i][j]) ;
			}
		}
	}

	cvNamedWindow( "DWT Image", 1 ); // creation of a visualisation window
	cvShowImage( "DWT Image", cvImg ); // image visualisation
	cvWaitKey();
	cvDestroyWindow("DWT Image");
	cvSaveImage("dwt.bmp",cvImg);

	 // Case 1 : Only 10% of the largest coefficients are considered

	          // Output is the 1D DWT vector

	          int n_coef1= int (output.size()/ 10);
	          cout << n_coef1 << endl;

	          // Finding Threshold Value corresponding to n_coef1

	          vector<double> temp1;
	             cout << "size: " << (int) temp1.size() << "\n";
	             cout << "capacity: " << (int) temp1.capacity() << "\n";
	             cout << "max_size: " << (int) temp1.max_size() << "\n";
	             for (unsigned int i =0; i < dwt_coef1.size(); i++) {
	           		  double tempval = abs(dwt_coef1[i]);
	           		  temp1.push_back(tempval);

	             }

	             double thresh1= 0.0;
	             findthresh(temp1,n_coef1,thresh1);
	             cout << "thresh" << thresh1 << endl;

	             ofstream temp("temp.txt");
	                       for (unsigned int i =0; i < temp1.size(); i++){
	                     	  temp << temp1[i] << " " ;
	                       }

	             // Reset coeffficients value depending on threshold value


	             for (unsigned int i =0; i < dwt_coef1.size(); i++) {
	       		    	double temp = abs(dwt_coef1[i]);

	       		    	if (temp < thresh1){
	       		    		dwt_coef1.at(i)= 0.0;

	       		    	}

	           	  }


	// Finding IDWT

	vector<vector<double> > idwt_output(rows, vector<double>(cols));

	idwt_2d_sym( dwt_coef1,flag, nm, idwt_output,length);



	//Displaying Reconstructed Image

	IplImage *dvImg;
	CvSize dvSize; // size of output image

	dvSize.width = idwt_output[0].size();
	dvSize.height = idwt_output.size();

	cout << idwt_output.size() << idwt_output[0].size() << endl;
	dvImg = cvCreateImage( dvSize, 8, 1 );

	for (int i = 0; i < dvSize.height; i++ )
		for (int j = 0; j < dvSize.width; j++ )
			((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
					(char) (idwt_output[i][j])  ;

	cvNamedWindow( "10% Coeff Reconstructed Image", 1 ); // creation of a visualisation window
	cvShowImage( "10% Coeff Reconstructed Image", dvImg ); // image visualisation
	cvWaitKey();
	cvDestroyWindow("10% Coeff Reconstructed Image");
	cvSaveImage("recon.bmp",dvImg);


	 // Case 2 : Only 2% of the largest coefficients are considered

		          // Output is the 1D DWT vector

		          int n_coef2= int (output.size()/ 50);
		          cout << n_coef2 << endl;

		          // Finding Threshold Value corresponding to n_coef1

		          vector<double> temp2;

		             for (unsigned int i =0; i < dwt_coef2.size(); i++) {
		           		  double tempval = abs(dwt_coef2[i]);
		           		  temp2.push_back(tempval);

		             }

		             double thresh2= 0.0;
		             findthresh(temp2,n_coef2,thresh2);
		             cout << "thresh" << thresh2 << endl;


		             // Reset coeffficients value depending on threshold value


		             for (unsigned int i =0; i < dwt_coef2.size(); i++) {
		       		    	double temp = abs(dwt_coef2[i]);

		       		    	if (temp < thresh2){
		       		    		dwt_coef2.at(i)= 0.0;

		       		    	}

		           	  }


		// Finding IDWT

		vector<vector<double> > idwt_output2(rows, vector<double>(cols));

		idwt_2d_sym( dwt_coef2,flag, nm, idwt_output2,length);



		//Displaying Reconstructed Image

		IplImage *dvImg2;
		CvSize dvSize2; // size of output image

		dvSize2.width = idwt_output2[0].size();
		dvSize2.height = idwt_output2.size();

		cout << idwt_output2.size() << idwt_output2[0].size() << endl;
		dvImg2 = cvCreateImage( dvSize2, 8, 1 );

		for (int i = 0; i < dvSize2.height; i++ )
			for (int j = 0; j < dvSize2.width; j++ )
				((uchar*)(dvImg2->imageData + dvImg2->widthStep*i))[j] =
						(char) (idwt_output2[i][j])  ;

		cvNamedWindow( "2% Coeff Reconstructed Image", 1 ); // creation of a visualisation window
		cvShowImage( "2% Coeff Reconstructed Image", dvImg2 ); // image visualisation
		cvWaitKey();
		cvDestroyWindow("2% Coeff Reconstructed Image");
		cvSaveImage("recon2.bmp",dvImg2);


	return 0;
}
