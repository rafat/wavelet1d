//============================================================================
// Name        : imagedemo1.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : 2D DWT and IDWT using OPENCV and a 512 by 512 grayscale image.
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "wavelet.h"

using namespace std;
using namespace cv;

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

void* minval(vector<vector<double> > &arr, double &min){
	min = 10000;
	for (unsigned int i =0; i < arr.size(); i++) {
		for (unsigned int j =0; j < arr[0].size(); j++) {
			if (min >= arr[i][j]){
				min = arr[i][j];
			}
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

     int rr1,cc1;
 	 string nm = "db4";
 	 // Finding DWT output dimensions as the DWT output is zeropadded
     dwt_output_dim(vec1, rr1, cc1 );
 	int J = 2;
 	vector<double> flag;
    vector<vector<double> >  dwt_output(rr1, vector<double>(cc1));
    cout << rr1 << cc1 << "size of op" << endl;
    // Computing 2d DWT ( vec1 is the signal while dwt_output is the DWT output)
 	dwt_2d(vec1, J, nm, dwt_output,flag );

   cout << "dwt size" << dwt_output.size() << dwt_output[0].size() << endl;

   double max,min;

   maxval(dwt_output,max);
   minval(dwt_output,min);
   cout << "maxval" << max << "  minval" << min << endl;

   //Displaying Scaled Image
      // Creating Image in OPENCV
   	  IplImage *cvImg; // image used for output
      CvSize imgSize; // size of output image

      imgSize.width = cc1;
      imgSize.height = rr1;

      cvImg = cvCreateImage( imgSize, 8, 1 );
      // dwt_hold is created to hold the dwt output as further operations need to be
      // carried out on dwt_output in order to display scaled images.
      vector<vector<double> > dwt_hold(rr1, vector<double>(cc1));
      dwt_hold = dwt_output;

     // Setting coefficients of created image to the scaled DWT output values
      for (int i = 0; i < imgSize.height; i++ ) {
      for (int j = 0; j < imgSize.width; j++ ){
		  if ( dwt_output[i][j] <= 0.0){
			  dwt_output[i][j] = 0.0;
		  }
		  if ( i <= (imgSize.height/pow(2.0,double(J))) && j <= (imgSize.width/pow(2.0,double(J))) ) {
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
	// Finding IDWT
	   vector<vector<double>  > final(rr1, vector<double>(cc1));
	   idwt_2d(dwt_hold,flag, nm ,final);

	   // Removing Zeropadding

	   zero_remove(vec1,final);

	   //Displaying Reconstructed Image

	   IplImage *dvImg;
	   CvSize dvSize; // size of output image

	   dvSize.width = final[0].size();
	   dvSize.height = final.size();

	   dvImg = cvCreateImage( dvSize, 8, 1 );

	      for (int i = 0; i < dvSize.height; i++ )
	      for (int j = 0; j < dvSize.width; j++ )
	          ((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
	          (char) (final[i][j]) ;

		  cvNamedWindow( "Reconstructed Image", 1 ); // creation of a visualisation window
	       cvShowImage( "Reconstructed Image", dvImg ); // image visualisation
		   cvWaitKey();
		   cvDestroyWindow("Reconstructed Image");
		   cvSaveImage("recon.bmp",dvImg);
       ofstream diff("diff.dat");
       for (unsigned int i=0; i < final.size(); i++) {
    	for (unsigned int j = 0; j < final[0].size(); j++) {
    		diff << final[i][j]-vec1[i][j] << " " ;
    	}
    	diff << endl;
        }
       ofstream recon("recon.dat");
       for (unsigned int i=0; i < final.size(); i++) {
    	for (unsigned int j = 0; j < final[0].size(); j++) {
    		recon << (int)final[i][j] << " " ;
    	}
    	recon << endl;
        }
       ofstream orig("orig.dat");
              for (unsigned int i=0; i < vec1.size(); i++) {
           	for (unsigned int j = 0; j < vec1[0].size(); j++) {
           		orig << (int)vec1[i][j] << " " ;
           	}
           	orig << endl;
               }
	return 0;
}
