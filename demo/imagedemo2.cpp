//============================================================================
// Name        : imagedemo2.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : Image Approximation
//============================================================================

// IMPORTANT - Algorithm used to display Image is imprecise because of int 8 overflow issues
// and it shouldn't be used to judge the performance of the DWT. The DWT and IDWT outputs
// should be used for performance measurements. I have used maximum value rescaling to
// solve overflow issues and , obviously, it is going to result in suboptimal performance but
// it is good enough for demonstration purposes.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "wavelet.h"

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
 	 string nm = "db3";
 	 // Finding DWT output dimensions as the DWT output is zeropadded
     dwt_output_dim(vec1, rr1, cc1 );
 	int J = 5;
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



   vector<vector<double> > dwt_coef1(rr1, vector<double>(cc1));
   vector<vector<double> > dwt_coef2(rr1, vector<double>(cc1));
   // Storing dwt_output for two different computations based on
   // different approximation values
   dwt_coef1 = dwt_output;
   dwt_coef2 = dwt_output;



   //Displaying Scaled Image
      // Creating Image in OPENCV
   	  IplImage *cvImg; // image used for output
      CvSize imgSize; // size of output image

      imgSize.width = cc1;
      imgSize.height = rr1;

      cvImg = cvCreateImage( imgSize, 8, 1 );

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
         // IMPORTANT -- dwt_output value has been modified above.

   	  cvNamedWindow( "DWT Image", 1 ); // creation of a visualisation window
          cvShowImage( "DWT Image", cvImg ); // image visualisation
   	   cvWaitKey();
   	   cvDestroyWindow("DWT Image");


      // Case 1 : Only 10% of the largest coefficients are considered

          // Total elements in the original image is rows * cols
          int n_coef1= int (( rows * cols)/ 10);

          // Finding Threshold Value corresponding to n_coef1

          vector<double> temp1;
          cout << "size: " << (int) temp1.size() << "\n";
          cout << "capacity: " << (int) temp1.capacity() << "\n";
          cout << "max_size: " << (int) temp1.max_size() << "\n";
          for (int i =0; i < rr1; i++) {
        	  for (int j = 0; j < cc1; j++){
        		  double tempval = abs(dwt_coef1[i][j]);
        		  temp1.push_back(tempval);
        	  }
          }

          double thresh1= 0.0;
          findthresh(temp1,n_coef1,thresh1);
          ofstream temp("temp.txt");
          for (int i =0; i < rows * cols; i++){
        	  temp << temp1[i] << " " ;
          }
          // Reset coeffficients value depending on threshold value


          for (int i =0; i < rr1; i++) {
        	  for (int j = 0; j < cc1; j++){
    		    	double temp = abs(dwt_coef1[i][j]);

    		    	if (temp < thresh1){
    		    		dwt_coef1[i][j] = 0.0;

    		    	}
    		    }
        	  }


    	// Finding IDWT
    	   vector<vector<double>  > final(rr1, vector<double>(cc1));
    	   idwt_2d(dwt_coef1,flag, nm ,final);

    	   // Removing Zeropadding

    	   zero_remove(vec1,final);

    	    double max1;
    	    maxval(final,max1);

    	   //Displaying Reconstructed Image

    	   IplImage *dvImg;
    	   CvSize dvSize; // size of output image

    	   dvSize.width = final[0].size();
    	   dvSize.height = final.size();

    	   dvImg = cvCreateImage( dvSize, 8, 1 );

    	      for (int i = 0; i < dvSize.height; i++ ) {
    	      for (int j = 0; j < dvSize.width; j++ ){
    	   		  if ( final[i][j] <= 0.0){
    	   			  final[i][j] = 0.0;
    	   		  }
    	          ((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
    	          (char) ((final[i][j]/max1) * 255) ;
    	      }
    	      }

    		  cvNamedWindow( "10% Coeff Reconstructed Image", 1 ); // creation of a visualisation window
    	       cvShowImage( "10% Coeff Reconstructed Image", dvImg ); // image visualisation
    		   cvWaitKey();
    		   cvDestroyWindow("10% Coeff Reconstructed Image");

    		   // Case 2 : Only 2% of the largest coefficients are considered

    		       // Total elements in the original image is rows * cols
    		       int n_coef2= int (( rows * cols) / 50);

    		       // Finding Threshold Value corresponding to n_coef1

    		       vector<double> temp2;
    		       for (int i =0; i < rr1; i++) {
    		     	  for (int j = 0; j < cc1; j++){
    		     		  double tempval = abs(dwt_coef2[i][j]);
    		     		  temp2.push_back(tempval);
    		     	  }
    		       }
    		       double thresh2= 0.0;
    		       findthresh(temp2,n_coef2,thresh2);

    		       // Reset coeffficients value depending on threshold value

    		       for (int i =0; i < rr1; i++) {
    		     	  for (int j = 0; j < cc1; j++){
    		 		    	double temp = abs(dwt_coef2[i][j]);

    		 		    	if (temp < thresh2){
    		 		    		dwt_coef2[i][j] = 0.0;

    		 		    	}
    		 		    }
    		     	  }





    		 	// Finding IDWT
    		 	   vector<vector<double>  > final2(rr1, vector<double>(cc1));
    		 	   idwt_2d(dwt_coef2,flag, nm ,final2);

    		 	   // Removing Zeropadding

    		 	   zero_remove(vec1,final2);
    		 	   double max2;
    		 	   maxval(final2,max2);

    		 	   //Displaying Reconstructed Image

    		 	   IplImage *dvImg2;
    		 	   CvSize dvSize2; // size of output image

    		 	   dvSize2.width = final2[0].size();
    		 	   dvSize2.height = final2.size();

    		 	   dvImg2 = cvCreateImage( dvSize2, 8, 1 );

    		 	      for (int i = 0; i < dvSize2.height; i++ ) {
    		 	      for (int j = 0; j < dvSize2.width; j++ ){
    			   		  if ( final2[i][j] <= 0.0){
    			   			  final2[i][j] = 0.0;
    			   		  }
    		 	          ((uchar*)(dvImg2->imageData + dvImg2->widthStep*i))[j] =
    		 	          (char) ((final2[i][j]/ max2)* 255) ;
    		 	      }
    		 	      }

    		 		  cvNamedWindow( "2% Coeff Reconstructed Image", 1 ); // creation of a visualisation window
    		 	       cvShowImage( "2% Coeff Reconstructed Image", dvImg2 ); // image visualisation
    		 		   cvWaitKey();
    		 		   cvDestroyWindow("2% Coeff Reconstructed Image");

    		 		   cout << thresh1 << " " << thresh2 << endl;


	return 0;
}
