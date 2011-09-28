//============================================================================
// Name        : swt2Ddemo.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : 2D SWT Demo using OPENCV
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <algorithm>
#include "wavelet2d.h"
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"

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
//   vector<double> l1,h1,l2,h2;
//   filtcoef(nm,l1,h1,l2,h2);


   vector<double> output;
   int J =3;
   swt_2d(vec1,J,nm,output);
   cout << "OUTPUT size" << output.size() << endl;
   cout << "LOOP OK" << endl;

   int row,col;
   row=vec1.size();
   col=vec1[0].size();

   // Extract and Display Low Pass Image at the Jth stage

   vector<vector<double> > blur(row,vector<double>(col));

   for (int i=0;i < row; i++){
       for (int j=0; j < col;j++){
        double temp = output[i*col + j];
        blur[i][j]= temp;
       }
   }

  double max;
  maxval(blur,max);

        // Creating Image in OPENCV
           IplImage *cvImg; // image used for output
        CvSize imgSize; // size of output image

        imgSize.width = col;
        imgSize.height = row;

        cvImg = cvCreateImage( imgSize, 8, 1 );

        for (int i = 0; i < imgSize.height; i++ ) {
        for (int j = 0; j < imgSize.width; j++ ){
            if ( blur[i][j] <= 0.0){
                blur[i][j] = 0.0;
            }

        ((uchar*)(cvImg->imageData + cvImg->widthStep*i))[j] =
         (char) ( (blur[i][j] / max) * 255.0);

        }
        }

        cvNamedWindow( "Low Pass Image", 1 ); // creation of a visualisation window
         cvShowImage( "Low Pass Image", cvImg ); // image visualisation
         cvWaitKey();
         cvDestroyWindow("Low Pass Image");
           cvSaveImage("blur.bmp",cvImg);

         // Displaying BandPass Images

    vector<vector<double> >  detail(3*row,vector<double>(J * col));

    for (int k=0; k < J; k++) {
    for (int i=0; i < 3*row; i++) {
        for(int j=0+ k*col; j < (k+1)*col; j++) {
            double temp = output[(3*k+1)*row*col+ i * col +j - k*col];
            detail[i][j]= temp;
        }
    }
    }
      IplImage *dvImg; // image used for output
      CvSize imgSz; // size of output image

      imgSz.width = J*col;
      imgSz.height = 3*row;

      dvImg = cvCreateImage( imgSz, 8, 1 );

      for (int i = 0; i < imgSz.height; i++ ) {
      for (int j = 0; j < imgSz.width; j++ ){
          if ( detail[i][j] <= 0.0){
              detail[i][j] = 0.0;
          }

      ((uchar*)(dvImg->imageData + dvImg->widthStep*i))[j] =
       (char) detail[i][j];

      }
      }

      cvNamedWindow( "Band Pass Image", 1 ); // creation of a visualisation window
       cvShowImage( "Band Pass Image", dvImg ); // image visualisation
       cvWaitKey();
       cvDestroyWindow("Band Pass Image");
       cvSaveImage("detail.bmp",dvImg);


        return 0;
}