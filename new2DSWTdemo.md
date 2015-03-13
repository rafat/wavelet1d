# 2D Stationary Wavelet Transform #

_Notes_

1. dwt\_2d and dwt\_2d\_sym functions are interchangeable in the following code as they take the same arguments.

2. Please make sure that right header file is included. Following sample program links to wavelet2d.h which is the header for shared libraries. wavelet2s.h must be included for static libraries.

3. I'm using OPENCV to handle images as i already use [OPENCV](http://opencv.willowgarage.com/wiki/) for other image processing work. You may want to use some simpler image libraries as OPENCV is a full image processing suite and is very bulky or you can just use 2D matrices/build your own image classes. Regardless, DWT/IDWT operations are more important than the choice of libraries.

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/wavelet2d/demo/swt2Ddemo.cpp)

```

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

```

Three level Stationary Wavelet Transform is computed using db2 wavelet. Approximation coefficients are stored only for the final (J=3) stage while the three detail coefficients( Horizontal, Vertical and Diagonal) are stored for each value. All 10 sets of coefficients are 512X512.

Input Image

![https://lh5.googleusercontent.com/-Y78s7sxlRcM/TlAItad_y4I/AAAAAAAAANM/fKDn4Um6Awo/lena512.jpg](https://lh5.googleusercontent.com/-Y78s7sxlRcM/TlAItad_y4I/AAAAAAAAANM/fKDn4Um6Awo/lena512.jpg)

Approximation Coefficient at level J = 3

![https://lh5.googleusercontent.com/-4E7Gtqj4MCs/TlAO_sX-kXI/AAAAAAAAANc/aXIMxdHv1Mw/blur.jpg](https://lh5.googleusercontent.com/-4E7Gtqj4MCs/TlAO_sX-kXI/AAAAAAAAANc/aXIMxdHv1Mw/blur.jpg)

Detail Coefficients at Levels 3,2,1(L-R). Coefficients are arranged Horizontally,Vertically and Diagonally from top to bottom.

![https://lh5.googleusercontent.com/-vqsGaPY5Pf8/TlAPFwwtzcI/AAAAAAAAANg/lV8npKIy_qc/s512/detail.jpg](https://lh5.googleusercontent.com/-vqsGaPY5Pf8/TlAPFwwtzcI/AAAAAAAAANg/lV8npKIy_qc/s512/detail.jpg)

_Resized Image (Actual Dimensions 1536X1536)_