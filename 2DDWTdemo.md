# 2D DWT Demo using dwt\_2d #

_Note- I am using OPENCV to handle image I/O and other basic functionalities so most of the code will be different if you are using different libraries or software. Even if you are using OPENCV, you may not want to use the suboptimal algorithm I have used to tackle the UINT8 image coefficient overflow issues. This isn't something I have paid much attention to and this code is just to demonstrate DWT/IDWT functionality of the library._

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/imagedemo1.cpp) to compute DWT/IDWT of a 512X512 grayscale image is as following

```
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
```

Input Image in this example is 512X512 "Lena" grayscale image.

![https://lh6.googleusercontent.com/_TwtGvT0Ma-0/TbqrRs2Ze1I/AAAAAAAAADY/7TQUtL2FFHk/orig.jpg](https://lh6.googleusercontent.com/_TwtGvT0Ma-0/TbqrRs2Ze1I/AAAAAAAAADY/7TQUtL2FFHk/orig.jpg)

Two Level DWT of this image is computed with Daubechies' db4 wavelet.

![https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqrZCUoi2I/AAAAAAAAADc/l_jWn24PQcs/dwt.jpg](https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqrZCUoi2I/AAAAAAAAADc/l_jWn24PQcs/dwt.jpg)

Reconstructed Image is computed using idwt\_2d function.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqrftMQhCI/AAAAAAAAADg/58m4OlfoxvM/recon.jpg](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqrftMQhCI/AAAAAAAAADg/58m4OlfoxvM/recon.jpg)

_Additional Note About dwt\_2d- This implementation is a NXN input and NXN output system designed to work with dyadic lengths N. If the image dimensions are not dyadic, dwt\_2d will zeropad to make it NXN and you will end up with lots of zeros depending on the original size. The algorithm below (dwt\_2d\_sym) is designed to handle any image size but there are redundancies in the DWT stage as it is not a NXN I/O system._

# 2D DWT Demo using dwt\_2d\_sym #

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/imagedemo_sym.cpp) to compute DWT/IDWT  using symmetric extension dwt\_2d\_sym is as follows

```
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
        IplImage* img = cvLoadImage("Fig10.04(a).jpg");
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

        string nm = "db3";
        vector<double> l1,h1,l2,h2;
        filtcoef(nm,l1,h1,l2,h2);
        // unsigned int lf=l1.size();
        //  int rows_n =(int) (rows+ J*(lf-1));
        //  int cols_n =(int)  (cols + J * ( lf -1));

        // Finding 2D DWT Transform of the image using symetric extension algorithm
        // Extension is set to 3 (eg., int e = 3)

        vector<int> length;
        vector<double> output,flag;
        int J =3;
        int e=3;
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

        // Finding IDWT

        vector<vector<double> > idwt_output(rows, vector<double>(cols));

        idwt_2d_sym( output,flag, nm, idwt_output,length);



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

        cvNamedWindow( "Reconstructed Image", 1 ); // creation of a visualisation window
        cvShowImage( "Reconstructed Image", dvImg ); // image visualisation
        cvWaitKey();
        cvDestroyWindow("Reconstructed Image");
        cvSaveImage("recon.bmp",dvImg);







        return 0;
   }
```

Input Image is a 486X486 grayscale image taken from [Image Processing Place Database](http://www.imageprocessingplace.com/root_files_V3/image_databases.htm).

![https://lh4.googleusercontent.com/_TwtGvT0Ma-0/Tbqr0N7BjFI/AAAAAAAAADo/Sq-0hyQlhvE/orig.jpg](https://lh4.googleusercontent.com/_TwtGvT0Ma-0/Tbqr0N7BjFI/AAAAAAAAADo/Sq-0hyQlhvE/orig.jpg)

3-Level DWT is computed using db4 wavelet.

![https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbqr57dIOrI/AAAAAAAAADs/TiNqMpscsD8/dwt.jpg](https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbqr57dIOrI/AAAAAAAAADs/TiNqMpscsD8/dwt.jpg)

Perfect Reconstruction is achieved using idwt\_2d\_sym function.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/Tbqr9ku8TJI/AAAAAAAAADw/yqE3StuDNmE/recon.jpg](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/Tbqr9ku8TJI/AAAAAAAAADw/yqE3StuDNmE/recon.jpg)