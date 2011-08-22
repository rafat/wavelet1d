#ifndef WAVELET2D_H
#define WAVELET2D_H
#include <vector>
#include <complex>
using namespace std;


     // the dll exports
#if defined WAVE_EXPORT
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __declspec(dllimport)
#endif


// 1D Functions


EXPORT void* dwt1(string, vector<double> &, vector<double> &, vector<double> &);

EXPORT void* dyadic_zpad_1d(vector<double> &);

EXPORT double convol(vector<double> &, vector<double> &, vector<double> &);

EXPORT int filtcoef(string , vector<double> &, vector<double> &, vector<double> &,
                vector<double> &);

EXPORT void downsamp(vector<double> &, int , vector<double> &);

EXPORT void upsamp(vector<double> &, int, vector<double> &);

EXPORT void circshift(vector<double> &, int );

EXPORT int sign(int);

EXPORT void* idwt1(string wname, vector<double> &, vector<double> &, vector<double> &);

EXPORT int vecsum(vector<double> &, vector<double> &, vector<double> &);



// 1D Symmetric Extension DWT Functions



EXPORT void* dwt_sym(vector<double> &, int ,string , vector<double> &,vector<double> &,
                vector<int> &);

EXPORT void* dwt1_sym(string , vector<double> &, vector<double> &, vector<double> &);

EXPORT void* idwt_sym(vector<double> &,vector<double> &, string,vector<double> &, vector<int> &);

EXPORT void* symm_ext(vector<double> &, int );

EXPORT void* idwt1_sym(string, vector<double> &, vector<double> &, vector<double> &); // Not Tested

// 1D Stationary Wavelet Transform

EXPORT void* swt(vector<double> &, int , string , vector<double> &, int &) ;

EXPORT void* iswt(vector<double> &,int , string, vector<double> &);

EXPORT void* per_ext(vector<double> &, int );




// 2D Functions

EXPORT void* branch_lp_dn(string , vector<double> &, vector<double> &);

EXPORT void* branch_hp_dn(string , vector<double> &, vector<double> &);

EXPORT void* branch_lp_hp_up(string ,vector<double> &, vector<double> &, vector<double> &);

//EXPORT void* dwt_2d(vector<vector<double> > &, int , string , vector<vector<double> > &
  //              , vector<double> &) ;

//EXPORT void* idwt_2d(vector<vector<double> > &,vector<double> &, string ,vector<vector<double> > &);

EXPORT void* dyadic_zpad_2d(vector<vector<double> > &,vector<vector<double> > &);

EXPORT void* dwt_output_dim(vector<vector<double> >&, int &, int & );

EXPORT void* zero_remove(vector<vector<double> > &,vector<vector<double> > &) ;

EXPORT void* getcoeff2d(vector<vector<double> > &, vector<vector<double> > &,
                vector<vector<double> > &,vector<vector<double> > &,vector<double> &, int &);

EXPORT void* idwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> >  &, vector<vector<double> > &);

EXPORT void* dwt2(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> > &, vector<vector<double> > &);

EXPORT void* downsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int);

EXPORT void* upsamp2(vector<vector<double> > &,vector<vector<double> > &, int, int);

// 2D DWT (Symmetric Extension) Functions

EXPORT void* dwt_2d_sym(vector<vector<double> > &, int , string , vector<double> &, vector<double> & ,
                vector<int> &);

EXPORT void* dwt2_sym(string ,vector<vector<double> > &, vector<vector<double> >  &,
                vector<vector<double> >  &, vector<vector<double> >  &, vector<vector<double> > &);

EXPORT void* idwt_2d_sym(vector<double>  &,vector<double> &, string ,vector<vector<double> > &,
                vector<int> &);

EXPORT void* circshift2d(vector<vector<double> > &, int , int );

EXPORT void symm_ext2d(vector<vector<double> > &,vector<vector<double> > &, int );

EXPORT void* dispDWT(vector<double> &,vector<vector<double> > &, vector<int> &, vector<int> &, int ) ;

EXPORT void* dwt_output_dim_sym(vector<int> &,vector<int> &, int );

//2D Stationary Wavelet Transform

EXPORT void* swt_2d(vector<vector<double> > &,int , string , vector<double> &);

EXPORT void* per_ext2d(vector<vector<double> > &,vector<vector<double> > &, int );

// FFT functions


EXPORT double convfft(vector<double> &, vector<double> &, vector<double> &);

EXPORT double convfftm(vector<double> &, vector<double> &, vector<double> &);

EXPORT void* fft(vector<complex<double> > &,int ,unsigned int);

EXPORT void* bitreverse(vector<complex<double> > &);

EXPORT void* freq(vector<double> &, vector<double> &);

//New


EXPORT void* dwt1_sym_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD);//FFTW3 for 2D

EXPORT void* idwt1_sym_m(string wname, vector<double> &X, vector<double> &app, vector<double> &detail);

EXPORT void* dwt(vector<double> &sig, int J, string nm, vector<double> &dwt_output
                , vector<double> &flag, vector<int> &length );

EXPORT void* idwt(vector<double> &,vector<double> &, string,vector<double> &, vector<int> &);

EXPORT void* dwt_2d(vector<vector<double> > &, int , string , vector<double> &, vector<double> & ,
                vector<int> &);
EXPORT void* dwt1_m(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) ;

EXPORT void* idwt_2d(vector<double>  &dwtop,vector<double> &flag, string nm,
                vector<vector<double> > &idwt_output, vector<int> &length);

EXPORT void* idwt1_m(string wname, vector<double> &X, vector<double> &cA, vector<double> &cD);

EXPORT void* dwt_output_dim2(vector<int> &length, vector<int> &length2, int J);


#endif/* WAVELET2D_H */
