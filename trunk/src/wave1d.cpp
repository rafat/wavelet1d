//============================================================================
// Name        : dwt_J_level.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   :
// Description : J-Level Wavelet Decomposition
//============================================================================

#include <iostream>
#include <fstream>
#include "wave1d.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std;

void* dwt(vector<double> &sig, int J, string nm, vector<double> &dwt_output
		, vector<double> &flag ) {

	int Max_Iter;
		    Max_Iter = (int) ceil(log( double(sig.size()))/log (2.0));
		    if ( Max_Iter < J) {
		    	cout << J << " Iterations are not possible with signal of length " << sig.size() << endl;
		    	exit(1);
		    }

    vector<double> original_copy, appx_sig, det_sig;
    original_copy = sig;
    cout <<  " original signal" << sig.size()<< endl;
    // Zero Pad the Signal to nearest 2^ M value ,where M is an integer.
    unsigned int n_size = sig.size();
    flag.push_back(0);
    double int_n_size = log10 (double (n_size)) / log10(2.0);
    if ( (pow(2.0, double(int(ceil(int_n_size)))) - pow(2.0, int_n_size)) != 0) {
    dyadic_zpad_1d(sig);
    flag.pop_back();
    int flag_1 = (int) (pow(2.0, double(int(ceil(int_n_size)))) - pow(2.0, int_n_size)) ;
    flag.push_back(flag_1);
    }

    flag.push_back(J);

    ofstream appx("appx.txt", ios::trunc);
    appx.close();
    ofstream det("det.txt", ios::trunc);
    det.close();
    	//<------------------------------------------------------------------->
        // Preparing for GnuPlot
    	// Storing Signal for Gnuplot
        unsigned int count = sig.size();
    	ofstream gnusig("gnusig.dat");
    	for (unsigned int i = 0;i < count; i++) {
    		gnusig << i << " " << sig[i] << endl;
    	}
    	gnusig.close();

    	//  Storing Filter Values for GnuPlot
    	     vector<double> lp1,hp1,lp2,hp2;
    	     filtcoef(nm,lp1,hp1,lp2,hp2);
        ofstream gnufilt("gnufilt.dat");
        unsigned int len_filt = lp1.size();
        for (unsigned int i = 0; i < len_filt; i++) {
        	gnufilt << i << " " << lp1[i] << " " << hp1[i] << " " << lp2[i] << " " << hp2[i] << endl;
        }
        gnufilt.close();
    	// <------------------------------------------------------------------->


    for (int iter = 0; iter < J; iter++) {
        dwt1(nm,sig, appx_sig, det_sig);
    	dwt_output.insert(dwt_output.begin(),det_sig.begin(),det_sig.end());



        if (iter == J-1 ) {
        	dwt_output.insert(dwt_output.begin(),appx_sig.begin(),appx_sig.end());

        }
  //      vector<double> lp1,hp1,lp2,hp2;
  //      filtcoef(nm,lp1,hp1,lp2,hp2);
  //      gnudwtplot(sig, appx_sig, det_sig,lp1,hp1,lp2,hp2);
        sig.clear();
        sig = appx_sig;
        appx_sig.clear();
        det_sig.clear();

    }

    // Outputting Dwt Output and Flag as text files

    ofstream dwtout ("dwtout.txt");
   	        for (unsigned int i = 0; i < dwt_output.size(); i++){
   	        	dwtout << dwt_output[i] << endl;
   	        }
   	    cout << "output size" << dwt_output.size() <<endl;
   	  //  cout << flag[0]<< " " << flag[1] << endl;

   	    ofstream flagged("flag.txt");
   	    for (int i =0; i < 2; i++){
   	    	flagged << flag[i] << endl;
   	    }


 // Processing Approximation and Detail Coefficients at each level for GnuPlot

    vector<double> appx_c , det_c;
    ifstream sig_inp("appx.txt");
    int appx_len = 0;

    int det_len = 0;
    while (sig_inp) {
    	double temp;
    	sig_inp >> temp;
    	appx_c.push_back(temp);
    	appx_len++;
    }
    appx_c.pop_back();
    appx_len--;

    ifstream sig_inp2("det.txt");
    while (sig_inp2) {
    	double temp;
    	sig_inp2 >> temp;
    	det_c.push_back(temp);
    	det_len++;
    }
    det_c.pop_back();
    det_len--;
  //  cout << det_len + 1 << endl;

    vector<vector<double>  > coeff(appx_len, vector<double>(2*J + 1,0));
     for (int i = 0; i < appx_len; i++) {
    	 coeff[i][0] = i;
   // 	 cout << coeff[i][0] << endl;
     }
     cout << "appx_len" << appx_len << endl;

     int val = 0;
     for (int j= 0; j < J  ; j++ ) {
    	 int temp  = count / (int) pow (2.0, (double) j+1);
    	 for (int i = 0; i <  temp; i++) {
    		 coeff[i][2 * (j + 1 ) - 1]=appx_c[val + i];
    		 coeff[i][2 * (j + 1)]=det_c[val + i];

    	 }
		 val+=temp;

     }

     ofstream gnuout("gnuout.dat");
     for (int i = 0; i < appx_len; i++){
    	 for (int j =0; j < 2* J +1; j++) {
    		 gnuout << coeff[i][j] << " ";
    	 }
    	 gnuout << endl;
     }

     sig = original_copy;
	return 0;
}

void circshift(vector<double> &sig_cir, int L){
	if ( abs(L) > sig_cir.size()) {
		L = sign(L) * (abs(L) % sig_cir.size());
	}

	if ( L < 0 ){
		L = (sig_cir.size() + L) % sig_cir.size();
	//	cout << "L" << L << endl;
	}
		for (int i = 0; i < L; i++){
			sig_cir.push_back(sig_cir[0]);
			sig_cir.erase(sig_cir.begin());
		}

}

double convol(vector<double> &a, vector<double> &b, vector<double> &c) {
     unsigned int len_c = a.size() + b.size() - 1;
     double*  oup= NULL;

     oup = new double[len_c];
     vector<double>::iterator a_it;
     a_it = a.end();
     signed int al = len_c - a.size();
     a.insert(a_it,al,0);


     vector<double>::iterator b_it;
     b_it = b.end();
     signed int bl = len_c - b.size();
     b.insert(b_it,bl, 0);


     for (unsigned int ini = 0; ini < len_c ; ini++){
    	 double ou1 = 0;
    	 oup[ini] = 0;
    	 double temp = 0;
    	 for (unsigned int jni = 0; jni <= ini; jni++) {
    		 ou1 = a[jni] * b[ini - jni];
    		 oup[ini]+= ou1;
    	 }
    	 temp = oup[ini];
    	 c.push_back(temp);
     }
     delete [] oup;
     oup = NULL;
     return 0;
}

void downsamp(vector<double> &sig, int M, vector<double> &sig_d){
	int len = sig.size();
	double len_n = ceil( (double) len / (double) M);
	for (int i = 0; i < (int) len_n; i++) {
		double temp = sig[i*M];
		sig_d.push_back(temp);
	}
}

// Function Declaration

// Global Variables

vector<double> lpd, hpd, lpr, hpr;



void* dwt1(string wname, vector<double> &signal, vector<double> &cA, vector<double> &cD) {

		filtcoef(wname,lpd,hpd,lpr,hpr);
//		for (unsigned int i = 0;  i < signal.size(); i++) {
//		cout << signal[i] <<  endl;
//		out2 << signal[i] <<endl;
//		}

		unsigned int temp_len = signal.size();

		if ( (temp_len % 2) != 0) {
			double temp =signal[temp_len - 1];
			signal.push_back(temp);
		}

		int len_sig = signal.size();
		int len_lpfilt = lpd.size();
		int len_hpfilt = hpd.size();
		int len_avg = (len_lpfilt + len_hpfilt) / 2;
		// cout << len_lpfilt << "Filter" << endl;
		circshift(signal,-len_avg / 2); // Signal is shifted circularly in order to perform
		// computations designed to deal with boundary distortions

		// Low Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));
		vector<double> cA_undec;
		// Convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
		convol(signal,lpd,cA_undec);
		int D = 2; // Downsampling Factor is 2

		// Downsampling by 2 gives cA
		downsamp(cA_undec, D, cA);
		int L = len_lpfilt / 2;
		int N = len_sig / 2;

		for (int i = 0; i < L; i++) {
			cA[i] = cA[N+i]+cA[i];
		}
		cA.resize(N);
	    circshift(cA,len_avg/2 );

		// High Pass Filtering Operations in the Analysis Filter Bank Section
//		int len_cA =(int)  floor(double (len_sig + len_lpfilt -1) / double (2));

		vector<double> cD_undec;
		// Convolving signal with lpd, Low Pass Filter, and O/P is stored in cA_undec
		convol(signal,hpd,cD_undec);
		 // Downsampling Factor is 2

		// Downsampling by 2 gives cA
		downsamp(cD_undec, D, cD);

		int Lh = len_hpfilt / 2;

		for (int i = 0; i < Lh; i++) {
			cD[i] = cD[N+i]+cD[i];
		}
		cD.resize(N);
	    circshift(cD,len_avg/2 );


		// File Outputs
		filtcoef(wname,lpd,hpd,lpr,hpr);

		fstream appx_s ("appx.txt", fstream::out | fstream::app);
		fstream det_s ("det.txt", fstream::out | fstream::app);
		for (int i = 0; i < N; i++){
			appx_s << cA[i] << endl;
			det_s << cD[i] << endl;
		}

//		cout << lpd.size() << "filter size" << endl;

  return 0;
}

void* dyadic_zpad_1d(vector<double> &signal) {
	unsigned int N = signal.size();
	double M = log10 (double (N)) / log10(2.0);
	int D = (int) ceil(M);
	double int_val = pow(2.0, double(D)) - pow(2.0, M);

	int zeros = (int) int_val;
    vector<double>::iterator a_it;
    a_it = signal.end();
 //   double val = signal[N-1];
    signal.insert(a_it,zeros,0);
    return 0;

}

void* gnudwtplot(int J) {

    ofstream gnudwt("gnudwt.gnu.txt");
    gnudwt << "# This Script will Plot Signal, DWT Coefficients and Filters Used" << endl;
    gnudwt << "set terminal wxt 1" << endl;
    gnudwt <<"set title \" Input Signal \" " << endl;
    gnudwt << "plot \"gnusig.dat\" using 1:2 with lines " << endl;
    gnudwt << "set terminal wxt 2" << endl;
    gnudwt << "set multiplot layout "<< J << "," << 2 << endl;
    for (int i = 0;i < J ; i++){
    gnudwt << "# Plot " << 2*(i + 1) -1 << endl;
    gnudwt << "set title \"Approximation Coefficients Level  " << i+ 1 << " \"" << endl;
    gnudwt << "plot \"gnuout.dat\" using 1:"<< 2*(i+1) << " with lines" << endl;
    gnudwt << "# Plot "<< 2*(i + 1)<< endl;
    gnudwt << "set title \"Detail Coefficients Level  " << i + 1 <<"\" " << endl;
    gnudwt << "plot \"gnuout.dat\" using 1:"<< 2*(i+1)+1<< " with lines" << endl;
    }
    gnudwt << "unset multiplot" <<endl;
    gnudwt << "set terminal wxt 3" << endl;
    gnudwt << "set multiplot layout 2,2" << endl;
    gnudwt << "set title \"Filter Coefficients\" " << endl;
    gnudwt << "# Plot 1" << endl;
    gnudwt << "set title \"Low Pass Decomposition Filter \" " << endl;
    gnudwt << "plot \"gnufilt.dat\" using 1:2 with impulse" << endl;
    gnudwt << "# Plot 2" << endl;
    gnudwt << "set title \"High Pass Decomposition Filter \" " << endl;
    gnudwt << "plot \"gnufilt.dat\" using 1:3 with impulse" << endl;
    gnudwt << "# Plot 3" << endl;
    gnudwt << "set title \"Low Pass Reconstruction Filter \" " << endl;
    gnudwt << "plot \"gnufilt.dat\" using 1:4 with impulse" << endl;
    gnudwt << "# Plot 4" << endl;
    gnudwt << "set title \"High Pass Reconstruction Filter \" " << endl;
    gnudwt << "plot \"gnufilt.dat\" using 1:5 with impulse" << endl;
    gnudwt << "unset multiplot" <<endl;
    gnudwt << "set terminal wxt 4" << endl;
    gnudwt <<"set title \" Reconstructed Signal \" " << endl;
    gnudwt << "plot \"gnurecon.dat\" using 1:2 with lines " << endl;




	return 0;

}

void* idwt(vector<double> &dwtop,vector<double> &flag, string nm,
		vector<double> &idwt_output) {

        cout << flag.size() << endl;
        int J =(int) flag[1];
 //       int zpad =(int) flag[0];

	    vector<double> app;
	    vector<double> detail;
	    unsigned int app_len = dwtop.size() / int(pow(2.0,J));
        vector<double>::iterator dwt;
        dwt = dwtop.begin();
        cout << app_len << endl;
        app.assign(dwt,dwtop.begin()+app_len);
        detail.assign(dwtop.begin()+app_len, dwtop.begin()+ 2* app_len);

	    for (int i = 0; i < J; i++) {

	    	idwt1(nm,idwt_output, app,detail);
	    	app_len = 2 * app_len;
	    	app.clear();
	    	detail.clear();
			if ( i < J - 1 ) {
	        detail.assign(dwtop.begin()+app_len, dwtop.begin()+ 2* app_len);
			}
            app = idwt_output;
	    }
 //       app.clear();
 //       detail.clear();

//	    vector<double> signal;
//	    idwt1(nm,signal, appx, det);

	    // Remove ZeroPadding

	    int zerop =(int) flag[0];
	    idwt_output.erase(idwt_output.end()- zerop,idwt_output.end());


        unsigned int count = idwt_output.size();
    	ofstream gnurecon("gnurecon.dat");
    	for (unsigned int i = 0;i < count; i++) {
    		gnurecon << i << " " << idwt_output[i] << endl;
    	}
    	gnurecon.close();

	    return 0;
}


void* idwt1(string wname, vector<double> &X, vector<double> &cA, vector<double> &cD) {
	vector<double> lpd1,hpd1, lpr1, hpr1;

	filtcoef(wname,lpd1,hpd1,lpr1,hpr1);
	int len_lpfilt = lpr1.size();
	int len_hpfilt = hpr1.size();
	int len_avg = (len_lpfilt + len_hpfilt) / 2;
	unsigned int N = 2 * cA.size();
	int U = 2; // Upsampling Factor

	// Operations in the Low Frequency branch of the Synthesis Filter Bank

	vector<double> cA_up;
	vector<double> X_lp;
    circshift(cA,-len_avg/2);

	upsamp(cA, U, cA_up);
	convol(cA_up, lpr1, X_lp);


	// Operations in the High Frequency branch of the Synthesis Filter Bank

	vector<double> cD_up;
	vector<double> X_hp;
    circshift(cD,-len_avg/2);
	upsamp(cD, U, cD_up);
	convol(cD_up, hpr1, X_hp);


	// Operations to obtain reconstructed signal by folding back and circular shifting

	cout << "size of N" << N << endl;
	cout << "size of filter" << len_lpfilt << endl;
	cout << "size of X_lp" << X_lp.size() << endl;
	for (unsigned int i = 0 ; i < X_lp.size() - N  ; i++){
		X_lp[i] = X_lp[N+i] + X_lp[i];
	}

	for (unsigned int i = 0 ; i < X_hp.size() - N; i++){
			X_hp[i] = X_hp[N+i] + X_hp[i];
		}

	X_lp.resize(N);
	X_hp.resize(N);

//	for (unsigned int i =0 ; i < N ; i++){
//    X[i] = X_lp[i] + X_hp[i];
//	}
	vecsum(X_lp,X_hp,X);
   // Work on circular shift
   circshift(X,len_avg +len_avg / 2 -1);

//    ofstream sig("sig.txt");
    cout << X.size() << endl;
	cout << "success " << endl;

	/*
    for (unsigned int i = 0; i < N; i++){
    	sig << X[i] << endl;
    }
    */

    return 0;
}

int sign(int X) {
	if (X >= 0)
		return 1;
	else
		return -1;
}

void upsamp(vector<double> &sig, int M, vector<double> &sig_u) {
	int len = sig.size();
	double len_n = ceil( (double) len * (double) M);

	for (int i = 0; i < (int) len_n; i++) {
		if ( i % M == 0) {
			        double temp = sig[i / M];
					sig_u.push_back(temp);

				}
				else
				{
					 sig_u.push_back(0);
				}

	}



}

double op_sum(double i, double j) {
	return (i+j);
}

int vecsum(vector<double> &a, vector<double> &b, vector<double> &c){

 cout << "Size a" << a.size() << endl;
 cout << "Size b" << b.size() << endl;

    c.resize(a.size());
	transform (a.begin(), a.end(), b.begin(), b.begin(), op_sum);
	c = b;
	cout << "Exit Ok" << endl;
		return 0;
}

int filtcoef(string name, vector<double> &lp1, vector<double> &hp1, vector<double> &lp2,
		vector<double> &hp2){
    if (name == "haar" || name == "db1" ) {
    	lp1.push_back(0.7071);lp1.push_back(0.7071);
    	hp1.push_back(-0.7071);hp1.push_back(0.7071);
    	lp2.push_back(0.7071);lp2.push_back(0.7071);
    	hp2.push_back(0.7071);hp2.push_back(-0.7071);
  //  	cout << lp2[1] << endl;
//    	hpd = {-0.7071, 0.7071};
//    	lpr = {0.7071, 0.7071};
//    	hpr = {0.7071, -0.7071};
     return 0;
    }
    else if ( name == "db2"){
    	double lp1_a[] = {-0.12940952255092145, 0.22414386804185735, 0.83651630373746899,
    			0.48296291314469025};
    	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

    	double hp1_a[] = {-0.48296291314469025, 0.83651630373746899, -0.22414386804185735,
    			-0.12940952255092145};
    	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

    	double lp2_a[] = {0.48296291314469025, 0.83651630373746899, 0.22414386804185735,
    			-0.12940952255092145};
    	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

    	 double hp2_a[] = {-0.12940952255092145, -0.22414386804185735, 0.83651630373746899,
    			 -0.48296291314469025};
    	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
    	 return 0;
    }

    else if ( name == "db3"){
    	double lp1_a[] = {0.035226291882100656, -0.085441273882241486, -0.13501102001039084,
    			0.45987750211933132, 0.80689150931333875, 0.33267055295095688};
    	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

    	double hp1_a[] = {-0.33267055295095688, 0.80689150931333875, -0.45987750211933132,
    			-0.13501102001039084, 0.085441273882241486, 0.035226291882100656 };
    	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

    	double lp2_a[] = {0.33267055295095688, 0.80689150931333875, 0.45987750211933132,
    			-0.13501102001039084, -0.085441273882241486, 0.035226291882100656 };
    	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

    	 double hp2_a[] = {0.035226291882100656, 0.085441273882241486, -0.13501102001039084,
    			 -0.45987750211933132, 0.80689150931333875, -0.33267055295095688 };
    	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
    	 return 0;
    }

    else if ( name == "db4"){
    	double lp1_a[] = {-0.010597401784997278, 0.032883011666982945, 0.030841381835986965,
    			-0.18703481171888114, -0.027983769416983849, 0.63088076792959036,
    			0.71484657055254153, 0.23037781330885523 };
    	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

    	double hp1_a[] = {-0.23037781330885523, 0.71484657055254153, -0.63088076792959036,
    			-0.027983769416983849, 0.18703481171888114, 0.030841381835986965,
    			-0.032883011666982945, -0.010597401784997278 };
    	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

    	double lp2_a[] = {0.23037781330885523, 0.71484657055254153, 0.63088076792959036,
    			-0.027983769416983849, -0.18703481171888114, 0.030841381835986965,
    			0.032883011666982945, -0.010597401784997278 };
    	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

    	 double hp2_a[] = {-0.010597401784997278, -0.032883011666982945, 0.030841381835986965,
    			 0.18703481171888114, -0.027983769416983849, -0.63088076792959036,
    			 0.71484657055254153, -0.23037781330885523 };
    	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
    	 return 0;
    }

    else if ( name == "db5"){
    	double lp1_a[] = {0.0033357252850015492, -0.012580751999015526, -0.0062414902130117052,
    			0.077571493840065148, -0.03224486958502952, -0.24229488706619015,
    			0.13842814590110342, 0.72430852843857441, 0.60382926979747287,
                0.16010239797412501 };
    	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

    	double hp1_a[] = {-0.16010239797412501, 0.60382926979747287, -0.72430852843857441,
    			0.13842814590110342, 0.24229488706619015, -0.03224486958502952,
    			-0.077571493840065148, -0.0062414902130117052, 0.012580751999015526,
    			0.0033357252850015492 };
    	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

    	double lp2_a[] = {0.16010239797412501, 0.60382926979747287, 0.72430852843857441,
    			0.13842814590110342, -0.24229488706619015, -0.03224486958502952,
    			0.077571493840065148, -0.0062414902130117052, -0.012580751999015526,
    			0.0033357252850015492 };
    	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

    	 double hp2_a[] = {0.0033357252850015492, 0.012580751999015526, -0.0062414902130117052,
    			 -0.077571493840065148, -0.03224486958502952, 0.24229488706619015,
    			 0.13842814590110342, -0.72430852843857441, 0.60382926979747287,
    			 -0.16010239797412501 };
    	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
    	 return 0;
    }

    else if ( name == "db6"){
    	double lp1_a[] = {-0.0010773010849955799,
    			0.0047772575110106514,
    			0.0005538422009938016,
    			-0.031582039318031156,
    			0.027522865530016288,
    			0.097501605587079362,
    			-0.12976686756709563,
    			-0.22626469396516913,
    			0.3152503517092432,
    			0.75113390802157753,
    			0.49462389039838539,
    			0.11154074335008017
};
    	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

    	double hp1_a[] = {-0.11154074335008017,
    			0.49462389039838539,
    			-0.75113390802157753,
    			0.3152503517092432,
    			0.22626469396516913,
    			-0.12976686756709563,
    			-0.097501605587079362,
    			0.027522865530016288,
    			0.031582039318031156,
    			0.0005538422009938016,
    			-0.0047772575110106514,
    			-0.0010773010849955799
};
    	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

    	double lp2_a[] = {0.11154074335008017,
    			0.49462389039838539,
    			0.75113390802157753,
    			0.3152503517092432,
    			-0.22626469396516913,
    			-0.12976686756709563,
    			0.097501605587079362,
    			0.027522865530016288,
    			-0.031582039318031156,
    			0.0005538422009938016,
    			0.0047772575110106514,
    			-0.0010773010849955799
};
    	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

    	 double hp2_a[] = {-0.0010773010849955799,
    			 -0.0047772575110106514,
    			 0.0005538422009938016,
    			 0.031582039318031156,
    			 0.027522865530016288,
    			 -0.097501605587079362,
    			 -0.12976686756709563,
    			 0.22626469396516913,
    			 0.3152503517092432,
    			 -0.75113390802157753,
    			 0.49462389039838539,
    			 -0.11154074335008017
};
    	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
    	 return 0;
    }

    else if ( name == "db7"){
        	double lp1_a[] = {0.00035371380000103988,
        			-0.0018016407039998328,
        			0.00042957797300470274,
        			0.012550998556013784,
        			-0.01657454163101562,
        			-0.038029936935034633,
        			0.080612609151065898,
        			0.071309219267050042,
        			-0.22403618499416572,
        			-0.14390600392910627,
        			0.4697822874053586,
        			0.72913209084655506,
        			0.39653931948230575,
        			0.077852054085062364
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {-0.077852054085062364,
        			0.39653931948230575,
        			-0.72913209084655506,
        			0.4697822874053586,
        			0.14390600392910627,
        			-0.22403618499416572,
        			-0.071309219267050042,
        			0.080612609151065898,
        			0.038029936935034633,
        			-0.01657454163101562,
        			-0.012550998556013784,
        			0.0004295779730047027,
        			0.0018016407039998328,
        			0.00035371380000103988
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.077852054085062364,
        			0.39653931948230575,
        			0.72913209084655506,
        			0.4697822874053586,
        			-0.14390600392910627,
        			-0.22403618499416572,
        			0.071309219267050042,
        			0.080612609151065898,
        			-0.038029936935034633,
        			-0.01657454163101562,
        			0.012550998556013784,
        			0.00042957797300470274,
        			-0.0018016407039998328,
        			0.00035371380000103988
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {0.00035371380000103988,
        			 0.0018016407039998328,
        			 0.00042957797300470274,
        			 -0.01255099855601378,
        			 -0.01657454163101562,
        			 0.038029936935034633,
        			 0.080612609151065898,
        			 -0.071309219267050042,
        			 -0.22403618499416572,
        			 0.14390600392910627,
        			 0.4697822874053586,
        			 -0.72913209084655506,
        			 0.39653931948230575,
        			 -0.077852054085062364
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "db8"){
        	double lp1_a[] = {-0.00011747678400228192,
        			0.00067544940599855677,
        			-0.00039174037299597711,
        			-0.0048703529930106603,
        			0.0087460940470156547,
        			0.013981027917015516,
        			-0.044088253931064719,
        			-0.017369301002022108,
        			0.12874742662018601,
        			0.00047248457399797254,
        			-0.28401554296242809,
        			-0.015829105256023893,
        			0.58535468365486909,
        			0.67563073629801285,
        			0.31287159091446592,
        			0.054415842243081609
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {-0.054415842243081609,
        			0.31287159091446592,
        			-0.67563073629801285,
        			0.58535468365486909,
        			0.015829105256023893,
        			-0.28401554296242809,
        			-0.00047248457399797254,
        			0.12874742662018601,
        			0.017369301002022108,
        			-0.044088253931064719,
        			-0.013981027917015516,
        			0.0087460940470156547,
        			0.0048703529930106603,
        			-0.00039174037299597711,
        			-0.00067544940599855677,
        			-0.00011747678400228192
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.054415842243081609,
        			0.31287159091446592,
        			0.67563073629801285,
        			0.58535468365486909,
        			-0.015829105256023893,
        			-0.28401554296242809,
        			0.00047248457399797254,
        			0.12874742662018601,
        			-0.017369301002022108,
        			-0.044088253931064719,
        			0.013981027917015516,
        			0.0087460940470156547,
        			-0.0048703529930106603,
        			-0.00039174037299597711,
        			0.00067544940599855677,
        			-0.00011747678400228192
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {-0.00011747678400228192,
        			 -0.00067544940599855677,
        			 -0.00039174037299597711,
        			 0.0048703529930106603,
        			 0.0087460940470156547,
        			 -0.013981027917015516,
        			 -0.044088253931064719,
        			 0.017369301002022108,
        			 0.12874742662018601,
        			 -0.00047248457399797254,
        			 -0.28401554296242809,
        			 0.015829105256023893,
        			 0.58535468365486909,
        			 -0.67563073629801285,
        			 0.31287159091446592,
        			 -0.054415842243081609
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "db9"){
        	double lp1_a[] = {3.9347319995026124e-05,
        			-0.00025196318899817888,
        			0.00023038576399541288,
        			0.0018476468829611268,
        			-0.0042815036819047227,
        			-0.004723204757894831,
        			0.022361662123515244,
        			0.00025094711499193845,
        			-0.067632829059523988,
        			0.030725681478322865,
        			0.14854074933476008,
        			-0.096840783220879037,
        			-0.29327378327258685,
        			0.13319738582208895,
        			0.65728807803663891,
        			0.6048231236767786,
        			0.24383467463766728,
        			0.038077947363167282
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {-0.038077947363167282,
        			0.24383467463766728,
        			-0.6048231236767786,
        			0.65728807803663891,
        			-0.13319738582208895,
        			-0.29327378327258685,
        			0.096840783220879037,
        			0.14854074933476008,
        			-0.030725681478322865,
        			-0.067632829059523988,
        			-0.00025094711499193845,
        			0.022361662123515244,
        			0.004723204757894831,
        			-0.0042815036819047227,
        			-0.0018476468829611268,
        			0.00023038576399541288,
        			0.00025196318899817888,
        			3.9347319995026124e-05
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.038077947363167282,
        			0.24383467463766728,
        			0.6048231236767786,
        			0.65728807803663891,
        			0.13319738582208895,
        			-0.29327378327258685,
        			-0.096840783220879037,
        			0.14854074933476008,
        			0.030725681478322865,
        			-0.067632829059523988,
        			0.00025094711499193845,
        			0.022361662123515244,
        			-0.004723204757894831,
        			-0.0042815036819047227,
        			0.0018476468829611268,
        			0.00023038576399541288,
        			-0.00025196318899817888,
        			3.9347319995026124e-05
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {3.9347319995026124e-05,
        			 0.00025196318899817888,
        			 0.00023038576399541288,
        			 -0.0018476468829611268,
        			 -0.0042815036819047227,
        			 0.004723204757894831,
        			 0.022361662123515244,
        			 -0.00025094711499193845,
        			 -0.067632829059523988,
        			 -0.030725681478322865,
        			 0.14854074933476008,
        			 0.096840783220879037,
        			 -0.29327378327258685,
        			 -0.13319738582208895,
        			 0.65728807803663891,
        			 -0.6048231236767786,
        			 0.24383467463766728,
        			 -0.038077947363167282
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "db10"){
        	double lp1_a[] = {-1.3264203002354869e-05,
        			9.3588670001089845e-05,
        			-0.0001164668549943862,
        			-0.00068585669500468248,
        			0.0019924052949908499,
        			0.0013953517469940798,
        			-0.010733175482979604,
        			0.0036065535669883944,
        			0.033212674058933238,
        			-0.029457536821945671,
        			-0.071394147165860775,
        			0.093057364603806592,
        			0.12736934033574265,
        			-0.19594627437659665,
        			-0.24984642432648865,
        			0.28117234366042648,
        			0.68845903945259213,
        			0.52720118893091983,
        			0.18817680007762133,
        			0.026670057900950818
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {-0.026670057900950818,
        			0.18817680007762133,
        			-0.52720118893091983,
        			0.68845903945259213,
        			-0.28117234366042648,
        			-0.24984642432648865,
        			0.19594627437659665,
        			0.12736934033574265,
        			-0.093057364603806592,
        			-0.071394147165860775,
        			0.029457536821945671,
        			0.033212674058933238,
        			-0.0036065535669883944,
        			-0.010733175482979604,
        			-0.0013953517469940798,
        			0.0019924052949908499,
        			0.00068585669500468248,
        			-0.0001164668549943862,
        			-9.3588670001089845e-05,
        			-1.3264203002354869e-05
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.026670057900950818,
        			0.18817680007762133,
        			0.52720118893091983,
        			0.68845903945259213,
        			0.28117234366042648,
        			-0.24984642432648865,
        			-0.19594627437659665,
        			0.12736934033574265,
        			0.093057364603806592,
        			-0.071394147165860775,
        			-0.029457536821945671,
        			0.033212674058933238,
        			0.0036065535669883944,
        			-0.010733175482979604,
        			0.0013953517469940798,
        			0.0019924052949908499,
        			-0.00068585669500468248,
        			-0.0001164668549943862,
        			9.3588670001089845e-05,
        			-1.3264203002354869e-05
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {-1.3264203002354869e-05,
        			 -9.3588670001089845e-05,
        			 -0.0001164668549943862,
        			 0.00068585669500468248,
        			 0.0019924052949908499,
        			 -0.0013953517469940798,
        			 -0.010733175482979604,
        			 -0.0036065535669883944,
        			 0.033212674058933238,
        			 0.029457536821945671,
        			 -0.071394147165860775,
        			 -0.093057364603806592,
        			 0.12736934033574265,
        			 0.19594627437659665,
        			 -0.24984642432648865,
        			 -0.28117234366042648,
        			 0.68845903945259213,
        			 -0.52720118893091983,
        			 0.18817680007762133,
        			 -0.026670057900950818
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "db12"){
            	double lp1_a[] = {-1.5290717580684923e-06,
            			1.2776952219379579e-05,
            			-2.4241545757030318e-05,
            			-8.8504109208203182e-05,
            			0.00038865306282092672,
            			6.5451282125215034e-06,
            			-0.0021795036186277044,
            			0.0022486072409952287,
            			0.0067114990087955486,
            			-0.012840825198299882,
            			-0.01221864906974642,
            			0.041546277495087637,
            			0.010849130255828966,
            			-0.09643212009649671,
            			0.0053595696743599965,
            			0.18247860592758275,
            			-0.023779257256064865,
            			-0.31617845375277914,
            			-0.044763885653777619,
            			0.51588647842780067,
            			0.65719872257929113,
            			0.37735513521420411,
            			0.10956627282118277,
            			0.013112257957229239
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.013112257957229239,
            			0.10956627282118277,
            			-0.37735513521420411,
            			0.65719872257929113,
            			-0.51588647842780067,
            			-0.044763885653777619,
            			0.31617845375277914,
            			-0.023779257256064865,
            			-0.18247860592758275,
            			0.0053595696743599965,
            			0.09643212009649671,
            			0.010849130255828966,
            			-0.041546277495087637,
            			-0.01221864906974642,
            			0.012840825198299882,
            			0.0067114990087955486,
            			-0.0022486072409952287,
            			-0.0021795036186277044,
            			-6.5451282125215034e-06,
            			0.00038865306282092672,
            			8.8504109208203182e-05,
            			-2.4241545757030318e-05,
            			-1.2776952219379579e-05,
            			-1.5290717580684923e-06
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.013112257957229239,
            			0.10956627282118277,
            			0.37735513521420411,
            			0.65719872257929113,
            			0.51588647842780067,
            			-0.044763885653777619,
            			-0.31617845375277914,
            			-0.023779257256064865,
            			0.18247860592758275,
            			0.0053595696743599965,
            			-0.09643212009649671,
            			0.010849130255828966,
            			0.041546277495087637,
            			-0.01221864906974642,
            			-0.012840825198299882,
            			0.0067114990087955486,
            			0.0022486072409952287,
            			-0.0021795036186277044,
            			6.5451282125215034e-06,
            			0.00038865306282092672,
            			-8.8504109208203182e-05,
            			-2.4241545757030318e-05,
            			1.2776952219379579e-05,
            			-1.5290717580684923e-06
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-1.5290717580684923e-06,
            			 -1.2776952219379579e-05,
            			 -2.4241545757030318e-05,
            			 8.8504109208203182e-05,
            			 0.00038865306282092672,
            			 -6.5451282125215034e-06,
            			 -0.0021795036186277044,
            			 -0.0022486072409952287,
            			 0.0067114990087955486,
            			 0.012840825198299882,
            			 -0.01221864906974642,
            			 -0.041546277495087637,
            			 0.010849130255828966,
            			 0.09643212009649671,
            			 0.0053595696743599965,
            			 -0.18247860592758275,
            			 -0.023779257256064865,
            			 0.31617845375277914,
            			 -0.044763885653777619,
            			 -0.51588647842780067,
            			 0.65719872257929113,
            			 -0.37735513521420411,
            			 0.10956627282118277,
            			 -0.013112257957229239
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "db13"){
            	double lp1_a[] = {5.2200350984547998e-07,
            			-4.7004164793608082e-06,
            			1.0441930571407941e-05,
            			3.0678537579324358e-05,
            			-0.00016512898855650571,
            			4.9251525126285676e-05,
            			0.00093232613086724904,
            			-0.0013156739118922766,
            			-0.002761911234656831,
            			0.0072555894016171187,
            			0.0039239414487955773,
            			-0.023831420710327809,
            			0.0023799722540522269,
            			0.056139477100276156,
            			-0.026488406475345658,
            			-0.10580761818792761,
            			0.072948933656788742,
            			0.17947607942935084,
            			-0.12457673075080665,
            			-0.31497290771138414,
            			0.086985726179645007,
            			0.58888957043121193,
            			0.61105585115878114,
            			0.31199632216043488,
            			0.082861243872901946,
            			0.0092021335389622788
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.0092021335389622788,
            			0.082861243872901946,
            			-0.31199632216043488,
            			0.61105585115878114,
            			-0.58888957043121193,
            			0.086985726179645007,
            			0.31497290771138414,
            			-0.12457673075080665,
            			-0.17947607942935084,
            			0.072948933656788742,
            			0.10580761818792761,
            			-0.026488406475345658,
            			-0.056139477100276156,
            			0.0023799722540522269,
            			0.023831420710327809,
            			0.0039239414487955773,
            			-0.0072555894016171187,
            			-0.002761911234656831,
            			0.0013156739118922766,
            			0.00093232613086724904,
            			-4.9251525126285676e-05,
            			-0.00016512898855650571,
            			-3.0678537579324358e-05,
            			1.0441930571407941e-05,
            			4.7004164793608082e-06,
            			5.2200350984547998e-07
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0092021335389622788,
            			0.082861243872901946,
            			0.31199632216043488,
            			0.61105585115878114,
            			0.58888957043121193,
            			0.086985726179645007,
            			-0.31497290771138414,
            			-0.12457673075080665,
            			0.17947607942935084,
            			0.072948933656788742,
            			-0.10580761818792761,
            			-0.026488406475345658,
            			0.056139477100276156,
            			0.0023799722540522269,
            			-0.023831420710327809,
            			0.0039239414487955773,
            			0.0072555894016171187,
            			-0.002761911234656831,
            			-0.0013156739118922766,
            			0.00093232613086724904,
            			4.9251525126285676e-05,
            			-0.00016512898855650571,
            			3.0678537579324358e-05,
            			1.0441930571407941e-05,
            			-4.7004164793608082e-06,
            			5.2200350984547998e-07
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {5.2200350984547998e-07,
            			 4.7004164793608082e-06,
            			 1.0441930571407941e-05,
            			 -3.0678537579324358e-05,
            			 -0.00016512898855650571,
            			 -4.9251525126285676e-05,
            			 0.00093232613086724904,
            			 0.0013156739118922766,
            			 -0.002761911234656831,
            			 -0.0072555894016171187,
            			 0.0039239414487955773,
            			 0.023831420710327809,
            			 0.0023799722540522269,
            			 -0.056139477100276156,
            			 -0.026488406475345658,
            			 0.10580761818792761,
            			 0.072948933656788742,
            			 -0.17947607942935084,
            			 -0.12457673075080665,
            			 0.31497290771138414,
            			 0.086985726179645007,
            			 -0.58888957043121193,
            			 0.61105585115878114,
            			 -0.31199632216043488,
            			 0.082861243872901946,
            			 -0.0092021335389622788
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }

    else if ( name == "db11"){
              	double lp1_a[] = {4.4942742772363519e-06,
              			-3.4634984186983789e-05,
              			5.4439074699366381e-05,
              			0.00024915252355281426,
              			-0.00089302325066623663,
              			-0.00030859285881515924,
              			0.0049284176560587777,
              			-0.0033408588730145018,
              			-0.015364820906201324,
              			0.020840904360180039,
              			0.031335090219045313,
              			-0.066438785695020222,
              			-0.04647995511667613,
              			0.14981201246638268,
              			0.066043588196690886,
              			-0.27423084681792875,
              			-0.16227524502747828,
              			0.41196436894789695,
              			0.68568677491617847,
              			0.44989976435603013,
              			0.14406702115061959,
              			0.018694297761470441
  };
              	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

              	double hp1_a[] = {-0.018694297761470441,
              			0.14406702115061959,
              			-0.44989976435603013,
              			0.68568677491617847,
              			-0.41196436894789695,
              			-0.16227524502747828,
              			0.27423084681792875,
              			0.066043588196690886,
              			-0.14981201246638268,
              			-0.04647995511667613,
              			0.066438785695020222,
              			0.031335090219045313,
              			-0.020840904360180039,
              			-0.015364820906201324,
              			0.0033408588730145018,
              			0.0049284176560587777,
              			0.00030859285881515924,
              			-0.00089302325066623663,
              			-0.00024915252355281426,
              			5.4439074699366381e-05,
              			3.4634984186983789e-05,
              			4.4942742772363519e-06
  };
              	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

              	double lp2_a[] = {0.018694297761470441,
              			0.14406702115061959,
              			0.44989976435603013,
              			0.68568677491617847,
              			0.41196436894789695,
              			-0.16227524502747828,
              			-0.27423084681792875,
              			0.066043588196690886,
              			0.14981201246638268,
              			-0.04647995511667613,
              			-0.066438785695020222,
              			0.031335090219045313,
              			0.020840904360180039,
              			-0.015364820906201324,
              			-0.0033408588730145018,
              			0.0049284176560587777,
              			-0.00030859285881515924,
              			-0.00089302325066623663,
              			0.00024915252355281426,
              			5.4439074699366381e-05,
              			-3.4634984186983789e-05,
              			4.4942742772363519e-06
  };
              	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

              	 double hp2_a[] = {4.4942742772363519e-06,
              			 3.4634984186983789e-05,
              			 5.4439074699366381e-05,
              			 -0.00024915252355281426,
              			 -0.00089302325066623663,
              			 0.00030859285881515924,
              			 0.0049284176560587777,
              			 0.0033408588730145018,
              			 -0.015364820906201324,
              			 -0.020840904360180039,
              			 0.031335090219045313,
              			 0.066438785695020222,
              			 -0.04647995511667613,
              			 -0.14981201246638268,
              			 0.066043588196690886,
              			 0.27423084681792875,
              			 -0.16227524502747828,
              			 -0.41196436894789695,
              			 0.68568677491617847,
              			 -0.44989976435603013,
              			 0.14406702115061959,
              			 -0.018694297761470441
  };
              	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
              	 return 0;
              }

    else if ( name == "db14"){
                  	double lp1_a[] = {-1.7871399683109222e-07,
                  			1.7249946753674012e-06,
                  			-4.3897049017804176e-06,
                  			-1.0337209184568496e-05,
                  			6.875504252695734e-05,
                  			-4.1777245770370672e-05,
                  			-0.00038683194731287514,
                  			0.00070802115423540481,
                  			0.001061691085606874,
                  			-0.003849638868019787,
                  			-0.00074621898926387534,
                  			0.012789493266340071,
                  			-0.0056150495303375755,
                  			-0.030185351540353976,
                  			0.026981408307947971,
                  			0.05523712625925082,
                  			-0.071548955503983505,
                  			-0.086748411568110598,
                  			0.13998901658445695,
                  			0.13839521386479153,
                  			-0.21803352999321651,
                  			-0.27168855227867705,
                  			0.21867068775886594,
                  			0.63118784910471981,
                  			0.55430561794077093,
                  			0.25485026779256437,
                  			0.062364758849384874,
                  			0.0064611534600864905
};
                  	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

                  	double hp1_a[] = {-0.0064611534600864905,
                  			0.062364758849384874,
                  			-0.25485026779256437,
                  			0.55430561794077093,
                  			-0.63118784910471981,
                  			0.21867068775886594,
                  			0.27168855227867705,
                  			-0.21803352999321651,
                  			-0.13839521386479153,
                  			0.13998901658445695,
                  			0.086748411568110598,
                  			-0.071548955503983505,
                  			-0.05523712625925082,
                  			0.026981408307947971,
                  			0.030185351540353976,
                  			-0.0056150495303375755,
                  			-0.012789493266340071,
                  			-0.00074621898926387534,
                  			0.003849638868019787,
                  			0.001061691085606874,
                  			-0.00070802115423540481,
                  			-0.00038683194731287514,
                  			4.1777245770370672e-05,
                  			6.875504252695734e-05,
                  			1.0337209184568496e-05,
                  			-4.3897049017804176e-06,
                  			-1.7249946753674012e-06,
                  			-1.7871399683109222e-07
};
                  	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

                  	double lp2_a[] = {0.0064611534600864905,
                  			0.062364758849384874,
                  			0.25485026779256437,
                  			0.55430561794077093,
                  			0.63118784910471981,
                  			0.21867068775886594,
                  			-0.27168855227867705,
                  			-0.21803352999321651,
                  			0.13839521386479153,
                  			0.13998901658445695,
                  			-0.086748411568110598,
                  			-0.071548955503983505,
                  			0.05523712625925082,
                  			0.026981408307947971,
                  			-0.030185351540353976,
                  			-0.0056150495303375755,
                  			0.012789493266340071,
                  			-0.00074621898926387534,
                  			-0.003849638868019787,
                  			0.001061691085606874,
                  			0.00070802115423540481,
                  			-0.00038683194731287514,
                  			-4.1777245770370672e-05,
                  			6.875504252695734e-05,
                  			-1.0337209184568496e-05,
                  			-4.3897049017804176e-06,
                  			1.7249946753674012e-06,
                  			-1.7871399683109222e-07
};
                  	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

                  	 double hp2_a[] = {-1.7871399683109222e-07,
                  			-1.7249946753674012e-06,
                  			-4.3897049017804176e-06,
                  			1.0337209184568496e-05,
                  			6.875504252695734e-05,
                  			4.1777245770370672e-05,
                  			-0.00038683194731287514,
                  			-0.00070802115423540481,
                  			0.001061691085606874,
                  			0.003849638868019787,
                  			-0.00074621898926387534,
                  			-0.012789493266340071,
                  			-0.0056150495303375755,
                  			0.030185351540353976,
                  			0.026981408307947971,
                  			-0.05523712625925082,
                  			-0.071548955503983505,
                  			0.086748411568110598,
                  			0.13998901658445695,
                  			-0.13839521386479153,
                  			-0.21803352999321651,
                  			0.27168855227867705,
                  			0.21867068775886594,
                  			-0.63118784910471981,
                  			0.55430561794077093,
                  			-0.25485026779256437,
                  			0.062364758849384874,
                  			-0.0064611534600864905
};
                  	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
                  	 return 0;
                  }
    else if ( name == "db15"){
            	double lp1_a[] = {6.1333599133037138e-08,
            			-6.3168823258794506e-07,
            			1.8112704079399406e-06,
            			3.3629871817363823e-06,
            			-2.8133296266037558e-05,
            			2.579269915531323e-05,
            			0.00015589648992055726,
            			-0.00035956524436229364,
            			-0.00037348235413726472,
            			0.0019433239803823459,
            			-0.00024175649075894543,
            			-0.0064877345603061454,
            			0.0051010003604228726,
            			0.015083918027862582,
            			-0.020810050169636805,
            			-0.025767007328366939,
            			0.054780550584559995,
            			0.033877143923563204,
            			-0.11112093603713753,
            			-0.039666176555733602,
            			0.19014671400708816,
            			0.065282952848765688,
            			-0.28888259656686216,
            			-0.19320413960907623,
            			0.33900253545462167,
            			0.64581314035721027,
            			0.49263177170797529,
            			0.20602386398692688,
            			0.046743394892750617,
            			0.0045385373615773762
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.0045385373615773762,
            			0.046743394892750617,
            			-0.20602386398692688,
            			0.49263177170797529,
            			-0.64581314035721027,
            			0.33900253545462167,
            			0.19320413960907623,
            			-0.28888259656686216,
            			-0.065282952848765688,
            			0.19014671400708816,
            			0.039666176555733602,
            			-0.11112093603713753,
            			-0.033877143923563204,
            			0.054780550584559995,
            			0.025767007328366939,
            			-0.020810050169636805,
            			-0.015083918027862582,
            			0.0051010003604228726,
            			0.0064877345603061454,
            			-0.00024175649075894543,
            			-0.0019433239803823459,
            			-0.00037348235413726472,
            			0.00035956524436229364,
            			0.00015589648992055726,
            			-2.579269915531323e-05,
            			-2.8133296266037558e-05,
            			-3.3629871817363823e-06,
            			1.8112704079399406e-06,
            			6.3168823258794506e-07,
            			6.1333599133037138e-08
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0045385373615773762,
            			0.046743394892750617,
            			0.20602386398692688,
            			0.49263177170797529,
            			0.64581314035721027,
            			0.33900253545462167,
            			-0.19320413960907623,
            			-0.28888259656686216,
            			0.065282952848765688,
            			0.19014671400708816,
            			-0.039666176555733602,
            			-0.11112093603713753,
            			0.033877143923563204,
            			0.054780550584559995,
            			-0.025767007328366939,
            			-0.020810050169636805,
            			0.015083918027862582,
            			0.0051010003604228726,
            			-0.0064877345603061454,
            			-0.00024175649075894543,
            			0.0019433239803823459,
            			-0.00037348235413726472,
            			-0.00035956524436229364,
            			0.00015589648992055726,
            			2.579269915531323e-05,
            			-2.8133296266037558e-05,
            			3.3629871817363823e-06,
            			1.8112704079399406e-06,
            			-6.3168823258794506e-07,
            			6.1333599133037138e-08
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {6.1333599133037138e-08,
            			 6.3168823258794506e-07,
            			 1.8112704079399406e-06,
            			 -3.3629871817363823e-06,
            			 -2.8133296266037558e-05,
            			 -2.579269915531323e-05,
            			 0.00015589648992055726,
            			 0.00035956524436229364,
            			 -0.00037348235413726472,
            			 -0.0019433239803823459,
            			 -0.00024175649075894543,
            			 0.0064877345603061454,
            			 0.0051010003604228726,
            			 -0.015083918027862582,
            			 -0.020810050169636805,
            			 0.025767007328366939,
            			 0.054780550584559995,
            			 -0.033877143923563204,
            			 -0.11112093603713753,
            			 0.039666176555733602,
            			 0.19014671400708816,
            			 -0.065282952848765688,
            			 -0.28888259656686216,
            			 0.19320413960907623,
            			 0.33900253545462167,
            			 -0.64581314035721027,
            			 0.49263177170797529,
            			 -0.20602386398692688,
            			 0.046743394892750617,
            			 -0.0045385373615773762
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior1.1"){
        	double lp1_a[] = {0.70710678118654757,
        			0.70710678118654757
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {-0.70710678118654757,
        			0.70710678118654757
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.70710678118654757,
        			0.70710678118654757
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {0.70710678118654757,
        			 -0.70710678118654757
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "bior1.3"){
            	double lp1_a[] = {-0.088388347648318447,
            			0.088388347648318447,
            			0.70710678118654757,
            			0.70710678118654757,
            			0.088388347648318447,
            			-0.088388347648318447,
    };
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			-0.70710678118654757,
            			0.70710678118654757,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
           			 0.0,
           			 0.70710678118654757,
           			 0.70710678118654757,
           			 0.0,
           			 0.0
    };
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.088388347648318447,
            			 -0.088388347648318447,
            			 0.70710678118654757,
            			 -0.70710678118654757,
            			 0.088388347648318447,
            			 0.088388347648318447
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }

    else if ( name == "bior1.5"){
        	double lp1_a[] = {0.01657281518405971,
        			-0.01657281518405971,
        			-0.12153397801643787,
        			0.12153397801643787,
        			0.70710678118654757,
        			0.70710678118654757,
        			0.12153397801643787,
        			-0.12153397801643787,
        			-0.01657281518405971,
        			0.01657281518405971
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {0.0,
        			0.0,
        			0.0,
        			0.0,
        			-0.70710678118654757,
        			0.70710678118654757,
        			0.0,
        			0.0,
        			0.0,
        			0.0
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.0,
        			0.0,
        			0.0,
        			0.0,
        			0.70710678118654757,
        			0.70710678118654757,
        			0.0,
        			0.0,
        			0.0,
        			0.0
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {0.01657281518405971,
        			 0.01657281518405971,
        			 -0.12153397801643787,
        			 -0.12153397801643787,
        			 0.70710678118654757,
        			 -0.70710678118654757,
        			 0.12153397801643787,
        			 0.12153397801643787,
        			 -0.01657281518405971,
        			 -0.01657281518405971
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "bior2.2"){
        	double lp1_a[] = {0.0,
        			-0.17677669529663689,
        			0.35355339059327379,
        			1.0606601717798214,
        			0.35355339059327379,
        			-0.17677669529663689
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {0.0,
        			0.35355339059327379,
        			-0.70710678118654757,
        			0.35355339059327379,
        			0.0,
        			0.0
};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.0,
        			0.35355339059327379,
        			0.70710678118654757,
        			0.35355339059327379,
        			0.0,
        			0.0
};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {0.0,
        			 0.17677669529663689,
        			 0.35355339059327379,
        			 -1.0606601717798214,
        			 0.35355339059327379,
        			 0.17677669529663689

};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "bior2.4"){
        	double lp1_a[] = {0.0,
        			0.033145630368119419,
        			-0.066291260736238838,
        			-0.17677669529663689,
        			0.4198446513295126,
        			0.99436891104358249,
        			0.4198446513295126,
        			-0.17677669529663689,
        			-0.066291260736238838,
        			0.033145630368119419
};
        	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

        	double hp1_a[] = {0.0,
        			0.0,
        			0.0,
        			0.35355339059327379,
        			-0.70710678118654757,
        			0.35355339059327379,
        			0.0,
        			0.0,
        			0.0,
        			0.0

};
        	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

        	double lp2_a[] = {0.0,
        			0.0,
        			0.0,
        			0.35355339059327379,
        			0.70710678118654757,
        			0.35355339059327379,
        			0.0,
        			0.0,
        			0.0,
        			0.0

};
        	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

        	 double hp2_a[] = {0.0,
        			 -0.033145630368119419,
        			 -0.066291260736238838,
        			 0.17677669529663689,
        			 0.4198446513295126,
        			 -0.99436891104358249,
        			 0.4198446513295126,
        			 0.17677669529663689,
        			 -0.066291260736238838,
        			 -0.033145630368119419
};
        	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
        	 return 0;
        }

    else if ( name == "bior2.6"){
            	double lp1_a[] = {0.0,
            			-0.0069053396600248784,
            			0.013810679320049757,
            			0.046956309688169176,
            			-0.10772329869638811,
            			-0.16987135563661201,
            			0.44746600996961211,
            			0.96674755240348298,
            			0.44746600996961211,
            			-0.16987135563661201,
            			-0.10772329869638811,
            			0.046956309688169176,
            			0.013810679320049757,
            			-0.0069053396600248784
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.35355339059327379,
            			-0.70710678118654757,
            			0.35355339059327379,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.35355339059327379,
            			0.70710678118654757,
            			0.35355339059327379,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0,
            			 0.0069053396600248784,
            			 0.013810679320049757,
            			 -0.046956309688169176,
            			 -0.10772329869638811,
            			 0.16987135563661201,
            			 0.44746600996961211,
            			 -0.96674755240348298,
            			 0.44746600996961211,
            			 0.16987135563661201,
            			 -0.10772329869638811,
            			 -0.046956309688169176,
            			 0.013810679320049757,
            			 0.0069053396600248784
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior2.8"){
            	double lp1_a[] = {0.0,
            			0.0015105430506304422,
            			-0.0030210861012608843,
            			-0.012947511862546647,
            			0.028916109826354178,
            			0.052998481890690945,
            			-0.13491307360773608,
            			-0.16382918343409025,
            			0.46257144047591658,
            			0.95164212189717856,
            			0.46257144047591658,
            			-0.16382918343409025,
            			-0.13491307360773608,
            			0.052998481890690945,
            			0.028916109826354178,
            			-0.012947511862546647,
            			-0.0030210861012608843,
            			0.0015105430506304422
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.35355339059327379,
            			-0.70710678118654757,
            			0.35355339059327379,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.35355339059327379,
            			0.70710678118654757,
            			0.35355339059327379,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0,
            			 -0.0015105430506304422,
            			 -0.0030210861012608843,
            			 0.012947511862546647,
            			 0.028916109826354178,
            			 -0.052998481890690945,
            			 -0.13491307360773608,
            			 0.16382918343409025,
            			 0.46257144047591658,
            			 -0.95164212189717856,
            			 0.46257144047591658,
            			 0.16382918343409025,
            			 -0.13491307360773608,
            			 -0.052998481890690945,
            			 0.028916109826354178,
            			 0.012947511862546647,
            			 -0.0030210861012608843,
            			 -0.0015105430506304422
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }

    else if ( name == "bior3.1"){
            	double lp1_a[] = {-0.35355339059327379,
            			1.0606601717798214,
            			1.0606601717798214,
            			-0.35355339059327379
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.17677669529663689,
            			0.53033008588991071,
            			-0.53033008588991071,
            			0.17677669529663689
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.17677669529663689,
            			0.53033008588991071,
            			0.53033008588991071,
            			0.17677669529663689
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.35355339059327379,
            			 -1.0606601717798214,
            			 1.0606601717798214,
            			 0.35355339059327379
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior3.3"){
            	double lp1_a[] = {0.066291260736238838,
            			-0.19887378220871652,
            			-0.15467960838455727,
            			0.99436891104358249,
            			0.99436891104358249,
            			-0.15467960838455727,
            			-0.19887378220871652,
            			0.066291260736238838
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			-0.17677669529663689,
            			0.53033008588991071,
            			-0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.17677669529663689,
            			0.53033008588991071,
            			0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.066291260736238838,
            			 0.19887378220871652,
            			 -0.15467960838455727,
            			 -0.99436891104358249,
            			 0.99436891104358249,
            			 0.15467960838455727,
            			 -0.19887378220871652,
            			 -0.066291260736238838
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior3.5"){
            	double lp1_a[] = {-0.013810679320049757,
            			0.041432037960149271,
            			0.052480581416189075,
            			-0.26792717880896527,
            			-0.071815532464258744,
            			0.96674755240348298,
            			0.96674755240348298,
            			-0.071815532464258744,
            			-0.26792717880896527,
            			0.052480581416189075,
            			0.041432037960149271,
            			-0.013810679320049757
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			-0.17677669529663689,
            			0.53033008588991071,
            			-0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.17677669529663689,
            			0.53033008588991071,
            			0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.013810679320049757,
            			 -0.041432037960149271,
            			 0.052480581416189075,
            			 0.26792717880896527,
            			 -0.071815532464258744,
            			 -0.96674755240348298,
            			 0.96674755240348298,
            			 0.071815532464258744,
            			 -0.26792717880896527,
            			 -0.052480581416189075,
            			 0.041432037960149271,
            			 0.013810679320049757
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }

    else if ( name == "bior3.7"){
            	double lp1_a[] = {0.0030210861012608843,
            			-0.0090632583037826529,
            			-0.016831765421310641,
            			0.074663985074019001,
            			0.031332978707362888,
            			-0.301159125922835,
            			-0.026499240945345472,
            			0.95164212189717856,
            			0.95164212189717856,
            			-0.026499240945345472,
            			-0.301159125922835,
            			0.031332978707362888,
            			0.074663985074019001,
            			-0.016831765421310641,
            			-0.0090632583037826529,
            			0.0030210861012608843
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			-0.17677669529663689,
            			0.53033008588991071,
            			-0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.17677669529663689,
            			0.53033008588991071,
            			0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0030210861012608843,
            			 0.0090632583037826529,
            			 -0.016831765421310641,
            			 -0.074663985074019001,
            			 0.031332978707362888,
            			 0.301159125922835,
            			 -0.026499240945345472,
            			 -0.95164212189717856,
            			 0.95164212189717856,
            			 0.026499240945345472,
            			 -0.301159125922835,
            			 -0.031332978707362888,
            			 0.074663985074019001,
            			 0.016831765421310641,
            			 -0.0090632583037826529,
            			 -0.0030210861012608843
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior3.9"){
            	double lp1_a[] = {-0.00067974437278369901,
            			0.0020392331183510968,
            			0.0050603192196119811,
            			-0.020618912641105536,
            			-0.014112787930175846,
            			0.09913478249423216,
            			0.012300136269419315,
            			-0.32019196836077857,
            			0.0020500227115698858,
            			0.94212570067820678,
            			0.94212570067820678,
            			0.0020500227115698858,
            			-0.32019196836077857,
            			0.012300136269419315,
            			0.09913478249423216,
            			-0.014112787930175846,
            			-0.020618912641105536,
            			0.0050603192196119811,
            			0.0020392331183510968,
            			-0.00067974437278369901
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			-0.17677669529663689,
            			0.53033008588991071,
            			-0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.17677669529663689,
            			0.53033008588991071,
            			0.53033008588991071,
            			0.17677669529663689,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.00067974437278369901,
            			 -0.0020392331183510968,
            			 0.0050603192196119811,
            			 0.020618912641105536,
            			 -0.014112787930175846,
            			 -0.09913478249423216,
            			 0.012300136269419315,
            			 0.32019196836077857,
            			 0.0020500227115698858,
            			 -0.94212570067820678,
            			 0.94212570067820678,
            			 -0.0020500227115698858,
            			 -0.32019196836077857,
            			 -0.012300136269419315,
            			 0.09913478249423216,
            			 0.014112787930175846,
            			 -0.020618912641105536,
            			 -0.0050603192196119811,
            			 0.0020392331183510968,
            			 0.00067974437278369901
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior4.4"){
            	double lp1_a[] = {0.0,
            			0.03782845550726404,
            			-0.023849465019556843,
            			-0.11062440441843718,
            			0.37740285561283066,
            			0.85269867900889385,
            			0.37740285561283066,
            			-0.11062440441843718,
            			-0.023849465019556843,
            			0.03782845550726404
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			-0.064538882628697058,
            			0.040689417609164058,
            			0.41809227322161724,
            			-0.7884856164055829,
            			0.41809227322161724,
            			0.040689417609164058,
            			-0.064538882628697058,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			-0.064538882628697058,
            			-0.040689417609164058,
            			0.41809227322161724,
            			0.7884856164055829,
            			0.41809227322161724,
            			-0.040689417609164058,
            			-0.064538882628697058,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0,
            			 -0.03782845550726404,
            			 -0.023849465019556843,
            			 0.11062440441843718,
            			 0.37740285561283066,
            			 -0.85269867900889385,
            			 0.37740285561283066,
            			 0.11062440441843718,
            			 -0.023849465019556843,
            			 -0.03782845550726404
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior5.5"){
            	double lp1_a[] = {0.0,
            			0.0,
            			0.03968708834740544,
            			0.0079481086372403219,
            			-0.054463788468236907,
            			0.34560528195603346,
            			0.73666018142821055,
            			0.34560528195603346,
            			-0.054463788468236907,
            			0.0079481086372403219,
            			0.03968708834740544,
            			0.0
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.013456709459118716,
            			-0.0026949668801115071,
            			0.13670658466432914,
            			-0.093504697400938863,
            			-0.47680326579848425,
            			0.89950610974864842,
            			-0.47680326579848425,
            			-0.093504697400938863,
            			0.13670658466432914,
            			-0.0026949668801115071,
            			-0.013456709459118716,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.013456709459118716,
            			-0.0026949668801115071,
            			-0.13670658466432914,
            			-0.093504697400938863,
            			0.47680326579848425,
            			0.89950610974864842,
            			0.47680326579848425,
            			-0.093504697400938863,
            			-0.13670658466432914,
            			-0.0026949668801115071,
            			0.013456709459118716,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0,
            			 0.0,
            			 0.03968708834740544,
            			 -0.0079481086372403219,
            			 -0.054463788468236907,
            			 -0.34560528195603346,
            			 0.73666018142821055,
            			 -0.34560528195603346,
            			 -0.054463788468236907,
            			 -0.0079481086372403219,
            			 0.03968708834740544,
            			 0.0
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "bior6.8"){
            	double lp1_a[] = {0.0,
            			0.0019088317364812906,
            			-0.0019142861290887667,
            			-0.016990639867602342,
            			0.01193456527972926,
            			0.04973290349094079,
            			-0.077263173167204144,
            			-0.09405920349573646,
            			0.42079628460982682,
            			0.82592299745840225,
            			0.42079628460982682,
            			-0.09405920349573646,
            			-0.077263173167204144,
            			0.04973290349094079,
            			0.01193456527972926,
            			-0.016990639867602342,
            			-0.0019142861290887667,
            			0.0019088317364812906
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0,
            			0.0,
            			0.0,
            			0.014426282505624435,
            			-0.014467504896790148,
            			-0.078722001062628819,
            			0.040367979030339923,
            			0.41784910915027457,
            			-0.75890772945365415,
            			0.41784910915027457,
            			0.040367979030339923,
            			-0.078722001062628819,
            			-0.014467504896790148,
            			0.014426282505624435,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.0,
            			0.0,
            			0.0,
            			0.014426282505624435,
            			0.014467504896790148,
            			-0.078722001062628819,
            			-0.040367979030339923,
            			0.41784910915027457,
            			0.75890772945365415,
            			0.41784910915027457,
            			-0.040367979030339923,
            			-0.078722001062628819,
            			0.014467504896790148,
            			0.014426282505624435,
            			0.0,
            			0.0,
            			0.0,
            			0.0
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {0.0,
            			 -0.0019088317364812906,
            			 -0.0019142861290887667,
            			 0.016990639867602342,
            			 0.01193456527972926,
            			 -0.04973290349094079,
            			 -0.077263173167204144,
            			 0.09405920349573646,
            			 0.42079628460982682,
            			 -0.82592299745840225,
            			 0.42079628460982682,
            			 0.09405920349573646,
            			 -0.077263173167204144,
            			 -0.04973290349094079,
            			 0.01193456527972926,
            			 0.016990639867602342,
            			 -0.0019142861290887667,
            			 -0.0019088317364812906
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }

    else if ( name == "coif1"){
            	double lp1_a[] = {-0.01565572813546454,
            			-0.072732619512853897,
            			0.38486484686420286,
            			0.85257202021225542,
            			0.33789766245780922,
            			-0.072732619512853897
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.072732619512853897,
            			0.33789766245780922,
            			-0.85257202021225542,
            			0.38486484686420286,
            			0.072732619512853897,
            			-0.01565572813546454
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {-0.072732619512853897,
            			0.33789766245780922,
            			0.85257202021225542,
            			0.38486484686420286,
            			-0.072732619512853897,
            			-0.01565572813546454
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.01565572813546454,
            			 0.072732619512853897,
            			 0.38486484686420286,
            			 -0.85257202021225542,
            			 0.33789766245780922,
            			 0.072732619512853897
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "coif2"){
            	double lp1_a[] = {-0.00072054944536451221,
            			-0.0018232088707029932,
            			0.0056114348193944995,
            			0.023680171946334084,
            			-0.059434418646456898,
            			-0.076488599078306393,
            			0.41700518442169254,
            			0.81272363544554227,
            			0.38611006682116222,
            			-0.067372554721963018,
            			-0.041464936781759151,
            			0.016387336463522112
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.016387336463522112,
            			-0.041464936781759151,
            			0.067372554721963018,
            			0.38611006682116222,
            			-0.81272363544554227,
            			0.41700518442169254,
            			0.076488599078306393,
            			-0.059434418646456898,
            			-0.023680171946334084,
            			0.0056114348193944995,
            			0.0018232088707029932,
            			-0.00072054944536451221
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.016387336463522112,
            			-0.041464936781759151,
            			-0.067372554721963018,
            			0.38611006682116222,
            			0.81272363544554227,
            			0.41700518442169254,
            			-0.076488599078306393,
            			-0.059434418646456898,
            			0.023680171946334084,
            			0.0056114348193944995,
            			-0.0018232088707029932,
            			-0.00072054944536451221
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-0.00072054944536451221,
            			 0.0018232088707029932,
            			 0.0056114348193944995,
            			 -0.023680171946334084,
            			 -0.059434418646456898,
            			 0.076488599078306393,
            			 0.41700518442169254,
            			 -0.81272363544554227,
            			 0.38611006682116222,
            			 0.067372554721963018,
            			 -0.041464936781759151,
            			 -0.016387336463522112
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "coif3"){
            	double lp1_a[] = {-3.4599772836212559e-05,
            			-7.0983303138141252e-05,
            			0.00046621696011288631,
            			0.0011175187708906016,
            			-0.0025745176887502236,
            			-0.0090079761366615805,
            			0.015880544863615904,
            			0.034555027573061628,
            			-0.082301927106885983,
            			-0.071799821619312018,
            			0.42848347637761874,
            			0.79377722262562056,
            			0.4051769024096169,
            			-0.061123390002672869,
            			-0.0657719112818555,
            			0.023452696141836267,
            			0.0077825964273254182,
            			-0.0037935128644910141
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.0037935128644910141,
            			0.0077825964273254182,
            			-0.023452696141836267,
            			-0.0657719112818555,
            			0.061123390002672869,
            			0.4051769024096169,
            			-0.79377722262562056,
            			0.42848347637761874,
            			0.071799821619312018,
            			-0.082301927106885983,
            			-0.034555027573061628,
            			0.015880544863615904,
            			0.0090079761366615805,
            			-0.0025745176887502236,
            			-0.0011175187708906016,
            			0.00046621696011288631,
            			7.0983303138141252e-05,
            			-3.4599772836212559e-05
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {-0.0037935128644910141,
            			0.0077825964273254182,
            			0.023452696141836267,
            			-0.0657719112818555,
            			-0.061123390002672869,
            			0.4051769024096169,
            			0.79377722262562056,
            			0.42848347637761874,
            			-0.071799821619312018,
            			-0.082301927106885983,
            			0.034555027573061628,
            			0.015880544863615904,
            			-0.0090079761366615805,
            			-0.0025745176887502236,
            			0.0011175187708906016,
            			0.00046621696011288631,
            			-7.0983303138141252e-05,
            			-3.4599772836212559e-05
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-3.4599772836212559e-05,
            			 7.0983303138141252e-05,
            			 0.00046621696011288631,
            			 -0.0011175187708906016,
            			 -0.0025745176887502236,
            			 0.0090079761366615805,
            			 0.015880544863615904,
            			 -0.034555027573061628,
            			 -0.082301927106885983,
            			 0.071799821619312018,
            			 0.42848347637761874,
            			 -0.79377722262562056,
            			 0.4051769024096169,
            			 0.061123390002672869,
            			 -0.0657719112818555,
            			 -0.023452696141836267,
            			 0.0077825964273254182,
            			 0.0037935128644910141
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "coif4"){
            	double lp1_a[] = {-1.7849850030882614e-06,
            			-3.2596802368833675e-06,
            			3.1229875865345646e-05,
            			6.2339034461007128e-05,
            			-0.00025997455248771324,
            			-0.00058902075624433831,
            			0.0012665619292989445,
            			0.0037514361572784571,
            			-0.0056582866866107199,
            			-0.015211731527946259,
            			0.025082261844864097,
            			0.039334427123337491,
            			-0.096220442033987982,
            			-0.066627474263425038,
            			0.4343860564914685,
            			0.78223893092049901,
            			0.41530840703043026,
            			-0.056077313316754807,
            			-0.081266699680878754,
            			0.026682300156053072,
            			0.016068943964776348,
            			-0.0073461663276420935,
            			-0.0016294920126017326,
            			0.00089231366858231456
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {-0.00089231366858231456,
            			-0.0016294920126017326,
            			0.0073461663276420935,
            			0.016068943964776348,
            			-0.026682300156053072,
            			-0.081266699680878754,
            			0.056077313316754807,
            			0.41530840703043026,
            			-0.78223893092049901,
            			0.4343860564914685,
            			0.066627474263425038,
            			-0.096220442033987982,
            			-0.039334427123337491,
            			0.025082261844864097,
            			0.015211731527946259,
            			-0.0056582866866107199,
            			-0.0037514361572784571,
            			0.0012665619292989445,
            			0.00058902075624433831,
            			-0.00025997455248771324,
            			-6.2339034461007128e-05,
            			3.1229875865345646e-05,
            			3.2596802368833675e-06,
            			-1.7849850030882614e-06
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {0.00089231366858231456,
            			-0.0016294920126017326,
            			-0.0073461663276420935,
            			0.016068943964776348,
            			0.026682300156053072,
            			-0.081266699680878754,
            			-0.056077313316754807,
            			0.41530840703043026,
            			0.78223893092049901,
            			0.4343860564914685,
            			-0.066627474263425038,
            			-0.096220442033987982,
            			0.039334427123337491,
            			0.025082261844864097,
            			-0.015211731527946259,
            			-0.0056582866866107199,
            			0.0037514361572784571,
            			0.0012665619292989445,
            			-0.00058902075624433831,
            			-0.00025997455248771324,
            			6.2339034461007128e-05,
            			3.1229875865345646e-05,
            			-3.2596802368833675e-06,
            			-1.7849850030882614e-06
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-1.7849850030882614e-06,
            			 3.2596802368833675e-06,
            			 3.1229875865345646e-05,
            			 -6.2339034461007128e-05,
            			 -0.00025997455248771324,
            			 0.00058902075624433831,
            			 0.0012665619292989445,
            			 -0.0037514361572784571,
            			 -0.0056582866866107199,
            			 0.015211731527946259,
            			 0.025082261844864097,
            			 -0.039334427123337491,
            			 -0.096220442033987982,
            			 0.066627474263425038,
            			 0.4343860564914685,
            			 -0.78223893092049901,
            			 0.41530840703043026,
            			 0.056077313316754807,
            			 -0.081266699680878754,
            			 -0.026682300156053072,
            			 0.016068943964776348,
            			 0.0073461663276420935,
            			 -0.0016294920126017326,
            			 -0.00089231366858231456
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else if ( name == "coif5"){
            	double lp1_a[] = {-9.517657273819165e-08,
            			-1.6744288576823017e-07,
            			2.0637618513646814e-06,
            			3.7346551751414047e-06,
            			-2.1315026809955787e-05,
            			-4.1340432272512511e-05,
            			0.00014054114970203437,
            			0.00030225958181306315,
            			-0.00063813134304511142,
            			-0.0016628637020130838,
            			0.0024333732126576722,
            			0.0067641854480530832,
            			-0.0091642311624818458,
            			-0.019761778942572639,
            			0.032683574267111833,
            			0.041289208750181702,
            			-0.10557420870333893,
            			-0.062035963962903569,
            			0.43799162617183712,
            			0.77428960365295618,
            			0.42156620669085149,
            			-0.052043163176243773,
            			-0.091920010559696244,
            			0.02816802897093635,
            			0.023408156785839195,
            			-0.010131117519849788,
            			-0.004159358781386048,
            			0.0021782363581090178,
            			0.00035858968789573785,
            			-0.00021208083980379827
};
            	lp1.assign (lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));

            	double hp1_a[] = {0.00021208083980379827,
            			0.00035858968789573785,
            			-0.0021782363581090178,
            			-0.004159358781386048,
            			0.010131117519849788,
            			0.023408156785839195,
            			-0.02816802897093635,
            			-0.091920010559696244,
            			0.052043163176243773,
            			0.42156620669085149,
            			-0.77428960365295618,
            			0.43799162617183712,
            			0.062035963962903569,
            			-0.10557420870333893,
            			-0.041289208750181702,
            			0.032683574267111833,
            			0.019761778942572639,
            			-0.0091642311624818458,
            			-0.0067641854480530832,
            			0.0024333732126576722,
            			0.0016628637020130838,
            			-0.00063813134304511142,
            			-0.00030225958181306315,
            			0.00014054114970203437,
            			4.1340432272512511e-05,
            			-2.1315026809955787e-05,
            			-3.7346551751414047e-06,
            			2.0637618513646814e-06,
            			1.6744288576823017e-07,
            			-9.517657273819165e-08
};
            	 hp1.assign (hp1_a,hp1_a + sizeof(hp1_a)/sizeof(double));

            	double lp2_a[] = {-0.00021208083980379827,
            			0.00035858968789573785,
            			0.0021782363581090178,
            			-0.004159358781386048,
            			-0.010131117519849788,
            			0.023408156785839195,
            			0.02816802897093635,
            			-0.091920010559696244,
            			-0.052043163176243773,
            			0.42156620669085149,
            			0.77428960365295618,
            			0.43799162617183712,
            			-0.062035963962903569,
            			-0.10557420870333893,
            			0.041289208750181702,
            			0.032683574267111833,
            			-0.019761778942572639,
            			-0.0091642311624818458,
            			0.0067641854480530832,
            			0.0024333732126576722,
            			-0.0016628637020130838,
            			-0.00063813134304511142,
            			0.00030225958181306315,
            			0.00014054114970203437,
            			-4.1340432272512511e-05,
            			-2.1315026809955787e-05,
            			3.7346551751414047e-06,
            			2.0637618513646814e-06,
            			-1.6744288576823017e-07,
            			-9.517657273819165e-08
};
            	 lp2.assign (lp2_a,lp2_a + sizeof(lp2_a)/sizeof(double));

            	 double hp2_a[] = {-9.517657273819165e-08,
            			 1.6744288576823017e-07,
            			 2.0637618513646814e-06,
            			 -3.7346551751414047e-06,
            			 -2.1315026809955787e-05,
            			 4.1340432272512511e-05,
            			 0.00014054114970203437,
            			 -0.00030225958181306315,
            			 -0.00063813134304511142,
            			 0.0016628637020130838,
            			 0.0024333732126576722,
            			 -0.0067641854480530832,
            			 -0.0091642311624818458,
            			 0.019761778942572639,
            			 0.032683574267111833,
            			 -0.041289208750181702,
            			 -0.10557420870333893,
            			 0.062035963962903569,
            			 0.43799162617183712,
            			 -0.77428960365295618,
            			 0.42156620669085149,
            			 0.052043163176243773,
            			 -0.091920010559696244,
            			 -0.02816802897093635,
            			 0.023408156785839195,
            			 0.010131117519849788,
            			 -0.004159358781386048,
            			 -0.0021782363581090178,
            			 0.00035858968789573785,
            			 0.00021208083980379827
};
            	 hp2.assign (hp2_a,hp2_a + sizeof(hp2_a)/sizeof(double));
            	 return 0;
            }
    else {
    	cout << "Filter Not in Database" << endl;
    	return -1;
    }

}


