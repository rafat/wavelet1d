# 1D Linear vs Nonlinear Approximation #

_Note_

1. dwt and dwt\_sym functions are interchangeable in the following code as they take the same arguments. You may have to make some trivial modifications due to different sized output coeffficients.

2. Please make sure that right header file is included. Following sample program links to wavelet2s.h which is the header for static libraries. wavelet2d.h must be included for shared libraries.

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/wavelet2d/demo/wavedemo2.cpp)

```
//============================================================================
// Name        : wavedemo2.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : Wavelet Demo comparing linear and non-linear approximation properties
//             : of a given wavelet. Implemented using dwt and idwt. Interchangeable with
//             : dwt_sym and idwt_sym.
//============================================================================

#include <iostream>
#include <fstream>
#include "wavelet2s.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

void findthresh(vector<double> vector1, int N, double& t){
	sort(vector1.begin(), vector1.end(), greater<double>());
	t = vector1.at(N-1);
}

int main() {
	cout << "********J- LEVEL DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********" << endl; // prints
	    cout << "This program accepts signal from the user in a file format " << endl;
	    cout << "and performs Discrete Wavelet Transform with specified   " << endl;
	    cout << "wavelet. "                                               << endl;
	    cout << "                                                             " << endl;
	    cout << " The Following Wavelets are in the Database:                 " << endl;
	    cout << " haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10,  "   << endl;
	    cout << " db11, db12, db13, db14, db15.                               "  << endl;
	    cout << " bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,bior2.6,bior2.8, " << endl;
	    cout << " bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4,"        << endl;
	    cout << " bior5.5, bior6.8."                                            << endl;
	    cout << " coif1, coif2, coif3, coif4, coif5."                           << endl;
	    cout << "Please Enter the Wavelet Name at the Prompt( No quotes)     :" << endl;

	    string nm; // nm will store the name of Wavelet Family
	    cin >> nm;
	    cout << "Enter the name of signal file at the Prompt eg., signal.txt :" << endl;
	    char inp[50];
	    cin >> inp;
	    vector<double> sig;
	    ifstream sig_inp(inp);
        if ( !sig_inp.good()){
        	cout << "The File doesn't exist"<< endl;
        	exit(1);
        }
	    while (sig_inp) {
	    	double temp;
	    	sig_inp >> temp;
	    	sig.push_back(temp);
	    }
	    sig.pop_back();
	    vector<double> original;
	    original = sig;
	    cout << "Please Enter the Number of DWT Stages J             :" << endl;

	    int J;
	    cin >> J ;

        vector<double> dwt_output, flag;
        vector<int> length1;

        // perform J-Level DWT
	    dwt(sig, J, nm, dwt_output,flag, length1);

	    // Performing Linear Approximation by using only first 100 coefficients
	    // Coefficients in dwt_output are stored as following
	    // dwt_output =[ Appx(J-1) Detail(J-1) Detail(J-2) .... Detail(0)]

	     int n_coef = 100; // Number of significant coefficients

	     int n_non_sig= dwt_output.size() - n_coef; // Number of Coefficients that will
	     // be set to zero
	     dwt_output.erase(dwt_output.end()- n_non_sig,dwt_output.end());
	     // Deleting last n_non_sig coefficients and replacing them with zeros
	     dwt_output.insert(dwt_output.end(),n_non_sig,0);

	     ofstream linearsig("linsig.txt");
	     for (unsigned int i = 0; i < dwt_output.size(); i++) {
	    	 linearsig << dwt_output[i] << endl;

	     }

	     // Finding IDWT with approximated coefficients

		    vector<double> output;
		    idwt(dwt_output, flag,nm,output, length1);

	        unsigned int count = output.size();
	    	ofstream gnulinappx("gnulinappx.dat");
	    	for (unsigned int i = 0;i < count; i++) {
	    		gnulinappx << output[i] << endl;
	    	}
	    	gnulinappx.close();

	    	// Performing Non Linear Approximation by using only most
	    	// significant coefficients

	        vector<double> dwt_output2, flag2;
	        vector<int> length2;

	        // perform J-Level DWT
		    dwt(sig, J, nm, dwt_output2,flag2,length2);

		    double thresh = 0.0;

		    vector<double> temp_dwtoutput;
		    for (unsigned int i =0; i < dwt_output2.size();i++){
		    	double temp = abs(dwt_output2[i]);
		    	temp_dwtoutput.push_back(temp);
		    }
		    /*
		    for (unsigned int i =0; i < temp_dwtoutput.size(); i++){
		    	cout << temp_dwtoutput[i] << endl;
		    }
		    */

		    findthresh(temp_dwtoutput,n_coef, thresh);

		    for (unsigned int i = 0; i < dwt_output2.size();i++){
		    	double temp = abs(dwt_output2[i]);
		    	if (temp < thresh){
		    		dwt_output2.at(i) = 0.0;

		    	}
		    }
		    /*
		    for (unsigned int i =0; i < dwt_output2.size(); i++){
		    	cout << dwt_output2[i] << endl;
		    }
		    */
		    
		         ofstream nonlinsig("nonlinsig.txt");
	     for (unsigned int i = 0; i < dwt_output2.size(); i++) {
	    	 nonlinsig << dwt_output2[i] << endl;

	     }


		     // Finding IDWT with approximated coefficients

			    vector<double> output2;
			    idwt(dwt_output2, flag2,nm,output2, length2);

		        unsigned int count2 = output2.size();
		        cout << count2 << endl;
		    	ofstream gnunlappx("gnunlappx.dat");
		    	for (unsigned int i = 0;i < count2; i++) {
		    		gnunlappx << output2[i] << endl;
		    	}
		    	gnunlappx.close();


	return 0;
}

```

This program computes DWT of the 256-point signal and then uses only 100 coefficients to reconstruct the signal. In the first case, first 100 coefficients are are used to compute linear approximation. In the second case, 100 largest coefficients are used by calculating a threshold that separates largest 100 coefficients from the others.The wavelet used is db3 and the signal is decomposed and reconstructed over 5 stages of DWT/IDWT. The results are plotted below.

![https://lh3.googleusercontent.com/-786mdxcFiFw/Tk-HwRh25GI/AAAAAAAAAMM/rE5xc6ngnv0/s800/1dappxscreen.png](https://lh3.googleusercontent.com/-786mdxcFiFw/Tk-HwRh25GI/AAAAAAAAAMM/rE5xc6ngnv0/s800/1dappxscreen.png)
_Console showing db3 5 level computation_

![https://lh5.googleusercontent.com/-wAxR0xM5I38/Tk9_iWhol_I/AAAAAAAAAL4/C9Sck_VZ7Ms/s912/1ddwtsig.png](https://lh5.googleusercontent.com/-wAxR0xM5I38/Tk9_iWhol_I/AAAAAAAAAL4/C9Sck_VZ7Ms/s912/1ddwtsig.png)
_Length 256 Piecewise Regular signal from WaveLab_

In Linear Approximation case, only first 100 coefficients are retained. Others are set to zero.

![https://lh3.googleusercontent.com/-WIIVJSzPjJI/Tk-HwGSTXxI/AAAAAAAAAMI/kzHzVRuIqeY/s912/linsigcoeff.png](https://lh3.googleusercontent.com/-WIIVJSzPjJI/Tk-HwGSTXxI/AAAAAAAAAMI/kzHzVRuIqeY/s912/linsigcoeff.png)
_Retained Coefficients for Linear Approximation case_

![https://lh5.googleusercontent.com/-vo62gHChC0U/Tk-HvyiCgyI/AAAAAAAAAMA/NrpX9FBIOXU/s912/1dlinappx.png](https://lh5.googleusercontent.com/-vo62gHChC0U/Tk-HvyiCgyI/AAAAAAAAAMA/NrpX9FBIOXU/s912/1dlinappx.png)
_Reconstructed Signal using Linear Approximation_

In Non-Linear case, largest 100 coefficients are retained.

![https://lh3.googleusercontent.com/-2hvwAm8QTyU/Tk-HvuIZwjI/AAAAAAAAAL8/LxS0wb_IFOk/s912/nonlinsigcoeff.png](https://lh3.googleusercontent.com/-2hvwAm8QTyU/Tk-HvuIZwjI/AAAAAAAAAL8/LxS0wb_IFOk/s912/nonlinsigcoeff.png)
_Retained Coefficients(Non-Linear Approximation case_

![https://lh4.googleusercontent.com/-OCLIjzcaO7U/Tk-Hv2ClRaI/AAAAAAAAAME/C-0CzXX8moI/s912/1dnlappx.png](https://lh4.googleusercontent.com/-OCLIjzcaO7U/Tk-Hv2ClRaI/AAAAAAAAAME/C-0CzXX8moI/s912/1dnlappx.png)
_Reconstructed Signal(Non-Linear Approximation Case)_