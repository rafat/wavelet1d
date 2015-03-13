# 1D DWT/IDWT Demo using dwt function #

_Notes_

1. dwt and dwt\_sym functions are interchangeable in the following code as they take the same arguments.

2. Please make sure that right header file is included. Following sample program links to wavelet2d.h which is the header for shared libraries. wavelet2s.h must be included for static libraries.

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/wavelet2d/demo/wavedemo1.cpp)

```

//============================================================================
// Name        : wavedemo1.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   :
// Description : 1D DWT Demo
//============================================================================

#include <iostream>
#include <fstream>
#include "wavelet2d.h"
#include <vector>
#include <string>
#include <cmath>
using namespace std;

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

        // perform J-Level DWT
		vector<int> length;

	    dwt_sym(sig, J, nm, dwt_output,flag,length);
            ofstream dwtout("dwtout.txt");
            for (unsigned int i = 0; i < dwt_output.size(); i++){
	    	 dwtout << dwt_output[i] << endl;

	    }


	    //Perform J-Level IDWT
	    vector<double> output;
	    idwt_sym(dwt_output, flag,nm,output,length);

	    ofstream sig1("recon.txt");
	    ofstream diff("diff.txt");

	    cout <<" Recon signal size" << output.size() << endl;
	    for (unsigned int i = 0; i < output.size(); i++){
	    	sig1 << output[i] << endl;
	    	diff << output[i] - original[i] << endl;

	    }
	    //gnudwtplot(J);
	    return 0;
}

```

This sample program asks user for wavelet name, signal file and number of DWT stages. It then computes DWT of the input signal and then the IDWT from the wavelet coefficients so obtained.
![https://lh3.googleusercontent.com/-wIrGfF_x4ug/Tk9_idI6YiI/AAAAAAAAAL0/jIFwli0Tndk/s912/1ddwtscreen.png](https://lh3.googleusercontent.com/-wIrGfF_x4ug/Tk9_idI6YiI/AAAAAAAAAL0/jIFwli0Tndk/s912/1ddwtscreen.png)

As can be seen, outputs are in _.txt_ format but they can be modified and plotted with any software of your choice.

![https://lh5.googleusercontent.com/-wAxR0xM5I38/Tk9_iWhol_I/AAAAAAAAAL4/C9Sck_VZ7Ms/s912/1ddwtsig.png](https://lh5.googleusercontent.com/-wAxR0xM5I38/Tk9_iWhol_I/AAAAAAAAAL4/C9Sck_VZ7Ms/s912/1ddwtsig.png)
_Length 256 Piecewise Regular signal from WaveLab_

![https://lh4.googleusercontent.com/-xG9LotALuPs/Tk9_iFkYsiI/AAAAAAAAALw/G0S0ipsAOKI/s912/1ddwtout.png](https://lh4.googleusercontent.com/-xG9LotALuPs/Tk9_iFkYsiI/AAAAAAAAALw/G0S0ipsAOKI/s912/1ddwtout.png)
_2 Level DWT Decomposition Coefficients obtained using db2 wavelet_

![https://lh6.googleusercontent.com/-DscYhVPOEPY/Tk9_iFvaZTI/AAAAAAAAALs/3q4hfboBLLA/s912/1ddwtrecon.png](https://lh6.googleusercontent.com/-DscYhVPOEPY/Tk9_iFvaZTI/AAAAAAAAALs/3q4hfboBLLA/s912/1ddwtrecon.png)
_Reconstructed Signal After IDWT computations_