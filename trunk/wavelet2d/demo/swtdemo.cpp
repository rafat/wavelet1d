//============================================================================
// Name        : swtdemo.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : 1D Stationary Wavelet Transform Demo
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
	    original = sig; // Make a copy of the signal if you want to use original signal
       // later on. The other option is to use IDWT output as the SWt/ISWT system is
       // Perefect Reconstruction system.
	    cout << "Please Enter the Number of DWT Stages J             :" << endl;

	    int J;
	    cin >> J ;

        vector<double> swt_output;

        // perform J-Level DWT
        int length;// All coefficients are of same length. Variable "length" returns length
        // of coefficients. It is not required for ISWT computations.



	    swt(sig, J, nm, swt_output, length);

            ofstream swtcoeff("swtcoeff.txt");
            for (unsigned int i=0; i < swt_output.size(); i++) {
              swtcoeff << swt_output[i] << endl;
             } 
	    vector<double> iswt_output;
	    iswt(swt_output,J, nm,iswt_output);
            ofstream sig1("recon.txt");
            ofstream diff("diff.txt");

            cout <<" Recon signal size" << iswt_output.size() << endl;
            for (unsigned int i = 0; i < iswt_output.size(); i++){
                sig1 << iswt_output[i] << endl;
                diff << iswt_output[i] - original[i] << endl;

            }


	    return 0;
}
