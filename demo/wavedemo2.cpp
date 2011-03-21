//============================================================================
// Name        : wavedemo2.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : Wavelet Demo comparing linear and non-linear approximation properties
//             : of a given wavelet.
//============================================================================

#include <iostream>
#include <fstream>
#include "wave1d.h"
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

        // perform J-Level DWT
	    dwt(sig, J, nm, dwt_output,flag );

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
		    idwt(dwt_output, flag,nm,output);

	        unsigned int count = output.size();
	    	ofstream gnulinappx("gnulinappx.dat");
	    	for (unsigned int i = 0;i < count; i++) {
	    		gnulinappx << i << " " << output[i] << endl;
	    	}
	    	gnulinappx.close();

	    	// Performing Non Linear Approximation by using only most
	    	// significant coefficients

	        vector<double> dwt_output2, flag2;

	        // perform J-Level DWT
		    dwt(sig, J, nm, dwt_output2,flag2 );

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


		     // Finding IDWT with approximated coefficients

			    vector<double> output2;
			    idwt(dwt_output2, flag2,nm,output2);

		        unsigned int count2 = output2.size();
		        cout << count2 << endl;
		    	ofstream gnunlappx("gnunlappx.dat");
		    	for (unsigned int i = 0;i < count2; i++) {
		    		gnunlappx << i << " " << output2[i] << endl;
		    	}
		    	gnunlappx.close();


	return 0;
}
