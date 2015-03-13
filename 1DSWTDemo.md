# 1D Stationary Wavelet Transform Demo #

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/swtdemo.cpp)

```
 //============================================================================
// Name        : swtdemo.cpp
// Author      : Rafat Hussain
// Version     :
// Copyright   : 
// Description : 1D Stationary Wavelet Transform Demo
//============================================================================

#include <iostream>
#include <fstream>
#include "wavelet.h"
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

        vector<double> swt_output;

        // perform J-Level DWT
        int length;// All coefficients are of same length. Variable "length" returns length
        // of coefficients. It is not required for ISWT computations.

        // If signal is not of dyadic length it is zeropadded as the algorithm is designed to
        // work with signals of 2^M length. These extra zeros can be removed after reconstruction.


            swt(sig, J, nm, swt_output, length);
            vector<double> iswt_output;
            iswt(swt_output,J, nm,iswt_output);
            gnuswtplot(J);

            return 0;
   }
```

The signal is the same Piecewise Regular signal that is used in the [1D Decimated DWT Demo](http://code.google.com/p/wavelet1d/wiki/1DDWTdemo).

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png)

> The filter used is Daubechies db2 Orthogonal length=4 filter. The undecimated detail and approximation coefficients are plotted as following using GNUPLOT.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TchqBRDzvzI/AAAAAAAAAFQ/OUMcPMT554s/swtdisp.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TchqBRDzvzI/AAAAAAAAAFQ/OUMcPMT554s/swtdisp.png)

> There is a very obvious trade-off here. The Transform is good for signal analysis as details and approximation coefficients are of the same length as the original signal at every decomposition stage but the system is highly redundant.

> The Inverse SWT is perfect reconstruction.

