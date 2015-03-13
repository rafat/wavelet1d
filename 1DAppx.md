# 1D Linear vs Nonlinear Approximation using dwt #

One of the important property of wavelets is that it gives better approximation results when N largest coefficients are considered instead of first N coefficients.

[A sample code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/wavedemo2.cpp) to compare the performance of linear vs. Non-linear approximations using dwt function is given below. The signal used is "Piecewise-Regular" function from Stanford University's Matlab Wavelab Toolbox.
```
 #include <iostream>
#include <fstream>
#include "wavelet.h"
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

```

This program computes DWT of the 256-point signal and then uses only 100 coefficients to reconstruct the signal. In the first case, first 100 coefficients are are used to compute linear approximation. In the second case, 100 largest coefficents are used by calculating a threshold that separates largest 100 coefficients from the other 156.The wavelet used is db3 and the signal is decomposed and reconstructd over 4 stages of DWT/IDWT. The results are plotted below.

Input Signal

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png)

Reconstructed Signal using only first 100 coefficients.(Linear Approximation)

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbssWS6M0AI/AAAAAAAAAEg/itai-HfR4Pc/linappx.gif](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbssWS6M0AI/AAAAAAAAAEg/itai-HfR4Pc/linappx.gif)

Reconstructed Signal using largest 100 coefficients.(Non-Linear Approximation)

![https://lh4.googleusercontent.com/_TwtGvT0Ma-0/TbssbItgEoI/AAAAAAAAAEk/GkOTJ5cyWUA/nlappx.gif](https://lh4.googleusercontent.com/_TwtGvT0Ma-0/TbssbItgEoI/AAAAAAAAAEk/GkOTJ5cyWUA/nlappx.gif)


# 1D Linear vs Nonlinear Approximation using dwt\_sym #

_Note- Symmetric Extension DWT Implementation dwt\_sym is not a N input N output system so it has more DWT coefficients which usually means you'll need more coefficients to match dwt approximation performance for same wavelet. the trade-off is better performance at signal boundaries. In this example same number of approximation coefficients (100) are considered._

We use the same signal and wavelet family as above.[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/wavedemo_sym2.cpp) is available at the source page.

```

 #include <iostream>
#include <fstream>
#include "wavelet.h"
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
        int ext = 2;

        // perform J-Level DWT
            dwt_sym(sig, J, nm, dwt_output,flag, length1,ext );

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
                    idwt_sym(dwt_output, flag,nm,output, length1);

                unsigned int count = output.size();
                ofstream gnulinappx("gnulinappx.dat");
                for (unsigned int i = 0;i < count; i++) {
                        gnulinappx << i << " " << output[i] << endl;
                }
                gnulinappx.close();

                // Performing Non Linear Approximation by using only most
                // significant coefficients

                vector<double> dwt_output2, flag2;
                vector<int> length2;

                // perform J-Level DWT
                    dwt_sym(sig, J, nm, dwt_output2,flag2,length2, ext );

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
                            idwt_sym(dwt_output2, flag2,nm,output2, length2);

                        unsigned int count2 = output2.size();
                        cout << count2 << endl;
                        ofstream gnunlappx("gnunlappx.dat");
                        for (unsigned int i = 0;i < count2; i++) {
                                gnunlappx << i << " " << output2[i] << endl;
                        }
                        gnunlappx.close();


        return 0;
   }
```

Reconstructed Signal using only first 100 coefficients.(Linear Approximation)

![https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbsre73jQxI/AAAAAAAAAEY/ffRAs9icfhM/symlinappx.gif](https://lh6.googleusercontent.com/_TwtGvT0Ma-0/Tbsre73jQxI/AAAAAAAAAEY/ffRAs9icfhM/symlinappx.gif)

Reconstructed Signal using largest 100 coefficients.(Non-Linear Approximation)

![https://lh4.googleusercontent.com/_TwtGvT0Ma-0/TbsrfWkT2sI/AAAAAAAAAEc/4vf49wiYHtM/symnlappx.gif](https://lh4.googleusercontent.com/_TwtGvT0Ma-0/TbsrfWkT2sI/AAAAAAAAAEc/4vf49wiYHtM/symnlappx.gif)