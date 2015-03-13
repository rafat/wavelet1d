# 1D DWT/IDWT Demo using dwt function #

[A sample code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/wavedemo1.cpp) for DWT/IDWT for a given signal is given below. The signal used is "Piecewise-Regular" function from Stanford University's Matlab Wavelab Toolbox.

```
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

        vector<double> dwt_output, flag;

        // perform J-Level DWT
	    dwt(sig, J, nm, dwt_output,flag );


	    //Perform J-Level IDWT
	    vector<double> output;
	    idwt(dwt_output, flag,nm,output);

	    ofstream sig1("recon.txt");
	    ofstream diff("diff.txt");

	    cout <<" Recon signal size" << output.size() << endl;
	    for (unsigned int i = 0; i < output.size(); i++){
	    	sig1 << output[i] << endl;
	    	diff << output[i] - original[i] << endl;

	    }
	    gnudwtplot(J);
	    return 0;
   }

```

The program asks user for wavelet name, signal file and number of DWT stages. It then computes DWT of the input signal and then the IDWT from the wavelet coefficients so obtained. Running the program from Eclipse console :

```
 ********J- LEVEL DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********
This program accepts signal from the user in a file format 
and performs Discrete Wavelet Transform with specified   
wavelet. 
                                                             
 The Following Wavelets are in the Database:                 
 haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10,  
 db11, db12, db13, db14, db15.                               
 bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,bior2.6,bior2.8, 
 bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4,
 bior5.5, bior6.8.
 coif1, coif2, coif3, coif4, coif5.
Please Enter the Wavelet Name at the Prompt( No quotes)     :
db3
Enter the name of signal file at the Prompt eg., signal.txt :
signal.txt
Please Enter the Number of DWT Stages J             :
3
```

gnudwtplot(J) is called which generates GNUPLOT script. GNUPLOT script is called from gnuplot as following.

gnuplot> cd '<Path to File Outputs>'

gnuplot> load "gnudwt.gnu.txt"

1. Input Signal : Piecewise-Regular signal from Stanford University Wavelab Toolbox.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqNlUktbI/AAAAAAAAAC8/YcwG3JqEw-s/wavedemo1-1.png)

2. Approximation and Detail coefficients at each iteration stage. Note- While the program is plotting approximation coefficients, the actual DWT doesn't store approximation values at each level for redundancy reasons. Only approximation coefficients at the Jth stage are stored.

![https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqqNoZC6fI/AAAAAAAAAC0/zWiH8iSn-2k/s912/wavedemo1-2.png](https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqqNoZC6fI/AAAAAAAAAC0/zWiH8iSn-2k/s912/wavedemo1-2.png)

3. Four filters associated with orthogonal db3 wavelet.

![https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqqOIJg76I/AAAAAAAAADA/EfV7aPfz6GY/wavedemo1-3.png](https://lh3.googleusercontent.com/_TwtGvT0Ma-0/TbqqOIJg76I/AAAAAAAAADA/EfV7aPfz6GY/wavedemo1-3.png)

4. Reconstructed Signal.

![https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqOaERebI/AAAAAAAAADI/EfxtbRPdYq4/wavedemo1-4.png](https://lh5.googleusercontent.com/_TwtGvT0Ma-0/TbqqOaERebI/AAAAAAAAADI/EfxtbRPdYq4/wavedemo1-4.png)


# 1D DWT/IDWT Demo using dwt\_sym function #

[Sample Code](http://code.google.com/p/wavelet1d/source/browse/trunk/demo/wavedemo_sym1.cpp) for DWT/IDWT using dwt\_sym function is also available.