/*
 * wave1d.h
 *
 *  Created on: Feb 26, 2011
 *      Author: Rafat Hussain
 */

#ifndef WAVE1D_H_
#define WAVE1D_H_
#include <vector>
using namespace std;

void* dwt(vector<double> &, int ,string , vector<double> &,  vector<double> &);
void* dwt1(string, vector<double> &, vector<double> &, vector<double> &);
void* dyadic_zpad_1d(vector<double> &);
double convol(vector<double> &, vector<double> &, vector<double> &);
int filtcoef(string , vector<double> &, vector<double> &, vector<double> &,
		vector<double> &);
void downsamp(vector<double> &, int , vector<double> &);
void upsamp(vector<double> &, int, vector<double> &);
void circshift(vector<double> &, int );

int sign(int);
void* gnudwtplot(int);
void* idwt(vector<double> &,vector<double> &, string , vector<double> &) ;

void* idwt1(string wname, vector<double> &, vector<double> &, vector<double> &);

int vecsum(vector<double> &, vector<double> &, vector<double> &);


#endif /* WAVE1D_H_ */
