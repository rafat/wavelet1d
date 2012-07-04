void* iswt_2d(vector<double> &swtop,int J, string nm, vector<vector<double> > &sig,int row, int col) {

// 2D Inverse Stationary Wavelet Transform 
// Works In Conjunction with C++ Wavelet Library http://code.google.com/p/wavelet1d/

// (C) Rafat Hussain under GNU GPL v3 
// See Project page http://code.google.com/p/wavelet1d/ 

// swtop [1D Output vector from swt_2d]
// J [Number of Levels]
// nm [Wavelet Name]
// sig [Output]
// row & col [Rows and cols of sig]

// I have not fully tested this code. 

// References :
 
// 1. Nason, G.P. and Silverman, B.W. (1995) The stationary wavelet transform and some statistical applications.
// In Antoniadis, A. and Oppenheim, G. (eds) Lecture Notes in Statistics, 103, 281--300. 

// 2. Pesquet, J.C.; H. Krim, H. Carfatan (1996), "Time-invariant orthonormal wavelet representations," 
// http://pdf.aminer.org/000/318/848/frequency_shift_invariant_orthonormal_wavelet_packet_representations.pdf

	int rows_n =row;
	int cols_n = col;
	vector<double> lp1,hp1,lp2,hp2;
	filtcoef(nm,lp1,hp1,lp2,hp2);
	vector<double> low_pass = lp2;
	vector<double> high_pass = hp2;
	int lf = low_pass.size();
    int sum =0;

    vector<vector<double> >  cLL(rows_n, vector<double>(cols_n));

	for (int iter=J; iter > 0; iter--) {

    	vector<vector<double> >  cLH(rows_n, vector<double>(cols_n));
    	vector<vector<double> >  cHL(rows_n, vector<double>(cols_n));
    	vector<vector<double> >  cHH(rows_n, vector<double>(cols_n));

        if (iter == J) {
            for (int i = 0 ; i < rows_n; i++ ){
                for (int j = 0; j < cols_n; j++){

                    cLL[i][j] = swtop[sum+ i * cols_n + j];

                    cHL[i][j] = swtop[sum+ rows_n * cols_n+ i * cols_n + j];

                    cLH[i][j] = swtop[sum+ 2 * rows_n * cols_n + i * cols_n + j];

                    cHH[i][j] = swtop[sum+ 3* rows_n * cols_n + i * cols_n + j];
                }
            }
            sum+= 4 * rows_n * cols_n;


        } else {
            for (int i = 0 ; i < rows_n; i++ ){
                for (int j = 0; j < cols_n; j++){

                    cHL[i][j] = swtop[sum+  i * cols_n + j];

                    cLH[i][j] = swtop[sum+ rows_n * cols_n + i * cols_n + j];

                    cHH[i][j] = swtop[sum+ 2* rows_n * cols_n + i * cols_n + j];
                    }
                }
            sum+= 3 * rows_n * cols_n;

        }


		int value = (int) pow(2.0,(double) (iter-1));

		for (int count = 0; count < value; count++) {

			vector<vector<double> >  cLL1(rows_n/value, vector<double>(cols_n/value));
			vector<vector<double> >  cLH1(rows_n/value, vector<double>(cols_n/value));
			vector<vector<double> >  cHL1(rows_n/value, vector<double>(cols_n/value));
			vector<vector<double> >  cHH1(rows_n/value, vector<double>(cols_n/value));

	    int row1 = 0;
		int col1 = 0;
			for (int index_r = count; index_r < rows_n; index_r+=value){
				for (int index_c = count; index_c < cols_n; index_c+=value){
					cLL1[row1][col1]=cLL[index_r][index_c];
					cLH1[row1][col1]=cLH[index_r][index_c];
					cHL1[row1][col1]=cHL[index_r][index_c];
					cHH1[row1][col1]=cHH[index_r][index_c];
					col1++;

				        }
				col1 = 0;
				row1++;
			          }
       
			//Implementing shift =0

	        unsigned int len_c = cLL1[0].size();
	        unsigned int len_r = cLL1.size();
	        vector<vector<double> > cL(rows_n,vector<double>(cols_n));
	        vector<vector<double> > cH(rows_n,vector<double>(cols_n));

			vector<vector<double> >  cLL2(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cLH2(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cHL2(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cHH2(len_r/2, vector<double>(len_c/2));

			row1=0;
			col1=0;

			for (int index_r = 0; index_r < len_r; index_r+=2){
							for (unsigned int index_c = 0; index_c < len_c; index_c+=2){
								cLL2[row1][col1]=cLL1[index_r][index_c];
								cLH2[row1][col1]=cLH1[index_r][index_c];
								cHL2[row1][col1]=cHL1[index_r][index_c];
								cHH2[row1][col1]=cHH1[index_r][index_c];
								col1++;

							        }
							col1 = 0;
							row1++;
						          }
            
			vector<vector<double> >  cLLu(len_r, vector<double>(len_c));
			vector<vector<double> >  cLHu(len_r, vector<double>(len_c));
			vector<vector<double> >  cHLu(len_r, vector<double>(len_c));
			vector<vector<double> >  cHHu(len_r, vector<double>(len_c));
           int row_u = cLL2.size();
           int col_u = cLL2[0].size();
			upsamp2(cLL2,cLLu,2,2);
			upsamp2(cLH2,cLHu,2,2);
			upsamp2(cHL2,cHLu,2,2);
			upsamp2(cHH2,cHHu,2,2);
			//cLL2.clear();cLH2.clear();cHL2.clear();cHH2.clear();



			

			vector<vector<double> >  cLL00(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cLH00(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cHL00(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cHH00(len_r+lf, vector<double>(len_c + lf));

			

			per_ext2d(cLLu,cLL00,lf/2);
			per_ext2d(cLHu,cLH00,lf/2);
			per_ext2d(cHLu,cHL00,lf/2);
			per_ext2d(cHHu,cHH00,lf/2);
			//cLLu.clear();cLHu.clear();cHLu.clear();cHHu.clear();

			//Shift By 1


			int U = 2; // Upsampling Factor

			

			vector<vector<double> >  cLL3(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cLH3(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cHL3(len_r/2, vector<double>(len_c/2));
			vector<vector<double> >  cHH3(len_r/2, vector<double>(len_c/2));
            col1 = 0;
            row1 = 0;


			for (int index_r = 1; index_r < len_r; index_r+=2){
				for (unsigned int index_c = 1; index_c < len_c; index_c+=2){
					cLL3[row1][col1]=cLL1[index_r][index_c];
					cLH3[row1][col1]=cLH1[index_r][index_c];
					cHL3[row1][col1]=cHL1[index_r][index_c];
					cHH3[row1][col1]=cHH1[index_r][index_c];
					col1++;

				}
				col1 = 0;
				row1++;
			}
			

			vector<vector<double> >  cLLu2(len_r, vector<double>(len_c));
			vector<vector<double> >  cLHu2(len_r, vector<double>(len_c));
			vector<vector<double> >  cHLu2(len_r, vector<double>(len_c));
			vector<vector<double> >  cHHu2(len_r, vector<double>(len_c));

			row_u = cLL3.size();
			col_u = cLL3[0].size();

			upsamp2(cLL3,cLLu2,2,2);
			upsamp2(cLH3,cLHu2,2,2);
			upsamp2(cHL3,cHLu2,2,2);
			upsamp2(cHH3,cHHu2,2,2);



			//cLL3.clear();cLH3.clear();cHL3.clear();cHH3.clear();

			cout << "STAGE 6" << endl;

			vector<vector<double> >  cLL01(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cLH01(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cHL01(len_r+lf, vector<double>(len_c + lf));
			vector<vector<double> >  cHH01(len_r+lf, vector<double>(len_c + lf));

			per_ext2d(cLLu2,cLL01,lf/2);
			per_ext2d(cLHu2,cLH01,lf/2);
			per_ext2d(cHLu2,cHL01,lf/2);
			per_ext2d(cHHu2,cHH01,lf/2);



			cLLu2.clear();cLHu2.clear();cHLu2.clear();cHHu2.clear();

			int rowsLL = cLL01.size();
			int colsLL = cLL01[0].size();

			vector<vector<double> >  cLLt00(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cLHt00(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cHLt00(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cHHt00(rowsLL + lf-1, vector<double>(colsLL));

			vector<vector<double> >  cLLt10(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cLHt10(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cHLt10(rowsLL + lf-1, vector<double>(colsLL));
			vector<vector<double> >  cHHt10(rowsLL + lf-1, vector<double>(colsLL));

			vector<vector<double> >  cLLt01(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cLHt01(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cHLt01(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cHHt01(rowsLL + lf-1, vector<double>(colsLL + lf-1));

			vector<vector<double> >  cLLt11(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cLHt11(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cHLt11(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cHHt11(rowsLL + lf-1, vector<double>(colsLL + lf-1));

			vector<vector<double> >  cLLt0(rowsLL + lf-1, vector<double>(colsLL + lf-1));
			vector<vector<double> >  cLLt1(rowsLL + lf-1, vector<double>(colsLL + lf-1));

			vector<vector<double> >  oupt0(len_r, vector<double>(len_c));
			vector<vector<double> >  oupt1(len_r, vector<double>(len_c));
			vector<vector<double> >  oupt(len_r, vector<double>(len_c));

			

			for (int j=0; j< colsLL;j++) {
			    vector<double> sigLL0;
        		vector<double> sigLH0;
        		vector<double> sigHL0;
        		vector<double> sigHH0;

				for (int i=0; i < rowsLL; i++) {
					double temp01 = cLL00[i][j];
					double temp02 = cLH00[i][j];
					double temp03 = cHL00[i][j];
					double temp04 = cHH00[i][j];
					sigLL0.push_back(temp01);
					sigLH0.push_back(temp02);
					sigHL0.push_back(temp03);
					sigHH0.push_back(temp04);

				}

				

				// First Convolution Step for LL,LH,HL,HH



    			vector<double> X0LL;
    			convfft(sigLL0, low_pass, X0LL);

				vector<double> X0LH;
    			convfft(sigLH0, low_pass, X0LH);

				vector<double> X0HL;
    			convfft(sigHL0, high_pass, X0HL);

				vector<double> X0HH;
    			convfft(sigHH0, high_pass, X0HH);

				int lent=X0LL.size();

				for (int i=0; i <  lent; i++) {

					cLLt00[i][j]=X0LL[i];
					cLHt00[i][j]=X0LH[i];
					cHLt00[i][j]=X0HL[i];
					cHHt00[i][j]=X0HH[i];


				}




			}

			


			for (int i=0; i< rowsLL+lf-1;i++) {
			    vector<double> sigLL0;
        		vector<double> sigLH0;
        		vector<double> sigHL0;
        		vector<double> sigHH0;

				for (int j=0; j < colsLL; j++) {
					double temp01 = cLLt00[i][j];
					double temp02 = cLHt00[i][j];
					double temp03 = cHLt00[i][j];
					double temp04 = cHHt00[i][j];
					sigLL0.push_back(temp01);
					sigLH0.push_back(temp02);
					sigHL0.push_back(temp03);
					sigHH0.push_back(temp04);

				}

				

				// Second Convolution Step for LL,LH,HL,HH

    			vector<double> X0LL;
    			convfft(sigLL0, low_pass, X0LL);

				vector<double> X0LH;
    			convfft(sigLH0, high_pass, X0LH);

				vector<double> X0HL;
    			convfft(sigHL0, low_pass, X0HL);

				vector<double> X0HH;
    			convfft(sigHH0, high_pass, X0HH);

				int lent=X0LL.size();
				

				for (int j=0; j <  lent; j++) {

					cLLt01[i][j]=X0LL[j];
					cLHt01[i][j]=X0LH[j];
					cHLt01[i][j]=X0HL[j];
					cHHt01[i][j]=X0HH[j];
					cLLt0[i][j]=X0LL[j]+X0LH[j]+X0HL[j]+X0HH[j];

				}





			}






			for (int j=0; j< colsLL;j++) {
			    vector<double> sigLL0;
        		vector<double> sigLH0;
        		vector<double> sigHL0;
        		vector<double> sigHH0;

				for (int i=0; i < rowsLL; i++) {
					double temp01 = cLL01[i][j];
					double temp02 = cLH01[i][j];
					double temp03 = cHL01[i][j];
					double temp04 = cHH01[i][j];
					sigLL0.push_back(temp01);
					sigLH0.push_back(temp02);
					sigHL0.push_back(temp03);
					sigHH0.push_back(temp04);

				}

				

				// First Convolution Step for LL,LH,HL,HH

    			vector<double> X0LL;
    			convfft(sigLL0, low_pass, X0LL);

				vector<double> X0LH;
    			convfft(sigLH0, low_pass, X0LH);

				vector<double> X0HL;
    			convfft(sigHL0, high_pass, X0HL);

				vector<double> X0HH;
    			convfft(sigHH0, high_pass, X0HH);

				int lent=X0LL.size();

				for (int i=0; i <  lent; i++) {

					cLLt10[i][j]=X0LL[i];
					cLHt10[i][j]=X0LH[i];
					cHLt10[i][j]=X0HL[i];
					cHHt10[i][j]=X0HH[i];

				}




			}



			for (int i=0; i< rowsLL + lf-1;i++) {
			    vector<double> sigLL0;
        		vector<double> sigLH0;
        		vector<double> sigHL0;
        		vector<double> sigHH0;

				for (int j=0; j < colsLL; j++) {
					double temp01 = cLLt10[i][j];
					double temp02 = cLHt10[i][j];
					double temp03 = cHLt10[i][j];
					double temp04 = cHHt10[i][j];
					sigLL0.push_back(temp01);
					sigLH0.push_back(temp02);
					sigHL0.push_back(temp03);
					sigHH0.push_back(temp04);

				}
				

				// Second Convolution Step for LL,LH,HL,HH

    			vector<double> X0LL;
    			convfft(sigLL0, low_pass, X0LL);

				vector<double> X0LH;
    			convfft(sigLH0, high_pass, X0LH);

				vector<double> X0HL;
    			convfft(sigHL0, low_pass, X0HL);

				vector<double> X0HH;
    			convfft(sigHH0, high_pass, X0HH);

				int lent=X0LL.size();

				for (int j=0; j <  lent; j++) {

					cLLt11[i][j]=X0LL[j];
					cLHt11[i][j]=X0LH[j];
					cHLt11[i][j]=X0HL[j];
					cHHt11[i][j]=X0HH[j];
					cLLt1[i][j]=X0LL[j]+X0LH[j]+X0HL[j]+X0HH[j];

				}




			}




			for (int i=0; i < len_r;i++) {

				for (int j=0; j < len_c;j++) {

					oupt0[i][j]=cLLt0[lf+i-1][lf+j-1];
					oupt1[i][j]=cLLt1[lf+i-1][lf+j-1];



				}

			}


			circshift2d(oupt1,-1,-1);
			

			for (int i=0; i < len_r;i++) {

				for (int j=0; j < len_c;j++) {

					double temp=oupt0[i][j]+oupt1[i][j];
					temp=temp/2.0;
					oupt[i][j]=temp;



				}

			}

			 row1 = 0;
			 col1 = 0;

			 cout << rows_n << cols_n << " " << oupt.size() << oupt[0].size() << endl;
			for (int index_r = count; index_r < rows_n; index_r+=value){
				for (int index_c = count; index_c < cols_n; index_c+=value){
					cLL[index_r][index_c]=oupt[row1][col1];
					col1++;

						}
					col1 = 0;
					row1++;
				}

			cout << "STAGE 14" << endl;

			}

		if (iter == 1) {
			sig=cLL;
		}
		}

return 0;

}