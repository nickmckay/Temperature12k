// COPYRIGHT
//   2010, Sami Hanhijärvi <sami.hanhijarvi@iki.fi>
//
// LICENSE
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
// 
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "mex.h"
//#include <unistd.h>
//#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs < 2 || nrhs > 4) {
		mexErrMsgTxt("The number of inputs must be 2-4.");        
    }
    
    double *input = mxGetPr(prhs[0]);
    int numel = mxGetNumberOfElements(prhs[0]);
    
    double *add = mxGetPr(prhs[1]);
    int addNumel = mxGetNumberOfElements(prhs[1]);

    if (nrhs == 2) {
        if (numel != mxGetNumberOfElements(prhs[1])) {
            mexErrMsgTxt("Number of elements has to be the same for both inputs.");        
        }
        
        for (int i = 0; i < numel; ++i) {
            input[i] += add[i];
        }
        
    } else if (nrhs == 4) {
        // Add with sparse representation
        if ((mxGetClassID(prhs[2]) != mxINT32_CLASS) | (mxGetClassID(prhs[3]) != mxINT32_CLASS)) {
            mexErrMsgTxt("Only integer indexing is supported.");
        }
        int *rows = (int *)mxGetPr(prhs[2]);
        int *cols = (int *)mxGetPr(prhs[3]);            
        int inputRows = mxGetM(prhs[0]);
        int inputCols = mxGetN(prhs[0]);        
        
        
        if (mxGetNumberOfElements(prhs[1]) == 1) {
            // Repeate added constant
            if (mxGetNumberOfElements(prhs[2]) == 1) {
                // Repeate rows and add value
                for (int i = 0; i < mxGetNumberOfElements(prhs[3]); i++) {
                    input[(int)(rows[0]-1+(cols[i]-1)*inputRows)] += add[0];
                }                
            } else if (mxGetNumberOfElements(prhs[3]) == 1) {
                // Repeate columns and add value
                for (int i = 0; i < mxGetNumberOfElements(prhs[2]); i++) {
                    input[(int)(rows[i]-1+(cols[0]-1)*inputRows)] += add[0];
                }                
            } else {
                // Repeate add value
                for (int i = 0; i < mxGetNumberOfElements(prhs[2]); i++) {
                    input[(int)(rows[i]-1+(cols[i]-1)*inputRows)] += add[0];
                }                
            }
        } else {
            // Different value for each index
            if (mxGetNumberOfElements(prhs[2]) == 1) {
                // Repeate rows 
                for (int i = 0; i < mxGetNumberOfElements(prhs[1]); i++) {
                    input[(int)(rows[0]-1+(cols[i]-1)*inputRows)] += add[i];
                }                
            } else if (mxGetNumberOfElements(prhs[3]) == 1) {
                // Repeate columns and add value
                for (int i = 0; i < mxGetNumberOfElements(prhs[1]); i++) {
                    input[(int)(rows[i]-1+(cols[0]-1)*inputRows)] += add[i];
                }                
            } else {
                // Repeate add value
                for (int i = 0; i < mxGetNumberOfElements(prhs[1]); i++) {
                    input[(int)(rows[i]-1+(cols[i]-1)*inputRows)] += add[i];
                }                
            }
        }
    } else {
        if (mxIsLogical(prhs[2])) {
            mxLogical *mask = mxGetLogicals(prhs[2]);
            if (numel == mxGetNumberOfElements(prhs[2])) {
                // The mask is a matrix of the same size as input
                int j = 0; 
                for (int i = 0; i < numel; ++i) {
                    if (mask[i]) {
                        input[i] += add[j];            
                        ++j;
                    }
                }
            } else {
                // The mask is a column vector
                int maskNumel = mxGetNumberOfElements(prhs[2]);
                int inputRows = mxGetM(prhs[0]);
                int inputCols = mxGetN(prhs[0]);
    //            mexPrintf("%dx%d matrix, mask %d\n",inputRows, inputCols, maskNumel);

                if (maskNumel != inputRows || inputCols != inputRows) {
                    mexErrMsgTxt("The mask vector has to be a column vector and the input matrix has to be a 2-dimensional square matrix.");                  
                }

                int j = 0; 
                for (int col = 0; col < inputCols; ++col) {
                    if (mask[col]) {
                        int index = col*inputRows;
                        for (int row = 0; row < inputRows; ++row, ++index) {
                            if (mask[row]) {
    //                            mexPrintf("Added %f + %f = ", input[index], add[j]);
                                input[index] += add[j];
    //                            mexPrintf("%f\n", input[index]);
                                ++j;
                            }
                        }
                    }
                }
            }
        } else if (mxIsDouble(prhs[2])) {
            double *index = mxGetPr(prhs[2]);
            int indexNumel = mxGetNumberOfElements(prhs[2]);
            int inputRows = mxGetM(prhs[0]);
            int inputCols = mxGetN(prhs[0]);

            if (inputCols == 1) {
                mexErrMsgTxt("Unsupported.");                
            } else if (inputRows == 1) {
                mexErrMsgTxt("Unsupported.");                
            } else if (indexNumel*indexNumel == addNumel) {
                int addi = 0;
                for (int ci = 0; ci < indexNumel; ++ci) {
                    int offset = inputRows*(index[ci]-1);
                    for (int ri = 0; ri < indexNumel; ++ri, ++addi) {
                        input[offset + (int)index[ri]-1] += add[addi];
                    }
                }
            } else {
                mexErrMsgTxt("Unsupported.");                
            }
        } else if (mxGetClassID(prhs[2]) == mxINT32_CLASS) {
            int *index = (int *)mxGetPr(prhs[2]);
            int indexNumel = mxGetNumberOfElements(prhs[2]);
            int inputRows = mxGetM(prhs[0]);
            int inputCols = mxGetN(prhs[0]);

            if (inputCols == 1) {
                mexErrMsgTxt("Unsupported.");                
            } else if (inputRows == 1) {
                mexErrMsgTxt("Unsupported.");                
            } else if (indexNumel*indexNumel == addNumel) {
                int addi = 0;
                for (int ci = 0; ci < indexNumel; ++ci) {
                    int offset = inputRows*(index[ci]-1);
                    for (int ri = 0; ri < indexNumel; ++ri, ++addi) {
                        input[offset + index[ri]-1] += add[addi];
                    }
                }
            } else {
                mexErrMsgTxt("Unsupported.");                
            }
                        
        } else {
            mexErrMsgTxt("Logical input assumed.");
        }
        
    }   
}
