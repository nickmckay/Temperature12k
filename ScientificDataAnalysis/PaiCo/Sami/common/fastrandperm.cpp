// FASTRANDPERM produces uniform permutations of the input data.
// Usage:
//	Y = FASTRANDPERM(X) randomly reorders the values within X. X can be a
//	double or logical matrix. 
//	FASTRANDPERM(X,SEED)  sets the SEED to the Mersenne-Twister
//	and produces a random reordering of X with the given seed. 

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

#include "MersenneTwister.h"
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

// Fix redefinition conflict between Matlab R2010a and MS Visual Studio 2010
#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    static MTRand mt;
    
    if (nrhs != 1 && nrhs != 2) {
		mexErrMsgTxt("The number of inputs must be 1 or 2.");        
    }

    if (nlhs != 1 && nlhs != 0) {
        mexErrMsgTxt("The number of outputs must be 0 or 1.");
    }
    
    if (nrhs == 2) {
        mt.seed((long int)mxGetScalar(prhs[1]));
    }

    
    int N = mxGetN(prhs[0]);
    int M = mxGetM(prhs[0]);
    int numel = N*M;
    
    if (mxIsLogical(prhs[0])) {
        mxLogical *input = mxGetLogicals(prhs[0]);
        mxLogical *output;

        if (nlhs == 1) {        
            plhs[0] = mxCreateLogicalMatrix(M,N);
            output = mxGetLogicals(plhs[0]);
            
            memcpy(output, input, sizeof(mxLogical)*numel);
        } else if (nlhs == 0) {
            output = input;
        } 

        for (int n = numel-1;n > 0; --n) {
            int i = mt.randInt(n);
            mxLogical t = output[i];
            output[i] = output[n];
            output[n] = t;
        }
    } else if (mxIsDouble(prhs[0])) {
        double *input = mxGetPr(prhs[0]);
        double *output;
        
        if (nlhs == 1) {        
            plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
            output = mxGetPr(plhs[0]);

            memcpy(output, input, sizeof(double)*numel);
        } else if (nlhs == 0) {
            output = input;
        } 

        for (int n = numel-1;n > 0; --n) {
            int i = mt.randInt(n);
            double t = output[i];
            output[i] = output[n];
            output[n] = t;
        }
    } else {
        mexErrMsgTxt("Input data needs to be double or logical matrix.");
    }
       
}
