/*==========================================================
 * CountOccurance.c
 *
 * Inputs:
 * inList - The input list of vertices
 * inFaces - The input list of faces
 *
 * Task:
 * Counts how often all given vertices are in the list of faces
 *
 * Output:
 * counter - the result of the operation above
 *
 * The calling syntax is:
 *
 * counter = CountOccurance(inList, inFaces)
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void CountOccurance(int *inList, int *inFaces, double *result, mwSize N, mwSize M)
{    
    for (int i=0; i < N; i++) {
        for (int j=0; j<3*M; j++) {
            if (inList[i] == inFaces[j]) {
                *result = *result+1;
            }
        }        
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int *inList;                    /* 1xN input list */
    int *inFaces;                   /* 3xM input vertices */
    int *result;                    /* output counter */
    int N;
    int M;

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* make sure the first input argument is type integer */
    if( !mxIsNumeric(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input list 1 must be integer.");
    }
    
    /* make sure the second input argument is type integer */
    if( !mxIsNumeric(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input vertices 2 must be integer.");
    }
    
    /* check that number of rows in first input argument equals 1 */
    if(mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 1 must be a list.");
    }
    
    /* check that number of rows in second input argument equals 3 */
    if(mxGetM(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be 3xM.");
    }

    /* create a pointer to the integer data in the input list  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inList = mxGetInt32s(prhs[0]);
    #else
    inList = mxGetPr(prhs[0]);
    #endif
    
    /* create a pointer to the integer data in the input vertices  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    inFaces = mxGetInt32s(prhs[1]);
    #else
    inFaces = mxGetPr(prhs[1]);
    #endif
    
    /* get dimensions of the input matrix */
    N = mxGetN(prhs[0]);
    M = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0.0);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    result = mxGetDoubles(plhs[0]);
    #else
    result = mxGetPr(plhs[0]);
    #endif

    /* call the computational routine */
    CountOccurance(inList, inFaces, result, N, M);
}
