/*==========================================================
 * MakeFaces.c
 *
 * Inputs:
 * inList - The input list of vertices
 * inFaces - The input list of faces
 * length - Length of the result
 *
 * Task:
 * Fills the result with faces that the vertices are in
 *
 * Output:
 * outFaces - Faces with the vertices in them
 *
 * The calling syntax is:
 *
 * outFaces = MakeFaces(inList, inFaces, length)
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void MakeFaces(int *inList, int *inFaces, double *result, int length, mwSize N, mwSize M)
{
    int counter = 0;
    int* mask = (int*)malloc(M*sizeof(int));
    int k = 0;
    for (int i=0; i < N; i++) {
        
        for (int j=0; j<M; j++) {
            mask[j] = 0;
        } 
        
        for (int j=0; j<3*M; j++) {
            if (inList[i] == inFaces[j]) {
                k = j / 3;
                mask[k] = 1;
            }
        }
        
        for (int j=0; j<M; j++) {
            if (mask[j] == 1) {
                result[counter] = (double)inFaces[3*j];
                result[counter+1] = (double)inFaces[3*j+1];
                result[counter+2] = (double)inFaces[3*j+2];
                result[counter+3] = j + 1; // +1 for MATLAB/FORTRAN counting
                counter = counter + 4;
            }
        }
    }
    free(mask);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int *inList;                    /* 1xN input list */
    int *inFaces;                   /* 3xM input vertices */
    int *length;                    /* Length of result */
    double *result;                 /* output*/
    int N;
    int M;

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
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
    
    if( !mxIsNumeric(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input vertices 3 must be integer.");
    }
    
    /* check that number of rows in first input argument equals 1 */
    if(mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 1 must be a list.");
    }
    
    /* check that number of rows in second input argument equals 3 */
    if(mxGetM(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be 3xM.");
    }
    
    if(mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 3 must be 1x1.");
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
    
    /* create a pointer to the integer data in the input vertices  */
    length = mxGetData(prhs[2]);
       
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(4, *length,mxREAL);
        
    /* get a pointer to the real data in the output matrix */    
    #if MX_HAS_INTERLEAVED_COMPLEX
    result = mxGetDoubles(plhs[0]);
    #else
    result = mxGetPr(plhs[0]);
    #endif
    
    /* call the computational routine */
    MakeFaces(inList, inFaces, result, *length, N, M);
}