/* iradon_c.c: 
   sub-routine of a modified iradon.m, i.e., the time consuming loop 
   of this routine in C. 
 
   Compilation from Matlab:
   mex iradon_c.c
   maybe a tiny bit faster code is generated by 
   mex -O COPTIMFLAGS='-O2' LDOPTIMFLAGS='-O2' iradon_c.c
 
   Usage from Matlab:
   iradon_c( p, theta, x, y );
 */

#include "mex.h"

#include <math.h> 


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int dim1, dim2, dim1_data;
  int i, no_of_angles, ctrIdx;
  double *data, *theta, *xorg, *yorg, *imgorg;
  
  /* Check for proper number of arguments. */
  if (nrhs != 4) 
    mexErrMsgTxt("Four input arguments required: Data, theta, x and y.");
  else if (nlhs != 1) 
    mexErrMsgTxt("One output argument has to be specified.");

  /* Input must be double. */
  if (mxIsDouble(prhs[0]) != 1)
    mexErrMsgTxt("Input 1 (data) must be of double precision floating point type.");
  if (mxIsDouble(prhs[1]) != 1)
    mexErrMsgTxt("Input 2 (theta) must be of double precision floating point type.");
  if (mxIsDouble(prhs[2]) != 1)
    mexErrMsgTxt("Input 3 (x) must be of double precision floating point type.");
  if (mxIsDouble(prhs[3]) != 1)
    mexErrMsgTxt("Input 4 (y) must be of double precision floating point type.");


  /* get number of different angles */
  if (mxGetM(prhs[1]) == 1) {
    no_of_angles = mxGetN(prhs[1]);
  } else {
    if (mxGetN(prhs[1]) == 1) {
      no_of_angles = mxGetM(prhs[1]);
    } else {
      mexErrMsgTxt("Theta has to be a vector, not an array.");
    }
  }
  
  /* get dimensions and check that they are consistent */
  dim1 = mxGetM(prhs[2]);
  dim2 = mxGetN(prhs[2]);
  if ((dim1 != mxGetM(prhs[3])) || (dim1 != mxGetN(prhs[3]))) 
    mexErrMsgTxt("x and y must have the same dimensions.");
  if (no_of_angles > mxGetN(prhs[0]))
    mexErrMsgTxt("The second dimension of data must be at least as large as the number of theta angles.");
dim1_data = mxGetM(prhs[0]);
  
  /* allocate memory for image data, to be returned */
  plhs[0] = 
    mxCreateNumericMatrix(dim1, dim2, mxDOUBLE_CLASS, mxREAL);
  if (plhs[0] == NULL)
    mexErrMsgTxt("Could not allocate memory for return data.");

  /* get pointers to input and output data */
  data = mxGetPr(prhs[0]);
  theta = mxGetPr(prhs[1]);
  xorg = mxGetPr(prhs[2]);
  yorg = mxGetPr(prhs[3]);
  imgorg = mxGetPr(plhs[0]);

  /* index to image center */
  ctrIdx = ceil(mxGetM(prhs[0]) / 2);
  
  for (i=0; i < no_of_angles; i++) {
    double *x = xorg;
    double *y = yorg;
    double *img = imgorg;
    /* temporary variables */
    double costheta = cos(*theta);
    double sintheta = sin(*theta);
    double *proj = &data[i*dim1_data +1];
    int j;
    for (j=0; j < dim2; j++) {
      int k;
      for (k=0; k < dim1; k++) {
        double t = *x * costheta + *y * sintheta;
        int a = floor(t);
        *img += (t-a) * proj[a+ctrIdx] + (a+1-t) * proj[a+ctrIdx-1];
        x++;
        y++;
        img++;
      }
    }
    theta++;
  }
  return;
}
