/*=================================================================
 *This function computes points on upper and lower parts of the surface for
 *a given vector of points on (Vi, Vj) plane
 *=================================================================*/



#include "mex.h"


void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
    //initializing variables
    double *Vi, *Vj, *branch_data, *delta_lower, *delta_upper;
    double t, D, delta_1, delta_2, w1, w2;
    mwSize N;
    
    Vi=mxGetPr(prhs[0]);
    Vj=mxGetPr(prhs[1]);
    branch_data=mxGetPr(prhs[2]);
    N=mxGetScalar(prhs[3]);
  
    //initialize output vectors
    plhs[0]=mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
    delta_lower=mxGetPr(plhs[0]);
    plhs[1]=mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
    delta_upper=mxGetPr(plhs[1]);
    
    for (int i=0; i<N; i++) {
        //set some intermediate variables for ease of computation
        t = branch_data[0] * Vi[i] / Vj[i] + branch_data[1] * Vj[i] / Vi[i] - branch_data[4] / (Vi[i]*Vj[i]);
        D = branch_data[5] - t*t;
        
        //check the determinant and compute the sine of delta
        if (D < 0.0) {//based on how we select points, this can happen only due to small numerical error => we assume that D=0
            delta_lower[i] = 5.0;
            delta_upper[i] = 5.0;
        }
        else {
            //compute the angle for the lower part
            delta_1 = asin(branch_data[6] * t - branch_data[7] * sqrt(D));
            delta_2 = -3.14159265358979323846 - delta_1;
            //check if the equation is satisfied
            w1 = fabs(branch_data[0] * Vi[i]*Vi[i] + branch_data[1] * Vj[i]*Vj[i] + 2.0 * Vi[i]*Vj[i]*(branch_data[2] * sin(delta_1) - branch_data[3] * cos(delta_1)) - branch_data[4]);
            w2 = fabs(branch_data[0] * Vi[i]*Vi[i] + branch_data[1] * Vj[i]*Vj[i] + 2.0 * Vi[i]*Vj[i]*(branch_data[2] * sin(delta_2) - branch_data[3] * cos(delta_2)) - branch_data[4]);
            if (w1 < w2)
                delta_lower[i]=delta_1;
            else
                delta_lower[i]=delta_2;
            
            
            //compute the angle for the upper part
            delta_1 = asin(branch_data[6] * t + branch_data[7] * sqrt(D));
            delta_2 = 3.14159265358979323846 - delta_1;
            //check if the equation is satisfied
            w1 = fabs(branch_data[0] * Vi[i]*Vi[i] + branch_data[1] * Vj[i]*Vj[i] + 2.0 * Vi[i]*Vj[i]*(branch_data[2] * sin(delta_1) - branch_data[3] * cos(delta_1)) - branch_data[4]);
            w2 = fabs(branch_data[0] * Vi[i]*Vi[i] + branch_data[1] * Vj[i]*Vj[i] + 2.0 * Vi[i]*Vj[i]*(branch_data[2] * sin(delta_2) - branch_data[3] * cos(delta_2)) - branch_data[4]);
            if (w1 < w2)
                delta_upper[i]=delta_1;
            else
                delta_upper[i]=delta_2;
        }
    }
}