/*=================================================================
*this function produces a conservative or relaxed linear approximation to 
*line flow constraint for one line
*=================================================================*/

#include "mex.h"
#include "stdint.h"
#include "matrix.h"
#include "string.h"
#include "line_flow.h"


/*extract value from the field with a given name from the Matlab structure*/
double value_of_structure_field(mxArray *struct_ptr, const char *fieldname, int index) {
    const mxArray *ptr;
    double *tmp;
    double result=0;
    ptr = mxGetField(struct_ptr, 0, fieldname);
    if (ptr == NULL)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "At least one required field of input structure is empty.");
    else if (mxGetNumberOfElements(ptr)<=index)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "Desired index exceeds the number of branches in the list");
    else {
        tmp=mxGetPr(ptr);
        result=tmp[index];
    }
    return result;
}

/*create C branch structure from provided Matlab structure*/
void branch_mat_to_C_structure(mxArray *struct_ptr, int index, LF_Branch *branch) {
    double V_i_min, V_i_max, V_j_min, V_j_max, g, b, b_sh, t_ratio, t_shift, I_max;
    //extract all values from the MATLAB branch structure
    V_i_min=value_of_structure_field(struct_ptr, "V_i_min", index);
    V_i_max=value_of_structure_field(struct_ptr, "V_i_max", index);
    V_j_min=value_of_structure_field(struct_ptr, "V_j_min", index);
    V_j_max=value_of_structure_field(struct_ptr, "V_j_max", index);
    g=value_of_structure_field(struct_ptr, "g", index);
    b=value_of_structure_field(struct_ptr, "b", index);
    b_sh=value_of_structure_field(struct_ptr, "b_sh", index);
    t_ratio=value_of_structure_field(struct_ptr, "t_ratio", index);
    t_shift=value_of_structure_field(struct_ptr, "t_shift", index);
    I_max=value_of_structure_field(struct_ptr, "I_max", index);
    
    //create C branch structure
    LF_set_branch_parameters(V_i_min, V_i_max, V_j_min, V_j_max,
            g, b, b_sh, t_ratio, t_shift, I_max, branch);
}

/*override default options with options provided by the user*/
void override_default_options(mxArray *struct_ptr, LF_Options *options) {
    const mxArray *ptr;
    double *tmp;
    ptr = mxGetField(struct_ptr, 0, "computation_mode");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_mode((int)tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "iter_Max");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_iter_max((int)tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "N_constraints_max");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_N_constraints_max((int)tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "N_adjustments");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_N_adjustments((int)tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "tr_model_type");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_transformer_model_type((int)tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "delta_max_user");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_delta_max_user(tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "eps_tolerance");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_eps_tolerance(tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "error_max");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_error_max(tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "max_error_change");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_max_error_change(tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "ratio_threshold");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_ratio_threshold(tmp[0], options);
    }
    ptr = mxGetField(struct_ptr, 0, "approximation");
    if (ptr) {
        tmp=mxGetPr(ptr);
        LF_set_approximation_type((int)tmp[0], options);
    }
}




/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    //initialize variables
    const char *fieldnames[6];
    LF_Branch branch;
    LF_Results *results;
    LF_Options *options;
    mxArray *A, *c;
    int N_lin_con, flow_side, index;
    double *v1, *v2, *A_ptr, *c_ptr;

    /* check proper inputs*/
    if (nrhs<2)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Not enough input arguments. Must be at least 2.");
    else if (nrhs>4)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Too many input arguments. Must be at most 4.");
    else if(!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct", "First input must be a structure.");
    else if(nrhs>=3 && !mxIsStruct(prhs[2]) && (mxGetM(prhs[2]) != 0 && mxGetN(prhs[2]) != 0))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct", "Third input must be a structure.");
    else if (nrhs==4 && mxGetScalar(prhs[3])<1)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Index of branch in the list has to be positive.");
        
    if (nrhs==4)
        index=mxGetScalar(prhs[3])-1;
    else
        index=0;
    flow_side=mxGetScalar(prhs[1]);
    
    //create branch C structure from provided Matlab structure
    branch_mat_to_C_structure(prhs[0], index, &branch);
    
    //create default options structure
    options=LF_get_default_options();
    
    //override default options by the options provided by the user
    if (nrhs>=3 && (mxGetM(prhs[2]) != 0 && mxGetN(prhs[2]) != 0))
        override_default_options(prhs[2], options);
    
   //call c function that constructs the approximation
   results=LF_construct(&branch, flow_side, options);
   
   //get values of (or pointers to) the fields of the results structure
   v1=LF_get_A_matrix(results);
   v2=LF_get_b_vector(results);
   N_lin_con=LF_get_number_constraints(results);
   
   
   //Initialize output structure
   //assign field names
   for (int i=0; i<6; i++)
       fieldnames[i] = (char*)mxMalloc(10);
   memcpy(fieldnames[0],"A", sizeof("A"));
   memcpy(fieldnames[1],"c", sizeof("c"));
   memcpy(fieldnames[2],"Ncon", sizeof("Ncon"));
   memcpy(fieldnames[3],"error", sizeof("error"));
   memcpy(fieldnames[4],"flag", sizeof("flag"));
   memcpy(fieldnames[5],"message", sizeof("message"));
   
   //allocate memory for the structure
   plhs[0] = mxCreateStructMatrix(1,1,6,fieldnames);
   
   //deallocate memory for the fieldnames
   for (int i=0; i<6; i++)
       mxFree(fieldnames[i]);
   
   A=mxCreateNumericMatrix(N_lin_con, 3, mxDOUBLE_CLASS, mxREAL);
   A_ptr=mxGetPr(A);
   c=mxCreateNumericMatrix(N_lin_con, 1, mxDOUBLE_CLASS, mxREAL);
   c_ptr=mxGetPr(c);
   
   //fill out output arrays
   for (int i=0; i<N_lin_con; i++)
       c_ptr[i]=v2[i];
  
   for (int j=0; j<3; j++) {
       for (int i=0; i<N_lin_con; i++)
           A_ptr[j*N_lin_con+i]=v1[i*3+j];
   }  
   
   //record everything into the output structure
   mxSetFieldByNumber(plhs[0], 0, 0, A);
   mxSetFieldByNumber(plhs[0], 0, 1, c);
   mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar(N_lin_con));
   mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar(LF_get_error(results)));
   mxSetFieldByNumber(plhs[0], 0, 4, mxCreateDoubleScalar(LF_get_flag(results)));
   mxSetFieldByNumber(plhs[0], 0, 5, mxCreateString(LF_get_message(results)));
   
   //free memory
   LF_free_results(results);
   free(options);
}






