/*=================================================================
 *this function produces a linear approximation to line flow constraints for the given system
 *=================================================================*/


#include "mex.h"
#include "stdint.h"
#include "matrix.h"
#include "line_flow.h"


/*check if the number of elements in a certain structure fields is equal to the given number*/
void check_number_of_elements_in_struct_field(mxArray *struct_ptr, const char *fieldname, int N_elements) {
    const mxArray *ptr;
    double *tmp;
    ptr = mxGetField(struct_ptr, 0, fieldname);
    if (ptr == NULL)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "At least one required field of input structure is empty.");
    else if (mxGetNumberOfElements(ptr)!=N_elements)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "All fields of branch structure must have the same number of elements.");
}

/*get the number of branches in the system and check that all structure fields have the same number of elements*/
int number_of_branches(mxArray *struct_ptr) {
    int N_lines=0;
    const mxArray *ptr;
    ptr = mxGetField(struct_ptr, 0, "V_i_min");
    if (ptr == NULL)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "At least one required field of input structure is empty.");
    else
        N_lines=mxGetNumberOfElements(ptr);
    check_number_of_elements_in_struct_field(struct_ptr, "V_i_max", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "V_j_min", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "V_j_max", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "g", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "b", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "b_sh", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "t_ratio", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "t_shift", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "I_max", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "ind_bus1", N_lines);
    check_number_of_elements_in_struct_field(struct_ptr, "ind_bus2", N_lines);
    if (N_lines==0)
        mexErrMsgIdAndTxt( "MATLAB:fieldEmpty", "All fields of branch structure are empty.");
}


/*override default options with options provided by the user*/
void override_default_options(mxArray *struct_ptr, LF_Options* options) {
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



void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        ) {
    //initialize variables
    LF_Branch branch;
    LF_Results *results;
    LF_Options *options;
    LF_Workspace *workspace;
    const char *fieldnames[4];
    int N_lin_con, N_lin_con_max, flow_side, N_columns, N_lines, N_nnz_I_max, N_nonlin_total, 
            counter_constraints, counter_nonzeros, ind_1, ind_2;
    int *ind_bus_1, *ind_bus_2;
    double indices_var[4], *b_vector_ptr, *rows, *cols, *vals, *A_row, *A_col, *A_val, *b_val, *g_all, *b_all, *indices_V, *indices_delta, *B_sh_all, 
            *v1, *v2, *Kt_ratio_all, *Kt_shift_all, *V_i_min_all, *V_i_max_all, *V_j_min_all, *V_j_max_all, *I_max_user_all, *max_errors_ptr, *n_con_ptr;
    mxArray *lhs1[1], *rhs1[5], *b_vector, *max_errors, *n_con;
    LF_ResultFlag flag_result;
    
    /*check proper inputs*/
    if (nrhs!=6 )
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "There have to be six input elements.");
    else if(!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct", "First input must be a structure.");
    else if (mxGetScalar(prhs[1])<1 || mxGetScalar(prhs[1])>3)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Flow side value has to be 1, 2, or 3.");
    else if(!mxIsStruct(prhs[2]) && (mxGetM(prhs[2]) != 0 && mxGetN(prhs[2]) != 0))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct", "Third input must be a structure.");
    else if (mxGetScalar(prhs[5])<1)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Number of columns in a matrix of constraints has to be positive.");
      
    //get pointers to/value of input arguments
    flow_side=mxGetScalar(prhs[1]);
    indices_V=mxGetPr(prhs[3]);
    indices_delta=mxGetPr(prhs[4]);
    N_columns=mxGetScalar(prhs[5]);
    
    //check if vectors of indices of V and delta are correct
    mexCallMATLAB(1, lhs1, 1, prhs[3], "min");
    ind_1=mxGetScalar(lhs1[0]);
    mexCallMATLAB(1, lhs1, 1, prhs[4], "min");
    ind_2=mxGetScalar(lhs1[0]);
    if (ind_1<1 || ind_2<1)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Indices of delta and V must be positive.");
    mexCallMATLAB(1, lhs1, 1, prhs[3], "max");
    ind_1=mxGetScalar(lhs1[0]);
    mexCallMATLAB(1, lhs1, 1, prhs[4], "max");
    ind_2=mxGetScalar(lhs1[0]);
    if (ind_1>N_columns || ind_2>N_columns)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Maximum index of V or delta exceeds the number of columns in constraint matrix.");
    
    //create default options structure
    options=LF_get_default_options();
    
    //override default options by the options provided by the user
    if (mxGetM(prhs[2]) != 0 && mxGetN(prhs[2]) != 0)
        override_default_options(prhs[2], options);
    
    //initialize output structure
    results = LF_initialize_results(options, flow_side);
    
    //check algorithm's options
    if (!LF_check_options(options, results, 1))
        mexErrMsgIdAndTxt( "MATLAB:maxlhs", "Algorithm's options are set incorrectly.");
    
    //get the number of lines in the system and check that all fields of branch structure have the same number of elements
    N_lines=number_of_branches(prhs[0]);
    
    //extract pointers to all fields of MATLAB structure so that it is easier to work with them later on
    g_all=mxGetPr(mxGetField(prhs[0], 0, "g"));
    b_all=mxGetPr(mxGetField(prhs[0], 0, "b"));
    B_sh_all=mxGetPr(mxGetField(prhs[0], 0, "b_sh"));
    Kt_ratio_all=mxGetPr(mxGetField(prhs[0], 0, "t_ratio"));
    Kt_shift_all=mxGetPr(mxGetField(prhs[0], 0, "t_shift"));
    V_i_min_all=mxGetPr(mxGetField(prhs[0], 0, "V_i_min"));
    V_i_max_all=mxGetPr(mxGetField(prhs[0], 0, "V_i_max"));
    V_j_min_all=mxGetPr(mxGetField(prhs[0], 0, "V_j_min"));
    V_j_max_all=mxGetPr(mxGetField(prhs[0], 0, "V_j_max"));
    I_max_user_all=mxGetPr(mxGetField(prhs[0], 0, "I_max"));
    ind_bus_1=mxGetData(mxGetField(prhs[0], 0, "ind_bus1"));
    ind_bus_2=mxGetData(mxGetField(prhs[0], 0, "ind_bus2"));
    
    
    //Initialize output structure
    //assign field names
    for (int i=0; i<4; i++)
        fieldnames[i] = (char*)mxMalloc(10);
    memcpy(fieldnames[0],"A", sizeof("A"));
    memcpy(fieldnames[1],"c", sizeof("c"));
    memcpy(fieldnames[2],"Ncons", sizeof("Ncons"));
    memcpy(fieldnames[3],"errors", sizeof("errors"));

    //allocate memory for the structure
    plhs[0] = mxCreateStructMatrix(1,1,4,fieldnames);
    
    //deallocate memory for the fieldnames
    for (int i=0; i<4; i++)
        mxFree(fieldnames[i]);
    
    //compute the number of branches that have nonzero thermal limits
    rhs1[0]=mxGetField(prhs[0], 0, "I_max");
    mexCallMATLAB(1, lhs1, 1, rhs1, "nnz");
    N_nnz_I_max=mxGetScalar(lhs1[0]);
    if (N_nnz_I_max==0) {
        mxSetFieldByNumber(plhs[0], 0, 0, mxCreateSparse(0,0,0,mxREAL));
        mxSetFieldByNumber(plhs[0], 0, 1, mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL));
        mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar(0.0));
        mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar(0.0));
        return;
    }
    
    //initialize LF_Workspace variables
    workspace=LF_initizalize_workspace(options, flow_side);
    
    //initialize intermediate vectors that are used to store approximation results
    if (flow_side==3)
        N_nonlin_total=2*N_nnz_I_max;
    else
        N_nonlin_total=N_nnz_I_max;
    N_lin_con_max=LF_get_N_constraints_max(options);
    A_row=mxMalloc(8*N_nonlin_total*N_lin_con_max*sizeof(double));
    A_col=mxMalloc(8*N_nonlin_total*N_lin_con_max*sizeof(double));
    A_val=mxMalloc(8*N_nonlin_total*N_lin_con_max*sizeof(double));
    b_val=mxMalloc(2*N_nonlin_total*N_lin_con_max*sizeof(double));
    max_errors=mxCreateNumericMatrix(N_lines, 1, mxDOUBLE_CLASS, mxREAL);
    max_errors_ptr=mxGetPr(max_errors);
    n_con=mxCreateNumericMatrix(N_lines, 1, mxDOUBLE_CLASS, mxREAL);
    n_con_ptr=mxGetPr(n_con);

    //loop over all branches and construct approximation
    counter_constraints=0;
    counter_nonzeros=0;
    for (int i=0; i<N_lines; i++) {
        //check if approximation has to be constructed
        if (I_max_user_all[i]==0)
            continue;
        
        //record indices of nodes that this branch is connected to
        ind_1=ind_bus_1[i]-1;
        ind_2=ind_bus_2[i]-1;
        
        //fill out fields of C structure
        LF_set_branch_parameters(V_i_min_all[i], V_i_max_all[i], V_j_min_all[i], V_j_max_all[i],
                g_all[i], b_all[i], B_sh_all[i], Kt_ratio_all[i], Kt_shift_all[i], I_max_user_all[i], &branch);
        
        //check branch parameters
        if (!LF_check_branch_data(&branch, results))
            printf("Caution: parameters of branch with index %i are incorrect.", i);

        //construct linear approximation
        LF_linearize_one_line(&branch, flow_side, options, workspace, results);
        n_con_ptr[i]=LF_get_number_constraints(results);
        
        //check if the linearization was successful
        if (LF_get_flag(results)==infeasible)
            printf("At least one nonlinear constraint is infeasible.\n");
        N_lin_con=LF_get_number_constraints(results);
        if (N_lin_con>0) {
            //get pointers to the fields of the results structure
            v1=LF_get_A_matrix(results);
            v2=LF_get_b_vector(results);
            //loop through all created linear constraints
            for (int j=0; j<N_lin_con; j++) {
                //Vi
                if (v1[j*3]!=0) { //only nonzero elements for voltage magnitudes
                    A_row[counter_nonzeros]=counter_constraints+1;
                    A_col[counter_nonzeros]=indices_V[ind_1];
                    A_val[counter_nonzeros]=v1[j*3];
                    counter_nonzeros=counter_nonzeros+1;
                }
                //Vj
                if (v1[j*3+1]!=0) { //only nonzero elements for voltage magnitudes
                    A_row[counter_nonzeros]=counter_constraints+1;
                    A_col[counter_nonzeros]=indices_V[ind_2];
                    A_val[counter_nonzeros]=v1[j*3+1];
                    counter_nonzeros=counter_nonzeros+1;
                }
                //delta_i
                A_row[counter_nonzeros]=counter_constraints+1;
                A_col[counter_nonzeros]=indices_delta[ind_1];
                A_val[counter_nonzeros]=v1[j*3+2];
                counter_nonzeros=counter_nonzeros+1;
                //delta_j
                A_row[counter_nonzeros]=counter_constraints+1;
                A_col[counter_nonzeros]=indices_delta[ind_2];
                A_val[counter_nonzeros]=-v1[j*3+2];
                counter_nonzeros=counter_nonzeros+1;
                
                //filling out rhs info
                b_val[counter_constraints]=v2[j];
                counter_constraints=counter_constraints+1;
            }
            //record error estimate
            max_errors_ptr[i]=LF_get_error(results);
        }
    }
    
    //create and record sparse matrix
    rhs1[0]=mxCreateNumericMatrix(counter_nonzeros, 1, mxDOUBLE_CLASS, mxREAL);
    rows=mxGetPr(rhs1[0]);
    rhs1[1]=mxCreateNumericMatrix(counter_nonzeros, 1, mxDOUBLE_CLASS, mxREAL);
    cols=mxGetPr(rhs1[1]);
    rhs1[2]=mxCreateNumericMatrix(counter_nonzeros, 1, mxDOUBLE_CLASS, mxREAL);
    vals=mxGetPr(rhs1[2]);
    for (int i=0; i<counter_nonzeros; i++) {
        rows[i]=A_row[i];
        cols[i]=A_col[i];
        vals[i]=A_val[i];
    }
    rhs1[3]=mxCreateDoubleScalar(counter_constraints);
    rhs1[4]=mxCreateDoubleScalar(N_columns);
    mexCallMATLAB(1, &lhs1[0], 5, rhs1, "sparse");
    mxSetFieldByNumber(plhs[0], 0, 0, lhs1[0]);
    
    //record rhs vector
    b_vector = mxCreateNumericMatrix(counter_constraints, 1, mxDOUBLE_CLASS, mxREAL);
    b_vector_ptr = mxGetPr(b_vector);
    for (int i=0; i<counter_constraints; i++)
        b_vector_ptr[i]=b_val[i];
    mxSetFieldByNumber(plhs[0], 0, 1, b_vector);
    
    //record vectors with number of constraints and error estimates
    mxSetFieldByNumber(plhs[0], 0, 2, n_con);
    mxSetFieldByNumber(plhs[0], 0, 3, max_errors);
    
    //free memory from library-related variables
	LF_free_workspace(workspace);
    LF_free_results(results);
    free(options);
    //deallocate internal arrays
    mxFree(A_row);
    mxFree(A_col);
    mxFree(A_val);
    mxFree(b_val);
    //deallocate mex arrays
    mxDestroyArray(rhs1[0]);
    mxDestroyArray(rhs1[1]);
    mxDestroyArray(rhs1[2]);
}