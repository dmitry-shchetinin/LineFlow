/*OVERVIEW: 
This example shows how to construct a linear approximation of a single line 
flow constraint in C/C++ using the following library function:
Result = LF_construct(LF_Branch *branch, int flow_side, LF_Options *options)
The function description is provided below

INPUTS:
 branch - structure of branch parameters with the following fields:
    Name      Type     Description 
   .g        (double)  branch conductance in p.u.
   .b        (double)  branch susceptance in p.u.
   .b_sh     (double)  total branch shunt susceptance in p.u.
   .t_ratio  (double)  transformer's tap ratio (set 0 if not a transformer)
   .t_shift  (double)  transformer's phase shift (set 0 if not a transformer)
   .I_max    (double)  value of thermal limit in p.u. (set 0 if unlimited)
   .V_i_min  (double)  lower bound on V at the beginning of the line in p.u.
   .V_j_min  (double)  lower bound on V at the end of the line in p.u.
   .V_i_max  (double)  upper bound on V at the beginning of the line in p.u.
   .V_j_max  (double)  upper bound on V at the end of the line in p.u.

 flow_side - scalar value showing the line flow constraint at which end
             of the line must be approximated:
             1 - beginning of the line
             2 - end of the line
             3 - both beginning and end of the line

 options - structure of algorithm's options. To use default options, set 
           options=NULL. To override a particular option, use the 
		   corresponding provided setter function. The fields are:
    Name                Type     Description
   .approximation      (enum)    approximation type:
                                 conservative (0) - inner (default)
                                 relaxed (1)      - outer

   .N_constraints_max  (int)     maximum number of linear constraints for one
                                 part (upper/lower) of the approximation at 
								 one end of the line (default 15)

   .error_max          (double)  maximum error in current magnitude in 
                                 percent (default 5)

   .computation_mode   (int)     mode of computing the approximation:
                                 1 - algorithm constructs the approximation 
								     with number of linear constraints equal 
									 to N_constraints_max, regardless of the
									 resulting approximation error
                                 2 - algorithm iteratively increases the 
								     number of constructed linear constraints 
									 until error is smaller than error_max, 
									 number of linear constraints equals 
									 N_constraints_max, or the error change at 
									 two consecutive iterations is smaller than 
									 max_error_change (default)

   .N_adjustments      (int)     maximum number of adjustments carried out to
                                 equalize the approximation errors associated
								 with individual linear constraints (default 4)

   .ratio_threshold    (double)  threshold value of the error ratio, which is
                                 defined as the minimum error associated with 
								 an individual linear constraint divided by the
								 maximum error associated with an individual
								 constraint. If actual ratio is higher than the
								 threshold, the adjustments stop (default 0.9)

   .max_error_change   (double)  threshold value of the change of maximum error 
                                 at two consecutive iterations. If the actual 
								 value is smaller than the threshold, the 
								 algorithm stops (default 0.1)

   .delta_max_user     (double)  threshold value of phase angle difference in
                                 degrees. The approximation is only constructed
								 for the points on the boundary surface, at 
								 which the angle value does not exceed the 
								 threshold (default 85.0)

   .tr_model_type      (int)     type of shunt model for transformers:
                                 0 - the given shunt susceptance is equally 
								     split between two line's ends (default)
                                 1 - the given shunt susceptance is only put at
								     the beginning of the line. This helps 
									 model the reactive power loss associated 
									 with the magnetizing current

   .eps_tolerance      (double)  threshold value of the change of the step in 
                                 the bisection algorithm and Newton-Raphson 
								 algorithm. If the actual step is smaller than 
								 the threshold, algorithms stop (default 1e-4)

   .iter_Max           (int)     maximum number of iterations in the bisection
                                 algorithm and Newton-Raphson algorithm, which 
								 are used by low-level functions (default 25)

OUTPUTS:
 Result - structure with the following fields:
    Name     Type      Description
   .Ncon    (int)      number of constructed linear constraints

   .A       (*double)  vector of constraints' normals with (3xNcon) elements.
                       For constraint with index k (starting from 0), A[3*k] 
					   is the coefficient for V_i, A[3*k+1] is the coefficient
					   for V_j, A[3*k+2] is the coefficient for theta_ij

   .c       (*double)  vector of constraints' 'offsets' with (Ncon) elements. 
                       Thus, original nonlinear constraint is replaced by 
					   A*[V_i;V_j;theta_ij]<=c

   .error   (double)   estimate of the maximum approximation error

   .flag    (enum)     shows the result of the approximation algorithm:
                       non_binding (0)       - constraint is never binding
                       infeasible (1)        - constraint is infeasible
                       success (2)           - approximation was constructed
					   error_branch_data (3) - error in the branch parameters
					   error_options (4)     - error in algorithm's options
					   zero_limit (5)        - branch has no thermal limit
					   error_other (6)       - other (currently not used)

   .message  (*char)   short information message regarding the result of 
                       the approximation algorithm

Short descriptions of other library functions are available in both the 
header file and source code. */


/*EXAMPLE USAGE*/


/*Uncomment the lines below if C++ is used*/
#include "stdafx.h"  
#ifdef __cplusplus
extern "C" {
#endif
#include "line_flow.h"
#ifdef __cplusplus
}
#endif 
using namespace std;


int main() {
	//input branch parameters
	LF_Branch branch;
	branch.g = 1.0;
	branch.b = -3.0;
	branch.b_sh = 0.01;
	branch.t_ratio = 0;
	branch.t_shift = 0;
	branch.I_max = 1.0;
	branch.V_i_min = 0.9;
	branch.V_i_max = 1.1;
	branch.V_j_min = 0.9;
	branch.V_j_max = 1.1;
	//alternatively, one can use the following function:
	/*LF_set_branch_parameters(V_i_min, V_i_max, V_j_min, V_j_max,
	  g, b, b_sh, t_ratio, t_shift, I_max, &branch);*/

	//select end of the line for which the approximation should be constructed
    int flow_side = 1;

	//override default options if need be (for default options, set options=NULL)
	LF_Options *options = LF_get_default_options();
	LF_set_approximation_type(1, options);
	LF_set_error_max(2.0, options);
	LF_set_N_constraints_max(20, options);

	//construct approximation
	LF_Results *results = LF_construct(&branch, flow_side, options);

	//get fields of the results structure
	int N_constraints = LF_get_number_constraints(results);
	double error = LF_get_error(results);
	double *A = LF_get_A_matrix(results);
	double *b = LF_get_b_vector(results);
	char *message = LF_get_message(results);
	LF_ResultFlag return_flag = LF_get_flag(results);

	//print results to screen
	printf("%s\n", message);
	printf("--------------------------------------------------\n");
	printf("Numeric value of return flag is %d.\n", return_flag);
	printf("--------------------------------------------------\n");
	printf("Number of constructed linear constraints is %d.\n", N_constraints);
	printf("--------------------------------------------------\n");
	printf("Estimate of maximum approximation error is %5.2f.\n", error);
	printf("--------------------------------------------------\n");
	if (N_constraints > 0) {
		printf("Constraints' normals and offsets:\n  V_i     V_j   theta   offset\n");
		for (int i = 0; i < N_constraints; i++)
			printf("%5.2f   %5.2f   %5.2f   %5.2f\n", A[3 * i], A[3 * i + 1], A[3 * i + 2], b[i]);
	}

	//free memory for dynamically allocated structures
	LF_free_results(results);
	free(options);
	
	printf("\n-------------------------------\nPress ENTER to exit the program\n");
	getchar();
	return 0;
}

