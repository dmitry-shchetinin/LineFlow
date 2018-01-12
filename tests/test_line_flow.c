#include "line_flow.h"

int main() {
        LF_Branch branch;
	//input branch parameters
	LF_set_branch_parameters(0.9, 1.1, 0.9, 1.1,
	  1.0, -3.0, 0.01, 0.0, 0.0, 1.0, &branch);

	//construct approximation
	LF_Results *results = LF_construct(&branch, 3, NULL);

	//get fields of the results structure
	int N_constraints = LF_get_number_constraints(results);
	double error = LF_get_error(results);
	double *A = LF_get_A_matrix(results);
	double *b = LF_get_b_vector(results);
	char *message = LF_get_message(results);
	LF_ResultFlag return_flag = LF_get_flag(results);

	//print results to screen
        printf("--------------------------------------------------\n");
        printf("THIS IS A TEST.\n");
        printf("--------------------------------------------------\n");
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
	
	printf("\n-------------------------------\nPress ENTER to exit the test\n");
	getchar();
	return 0;
}

