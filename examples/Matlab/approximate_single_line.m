%% Description
% OVERVIEW: this example shows how to construct a linear approximation of a
% single line flow constraint using the following mex function:
% Result = LF_linearize_line(branch, flow_side, options)
% Below is the function description
% INPUTS: 
%  branch - structure of branch parameters, all fields are required:
%    .g       - branch conductance in p.u.
%    .b       - branch susceptance in p.u.
%    .b_sh    - total branch shunt susceptance in p.u.
%    .t_ratio - transformer's tap ratio (set 0 if not a transformer)
%    .t_shift - transformer's phase shift (set 0 if not a transformer)
%    .I_max   - value of thermal limit in p.u. (set 0 if unlimited)
%    .V_i_min - lower bound on V at the beginning of the line in p.u.
%    .V_j_min - lower bound on V at the end of the line in p.u.
%    .V_i_max - upper bound on V at the beginning of the line in p.u.
%    .V_j_max - upper bound on V at the end of the line in p.u.
%  flow_side - scalar value showing the line flow constraint at which end 
%              of the line must be approximated (required):
%              1 - beginning of the line
%              2 - end of the line
%              3 - both beginning and end of the line
%  options - structure of algorithm's options, all fields are optional. To
%            use default options, set options=[] :
%    .approximation     - approximation type:
%                         0 - inner (conservative) approximation (default)
%                         1 - outer (relaxed) approximation
%    .N_constraints_max - maximum number of linear constraints for one part 
%                         (upper/lower) of the approximation at one end of 
%                         the line (default 15)
%    .error_max         - maximum error in current magnitude in percent
%                         (default 5)
%    .computation_mode  - mode of computing the approximation:
%                         1 - algorithm constructs the approximation with 
%                             number of linear constraints equal to 
%                             N_constraints_max, regardless of the 
%                             resulting approximation error
%                         2 - algorithm iteratively increases the number of 
%                             constructed linear constraints until error is
%                             smaller than error_max, number of linear 
%                             constraints equals N_constraints_max, or the 
%                             error change at two consecutive iterations is 
%                             smaller than max_error_change (default)
%    .N_adjustments     - maximum number of adjustments carried out to 
%                         equalize the approximation errors associated 
%                         with individual linear constraints (default 4)
%    .ratio_threshold   - threshold value of the error ratio, which is 
%                         defined as the minimum error associated with an 
%                         individual linear constraint divided by the 
%                         maximum error associated with an individual 
%                         constraint. If actual ratio is higher than the 
%                         threshold, the adjustments stop (default 0.9)
%    .max_error_change  - threshold value of the change of maximum error at 
%                         two consecutive iterations. If the actual value 
%                         is smaller than the threshold, the algorithm 
%                         stops (default 0.1)
%    .delta_max_user    - threshold value of phase angle difference in 
%                         degrees. The approximation is only constructed 
%                         for the points on the boundary surface, at which 
%                         the angle value does not exceed the threshold
%                         (default 85.0)
%    .tr_model_type     - type of shunt model for transformers:
%                         0 - the given shunt susceptance is equally split 
%                             between two ends of the line (default)
%                         1 - the given shunt susceptance is only put at 
%                             the beginning of the line. This helps model 
%                             the reactive power loss associated with the 
%                             magnetizing current
%    .eps_tolerance     - threshold value of the change of the step in the 
%                         bisection algorithm and Newton-Raphson algorithm. 
%                         If the actual step is smaller than the threshold, 
%                         the algorithms stop (default 1e-4)
%    .iter_Max          - maximum number of iterations in the bisection 
%                         algorithm and Newton-Raphson algorithm, which are 
%                         used by low-level functions (default 25)
% OUTPUTS: 
%  Result - structure with the following fields:
%    .Ncon    - number of constructed linear constraints
%    .A       - (Ncon x 3) dense matrix of constraints' normals. The first 
%               column contains the coefficients for V_i, the second column
%               contains the coefficients for V_j, the last column contains
%               coefficients for phase angle difference theta_ij
%    .c       - (Ncon x 1) vector of constraints' 'offsets'. Thus, original
%               nonlinear constraint is replaced by A*[V_i;V_j;theta_ij]<=c
%    .error   - estimate of the maximum approximation error
%    .flag    - shows the result of the approximation algorithm:
%               0 - The nonlinear constraint cannot become binding
%               1 - The nonlinear constraint is infeasible
%               2 - The approximation was successfully constructed
%               3 - There was an error in the input branch parameters
%               4 - There was an error in the input algorithm's options
%               5 - There is no thermal limit for the given branch
%               6 - Other (currently not used)
%    .message - short information message regarding the algorithm's result

clear;

%% Define input parameters
%branch parameters
branch.g=1;
branch.b=-3;
branch.b_sh=0.01;
branch.t_ratio=0;
branch.t_shift=0;
branch.I_max=0.5;
branch.V_i_min=0.9;
branch.V_i_max=1.1;
branch.V_j_min=0.9;
branch.V_j_max=1.1;

%end of the line for which the approximation should be constructed
flow_side=1;

%algorithm's options
options.approximation=0;
options.N_constraints_max=15;
options.error_max=5;
options.computation_mode=2;
options.N_adjustments=4;
options.ratio_threshold=0.9;
options.max_error_change=0.1;
options.delta_max_user=85.0;
options.transformer_model_type=0;
options.eps_tolerance=1e-4;
options.iter_Max=25;
%to use default options, set options=[];


%% Construct approximation
Result = LF_linearize_line(branch, flow_side, options);





