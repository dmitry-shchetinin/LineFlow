%% Description
% OVERVIEW: this example shows how to construct linear approximation of all
% line flow constraints in the system using the following mex function:
% Result=LF_linearize_system(branches,flow_side,options,ind_V,ind_theta,Ncol)
% Below is the function description
% INPUTS: 
%  branches - structure of parameters of branches, all fields are required.
%             Each field is a (n x 1) vector, where n is the total number 
%             of branches in the system:
%    .g        - branch conductances in p.u.
%    .b        - branch susceptances in p.u.
%    .b_sh     - total branch shunt susceptances in p.u.
%    .t_ratio  - transformer's tap ratios (set 0 if not a transformer)
%    .t_shift  - transformer's phase shifts (set 0 if not a transformer)
%    .I_max    - values of thermal limits in p.u. (set 0 if unlimited)
%    .V_i_min  - lower bounds on V at the beginning of the line in p.u.
%    .V_j_min  - lower bounds on V at the end of the line in p.u.
%    .V_i_max  - upper bounds on V at the beginning of the line in p.u.
%    .V_j_max  - upper bounds on V at the end of the line in p.u.
%    .ind_bus1 - for a given branch, the corresponding element contains the
%                index of the bus that the beginning of this branch is 
%                connected to in the list of all buses (must have int type)
%    .ind_bus2 - for a given branch, the corresponding element contains the
%                index of the bus that the end of this branch is connected 
%                to in the list of all buses (must have int type)
%  flow_side - scalar value showing the line flow constraints at which end 
%              of the lines must be approximated (required):
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
%  ind_V - (m x 1) vector, where m is the number of buses in the system.
%          The i-th element contains the index of the voltage magnitude of
%          i-th bus in the vector of variables for the OPF problem
%  ind_theta - (m x 1) vector. The i-th element contains the index of the 
%              voltage phase angle of i-th bus in the vector of variables 
%              for the OPF problem
%  Ncol - scalar, whose value is equal to the total number of elements in 
%         the vector of variables for the OPF problem
% OUTPUTS: 
%  Result - structure with the following fields:
%    .Ncons  - (n x 1) vector, where n is the total number of branches in 
%              the system. i-th element contains the number of linear
%              constraints constructed for i-th branch
%    .errors - (n x 1) vector of estimates of maximum appproximation errors
%              for all branches
%    .A      - (sum(Result.Ncons) x Ncol) sparse matrix of the normals of
%              constructed constraints
%    .c      - (sum(Result.Ncons) x Ncol) vector of constraints' 'offsets'.
%              Thus, original nonlinear line flow constraints are replaced
%              by Ax<=c, where x is the vector of variables for OPF problem

clear;

%% Define input parameters
%parameters are taken from the 9-bus system ('case9' in Matpower)
branches.g=[0;1.942191248714727;1.282009138424115;0;1.155087480890097; ...
    1.617122473246136;0;1.187604379291148;1.365187713310580];
branches.b=[-17.361111111111110;-10.510682051867931;-5.588244962361526;...
    -17.064846416382252;-9.784270426363173;-13.697978596908444;-16;...
    -5.975134533308591;-11.604095563139930];
branches.b_sh=[0;0.158;0.358;0;0.209;0.149;0;0.306;0.176];
branches.t_ratio=[0;0;0;0;0;0;0;0;0];
branches.t_shift=[0;0;0;0;0;0;0;0;0];
branches.V_i_min=[0.9;0.9;0.9;0.9;0.9;0.9;0.9;0.9;0.9];
branches.V_j_min=[0.9;0.9;0.9;0.9;0.9;0.9;0.9;0.9;0.9];
branches.V_i_max=[1.1;1.1;1.1;1.1;1.1;1.1;1.1;1.1;1.1];
branches.V_j_max=[1.1;1.1;1.1;1.1;1.1;1.1;1.1;1.1;1.1];
branches.I_max=[2.5;2.5;1.5;3;1.5;2.5;2.5;2.5;2.5];
branches.ind_bus1=int32([1;4;5;3;6;7;8;8;9]);
branches.ind_bus2=int32([4;5;6;6;7;8;2;9;4]);

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

%indices of voltage magnitudes and angles
ind_V=[10;11;12;13;14;15;16;17;18];
ind_theta=[1;2;3;4;5;6;7;8;9];
%total number of variables in the optimization vector
Ncol=24;


%% Construct approximation
Result = LF_linearize_system(branches, flow_side, options, ind_V, ind_theta, Ncol);





