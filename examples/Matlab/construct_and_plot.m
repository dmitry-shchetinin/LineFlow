%% Description
% This example shows how to plot a linear approximation of a line flow 
% constraint. For plotting, Matlab function 'plot_constraints.m' is used. 
% Its example usage and description of specific inputs are given below.

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
%to use default options, set options=[];


%% Construct approximation
Result = LF_linearize_line(branch, flow_side, options);


%% Plot the approximation
%set plotting options (all are required)
plot_opts.plot_surface_upper=1; %if 1, upper part of the surface is plotted
plot_opts.plot_planes_upper=1; %if 1, upper planes are plotted
plot_opts.plot_surface_lower=0; %if 1, lower part of the surface is plotted
plot_opts.plot_planes_lower=0; %if 1, lower planes are plotted
plot_opts.N_samples_V=50; %number of sample points of V along each V axis
plot_opts.margin_V_axis=0.05; %margin from surface to plot's edge in p.u. along V axis
plot_opts.margin_delta_axis=0.05; %margin from surface to plot's edge in p.u. along theta axis
plot_opts.fontsize=16; %fontsize for axes
plot_opts.linewidth=1; %width of lines at 2D plots

%plot
plot_constraints (Result, branch, flow_side, options, plot_opts);





