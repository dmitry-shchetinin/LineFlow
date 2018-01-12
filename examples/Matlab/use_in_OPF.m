%% Description
% This example shows how to linearize the line flow constraints for a
% selected system and solve the AC OPF with linearized line flow
% constraints using Matpower.
clear;
mpopt = mpoption;

%% Select the system and solution options
%test system
system_name='case24_ieee_rts';

%OPF options
mpopt.opf.ac.solver='IPOPT';
mpopt.out.all=0; %do not print system information to screen

%choose end of the line for which the approximation should be constructed
flow_side=3;

%linearization algorithm options
options.approximation=0;
options.N_constraints_max=15; %number of linear constraints
options.error_max=5; %maximum error in current magnitude in percent

%% Load the system and preprocess data
%load the system
mpc=loadcase(system_name);

%record thermal limits and delete them from mpc structure
I_max=mpc.branch(:,6);
mpc.branch(:,6)=0;

%create structure with branch parameters
Nbus=size(mpc.bus,1);
e2i=sparse(mpc.bus(:,1),ones(Nbus,1),1:Nbus,max(mpc.bus(:,1)), 1);
ind_bus1=int32(full(e2i(mpc.branch(:,1))));
ind_bus2=int32(full(e2i(mpc.branch(:,2))));
Y=1./complex(mpc.branch(:,3),mpc.branch(:,4));
branches=struct('g', real(Y), 'b', imag(Y), 'b_sh', mpc.branch(:,5), ...
    't_ratio', mpc.branch(:,9), 't_shift', mpc.branch(:,10), 'V_i_min', ...
    mpc.bus(ind_bus1,13), 'V_j_min', mpc.bus(ind_bus2,13), 'V_i_max', ...
    mpc.bus(ind_bus1,12), 'V_j_max', mpc.bus(ind_bus2,12), 'I_max', ...
    I_max/mpc.baseMVA, 'ind_bus1', ind_bus1, 'ind_bus2', ind_bus2);

%fill out indices of voltage magnitudes and angles
ind_V=(Nbus+1:2*Nbus)';
ind_theta=(1:Nbus)';
    
%compute total number of variables in the optimization vector
Ngen=sum(mpc.gen(:,1)>-1);
Ncol=2*(Nbus+Ngen);

%% Construct linear approximation
Res_LF=LF_linearize_system(branches,flow_side,options,ind_V,ind_theta,Ncol);

%make ordering of columns of A consistent with internal Matpower orderting
Res_LF.A=e2i_data(ext2int(mpc), Res_LF.A, {'bus','bus'}, 2);

%compute the total number of created constraints
Ncon=length(Res_LF.c);

%% Run optimization
[results, success] = opf(mpc, Res_LF.A, -inf(Ncon,1), Res_LF.c, mpopt);
printpf(results, 1, mpopt);

%% Check violations of thermal limits
ind=I_max>0 & mpc.branch(:,11)==1; %indices of branches with nonzero limits
I_limits=repmat(I_max(ind),2,1);
I_dist_to_limit=I_limits-[sqrt(results.branch(ind,14).^2+...
    results.branch(ind,15).^2)./results.bus(ind_bus1(ind),8); ...
    sqrt(results.branch(ind,16).^2+...
    results.branch(ind,17).^2)./results.bus(ind_bus2(ind),8)];
%compute violations of thermal limits
Violations.values=abs(I_dist_to_limit(I_dist_to_limit<-1e-5));
Violations.percent=Violations.values*100./I_limits(I_dist_to_limit<-1e-5);
Violations.number=numel(Violations.values);


