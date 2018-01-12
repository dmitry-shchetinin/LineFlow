function [] = plot_constraints (Result, branch, flow_side, options, plot_options)
    %this function plots the nonlinear surface of I=I_max and its linear appoximation

    %unpack result structure
    A_matrix=Result.A;
    b_vector=Result.c;
    flag_result=Result.flag;
    
    %check if any plotting is needed
    if (plot_options.plot_surface_upper+plot_options.plot_surface_lower>0 || ...
        (plot_options.plot_planes_upper==1 && ~isempty(b_vector(A_matrix(:,3)==1))) || ...
        (plot_options.plot_planes_lower==1 && ~isempty(b_vector(A_matrix(:,3)==-1))))
        hold;
    else
        return;
    end
   
    %add default algorithm's options
    if (~isfield(options, 'delta_max_user'))
        options.delta_max_user=85;
    end
    
    %add some internal plotting options
    plot_options.margin_in_delta=0.0001; %in radian
    plot_options.coef_height=1; %keep at 1
    
    %% Preprocessing
    %check what kind of branch we have
    if (branch.V_i_min == branch.V_i_max && branch.V_j_min == branch.V_j_max)
        branch_type = 0;
    elseif (branch.V_i_min ~= branch.V_i_max && branch.V_j_min ~= branch.V_j_max)
        branch_type = 2;
    else
        branch_type = 1;
    end

    %turn angles into radians
    Kt_shift=branch.t_shift*pi/180;
    delta_max_user=options.delta_max_user*pi/180;
    
    %compute values used for plotting
    [branch_data_1, ratio, coef_Vij_1, feasible_1, offset_1]=compute_branch_parameters(branch, 1, options);
    [branch_data_2, ~, coef_Vij_2, feasible_2, offset_2]=compute_branch_parameters(branch, 2, options);
    
    %get min, max, and middle values of delta in the V box (for better plotting) as well as some info about feasible region
    if (flow_side==1)
        delta_limits=compute_delta_limits_for_plotting(branch, plot_options, branch_data_1, Kt_shift, feasible_1);
        offset=offset_1;
        coef_Vij_min=coef_Vij_1; 
        coef_Vij_max=coef_Vij_1;
    elseif (flow_side==2)
        delta_limits=compute_delta_limits_for_plotting(branch, plot_options, branch_data_2, Kt_shift, feasible_2);
        offset=offset_2;
        coef_Vij_min=coef_Vij_2; 
        coef_Vij_max=coef_Vij_2;
    else
        delta_temp1=compute_delta_limits_for_plotting(branch, plot_options, branch_data_1, Kt_shift, feasible_1);
        delta_temp2=compute_delta_limits_for_plotting(branch, plot_options, branch_data_2, Kt_shift, feasible_2);
        delta_limits=[min(delta_temp1(1),delta_temp2(1)), max(delta_temp1(2),delta_temp2(2))];
        if (coef_Vij_1>coef_Vij_2)
            offset.min=offset_2.min;
            offset.max=offset_1.max;
            coef_Vij_min=coef_Vij_2;
            coef_Vij_max=coef_Vij_1;
        else
            offset.min=offset_1.min;
            offset.max=offset_2.max;
            coef_Vij_min=coef_Vij_1;
            coef_Vij_max=coef_Vij_2;
        end
    end
    
    %extract information about upper and lower parts for plotting
    N_planes_lower1=0; N_planes_upper1=0; N_planes_lower2=0; N_planes_upper2=0;
    A_matrix_lower1=[]; A_matrix_upper1=[]; b_vector_lower1=[]; b_vector_upper1=[];
    A_matrix_lower2=[]; A_matrix_upper2=[]; b_vector_lower2=[]; b_vector_upper2=[];
    if (flag_result==2) %approximation was constructed
        if (flow_side==3 && sum(abs(diff(A_matrix(:,3))))>2)
            ind_switch_sides=find(diff(A_matrix(:,3))==2, 1, 'first');
            [A_matrix_lower1, A_matrix_upper1, b_vector_lower1, b_vector_upper1, N_planes_lower1, N_planes_upper1]=...
                extract_desired_approximation_part(A_matrix(1:ind_switch_sides,:), b_vector(1:ind_switch_sides));
            [A_matrix_lower2, A_matrix_upper2, b_vector_lower2, b_vector_upper2, N_planes_lower2, N_planes_upper2]=...
                extract_desired_approximation_part(A_matrix(ind_switch_sides+1:numel(b_vector),:), b_vector(ind_switch_sides+1:numel(b_vector)));
        else
            [A_matrix_lower1, A_matrix_upper1, b_vector_lower1, b_vector_upper1, N_planes_lower1, N_planes_upper1]=...
                extract_desired_approximation_part(A_matrix, b_vector);
        end
     end
    


    %% Plotting
    %plotting depends on the type of line
    if (branch_type==0) %line between two generators
        %plot max current line
        plot([-delta_max_user, delta_max_user]-Kt_shift,[branch.I_max;branch.I_max],'g', 'LineWidth',plot_options.linewidth); 
        
        %plot nonlinear function
        Vi=branch.V_i_min;
        Vj=branch.V_j_min;
        I_ij= @(delta_ij) sqrt(branch_data_1(1)*Vi.^2+branch_data_1(2)*Vj.^2+2*Vi.*Vj.*...
            (branch_data_1(3)*sin(delta_ij-Kt_shift)-branch_data_1(4)*cos(delta_ij-Kt_shift)))/ratio;
        I_ji= @(delta_ij) sqrt(branch_data_2(1)*Vi.^2+branch_data_2(2)*Vj.^2+2*Vi.*Vj.*...
            (branch_data_2(3)*sin(delta_ij-Kt_shift)-branch_data_2(4)*cos(delta_ij-Kt_shift)));
        
        if (plot_options.plot_surface_upper==1 || plot_options.plot_surface_lower==1)
            if (flow_side==1 || flow_side==3)
                fplot(I_ij,[-delta_max_user, delta_max_user]-Kt_shift,'b', 'LineWidth',plot_options.linewidth);
            end
            if (flow_side==3 || flow_side==3)
                fplot(I_ji,[-delta_max_user, delta_max_user]-Kt_shift,'b', 'LineWidth',plot_options.linewidth);
            end
        end

        %plot linear approximation
        if (plot_options.plot_planes_lower==1 && N_planes_lower1>0)
            plot(-[b_vector_lower1(1); b_vector_lower1(1)],[0, branch.I_max],'r', 'LineWidth',plot_options.linewidth);
        end
        if (plot_options.plot_planes_upper==1 && N_planes_upper1>0)
            plot([b_vector_upper1(1); b_vector_upper1(1)],[0, branch.I_max],'r', 'LineWidth',plot_options.linewidth);
        end
        if (plot_options.plot_planes_lower==1 && N_planes_lower2>0)
            plot(-[b_vector_lower2(1); b_vector_lower2(1)],[0, branch.I_max],'r', 'LineWidth',plot_options.linewidth);
        end
        if (plot_options.plot_planes_upper==1 && N_planes_upper2>0)
            plot([b_vector_upper2(1); b_vector_upper2(1)],[0, branch.I_max],'r', 'LineWidth',plot_options.linewidth);
        end
        title('I_{ij}=f(\theta_{ij}) and its linear approximation.');
        xlabel('\theta_{ij} (rad)');
        ylabel('I_{ij} (p.u.)');
        xlim([-delta_max_user-Kt_shift-2*delta_max_user*plot_options.margin_V_axis, delta_max_user-Kt_shift+2*delta_max_user*plot_options.margin_V_axis]);
        
        %fill out legend
        if (plot_options.plot_surface_upper==1 || plot_options.plot_surface_lower==1)
            if (plot_options.plot_planes_upper==1 || plot_options.plot_planes_lower==1)
                if (flow_side==3)
                    legend('I_{ij}^{max}','I_{ij}','I_{ji}','Approximation');
                else
                    legend('I_{ij}^{max}','I_{ij}','Approximation');
                end
            else
                legend('I_{ij}^{max}','I_{ij}');
            end
        elseif (plot_options.plot_planes_upper==1 || plot_options.plot_planes_lower==1)
            legend('I_{ij}^{max}','Approximation');
        end
    elseif (branch_type==1) %line between generator and load
        if (branch.V_i_min==branch.V_i_max)
            Vi=branch.V_i_min;
            if (flow_side==1)
                V_limits=[coef_Vij_1*Vi+offset_1.min; coef_Vij_1*Vi+offset_1.max];
            elseif (flow_side==2)
                V_limits=[coef_Vij_2*Vi+offset_2.min; coef_Vij_2*Vi+offset_2.max];
            else
                V_limits=[min(coef_Vij_1*Vi+offset_1.min, coef_Vij_2*Vi+offset_2.min); max(coef_Vij_1*Vi+offset_1.max, coef_Vij_2*Vi+offset_2.max)];
            end
            f_surface1= @(Vj, delta_ij) sqrt(branch_data_1(1)*Vi.^2+branch_data_1(2)*Vj.^2+2*Vi.*Vj.*...
                (branch_data_1(3)*sin(delta_ij-Kt_shift)-branch_data_1(4)*cos(delta_ij-Kt_shift)))/ratio-branch.I_max;
            f_surface2= @(Vj, delta_ij) sqrt(branch_data_2(1)*Vi.^2+branch_data_2(2)*Vj.^2+2*Vi.*Vj.*...
                (branch_data_2(3)*sin(delta_ij-Kt_shift)-branch_data_2(4)*cos(delta_ij-Kt_shift)))-branch.I_max;
            title('\theta_{ij}=f(V_j) and its linear approximation.');
                xlabel('V_j (p.u.)');
        else
            Vj=branch.V_j_min;
            if (flow_side==1)
                V_limits=[(Vj-offset_1.max)/coef_Vij_1; (Vj-offset_1.min)/coef_Vij_1];
            elseif (flow_side==2)
                V_limits=[(Vj-offset_2.max)/coef_Vij_2; (Vj-offset_2.min)/coef_Vij_2];
            else
                V_limits=[min((Vj-offset_1.max)/coef_Vij_1, (Vj-offset_2.max)/coef_Vij_2); max((Vj-offset_1.min)/coef_Vij_1, (Vj-offset_2.min)/coef_Vij_2)];
            end
            f_surface1= @(Vi, delta_ij) sqrt(branch_data_1(1)*Vi.^2+branch_data_1(2)*Vj.^2+2*Vi.*Vj.*...
                (branch_data_1(3)*sin(delta_ij-Kt_shift)-branch_data_1(4)*cos(delta_ij-Kt_shift)))/ratio-branch.I_max;
            f_surface2= @(Vi, delta_ij) sqrt(branch_data_2(1)*Vi.^2+branch_data_2(2)*Vj.^2+2*Vi.*Vj.*...
                (branch_data_2(3)*sin(delta_ij-Kt_shift)-branch_data_2(4)*cos(delta_ij-Kt_shift)))-branch.I_max;
            title('\theta_{ij}=f(V_i) and its linear approximation.');
                xlabel('V_i (p.u.)');
        end
        ylabel('\theta_{ij} (rad)');
        
        %plot nonlinear function
        if (plot_options.plot_surface_upper==1 || plot_options.plot_surface_lower==1)
            if (flow_side==1 || flow_side==3)
                fimplicit(f_surface1,[V_limits(1)*0.9999, V_limits(2)*1.0001, delta_limits],'b', 'LineWidth',plot_options.linewidth);
            end
            if (flow_side==2 || flow_side==3)
                fimplicit(f_surface2,[V_limits(1)*0.9999, V_limits(2)*1.0001, delta_limits],'b', 'LineWidth',plot_options.linewidth);
            end
        else
            V_limits=[];
            delta_limits=[];
        end
        
        %plot linear approximation
        if (plot_options.plot_planes_upper==1 && N_planes_upper1>0)
            Vertices_upper1=compute_2D_polyhedron_vertices(A_matrix_upper1, b_vector_upper1, ...
                N_planes_upper1, A_matrix_lower1, b_vector_lower1, N_planes_lower1, branch, feasible_1);
            plot(Vertices_upper1.V, Vertices_upper1.delta,'r', 'LineWidth',plot_options.linewidth);
        else
            Vertices_upper1=struct('V',[],'delta',[]);
        end
        if (plot_options.plot_planes_lower==1 && N_planes_lower1>0)
            Vertices_lower1=compute_2D_polyhedron_vertices(A_matrix_lower1, b_vector_lower1, ...
                N_planes_lower1, A_matrix_upper1, b_vector_upper1, N_planes_upper1, branch, feasible_1);
            plot(Vertices_lower1.V, Vertices_lower1.delta,'r', 'LineWidth',plot_options.linewidth);
        else
            Vertices_lower1=struct('V',[],'delta',[]);
        end
        if (plot_options.plot_planes_upper==1 && N_planes_upper2>0)
            Vertices_upper2=compute_2D_polyhedron_vertices(A_matrix_upper2, b_vector_upper2, ...
                N_planes_upper2, A_matrix_lower2, b_vector_lower2, N_planes_lower2, branch, feasible_2);
            plot(Vertices_upper2.V, Vertices_upper2.delta,'r', 'LineWidth',plot_options.linewidth);
        else
            Vertices_upper2=struct('V',[],'delta',[]);
        end
        if (plot_options.plot_planes_lower==1 && N_planes_lower2>0)
            Vertices_lower2=compute_2D_polyhedron_vertices(A_matrix_lower2, b_vector_lower2, ...
                N_planes_lower2, A_matrix_upper2, b_vector_upper2, N_planes_upper2, branch, feasible_2);
            plot(Vertices_lower2.V, Vertices_lower2.delta,'r', 'LineWidth',plot_options.linewidth);
        else
            Vertices_lower2=struct('V',[],'delta',[]);
        end
        
        %choose axes limits
        V_limits_plot=[min([V_limits; Vertices_lower1.V; Vertices_upper1.V; Vertices_lower2.V; Vertices_upper2.V]), ...
            max([V_limits; Vertices_lower1.V; Vertices_upper1.V; Vertices_lower2.V; Vertices_upper2.V])];
        delta_limits_plot=[min([delta_limits'; Vertices_lower1.delta; Vertices_upper1.delta; Vertices_lower2.delta; Vertices_upper2.delta]), ...
            max([delta_limits'; Vertices_lower1.delta; Vertices_upper1.delta; Vertices_lower2.delta; Vertices_upper2.delta])];
        xlim([V_limits_plot(1)-diff(V_limits_plot)*plot_options.margin_V_axis, V_limits_plot(2)+diff(V_limits_plot)*plot_options.margin_V_axis]);
        ylim([delta_limits_plot(1)-diff(delta_limits_plot)*plot_options.margin_delta_axis, delta_limits_plot(2)+diff(delta_limits_plot)*plot_options.margin_delta_axis]);
        
        %fill out legend
        if (plot_options.plot_surface_upper==1 || plot_options.plot_surface_lower==1)
            if (plot_options.plot_planes_upper==1 || plot_options.plot_planes_lower==1)
                if (flow_side==3)
                    legend('Boundary1', 'Boundary2', 'Approximation');
                else
                    legend('Boundary', 'Approximation');
                end
            else
                legend('Boundary');
            end
        else
            legend('Approximation');
        end
    else %line between two loads
        f_surface1= @(Vi, Vj, delta_ij) sqrt(branch_data_1(1)*Vi.^2+branch_data_1(2)*Vj.^2+2*Vi.*Vj.*...
            (branch_data_1(3)*sin(delta_ij-Kt_shift)-branch_data_1(4)*cos(delta_ij-Kt_shift)))/ratio-branch.I_max;
        f_surface2= @(Vi, Vj, delta_ij) sqrt(branch_data_2(1)*Vi.^2+branch_data_2(2)*Vj.^2+2*Vi.*Vj.*...
            (branch_data_2(3)*sin(delta_ij-Kt_shift)-branch_data_2(4)*cos(delta_ij-Kt_shift)))-branch.I_max;
        
        %plot nonlinear function
        if (plot_options.plot_surface_upper==1 || plot_options.plot_surface_lower==1)
            %compute V coordinates of end points of outermost intersection lines
            Line_points1=compute_intersections_of_line_with_Vbox(coef_Vij_min, offset.min, branch);
            Line_points2=compute_intersections_of_line_with_Vbox(coef_Vij_max, offset.max, branch);
            V_i_limits=[min(Line_points1.V_i_begin, Line_points2.V_i_begin), max(Line_points1.V_i_end, Line_points2.V_i_end)];
            V_j_limits=[min(Line_points1.V_j_begin, Line_points2.V_j_begin), max(Line_points1.V_j_end, Line_points2.V_j_end)];
            if (flow_side==1 || flow_side==3)
                fimplicit3(f_surface1, [V_i_limits, V_j_limits, delta_limits]);
            end
            if (flow_side==2 || flow_side==3)
                fimplicit3(f_surface2, [V_i_limits, V_j_limits, delta_limits]);
            end
            Limits_surface=struct('V_i', V_i_limits, 'V_j', V_j_limits, 'delta', delta_limits);
        else
            Limits_surface=struct('V_i', [], 'V_j', [], 'delta', []);
        end
        
        %plot linear approximation
        if (plot_options.plot_planes_upper==1 && N_planes_upper1>0)
            Vertices_upper1=compute_3D_polyhedron_vertices(A_matrix_upper1, b_vector_upper1, ...
                N_planes_upper1, A_matrix_lower1, b_vector_lower1, N_planes_lower1, branch, feasible_1, options);
            Limits_upper1=compute_min_and_max_vertices(Vertices_upper1, N_planes_upper1);
            for i=1:N_planes_upper1
                fill3(Vertices_upper1.V_i{i}, Vertices_upper1.V_j{i}, Vertices_upper1.delta{i}*plot_options.coef_height, Vertices_upper1.delta{i});
            end
        else
            Limits_upper1=struct('V_i', [], 'V_j', [], 'delta', []);
        end
        if (plot_options.plot_planes_lower==1 && N_planes_lower1>0)
            Vertices_lower1=compute_3D_polyhedron_vertices(A_matrix_lower1, b_vector_lower1, ...
                N_planes_lower1, A_matrix_upper1, b_vector_upper1, N_planes_upper1, branch, feasible_1, options);
            Limits_lower1=compute_min_and_max_vertices(Vertices_lower1, N_planes_lower1);
            for i=1:N_planes_lower1
                fill3(Vertices_lower1.V_i{i}, Vertices_lower1.V_j{i}, Vertices_lower1.delta{i}*plot_options.coef_height, Vertices_lower1.delta{i});
            end
        else
            Limits_lower1=struct('V_i', [], 'V_j', [], 'delta', []);
        end
        if (plot_options.plot_planes_upper==1 && N_planes_upper2>0)
            Vertices_upper2=compute_3D_polyhedron_vertices(A_matrix_upper2, b_vector_upper2, ...
                N_planes_upper2, A_matrix_lower2, b_vector_lower2, N_planes_lower2, branch, feasible_2, options);
            Limits_upper2=compute_min_and_max_vertices(Vertices_upper2, N_planes_upper2);
            for i=1:N_planes_upper2
                fill3(Vertices_upper2.V_i{i}, Vertices_upper2.V_j{i}, Vertices_upper2.delta{i}*plot_options.coef_height, Vertices_upper2.delta{i});
            end
        else
            Limits_upper2=struct('V_i', [], 'V_j', [], 'delta', []);
        end
        if (plot_options.plot_planes_lower==1 && N_planes_lower2>0)
            Vertices_lower2=compute_3D_polyhedron_vertices(A_matrix_lower2, b_vector_lower2, ...
                N_planes_lower2, A_matrix_upper2, b_vector_upper2, N_planes_upper2, branch, feasible_2, options);
            Limits_lower2=compute_min_and_max_vertices(Vertices_lower2, N_planes_lower2);
            for i=1:N_planes_lower1
                fill3(Vertices_lower2.V_i{i}, Vertices_lower2.V_j{i}, Vertices_lower2.delta{i}*plot_options.coef_height, Vertices_lower2.delta{i});
            end
        else
            Limits_lower2=struct('V_i', [], 'V_j', [], 'delta', []);
        end
        
        title('\theta_{ij}=f(V_i,V_j) and its linear approximation.');
        xlabel('V_i (p.u.)');
        ylabel('V_j (p.u.)');
        zlabel('\theta_{ij} (rad)');
        
        %choose axes limits
        V_i_limits_plot=[min([Limits_surface.V_i, Limits_lower1.V_i, Limits_upper1.V_i, Limits_lower2.V_i, Limits_upper2.V_i]), ...
            max([Limits_surface.V_i, Limits_lower1.V_i, Limits_upper1.V_i, Limits_lower2.V_i, Limits_upper2.V_i])];
        V_j_limits_plot=[min([Limits_surface.V_j, Limits_lower1.V_j, Limits_upper1.V_j, Limits_lower2.V_j, Limits_upper2.V_j]), ...
            max([Limits_surface.V_j, Limits_lower1.V_j, Limits_upper1.V_j, Limits_lower2.V_j, Limits_upper2.V_j])];
        delta_limits_plot=[min([Limits_surface.delta, Limits_lower1.delta, Limits_upper1.delta, Limits_lower2.delta, Limits_upper2.delta]), ...
            max([Limits_surface.delta, Limits_lower1.delta, Limits_upper1.delta, Limits_lower2.delta, Limits_upper2.delta])];
        xlim([V_i_limits_plot(1)-diff(V_i_limits_plot)*plot_options.margin_V_axis, V_i_limits_plot(2)+diff(V_i_limits_plot)*plot_options.margin_V_axis]);
        ylim([V_j_limits_plot(1)-diff(V_j_limits_plot)*plot_options.margin_V_axis, V_j_limits_plot(2)+diff(V_j_limits_plot)*plot_options.margin_V_axis]);
        zlim([delta_limits_plot(1)-diff(delta_limits_plot)*plot_options.margin_delta_axis, delta_limits_plot(2)+diff(delta_limits_plot)*plot_options.margin_delta_axis]);
    end
    
    %set fontsize
    set(gca,'fontsize',plot_options.fontsize);
end



%% Supporting Functions
%compute vertices of upper or lower part of approximation for line between two loads
function Vertices=compute_3D_polyhedron_vertices(A_this, b_this, N_this, A_opposite, b_opposite, N_opposite, branch, feasible, options)
    %initialize output
    Vertices.V_i=cell(N_this,1);
    Vertices.V_j=cell(N_this,1);
    Vertices.delta=cell(N_this,1);
    
    %compute intersections of intersection lines with Vbox
    %last line
    if (N_opposite==0 || feasible.top_left==1)
        Line_points(N_this+1)=struct('V_i_begin', branch.V_i_min, 'V_j_begin', branch.V_j_max, 'V_i_end', branch.V_i_min, 'V_j_end', branch.V_j_max);
    else
        Line_parameters=parameters_of_intersection_line(A_this(end,:), b_this(end), A_opposite(end,:), b_opposite(end));
        Line_points(N_this+1)=compute_intersections_of_line_with_Vbox(Line_parameters.slope, Line_parameters.offset, branch);
        if (Line_points(N_this+1).V_i_begin<branch.V_i_min*0.9999 || Line_points(N_this+1).V_i_begin>branch.V_i_max*1.0001 || ...
                Line_points(N_this+1).V_j_begin<branch.V_j_min*0.9999 || Line_points(N_this+1).V_j_begin>branch.V_j_max*1.0001 || ...
                Line_points(N_this+1).V_i_end<branch.V_i_min*0.9999 || Line_points(N_this+1).V_i_end>branch.V_i_max*1.0001 || ...
                Line_points(N_this+1).V_j_end<branch.V_j_min*0.9999 || Line_points(N_this+1).V_j_end>branch.V_j_max*1.0001)
            Line_points(N_this+1)=struct('V_i_begin', branch.V_i_min, 'V_j_begin', branch.V_j_max, 'V_i_end', branch.V_i_min, 'V_j_end', branch.V_j_max);
        end
    end

    %first line
    if (N_opposite==0 || feasible.bottom_right==1)
        Line_points(1)=struct('V_i_begin', branch.V_i_max, 'V_j_begin', branch.V_j_min, 'V_i_end', branch.V_i_max, 'V_j_end', branch.V_j_min);
    else
        Line_parameters=parameters_of_intersection_line(A_this(1,:), b_this(1), A_opposite(1,:), b_opposite(1));
        Line_points(1)=compute_intersections_of_line_with_Vbox(Line_parameters.slope, Line_parameters.offset, branch);
        if (Line_points(1).V_i_begin<branch.V_i_min*0.9999 || Line_points(1).V_i_begin>branch.V_i_max*1.0001 || ...
                Line_points(1).V_j_begin<branch.V_j_min*0.9999 || Line_points(1).V_j_begin>branch.V_j_max*1.0001 || ...
                Line_points(1).V_i_end<branch.V_i_min*0.9999 || Line_points(1).V_i_end>branch.V_i_max*1.0001 || ...
                Line_points(1).V_j_end<branch.V_j_min*0.9999 || Line_points(1).V_j_end>branch.V_j_max*1.0001)
            Line_points(1)=struct('V_i_begin', branch.V_i_max, 'V_j_begin', branch.V_j_min, 'V_i_end', branch.V_i_max, 'V_j_end', branch.V_j_min);
        end
    end
    
    %lines from 2nd to 2nd to last
    for i=2:N_this
        Line_parameters=parameters_of_intersection_line(A_this(i-1,:), b_this(i-1), A_this(i,:), b_this(i));
        Line_points(i)=compute_intersections_of_line_with_Vbox(Line_parameters.slope, Line_parameters.offset, branch);
    end
    
    %override in case of conservative approximation with one plane and both corners infeasible
    if (options.approximation==0 && N_this==1 && feasible.bottom_right==0 && feasible.top_left==0)
        Line_points(1)=struct('V_i_begin', branch.V_i_max, 'V_j_begin', branch.V_j_min, 'V_i_end', branch.V_i_max, 'V_j_end', branch.V_j_min);
        Line_points(N_this+1)=struct('V_i_begin', branch.V_i_min, 'V_j_begin', branch.V_j_max, 'V_i_end', branch.V_i_min, 'V_j_end', branch.V_j_max);
    end
    
    %compute vertices of planes (add vertices in [V_i_min, V_j_min], [V_i_max, V_j_max] corners)
    for i=1:N_this
        corner_min=[];
        corner_max=[];
        if (Line_points(i).V_j_begin==branch.V_j_min && Line_points(i+1).V_j_begin>branch.V_j_min) %add vertex in (V_i_min, V_j_min) corner
            corner_min=1;
        end
        if (Line_points(i).V_j_end<branch.V_j_max && Line_points(i+1).V_j_end==branch.V_j_max) %add vertex in (V_i_max, V_j_max) corner
            corner_max=1;
        end
        Vertices.V_i{i}=[Line_points(i).V_i_begin; branch.V_i_min*corner_min; Line_points(i+1).V_i_begin; Line_points(i+1).V_i_end; branch.V_i_max*corner_max; Line_points(i).V_i_end];
        Vertices.V_j{i}=[Line_points(i).V_j_begin; branch.V_j_min*corner_min; Line_points(i+1).V_j_begin; Line_points(i+1).V_j_end; branch.V_j_max*corner_max; Line_points(i).V_j_end];
        Vertices.delta{i}=(b_this(i)-A_this(i,1)*Vertices.V_i{i}-A_this(i,2)*Vertices.V_j{i})/A_this(i,3);
    end
end


%compute min and max values of vertices for Vi, Vj, delta for all planes
function Limits=compute_min_and_max_vertices(Vertices, N_planes)
    Limits=struct('V_i', [inf, -inf], 'V_j', [inf, -inf], 'delta', [inf, -inf]);
    for i=1:N_planes
        Limits.V_i=update_min_and_max_values(Limits.V_i, Vertices.V_i{i});
        Limits.V_j=update_min_and_max_values(Limits.V_j, Vertices.V_j{i});
        Limits.delta=update_min_and_max_values(Limits.delta, Vertices.delta{i});
    end
end


%update recorded min and max values given the new array
function Values=update_min_and_max_values(Values, Vector)
    if (Values(1)>min(Vector))
        Values(1)=min(Vector);
    end
    if (Values(2)<max(Vector))
        Values(2)=max(Vector);
    end
end

%compute offset and slope of projection of intersection line in (Vi, Vj) plane
function Line_parameters=parameters_of_intersection_line(Normal_plane_1, Offset_plane_1, Normal_plane_2, Offset_plane_2)
    %compute value of V_i when V_j=0
    V_i=intersection_of_2D_lines([Normal_plane_1(1); Normal_plane_2(1)], [Normal_plane_1(3); Normal_plane_2(3)], [Offset_plane_1; Offset_plane_2]);
    if (V_i==0) %offset is zero
        %compute value of V_j when V_i=1
        V_j=intersection_of_2D_lines([Normal_plane_1(2); Normal_plane_2(2)], [Normal_plane_1(3); Normal_plane_2(3)], [Offset_plane_1-Normal_plane_1(1); Offset_plane_2-Normal_plane_2(1)]);
        %the offset and slope are then given as follows:
        Line_parameters.offset=0;
        Line_parameters.slope=V_j;
    else
        %compute value of V_j when V_i=0
        V_j=intersection_of_2D_lines([Normal_plane_1(2); Normal_plane_2(2)], [Normal_plane_1(3); Normal_plane_2(3)], [Offset_plane_1; Offset_plane_2]);
        %the offset and slope are then given as follows:
        Line_parameters.offset=V_j;
        Line_parameters.slope=-V_j/V_i;
    end
end


%compute vertices of upper or lower part of approximation for line between load and generator
function Vertices=compute_2D_polyhedron_vertices(A_this, b_this, N_this, A_opposite, b_opposite, N_opposite, branch, feasible)
    %initialize output
    Vertices.V=zeros(N_this+1,1);
    Vertices.delta=zeros(N_this+1,1);
    
    %extract relevant parameters for a variable voltage
    if (branch.V_i_min==branch.V_i_max)
        V_min=branch.V_j_min;
        V_max=branch.V_j_max;
        normals_V=A_this(:,2);
        if (N_opposite>0)
            normals_V_opposite=A_opposite(:,2);
        end
        feas_V_min=feasible.bottom_right;
        feas_V_max=feasible.top_left;
    else
        %flip order to make sure we'll start plotting with V_min
        if (N_this>1)
            A_this=flip(A_this); b_this=flip(b_this); 
        end
        if (N_opposite>1)
            A_opposite=flip(A_opposite); b_opposite=flip(b_opposite);
        end
        V_min=branch.V_i_min;
        V_max=branch.V_i_max;
        normals_V=A_this(:,1);
        if (N_opposite>0)
            normals_V_opposite=A_opposite(:,1);
        end
        feas_V_min=feasible.top_left;
        feas_V_max=feasible.bottom_right;
    end
    
    %compute the first vertex
    Vertices.V(1)=V_min; %this is what we consider by default
    if (N_opposite>0 && feas_V_min==0)
        V_temp=intersection_of_2D_lines([normals_V(1); normals_V_opposite(1)], [A_this(1,3); A_opposite(1,3)], [b_this(1); b_opposite(1)]);
        if (V_temp>V_min && V_temp<V_max)
            Vertices.V(1)=V_temp;
        end
    end
    Vertices.delta(1)=(b_this(1)-Vertices.V(1)*normals_V(1))/A_this(1,3);
    
    %compute all vertices from 2nd to 2nd to last
    for i=2:N_this
        Vertices.V(i)=intersection_of_2D_lines(normals_V(i-1:i), A_this(i-1:i,3), b_this(i-1:i));
        Vertices.delta(i)=(b_this(i)-Vertices.V(i)*normals_V(i))/A_this(i,3);
    end
    
    %compute the last vertex
    Vertices.V(end)=V_max; %this is what we consider by default
    if (N_opposite>0 && feas_V_max==0)
        V_temp=intersection_of_2D_lines([normals_V(end); normals_V_opposite(end)], [A_this(end,3); A_opposite(end,3)], [b_this(end); b_opposite(end)]);
        if (V_temp>V_min && V_temp<V_max)
            Vertices.V(end)=V_temp;
        end
    end
    Vertices.delta(end)=(b_this(end)-Vertices.V(end)*normals_V(end))/A_this(end,3);
end


function x1=intersection_of_2D_lines(coef_x1, coef_x2, coef_b)
    x1=(coef_b(1)*coef_x2(2)-coef_b(2)*coef_x2(1))/(coef_x1(1)*coef_x2(2)-coef_x1(2)*coef_x2(1));
end


function Line_points=compute_intersections_of_line_with_Vbox(line_slope, line_offset, Vbox)
    %beginning of the line segment
    if (line_slope * Vbox.V_i_min + line_offset >= Vbox.V_j_min)
        Line_points.V_i_begin = Vbox.V_i_min;
        Line_points.V_j_begin = line_slope * Vbox.V_i_min + line_offset;
    else
        Line_points.V_i_begin = (Vbox.V_j_min - line_offset) / line_slope;
        Line_points.V_j_begin = Vbox.V_j_min;
    end
    %end of the line segment
    if (line_slope * Vbox.V_i_max + line_offset <= Vbox.V_j_max)
        Line_points.V_i_end = Vbox.V_i_max;
        Line_points.V_j_end = line_slope * Vbox.V_i_max + line_offset;
    else
        Line_points.V_i_end = (Vbox.V_j_max - line_offset) / line_slope;
        Line_points.V_j_end = Vbox.V_j_max;
    end
end


%compute limits on delta for plotting
function delta_limits=compute_delta_limits_for_plotting(branch, plot_options, branch_data, Kt_shift, feasible)
    N_V=plot_options.N_samples_V;
    N_all=N_V^2;
    [V_i, V_j]=create_uniformly_distributed_points(branch, N_V);
    V_i_all=reshape(repmat(V_i, 1, N_V)',N_all,1);
    V_j_all=repmat(V_j, N_V, 1);
    [delta_nonlin_lower, delta_nonlin_upper]=compute_surface_angles(V_i_all, V_j_all, branch_data, N_all);
    delta_nonlin_lower(delta_nonlin_lower==5)=[];
    delta_nonlin_upper(delta_nonlin_upper==5)=[];
    deltas.delta_min_lower=min(delta_nonlin_lower)+Kt_shift-plot_options.margin_in_delta;
    deltas.delta_max_upper=max(delta_nonlin_upper)+Kt_shift+plot_options.margin_in_delta;
    deltas.delta_min_upper=min(delta_nonlin_upper)+Kt_shift;
    deltas.delta_max_lower=max(delta_nonlin_lower)+Kt_shift;
    deltas.delta_middle=mean(delta_nonlin_lower+delta_nonlin_upper)/2+Kt_shift;
    
    if (plot_options.plot_surface_upper==1 && plot_options.plot_surface_lower==1)
        delta_limits=[deltas.delta_min_lower, deltas.delta_max_upper];
    elseif (plot_options.plot_surface_upper==1)
        if (feasible.bottom_right==1 && feasible.top_left==1)
            delta_limits=[max(deltas.delta_min_upper, deltas.delta_middle), deltas.delta_max_upper];
        else
            delta_limits=[deltas.delta_middle, deltas.delta_max_upper];
        end
    else
        if (feasible.bottom_right==1 && feasible.top_left==1)
            delta_limits=[deltas.delta_min_lower, min(deltas.delta_max_lower, deltas.delta_middle)];
        else
            delta_limits=[deltas.delta_min_lower, deltas.delta_middle];
        end
    end
end


function [V_i, V_j]=create_uniformly_distributed_points(branch, N_steps_V)
    step_V_i=(branch.V_i_max-branch.V_i_min)/(N_steps_V-1);
    step_V_j=(branch.V_j_max-branch.V_j_min)/(N_steps_V-1);
    %get vector of data points
    if (branch.V_i_min==branch.V_i_max)
        V_i=ones(N_steps_V,1)*branch.V_i_min;
    else
        V_i=(branch.V_i_min:step_V_i:branch.V_i_max)';
    end
    if (branch.V_j_min==branch.V_j_max)
        V_j=ones(N_steps_V,1)*branch.V_j_min;
    else
        V_j=(branch.V_j_min:step_V_j:branch.V_j_max)';
    end
end


function [branch_data, ratio, coef_Vij, feasible, offset]=compute_branch_parameters(branch, flow_side, options)
    %take transformer into account
    I_max_user = branch.I_max;
    if (branch.t_ratio == 0.0)
        b_sh = branch.b_sh / 2.0;
        ratio = 1.0;
    else
        %compute parameters depending on the transformer's model
        if (options.tr_model_type == 0) %if Matpower convention is used
            b_sh = branch.b_sh / 2.0;
        else %if Russian convention is used
            if (flow_side == 1)
                b_sh = branch.b_sh;
            else
                b_sh = 0.0;
            end
        end
        ratio = branch.t_ratio;
        if (flow_side == 1)
            I_max_user = I_max_user*ratio;
        end
    end
    buf=2*(branch.g^2+(branch.b+b_sh)^2)*(branch.g^2+branch.b^2)/ratio^2;
    if (flow_side==1)
        coef_Vij=sqrt(((branch.g^2+(branch.b+b_sh)^2)/ratio^2)/(branch.g^2+branch.b^2));
        branch_data=[(branch.g^2+(branch.b+b_sh)^2)/ratio^2; branch.g^2+branch.b^2; branch.g*b_sh/ratio; ...
            (branch.g^2+branch.b^2+branch.b*b_sh)/ratio; (I_max_user)^2; 2*buf; ...
            -(branch.g*b_sh/ratio)/buf; abs((branch.g^2+branch.b^2+branch.b*b_sh)/ratio)/buf];
    else
        coef_Vij=sqrt(((branch.g^2+branch.b^2)/ratio^2)/(branch.g^2+(branch.b+b_sh)^2));
        branch_data=[(branch.g^2+branch.b^2)/ratio^2; branch.g^2+(branch.b+b_sh)^2; -branch.g*b_sh/ratio; ...
            (branch.g^2+branch.b^2+branch.b*b_sh)/ratio; (I_max_user)^2; 2*buf; ...
            (branch.g*b_sh/ratio)/buf; abs((branch.g^2+branch.b^2+branch.b*b_sh)/ratio)/buf];
    end
    
    temp_offset = sqrt(branch_data(5) / branch_data(2)); %abs value of the offset of boundary lines
    if (branch.V_j_max - coef_Vij * branch.V_i_min > temp_offset)
        offset.max = temp_offset;
        feasible.top_left=0;
    else
        offset.max = branch.V_j_max - coef_Vij * branch.V_i_min;
        feasible.top_left=1;
    end
    if (branch.V_j_min - coef_Vij * branch.V_i_max < -temp_offset)
        offset.min = -temp_offset;
        feasible.bottom_right=0;
    else
        offset.min = branch.V_j_min - coef_Vij * branch.V_i_max;
        feasible.bottom_right=1;
    end
end


function [A_lower, A_upper, b_lower, b_upper, N_lower, N_upper]=extract_desired_approximation_part(A_matrix, b_vector)
    A_lower=A_matrix(A_matrix(:,3)==-1,:);
    A_upper=A_matrix(A_matrix(:,3)==1,:);
    b_lower=b_vector(A_matrix(:,3)==-1);
    b_upper=b_vector(A_matrix(:,3)==1);
    if (~isempty(b_upper))
        N_upper=length(b_upper);
    end
    if (~isempty(b_lower))
        N_lower=length(b_lower);
    end
end





