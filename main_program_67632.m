clear all; close all;



%% ******************************************************************************
% Course: Non Linear Finite Element Method, Summer 2023
% Author: Brijraj Jadeja
% ---------------  FEM Analysis of a elastic-plastic plate with-------------
% ---------------- a hole under uniform farfield tractions -----------------
% Supervisor: Dr.-Ing. Stefan Prüger
% ******************************************************************************
% main routine for nonlinear finite element program




%% ******************************************************************************
% material properties
Emod=240;         % Young's modulus in GPa
nu=0.2;          % Poisson's ratio
sigmaY0=0.35; %240 for linear;      % initial yield stress in GPa
hard_iso=0.7;     % linear, isotropic hardening modulus in GPa
deltaY=0;      % increase in yield stress stemming from exponential hardening in GPa
h_sh=15;          % shape parameter for exponential hardening| rate at which the yield stress increases with plastic deformation
hard_kin=0.0;     % linear kinematic hardening modulus in GPa
hard_funct='lin'; % specifies how the yield stress evolves as plastic deformation progresses.
%
props=[Emod,nu,sigmaY0,hard_iso,deltaY,h_sh,hard_kin];
n_sdv_p_GP=7;
%
visualize_mesh_data=true;
%
% input file
foi='Mesh_Tri_Coarse.mat'; % 'Mesh_Quad_Coarse.mat' 'Mesh_Tri_uniform.mat'
elem_type='Tri';           % 'Quad', 'Tri'
%
% read mesh information
load(foi);
[n_nodes,n_dim]=size(node_list); n_dim=n_dim-1;
[n_elem,n_node_p_elem]=size(elem_list); n_node_p_elem=n_node_p_elem-1;
%
% visualization of initial mesh data
if visualize_mesh_data
  figID=1;
  visualize_mesh (node_list,elem_list,n_elem,elem_type,figID);
end
%
% -------------------------------------------------------------------------
% information how work with element sets and node sets and to call material
% routine
% -------------------------------------------------------------------------
%
% name of element sets
for i=1:4
    disp(['element set: ',elem_sets{i,1}])
end
% name of node sets
for i=1:4
    disp(['node set: ',node_sets{i,1}])
end



sdv = zeros(1, n_sdv_p_GP);
tau = 0;
sigmainf = 0.335;
no_of_t_inc = 10; %10;
u_global = zeros(2 * n_nodes, 1);
bk_u_global = zeros(2 * n_nodes, 1); % backup global displacement
u_e = zeros(2*n_node_p_elem,1);
d_tau = sigmainf/no_of_t_inc;
tau_store = [];
store_n50 = [];
tau = tau+d_tau;
     sdv_values = cell((size(elem_list, 1)), 6);
     for o = 1:6
         for o_e = 1:size(elem_list, 1)
             sdv_values{o_e, o} = sdv;
         end
     end

     n_50 = (node_list(:, 2) == 50) & (node_list(:, 3) == 0);
     no_50 = node_list(n_50, 1);
     
  

            ray_0 = [];
            ray_45 = [];
            ray_90 = [];
            sigma_analytic_0 = [];
            sigma_analytic_45 = [];
            sigma_analytic_90 = [];
            stored_stress = [];  % ensures when load scaling starts the values are zero
            stored_p_stress = [];
            stored_sigma_0 = [];
            stored_sigma_45 = [];
            stored_sigma_90 = [];
while tau <= sigmainf     %load scaling while loop
    
    
    tau ;
    del_ug = zeros(2*n_nodes,1);
    bk_u_global = u_global; 
    prev_stored_stress = stored_stress;   %if code diverges we can start from the previous values.

    k=1;
    while true          %Newton raphson loop
    
        sigmainf;
        k;

        K_global = zeros(2*n_nodes, 2*n_nodes);
        Res_global = zeros(2*n_nodes, 1);

        stored_stress = [];  %for each iteration storing stress values starting from zero at first increment

        
        for elem = 1:n_elem


            [K_global, Res_global, ray_0, ray_45, ray_90, stored_stress, stored_p_stress, sdv_values] = elementroutine(elem, elem_list,elem_sets, node_list, node_sets, n_node_p_elem, n_nodes,  u_global, sdv_values, sdv, props ,hard_funct, stored_stress, stored_p_stress, ray_0, ray_45, ray_90, sigma_analytic_0, sigma_analytic_45, sigma_analytic_90, tau, K_global, Res_global, sigmainf);

            

       
        end

       
        norm(Res_global);

        first_set = node_sets(1,2);
        fourth_set = node_sets(4,2);

        % Convert the cell array elements to numeric arrays
        first_set = cell2mat(first_set);
        fourth_set = cell2mat(fourth_set);

        [K_red, Res_red, indices_to_delete, indices] = applyBoundaryConditions(K_global, Res_global, first_set, fourth_set, n_nodes);

        del_ug_red = -K_red\Res_red;
        
        
        del_ug(indices) = del_ug_red;
        u_global = u_global + del_ug;

        u_global ;
       
      
      norm(Res_red);

    %Convergence conditions  
   if ((norm(Res_red) < 1e-06) && norm(del_ug) < 1e-06)
       disp("convergence achieved")
       tau;
       sdv_up_values = cell((size(elem_list, 1)), 6);    %a loop like this will run below as well, this is for creating a cell array and storing current sdv values, and in the loop below we have to store the values of pvs sdv, it should not clash. 
       for o = 1:6
         for o_e = 1:size(elem_list, 1)
             sdv_up_values{o_e, o} = sdv_values;
         end
       end
       store_n50 = [store_n50; u_global(no_50*2 - 1)];
       tau_store = [tau_store; (tau/sigmainf)];
       tau = tau + d_tau;
       k;
       u_global;
       break
       
   elseif k>20
       disp("Iterations are exceeding desired value")
       tau = tau - d_tau;
       d_tau = d_tau/10;
        
       tau = tau + d_tau;
       u_global = bk_u_global;
       stored_stress =  prev_stored_stress;
       k;
       break
   else

       k = k + 1; 
       continue
       

   end     



    end
end


% Having the sorted stress values for each case: sorted_sigma_0, sorted_sigma_45, sorted_sigma_90
a=2.5;
[sorted_sigma_0, sorted_sigma_45, sorted_sigma_90, ray_0, ray_45, ray_90] = plotting(stored_stress, sigmainf,a);

% Plotting the first subplot for Case 0
subplot(1, 3, 1);
plot(sorted_sigma_0(:, 1), sorted_sigma_0(:, 3), 'r-', 'LineWidth', 1.5); % Plots stress_rr vs r in red
hold on;
plot(ray_0(:, 1), ray_0(:, 3), 'k--', 'LineWidth', 1.2); % Plots ray_0 stress
plot(ray_0(:, 1), ray_0(:, 4), 'm--', 'LineWidth', 1.2); % Plots ray_0 stress
plot(ray_0(:, 1), ray_0(:, 5), 'c--', 'LineWidth', 1.2); % Plots ray_0 stress
plot(sorted_sigma_0(:, 1), sorted_sigma_0(:, 4), 'g-', 'LineWidth', 1.5); % Plots stress_thth vs r in green
plot(sorted_sigma_0(:, 1), sorted_sigma_0(:, 5), 'b-', 'LineWidth', 1.5); % Plots stress_rth vs r in blue
hold off;
xlabel('r');
ylabel('Stress');
title('Stress vs r (Case 0)');
legend('Stress_{rr}', 'Ray_0 Stress_{rr}', 'Ray_0 Stress_{\theta\theta}', 'Ray_0 Stress_{r\theta}', 'Stress_{\theta\theta}', 'Stress_{r\theta}');
grid on;

% Plotting the second subplot for Case 45
subplot(1, 3, 2);
plot(sorted_sigma_45(:, 1), sorted_sigma_45(:, 3), 'r-', 'LineWidth', 1.5); % Plots stress_rr vs r in red
hold on;
plot(ray_45(:, 1), ray_45(:, 3), 'k--', 'LineWidth', 1.2); % Plots ray_45 stress
plot(ray_45(:, 1), ray_45(:, 4), 'm--', 'LineWidth', 1.2); % Plots ray_45 stress
plot(ray_45(:, 1), ray_45(:, 5), 'c--', 'LineWidth', 1.2); % Plots ray_45 stress
plot(sorted_sigma_45(:, 1), sorted_sigma_45(:, 4), 'g-', 'LineWidth', 1.5); % Plots stress_thth vs r in green
plot(sorted_sigma_45(:, 1), sorted_sigma_45(:, 5), 'b-', 'LineWidth', 1.5); % Plots stress_rth vs r in blue
hold off;
xlabel('r');
ylabel('Stress');
title('Stress vs r (Case 45)');
legend('Stress_{rr}', 'Ray_{45} Stress_{rr}', 'Ray_{45} Stress_{\theta\theta}', 'Ray_{45} Stress_{r\theta}', 'Stress_{\theta\theta}', 'Stress_{r\theta}');
grid on;

% Plotting the third subplot for Case 90
subplot(1, 3, 3);
plot(sorted_sigma_90(:, 1), sorted_sigma_90(:, 3), 'r-', 'LineWidth', 1.5); % Plots stress_rr vs r in red
hold on;
plot(ray_90(:, 1), ray_90(:, 3), 'k--', 'LineWidth', 1.2); % Plots ray_90 stress
plot(ray_90(:, 1), ray_90(:, 4), 'm--', 'LineWidth', 1.2); % Plots ray_90 stress
plot(ray_90(:, 1), ray_90(:, 5), 'c--', 'LineWidth', 1.2); % Plots ray_90 stress
plot(sorted_sigma_90(:, 1), sorted_sigma_90(:, 4), 'g-', 'LineWidth', 1.5); % Plots stress_thth vs r in green
plot(sorted_sigma_90(:, 1), sorted_sigma_90(:, 5), 'b-', 'LineWidth', 1.5); % Plots stress_rth vs r in blue
hold off;
xlabel('r');
ylabel('Stress');
title('Stress vs r (Case 90)');
legend('Stress_{rr}', 'Ray_{90} Stress_{rr}', 'Ray_{90} Stress_{\theta\theta}', 'Ray_{90} Stress_{r\theta}', 'Stress_{\theta\theta}', 'Stress_{r\theta}');
grid on;

% Adjust the subplot sizes for better visualization
set(gcf, 'Position', get(0, 'Screensize'));
% Create a new figure for displacement plot
figure;
plot(tau_store, store_n50, 'b', 'Marker', 'o', 'DisplayName', 'Displacement');
xlabel('load')
ylabel('displacement')
title('Displacement plot')
%%
