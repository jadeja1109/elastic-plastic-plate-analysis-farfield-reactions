function [K_red, Res_red, indices_to_delete, indices] = applyBoundaryConditions(K_global, Res_global, first_set, fourth_set, n_nodes)
    
    
    % Find the indices of the nodes in the first and fourth sets in the global displacement vector
    u_indices_x = (fourth_set - 1) * 2 + 1; % Indices for x-direction displacements
    u_indices_y = first_set * 2; % Indices for y-direction displacements
    
    indices_to_delete = [u_indices_x u_indices_y];
    %initializing a matrix and getting another set of indices which will be used later
    mat = 1: 2 * n_nodes ;
    indices = setdiff(mat, indices_to_delete);

    % Apply boundary conditions to the stiffness matrix and force vector
    K_red = K_global;
    K_red(indices_to_delete, :) = []; % Set rows corresponding to x-direction displacements to zero
    K_red(:, indices_to_delete) = []; % Set columns corresponding to x-direction displacements to zero

    % Entries in the force vector corresponding to y-direction displacements to zero
    Res_red = Res_global;
    Res_red(indices_to_delete) = []; % Set entries in the force vector corresponding to x-direction displacements to zero
end
