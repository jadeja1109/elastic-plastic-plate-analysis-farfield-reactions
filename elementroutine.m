function  [K_global, Res_global, ray_0, ray_45, ray_90, stored_stress, stored_p_stress, sdv_values] = elementroutine(elem, elem_list,elem_sets, node_list, node_sets, n_node_p_elem, n_nodes,  u_global, sdv_values, sdv, props ,hard_funct, stored_stress, stored_p_stress, ray_0, ray_45, ray_90, sigma_analytic_0, sigma_analytic_45, sigma_analytic_90, tau, K_global, Res_global, sigmainf);

 
            % Gauss quadrature points and weights
            L1 = 0.816847572980459;
            L2 = 0.091576213509771;
            L3 = 0.108103018168070;
            L4 = 0.445948490915965;
            W1 = 0.109951743655322 / 2.0;
            W2 = 0.223381589678011 / 2.0;
            
            si = [L1, L2, L2, L3, L4, L4];
            ti = [L2, L1, L2, L4, L3, L4];
            wi = [W1, W1, W1, W2, W2, W2];


            % Get the node numbers for the current element
            element_nodes = elem_list(elem, 2:end);
            
            % Get the corresponding node coordinates from node_list
            x_e = node_list(element_nodes, 2:3); % Assuming 2D coordinates

            
            
           
            %Assembly Matrix calculations
            Ae = sparse(2 * n_node_p_elem, n_nodes * 2);

    
            x = elem_list(elem, 2:end);

            for p = 1:n_node_p_elem

                Ae((p*2)-1, (x(1,p))*2-1) = 1;
                Ae(p*2, (x(1,p))*2) = 1;

            end

            u_e = Ae * u_global;
            u_e;
        
            % Initialize element stiffness matrix and external/internal force vectors
            Kt_e = zeros(2*n_node_p_elem, 2*n_node_p_elem);
            Fint_e = zeros(2*n_node_p_elem, 1);
            
        
            % Loop over Gauss quadrature points
            for q = 1:length(wi)
                s = si(q);
                t = ti(q);
                w = wi(q);
                
                
                N1 = (1-s-t)*(2*(1-s-t)-1);
                N2 = (2*s^2 - s);
                N3 = t*(2*t - 1);
                N4 = 4*s*(1 - s - t);
                N5 = 4*s*t;
                N6 = 4*t*(1 - s - t);
        
                N = [ N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0; 
                      0, N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6];

                N_matrix = [N1, N2, N3, N4, N5, N6];
        
                %dN_dz(:, :, q) = [-3 + 4*t + 4*s, 0, 4*s - 1, 0, 0, 0, 4 - 8*s - 4*t, 0, 4*t, 0, -4*t, 0;
                                   % 0, -3 + 4*t + 4*s, 0, 0, 0, 4*t - 1, 0, -4*s, 0, 4*s, 0, 4 - 4*s - 8*t];
                
                t_dN_dz = [-3 + 4*t + 4*s, 4*s - 1, 0, 4 - 8*s - 4*t, 4*t, -4*t;
                                  -3 + 4*t + 4*s, 0, 4*t - 1, -4*s, 4*s, 4 - 4*s - 8*t];
               
        
                % Compute the Jacobian matrix (J) for the current element
                J = t_dN_dz * x_e;
                
                % Compute the determinant of the Jacobian matrix
                detJ = det(J);
                
                % Calculate the inverse of the Jacobian matrix
                invJ = inv(J);
        
                dN_dz = invJ * t_dN_dz;
        
            
                %B matrix
                
                B=  [dN_dz(1,1),0,dN_dz(1,2),0,dN_dz(1,3),0,dN_dz(1,4),0,dN_dz(1,5),0,dN_dz(1,6),0;
                    0,dN_dz(2,1),0,dN_dz(2,2),0,dN_dz(2,3),0,dN_dz(2,4),0,dN_dz(2,5),0,dN_dz(2,6);
                    dN_dz(2,1),dN_dz(1,1),dN_dz(2,2),dN_dz(1,2),dN_dz(2,3),dN_dz(1,3),dN_dz(2,4),dN_dz(1,4),dN_dz(2,5),dN_dz(1,5),dN_dz(2,6),dN_dz(1,6)];
               
   
                %strain
                strain = B*u_e;

                
                sdv = sdv_values{elem, q};
   
                %calling material routine
                [stress,C_t,sdv,eps33] = elastic_plastic_von_mises_model_plane_stress(strain, sdv, props, hard_funct);
                %[stress, C_t, sdv, eps33] = elastic_plastic_von_mises_model_plane_stress(strain, sdv, props, hard_funct);
                
                sdv_values{elem , q} = sdv;
               
                
                % Compute the element stiffness matrix contribution at the current Gauss point
                Kt_e = Kt_e + B' * C_t * B * detJ * w;
                
               % Compute the internal force vector contribution at the current Gauss point
                Fint_e = Fint_e + B' * stress * detJ * w;
                %disp(Fint_e)


                %storing of stress values
                cd = N_matrix * x_e;  %cd stands for coordinates in the physical space
                st_stress = [cd stress'];  %stored stress in a 5x1 matrix
                stored_stress  = [stored_stress; st_stress];
                


            end 
         
        
            %Compute external force vector
        if ismember(elem_list(elem), cell2mat(elem_sets(2,2)))  %this is the condition to check that our elem lies on the right edge.
            Fext_e = Fextcalc(elem_sets, node_sets, elem_list, node_list, n_node_p_elem, x_e, tau,elem);

        else

            Fext_e = zeros(n_node_p_elem*2, 1);
            
        end
 
            % Compute the external force vector contribution at the current Gauss point
                
        
        
            %Assembly of global matrices
            K_global = K_global + Ae' * Kt_e * Ae;
            Res_global = Res_global + Ae' * (Fint_e - Fext_e);

end
