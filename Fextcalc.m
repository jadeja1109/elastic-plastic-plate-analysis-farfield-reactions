function Fext_e = Fextcalc(elem_sets, node_sets, elem_list, node_list, n_node_p_elem, x_e, tau,elem)


%finding the element on the right side
right_side_elem = cell2mat(elem_sets(2,2));


%finding the nodes on the right side
right_side_nodes = cell2mat(node_sets(2,2));


%finding the nodes of the elements on the right hand side
right_side_elem_nodes = elem_list(elem, 2:end);
% Display the nodes of the elements on the right side 


Fext_e = zeros(n_node_p_elem*2, 1);
w = 0.5;



for i = 1:(size(right_side_elem_nodes, 1))
    nodes = right_side_elem_nodes(i, :);
    a = node_list((right_side_elem_nodes(i, :)), 2:end);
    
    if ismember(nodes(1), right_side_nodes) && ismember(nodes(2), right_side_nodes)
        edge = 'bottom edge';
        
        sl = [ 0.211324865405187,  0.788675134594813];
        tl = [0,  0];

    elseif ismember(nodes(2), right_side_nodes) && ismember(nodes(3), right_side_nodes)
        edge = 'hypotenuse';
        tl = [0.211324865405187, 0.788675134594813];
        %[ 0.211324865405187,  0.788675134594813];
        sl = [ 1 - 0.211324865405187, 1 - 0.788675134594813];

    else 
        edge = 'side edge';
        sl = [ 0,  0];
        tl = [ 0.211324865405187,  0.788675134594813];
    end
        
        for i = 1:2
            i;
            s = sl(i);
            t = tl(i);

            N1 = (1-s-t)*(2*(1-s-t)-1);
            N2 = (2*s^2 - s);
            N3 = t*(2*t - 1);
            N4 = 4*s*(1 - s - t);
            N5 = 4*s*t;
            N6 = 4*t*(1 - s - t);
    
            N = [ N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0; 
                  0, N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6];
            

            
    
            t_dN_dz = [-3 + 4*t + 4*s, 4*s - 1, 0, 4 - 8*s - 4*t, 4*t, -4*t;
                       -3 + 4*t + 4*s, 0, 4*t - 1, -4*s, 4*s, 4 - 4*s - 8*t];
    

            
            % Compute the Jacobian matrix (J) for the current element
            J = t_dN_dz * a;
            
            J = J';
           

            if strcmp(edge, 'bottom edge')
                J_u = norm(J(:, 1));
            elseif strcmp(edge, 'side edge')
                J_u = norm(J(:, 2));
            else 
                J_u = norm(J * [-1; 1]);

            end
           
             

            trn = [tau; 0];

           
        
            Fext_e = Fext_e + N' * trn * (J_u) * w;

            
         end

end





end





















































































% function Fext = Fextval(elem_sets, elem_list, node_sets, node_list);