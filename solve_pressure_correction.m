function [P_new, U_new, V_new] = solve_pressure_correction(P, U_star, V_star, domain_ID, y, dy, dz, rho, alpha_P)
    % Solves the pressure correction equation and updates the P, U, and V fields
    % using the SIMPLEC method.
    
    [Ny, Nz] = size(P);
    N_total = Ny * Nz;
    
    I_p = []; J_p = []; V_p = [];
    b_p = zeros(N_total, 1);
    d = ones(Ny, Nz); % This would hold the 'd' terms from the momentum eq. For simplicity, we use 1.

    for j = 1:Nz
        for i = 1:Ny
            idx = (j-1)*Ny + i;

            % Boundary Conditions
            if domain_ID(i,j) == 2 % Solid (no pressure correction)
                I_p(end+1) = idx; J_p(end+1) = idx; V_p(end+1) = 1;
            elseif j == Nz % Outlet (fixed pressure)
                I_p(end+1) = idx; J_p(end+1) = idx; V_p(end+1) = 1;
            elseif j == 1 % Inlet (zero gradient)
                I_p(end+1:end+2) = [idx,idx]; J_p(end+1:end+2) = [idx,idx+Ny]; V_p(end+1:end+2) = [1,-1];
            elseif i == 1 || i == Ny % Top/Bottom Walls (zero gradient)
                I_p(end+1:end+2) = [idx,idx]; J_p(end+1:end+2) = [idx, idx+(i==1)*2-1]; V_p(end+1:end+2) = [1,-1];
            else % Internal Fluid Cell
                An=dy(i);As=dy(i);Ae=dz(j);Aw=dz(j);
                
                % Mass imbalance term
                m_dot_out = rho * (U_star(i,j+1)*An + V_star(i+1,j)*Ae);
                m_dot_in  = rho * (U_star(i,j)*As + V_star(i,j)*Aw);
                b_p(idx) = m_dot_in - m_dot_out;

                aE = rho*d(i,j)*Ae; aW = rho*d(i-1,j)*Aw; 
                aN = rho*d(i,j)*An; aS = rho*d(i,j-1)*As;
                aP = aE + aW + aN + aS;
                
                I_p(end+1:end+5) = [idx,idx,idx,idx,idx];
                J_p(end+1:end+5) = [idx,idx+1,idx-1,idx+Ny,idx-Ny];
                V_p(end+1:end+5) = [aP,-aE,-aW,-aN,-aS];
            end
        end
    end
    
    A_p = sparse(I_p, J_p, V_p, N_total, N_total);
    
    % Solve the linear system for P_prime
    P_prime_vec = A_p \ b_p;
    P_prime = reshape(P_prime_vec, [Ny, Nz]);
    
    % Update pressure and velocities
    P_new = P + alpha_P * P_prime;
    U_new = U_star;
    V_new = V_star;
    
    % Correct U-velocity
    for j = 2:Nz
        for i = 1:Ny
            if domain_ID(i,j) == 1
                U_new(i,j) = U_star(i,j) + d(i,j) * (P_prime(i,j-1) - P_prime(i,j));
            end
        end
    end
    
    % Correct V-velocity
    for j = 1:Nz
        for i = 2:Ny
            if domain_ID(i,j) == 1
                V_new(i,j) = V_star(i,j) + d(i,j) * (P_prime(i-1,j) - P_prime(i,j));
            end
        end
    end
end
