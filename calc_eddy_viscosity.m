function mu_eff = calc_eddy_viscosity(U, y, dy, rho, mu, domain_ID)
    % Calculates the effective viscosity using a mixing-length turbulence model.
    
    [Ny, Nz] = size(U);
    mu_eff = mu * ones(Ny, Nz);
    kappa = 0.41;     % von Karman constant
    A_plus = 26;      % Van Driest damping constant
    H_half = 0.025;   % Channel half-height

    for j = 1:Nz
        % Calculate distance to the nearest wall for each cell in the column
        y_wall = min(abs(y - H_half), abs(y + H_half)); 
        if any(domain_ID(:, j) == 2) % If catalyst is in this z-slice
             y_wall = min([y_wall; abs(y-0.010); abs(y-0.0035); abs(y+0.0035); abs(y+0.010)], [], 1);
        end
        
        for i = 2:Ny-1
            if domain_ID(i,j) == 1 % Only calculate for fluid cells
                dUdy = (U(i+1,j) - U(i-1,j)) / (dy(i-1) + dy(i));
                u_tau = 0.05 * U(i,j) + 1e-6; % Friction velocity approximation
                y_plus = (rho * u_tau * y_wall(i)) / mu;
                
                % Van Driest damping function
                Damping = 1 - exp(-y_plus / A_plus);
                
                % Mixing length, limited by a fraction of the channel height
                l_m = kappa * y_wall(i) * Damping;
                l_m = min(l_m, 0.09 * (2 * H_half));
                
                % Turbulent viscosity and effective viscosity
                mu_t = rho * (l_m^2) * abs(dUdy);
                mu_eff(i,j) = mu + mu_t;
            end
        end
    end
end
