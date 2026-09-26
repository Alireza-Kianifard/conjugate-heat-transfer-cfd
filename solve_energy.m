function T = solve_energy(T_old, U, V, domain_ID, y, dy, dz, dt, rho_f, cp_f, k_f, rho_s, cp_s, k_s, T_inlet, q_dot_solid)
    % Solves the transient energy equation for both fluid and solid domains.

    [Ny, Nz] = size(T_old);
    N_total = Ny * Nz;
    
    % Preallocate sparse matrix arrays
    max_entries = 5 * N_total;
    I = zeros(max_entries, 1);
    J = zeros(max_entries, 1);
    V_val = zeros(max_entries, 1);
    b = zeros(N_total, 1);
    counter = 0;
    
    for j = 1:Nz
        for i = 1:Ny
            idx = (j-1)*Ny + i;

            % Boundary Conditions
            if j == 1 % Inlet (Fixed Temperature)
                counter=counter+1; I(counter)=idx; J(counter)=idx; V_val(counter)=1; b(idx)=T_inlet;
                continue;
            elseif j == Nz % Outlet (Zero Gradient)
                counter=counter+1; I(counter)=idx; J(counter)=idx; V_val(counter)=1;
                counter=counter+1; I(counter)=idx; J(counter)=idx-Ny; V_val(counter)=-1;
                continue;
            elseif i == 1 % Bottom Wall (Adiabatic)
                counter=counter+1; I(counter)=idx; J(counter)=idx; V_val(counter)=1;
                counter=counter+1; I(counter)=idx; J(counter)=idx+1; V_val(counter)=-1;
                continue;
            elseif i == Ny % Top Wall (Adiabatic)
                counter=counter+1; I(counter)=idx; J(counter)=idx; V_val(counter)=1;
                counter=counter+1; I(counter)=idx; J(counter)=idx-1; V_val(counter)=-1;
                continue;
            end
            
            % Internal Cells (Fluid and Solid)
            is_fluid_P = (domain_ID(i,j) == 1);
            
            if is_fluid_P
                rho = rho_f; cp = cp_f; k_P = k_f;
            else
                rho = rho_s; cp = cp_s; k_P = k_s;
            end

            Vol = dy(i)*dz(j);
            ap0 = rho*cp*Vol/dt;
            
            % Face Areas
            Ae = dy(i); % East/West face area
            Aw = dy(i);
            An = dz(j); % North/South face area
            As = dz(j);
            
            % Distances to neighbors
            dist_E = dz(j); dist_W = dz(j);
            dist_N = dy(i); dist_S = dy(i);

            % Neighboring conductivities
            k_E_val = k_s; if domain_ID(i,j+1)==1, k_E_val=k_f; end
            k_W_val = k_s; if domain_ID(i,j-1)==1, k_W_val=k_f; end
            k_N_val = k_s; if domain_ID(i+1,j)==1, k_N_val=k_f; end
            k_S_val = k_s; if domain_ID(i-1,j)==1, k_S_val=k_f; end

            % Harmonic Mean for Interface Conductivity
            ke = 2 * k_P * k_E_val / (k_P + k_E_val + 1e-12);
            kw = 2 * k_P * k_W_val / (k_P + k_W_val + 1e-12);
            kn = 2 * k_P * k_N_val / (k_P + k_N_val + 1e-12);
            ks_eff = 2 * k_P * k_S_val / (k_P + k_S_val + 1e-12);

            % Diffusive terms (W/K)
            De = ke * Ae / dist_E;
            Dw = kw * Aw / dist_W;
            Dn = kn * An / dist_N;
            Ds = ks_eff * As / dist_S;
            
            % Convective mass fluxes at cell faces
            Fe = 0; Fw = 0; Fn = 0; Fs_conv = 0;
            if is_fluid_P
                % East face
                if domain_ID(i, j+1) == 1
                    u_e = 0.5 * (U(i,j) + U(i,j+1));
                else
                    u_e = 0; % Solid wall
                end
                % West face
                if domain_ID(i, j-1) == 1
                    u_w = 0.5 * (U(i,j) + U(i,j-1));
                else
                    u_w = 0; % Solid wall
                end
                % North face
                if domain_ID(i+1, j) == 1
                    v_n = 0.5 * (V(i,j) + V(i+1,j));
                else
                    v_n = 0; % Solid wall
                end
                % South face
                if domain_ID(i-1, j) == 1
                    v_s = 0.5 * (V(i,j) + V(i-1,j));
                else
                    v_s = 0; % Solid wall
                end

                Fe = rho_f * cp_f * u_e * Ae;
                Fw = rho_f * cp_f * u_w * Aw;
                Fn = rho_f * cp_f * v_n * An;
                Fs_conv = rho_f * cp_f * v_s * As;
            end
            
            % Heat source term
            S_h = (1-is_fluid_P) * q_dot_solid * Vol;

            % Upwind formulation for neighbor coefficients
            aW = Dw + max(Fw, 0);
            aE = De + max(-Fe, 0);
            aS_coeff = Ds + max(Fs_conv, 0);
            aN = Dn + max(-Fn, 0);
            
            % Diagonally dominant aP (enforcing conservative boundedness)
            aP = aW + aE + aS_coeff + aN + ap0;
            
            % Populate matrix
            counter=counter+1; I(counter)=idx; J(counter)=idx; V_val(counter)=aP;
            counter=counter+1; I(counter)=idx; J(counter)=idx-Ny; V_val(counter)=-aW;
            counter=counter+1; I(counter)=idx; J(counter)=idx+Ny; V_val(counter)=-aE;
            counter=counter+1; I(counter)=idx; J(counter)=idx-1; V_val(counter)=-aS_coeff;
            counter=counter+1; I(counter)=idx; J(counter)=idx+1; V_val(counter)=-aN;
            
            b(idx) = ap0*T_old(i,j) + S_h;
        end
    end
    
    % Trim arrays and solve
    I = I(1:counter); J = J(1:counter); V_val = V_val(1:counter);
    A = sparse(I, J, V_val, N_total, N_total);
    T_vec = A \ b;
    T = reshape(T_vec, [Ny, Nz]);
end
