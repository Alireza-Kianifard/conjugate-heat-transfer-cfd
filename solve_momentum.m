function [U, V] = solve_momentum(U_old, V_old, P, domain_ID, y, dy, dz, rho, mu_eff, U_inlet, alpha_U)
    % Solves the discretized U and V momentum equations using a finite volume method.
    
    [Ny, Nz] = size(U_old);
    N_total = Ny * Nz;

    % U-Momentum (Z-direction)
    I_u = []; J_u = []; V_u = [];
    b_u = zeros(N_total, 1);
    
    for j = 1:Nz
        for i = 1:Ny
            idx = (j-1)*Ny + i;
            
            % Boundary Conditions & Solid Domain
            if domain_ID(i,j) == 2 % Solid
                I_u(end+1) = idx; J_u(end+1) = idx; V_u(end+1) = 1; b_u(idx) = 0;
            elseif j == 1 % Inlet
                I_u(end+1) = idx; J_u(end+1) = idx; V_u(end+1) = 1; b_u(idx) = U_inlet;
            elseif j == Nz % Outlet
                I_u(end+1:end+2) = [idx, idx]; J_u(end+1:end+2) = [idx, idx-Ny]; V_u(end+1:end+2) = [1, -1];
            elseif i == 1 || i == Ny % Top/Bottom Walls (No-slip)
                I_u(end+1) = idx; J_u(end+1) = idx; V_u(end+1) = 1; b_u(idx) = 0;
            else % Internal Fluid Cell
                An=dy(i);As=dy(i);Ae=dz(j);Aw=dz(j);
                Fe=rho*V_old(i,j)*Ae;Fw=rho*V_old(i-1,j)*Aw;Fn=rho*U_old(i,j)*An;Fs=rho*U_old(i,j-1)*As;
                De=mu_eff(i,j)*Ae/dy(i);Dw=mu_eff(i,j)*Aw/dy(i);Dn=mu_eff(i,j)*An/dz(j);Ds=mu_eff(i,j)*As/dz(j);
                aW=Dw+max(Fw,0);aE=De+max(-Fe,0);aS=Ds+max(Fs,0);aN=Dn+max(-Fn,0);
                Sp=(P(i,j-1)-P(i,j))*dy(i);
                aP=aW+aE+aS+aN+(Fe-Fw)+(Fn-Fs);
                I_u(end+1:end+5)=[idx,idx,idx,idx,idx]; J_u(end+1:end+5)=[idx,idx-1,idx+1,idx-Ny,idx+Ny];
                V_u(end+1:end+5)=[aP/alpha_U,-aW,-aE,-aS,-aN];
                b_u(idx) = Sp + (1-alpha_U)*(aP/alpha_U)*U_old(i,j);
            end
        end
    end
    A_u = sparse(I_u, J_u, V_u, N_total, N_total);
    U_vec = A_u \ b_u;
    U = reshape(U_vec, [Ny, Nz]);

    % V-Momentum (Y-direction)
    I_v = []; J_v = []; V_v = [];
    b_v = zeros(N_total, 1);
    
    for j = 1:Nz
        for i = 1:Ny
            idx = (j-1)*Ny + i;
            
            % Boundary Conditions & Solid Domain
            if domain_ID(i,j)==2 || j==1 || i==1 || i==Ny % Solid, Inlet, Walls
                I_v(end+1)=idx; J_v(end+1)=idx; V_v(end+1)=1; b_v(idx) = 0;
            elseif j==Nz % Outlet
                I_v(end+1:end+2)=[idx,idx]; J_v(end+1:end+2)=[idx,idx-Ny]; V_v(end+1:end+2)=[1,-1];
            else % Internal Fluid Cell
                An=dy(i);As=dy(i);Ae=dz(j);Aw=dz(j);
                Fe=rho*V_old(i,j)*Ae;Fw=rho*V_old(i-1,j)*Aw;Fn=rho*U_old(i,j)*An;Fs=rho*U_old(i,j-1)*As;
                De=mu_eff(i,j)*Ae/dy(i);Dw=mu_eff(i,j)*Aw/dy(i);Dn=mu_eff(i,j)*An/dz(j);Ds=mu_eff(i,j)*As/dz(j);
                aW=Dw+max(Fw,0);aE=De+max(-Fe,0);aS=Ds+max(Fs,0);aN=Dn+max(-Fn,0);
                Sp=(P(i-1,j)-P(i,j))*dz(j);
                aP=aW+aE+aS+aN+(Fe-Fw)+(Fn-Fs);
                I_v(end+1:end+5)=[idx,idx,idx,idx,idx]; J_v(end+1:end+5)=[idx,idx-1,idx+1,idx-Ny,idx+Ny];
                V_v(end+1:end+5)=[aP/alpha_U,-aW,-aE,-aS,-aN];
                b_v(idx)=Sp+(1-alpha_U)*(aP/alpha_U)*V_old(i,j);
            end
        end
    end
    A_v = sparse(I_v, J_v, V_v, N_total, N_total);
    V_vec = A_v \ b_v;
    V = reshape(V_vec, [Ny, Nz]);
end
