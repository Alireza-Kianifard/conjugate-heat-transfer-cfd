function [y, z, dy, dz, Ny, Nz, domain_ID] = mesh_generator(Ny_input, Nz_input)
    % Defines the geometry and creates a stretched mesh for the 2D planar domain.

    % Geometric dimensions
    D_pipe = 0.050;       % Full height of the channel [m]
    L_pipe = 0.500;       % Full length of the channel [m]
    H_half = D_pipe / 2;  % Half-height of the channel
    
    % Catalyst dimensions
    y_inner = 0.0035; % Inner y-coordinate of the top catalyst block
    y_outer = 0.010;  % Outer y-coordinate of the top catalyst block
    L_cat = 0.020;    % Length of the catalyst blocks
    
    z_start_cat = (L_pipe / 2) - (L_cat / 2);
    z_end_cat   = (L_pipe / 2) + (L_cat / 2);

    % Create stretched mesh in Z-direction (finer around catalyst)
    Nz_1 = round(Nz_input * 0.4);
    Nz_2 = round(Nz_input * 0.2);
    Nz_3 = Nz_input - Nz_1 - Nz_2;
    z1 = linspace(0, z_start_cat, Nz_1);
    z2 = linspace(z_start_cat, z_end_cat, Nz_2);
    z3 = linspace(z_end_cat, L_pipe, Nz_3);
    z_faces = unique([z1, z2, z3]);
    Nz = length(z_faces) - 1;
    z = (z_faces(1:end-1) + z_faces(2:end)) / 2;
    dz = diff(z_faces)';

    % Create stretched mesh in Y-direction (finer around catalysts)
    Ny_outer_gap = round(Ny_input * 0.20); % Region between outer wall and catalyst
    Ny_solid     = round(Ny_input * 0.15); % Region of the catalyst block
    Ny_inner_gap = round(Ny_input * 0.30); % Region between the two catalysts
    
    y_faces = [linspace(-H_half, -y_outer, Ny_outer_gap), ...
               linspace(-y_outer, -y_inner, Ny_solid), ...
               linspace(-y_inner, y_inner, Ny_inner_gap), ...
               linspace(y_inner, y_outer, Ny_solid), ...
               linspace(y_outer, H_half, Ny_outer_gap)];
    y_faces = unique(y_faces);
    Ny = length(y_faces) - 1;
    y = (y_faces(1:end-1) + y_faces(2:end)) / 2;
    dy = diff(y_faces)';

    % Identify domain for each cell (1=Fluid, 2=Solid)
    domain_ID = ones(Ny, Nz);
    for i = 1:Ny
        for j = 1:Nz
            is_in_z_range = (z(j) >= z_start_cat && z(j) <= z_end_cat);
            is_in_y_range = ((y(i) >= y_inner && y(i) <= y_outer) || ...
                             (y(i) >= -y_outer && y(i) <= -y_inner));
            
            if is_in_z_range && is_in_y_range
                domain_ID(i,j) = 2; % Solid domain
            end
        end
    end
end
