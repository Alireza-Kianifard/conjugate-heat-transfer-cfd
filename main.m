%% Conjugate Heat Transfer (CHT) MATLAB Solver
% Description: 2D Planar CHT simulation for turbulent fluid flow over heated catalyst blocks in a channel using the SIMPLEC algorithm.
% Features:
%  - Geometry: 2D Planar Cartesian (Stretched Mesh)
%  - Flow: Steady-state, Incompressible, Turbulent (Mixing-Length)
%  - Heat Transfer: Transient conduction/convection
% Usage: Run this script directly. Ensure all helper functions are in the same directory.
clear; clc; close all;

%% 1. Parameters & Properties
fprintf('--- Initializing Simulation Parameters ---\n');
% Solver Controls
max_flow_iter = 5000;    % Max iterations for the flow solver
flow_tolerance = 5e-5;   % Convergence criterion for the flow solver
alpha_U = 0.5;           % Velocity under-relaxation factor
alpha_P = 0.2;           % Pressure under-relaxation factor
dt = 0.01;               % Time step for thermal simulation [s]
t_end = 100;             % Total simulation time [s]

% Fluid Properties (Air)
rho_f = 1.225;      % Density [kg/m^3]
mu_f = 1.7894e-5;   % Dynamic viscosity [kg/(m*s)]
cp_f = 1006.43;     % Specific heat [J/kg-K]
k_f = 0.0242;       % Thermal conductivity [W/m-K]

% Solid Properties (Nickel Catalyst)
rho_s = 8900;       % Density [kg/m^3]
cp_s = 460.6;       % Specific heat [J/kg-K]
k_s = 91.74;        % Thermal conductivity [W/m-K]

q_dot_solid = 0;  % Volumetric Heat Generation [W/m^3]
% A value of ~1.7e7 should result in a max temp of ~150 C.

% Boundary & Initial Conditions
U_inlet = 2.0;                 % Inlet velocity [m/s]
T_inlet = 100 + 273.15;        % Inlet air temperature [K]
T_initial_solid = 25 + 273.15; % Initial solid temperature [K]

%% 2. Meshing and Initialization
fprintf('--- Generating Stretched Mesh and Initializing Fields ---\n');
% A higher 'Ratio' creates a finer mesh for better accuracy.
% Recommended value: more than 5 for this geometry
Ratio = 4; 
Ny_input = 10 * Ratio; 
Nz_input = 100 * Ratio;
% Change 10 and 100 if you completely changed geometry
[y, z, dy, dz, Ny, Nz, domain_ID] = mesh_generator(Ny_input, Nz_input);

% Initialize fields
U = ones(Ny, Nz) * U_inlet; 
V = zeros(Ny, Nz); 
P = zeros(Ny, Nz);
T = T_inlet * ones(Ny, Nz);
T(domain_ID == 2) = T_initial_solid; % Set initial temp for solid
U_old = U;

%% 3. Steady-State Flow Solver (SIMPLEC Algorithm)
fprintf('--- Starting Steady-State Flow Solver (SIMPLEC) ---\n');
iter_history = []; U_res_history = [];

for iter = 1:max_flow_iter
    % Calculate effective viscosity (laminar + turbulent)
    mu_eff = calc_eddy_viscosity(U, y, dy, rho_f, mu_f, domain_ID);

    % Solve momentum and pressure correction equations
    [U_star, V_star] = solve_momentum(U, V, P, domain_ID, y, dy, dz, rho_f, mu_eff, U_inlet, alpha_U);
    [P, U, V] = solve_pressure_correction(P, U_star, V_star, domain_ID, y, dy, dz, rho_f, alpha_P);

    % This is critical to prevent non-physical vortex shedding oscillations
    % and ensure a stable, symmetric solution.
    U = (U + flipud(U)) / 2;
    V = (V - flipud(V)) / 2;
    P = (P + flipud(P)) / 2;

    % Check for convergence
    residual = sum(abs(U(:) - U_old(:))) / (Ny * Nz);
    iter_history = [iter_history, iter];
    U_res_history = [U_res_history, residual];

    if isnan(residual), error('Solver diverged. Check parameters or boundary conditions.'); end
    if residual < flow_tolerance, fprintf('Flow converged successfully after %d iterations.\n', iter); break; end

    U_old = U;
    if mod(iter, 50) == 0, fprintf('Iteration %d | U-residual = %.2e\n', iter, residual); end
end
if iter == max_flow_iter, fprintf('Warning: Flow solver did not converge within max iterations.\n'); end

%% 4. Unsteady Conjugate Heat Transfer (CHT) Solver
fprintf('--- Starting Unsteady Thermal Solver ---\n');
time = 0; 
time_history = []; 
T_max_history = [];
max_P_history = [];

while time < t_end
    T_old = T;
    
    % Solve energy equation for one time step
    T = solve_energy(T_old, U, V, domain_ID, y, dy, dz, dt, ...
                     rho_f, cp_f, k_f, rho_s, cp_s, k_s, T_inlet, q_dot_solid);
                 
    time = time + dt;
    
    % Store history for plotting
    if any(domain_ID(:) == 2), T_max_history = [T_max_history, max(T(domain_ID == 2))]; end
    time_history = [time_history, time];
    max_P_history = [max_P_history, max(P(:))];                     

    if mod(round(time), 10) == 0 && mod(time, 1) < dt
        fprintf('Time = %.1f s | Max Flow Temp = %.2f K (%.2f C)\n', time, max(T(:)), max(T(:))-273.15);
    end
end
fprintf('--- Simulation Complete ---\n');

%% 5. Visualization
fprintf('--- Generating Plots ---\n');
[Z, Y] = meshgrid(z, y);

% Plot 1: Flow Solver Convergence
figure('Name', 'Flow Convergence', 'Position', [100 100 600 400]);
semilogy(iter_history, U_res_history, 'r-', 'LineWidth', 2);
grid on; xlabel('Iteration Number'); ylabel('Residual'); title('Flow Convergence');

% Plot 2: Velocity Field
figure('Name', 'Velocity Field', 'Position', [700 100 800 400]);
contourf(Z, Y, sqrt(U.^2 + V.^2), 50, 'LineColor', 'none'); hold on;
rectangle('Position', [0.240, 0.0035, 0.020, 0.0065], 'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
rectangle('Position', [0.240, -0.010, 0.020, 0.0065], 'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
streamslice(Z, Y, U, V, 2); 
colorbar; xlabel('z [m]'); ylabel('y [m]'); title('Velocity Magnitude & Streamlines [m/s]'); axis equal tight;

% Plot 3: Temperature Field
figure('Name', 'Temperature Field', 'Position', [700 500 800 400]);
contourf(Z, Y, T - 273.15, 50, 'LineColor', 'none'); hold on;
rectangle('Position', [0.240, 0.0035, 0.020, 0.0065], 'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
rectangle('Position', [0.240, -0.010, 0.020, 0.0065], 'EdgeColor', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
colorbar; xlabel('z [m]'); ylabel('y [m]'); title('Static Temperature [°C]'); axis equal tight;

% Plot 4: Transient Metrics
figure('Name', 'Transient Metrics', 'Position', [100 100 600 400]);
semilogy(time_history, T_max_history - 273.15, 'b-', 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Max Temp [°C]'); title('Max Catalyst Temperature');
