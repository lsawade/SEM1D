% WAVE - 1D wave equation solver using Spectral Element Method (SOLUTION)
% Based on Fortran implementation from hw8A/wave.f90

function wave(bc_type_str)
% WAVE 1D wave equation solver
%
% Usage:
%   wave()              - runs with Dirichlet BC (default)
%   wave('dirichlet')   - runs with Dirichlet BC
%   wave('neumann')     - runs with Neumann BC

if nargin < 1
    bc_type_str = 'dirichlet';
end

% Add functions directory to path
script_path = fileparts(mfilename('fullpath'));
addpath(fullfile(script_path, 'functions'));

% Boundary condition types
BOTH_DIRICHLET = 0;
BOTH_NEUMANN = 1;

% Parse boundary condition type
if strcmpi(bc_type_str, 'dirichlet')
    bc_type = BOTH_DIRICHLET;
elseif strcmpi(bc_type_str, 'neumann')
    bc_type = BOTH_NEUMANN;
else
    error('Invalid BC type. Use ''dirichlet'' or ''neumann''');
end

if bc_type == BOTH_DIRICHLET
    fprintf('--------------------------------------------------------\n');
    fprintf('BC type: Dirichlet\n');
    fprintf('--------------------------------------------------------\n');
else
    fprintf('--------------------------------------------------------\n');
    fprintf('BC type: Neumann\n');
    fprintf('--------------------------------------------------------\n');
end

%% Parameters

% Number of spectral elements
NSPEC = 20;

% Number of GLL points (polynomial degree plus one)
NGLL = 7;

% Number of global points
NGLOB = (NGLL-1)*NSPEC + 1;

% Number of timesteps
NSTEP = 2000;

% Time step in seconds
DT_REF = 0.25;

% Fixed boundary conditions
FIXED_BC = true;

% Model parameters (SI)
LENGTH = 100;

%% Define material properties as functions of position
% Users can modify these functions to create heterogeneous models

% Density as a function of x
get_density = @(x) ones(size(x));  % Default: homogeneous with rho = 1

% Shear modulus as a function of x
get_shear = @(x) ones(size(x));    % Default: homogeneous with mu = 1

% Example heterogeneous models (uncomment to use):
%
% % Two-layer model (interface at x = 50)
% get_density = @(x) 1.0 + 0.5 * (x > 50);
% get_shear = @(x) 1.0 + 1.0 * (x > 50);
%
% % Gradient model
% get_density = @(x) 1.0 + 0.01 * x;
% get_shear = @(x) 1.0 + 0.02 * x;
%
% % Low-velocity zone
% get_density = @(x) 1.0 - 0.3 * exp(-0.01 * (x - 50).^2);
% get_shear = @(x) 1.0 - 0.5 * exp(-0.01 * (x - 50).^2);

%% Setup derivative matrix and GLL points

[xigll, wgll, hprime] = define_derivative_matrix(NGLL);

%% Setup mesh

% Evenly spaced anchors between 0 and LENGTH
x1 = zeros(NSPEC, 1);
x2 = zeros(NSPEC, 1);
for ispec = 1:NSPEC
    x1(ispec) = LENGTH*(ispec-1)/NSPEC;
    x2(ispec) = LENGTH*ispec/NSPEC;
end

% Set up the mesh properties
rho = zeros(NGLL, NSPEC);
shear = zeros(NGLL, NSPEC);
inverse_jacobian = zeros(NGLL, NSPEC);
jacobian = zeros(NGLL, NSPEC);

for ispec = 1:NSPEC
    for i = 1:NGLL
        % Compute physical coordinate at this GLL point
        x_gll = 0.5*(1-xigll(i))*x1(ispec) + 0.5*(1+xigll(i))*x2(ispec);

        % Evaluate material properties at this position
        rho(i, ispec) = get_density(x_gll);
        shear(i, ispec) = get_shear(x_gll);

        % These Jacobians only work because the mesh is not deformed
        inverse_jacobian(i, ispec) = 2.0 / (x2(ispec) - x1(ispec));
        jacobian(i, ispec) = (x2(ispec) - x1(ispec)) / 2.0;
    end
end

% Set up local to global numbering
ibool = zeros(NGLL, NSPEC);
iglob = 1;
for ispec = 1:NSPEC
    for i = 1:NGLL
        if i > 1
            iglob = iglob + 1;
        end
        ibool(i, ispec) = iglob;
    end
end

% Get the global grid points and material properties at global points
x = zeros(NGLOB, 1);
rho_global = zeros(NGLOB, 1);
shear_global = zeros(NGLOB, 1);
for ispec = 1:NSPEC
    for i = 1:NGLL
        iglob = ibool(i, ispec);
        x(iglob) = 0.5*(1-xigll(i))*x1(ispec) + 0.5*(1+xigll(i))*x2(ispec);
        rho_global(iglob) = rho(i, ispec);
        shear_global(iglob) = shear(i, ispec);
    end
end

%% Visualize material properties
figure('Position', [100, 600, 900, 400]);
subplot(2,1,1);
plot(x, rho_global, '-b', 'LineWidth', 2);
xlabel('Position (x)');
ylabel('Density (\rho)');
title('Material Properties');
grid on;

subplot(2,1,2);
plot(x, shear_global, '-r', 'LineWidth', 2);
xlabel('Position (x)');
ylabel('Shear Modulus (\mu)');
grid on;
drawnow;

%% Calculate the global mass matrix 'mass_global'
% SOLUTION: Assemble global mass matrix from local contributions
mass_global = zeros(NGLOB, 1);
for ispec = 1:NSPEC
    for i = 1:NGLL
        iglob = ibool(i, ispec);
        mass_local = wgll(i) * rho(i, ispec) * jacobian(i, ispec);
        mass_global(iglob) = mass_global(iglob) + mass_local;
    end
end

%% Estimate the time step
% SOLUTION: Set up time step parameters
dt = DT_REF;
halfdt = 0.5 * dt;
halfdt2 = 0.5 * dt * dt;

fprintf('Time step estimate: %.2f seconds\n', dt);

%% Set up the initial and boundary conditions
% SOLUTION: Initialize Source as a Gaussian pulse
displ_1 = 0.0;
displ_NGLOB = 0.0;
displ = exp(-0.1 * (x - 50.0).^2);
veloc = zeros(NGLOB, 1);
accel = zeros(NGLOB, 1);

% Apply boundary conditions
if bc_type == BOTH_DIRICHLET
    displ(1) = displ_1;
    displ(NGLOB) = displ_NGLOB;
elseif bc_type == BOTH_NEUMANN
    % Nothing to do for traction free BC
else
    error('Unsupported BC type: %d', bc_type);
end


%% Setup real-time visualization
fig = figure('Position', [100, 100, 900, 500]);
h_line = plot(x, displ, '-b', 'LineWidth', 2);
ylim([-1 1]);
xlim([0 LENGTH]);
xlabel('Position (x)');
ylabel('Displacement (u)');
grid on;
if bc_type == BOTH_DIRICHLET
    title_str = 'Dirichlet BC - t = 0.00 s';
else
    title_str = 'Neumann BC - t = 0.00 s';
end
h_title = title(title_str);
drawnow;

%% Time stepping loop

F_global = zeros(NGLOB, 1);
daccel = zeros(NGLOB, 1);

for itime = 1:NSTEP

    %% "Predictor" displacement, velocity, initialize acceleration
    % SOLUTION: Predictor step of Newmark scheme
    displ = displ + dt * veloc + halfdt2 * accel;
    veloc = veloc + halfdt * accel;
    accel = zeros(NGLOB, 1);

    %% Boundary conditions
    % SOLUTION: Apply boundary conditions
    if bc_type == BOTH_DIRICHLET
        displ(1) = displ_1;
        displ(NGLOB) = displ_NGLOB;
    elseif bc_type == BOTH_NEUMANN
        % Nothing to do for traction free BC
    else
        error('Unsupported BC type: %d', bc_type);
    end

    %% Solver
    % SOLUTION: Assemble global force vector and solve for acceleration
    F_global(:) = 0.0;
    for ispec = 1:NSPEC
        for i = 1:NGLL
            iglob = ibool(i, ispec);
            f_local = 0.0;
            for j = 1:NGLL
                iglobj = ibool(j, ispec);
                % Compute stiffness matrix entry
                stiffness_local = -sum(wgll .* shear(:, ispec) .* ...
                                     inverse_jacobian(:, ispec) .* ...
                                     hprime(i, :)' .* hprime(j, :)');
                f_local = f_local + stiffness_local * displ(iglobj);
            end
            F_global(iglob) = F_global(iglob) + f_local;
        end
    end

    % Solve for acceleration
    daccel = F_global ./ mass_global;

    %% "Corrector" acceleration, velocity, displacement
    % SOLUTION: Corrector step of Newmark scheme
    accel = accel + daccel;
    veloc = veloc + halfdt * accel;
    % displ remains unchanged

    % Update plot
    set(h_line, 'YData', displ);
    current_time = dt * (itime - 1);
    if bc_type == BOTH_DIRICHLET
        title_str = sprintf('Dirichlet BC - t = %.2f s', current_time);
    else
        title_str = sprintf('Neumann BC - t = %.2f s', current_time);
    end
    set(h_title, 'String', title_str);
    drawnow;

    % Progress indicator
    fprintf('\rTime step: %4d / %4d', itime, NSTEP);
end % end time loop

fprintf('\nSimulation complete!\n');
end
