function result= spm_simulation(t0,tf,time_step,SOC,initialState,input_density)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_SIMULATION
% Solves the Single Particle Model (SPM) over a specified time interval.
%
%       result = spm_simulation(t0, tf, time_step, SOC, initialState, input_density)
%
% INPUT:
%   t0                  - Initial time [s]
%   tf                  - Final time [s]
%   time_step           - Sime step for the simulation [s]
%   SOC                 - Initial state of charge [fraction]
%   initialState        - Structure containing the initial states of the system
%   input_density       - Applied current density to single celle per unit of area [A/m²]
% 
% OUTPUT:
%   result              - Structure containing the simulation results
%
%   This file is part of the POLIDEMO software.
%
%   Official GitHub:   https://github.com/Poli-ON/POLIDEMO
%
%   Reference:
%   F. Pistorio, D. Clerici, A. Somà "POLIDEMO: An electrochemical-mechanical 
%   framework for modeling lithium-ion batteries degradation.", 
%   Applied Energy 404 (2026): 126744.
%   Please cite this work if you use POLIDEMO in your research.
%
%   Copyright (c) 2025:
%       Francesca Pistorio, Davide Clerici
%       Department of Mechanical and Aerospace Engineering,
%       Politecnico di Torino, Corso Duca degli Abruzzi 24,
%       10129, Torino, Italy
%
%   POLIDEMO is free MATLAB-based software distributed under the BSD 3-Clause
%   License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param

%% Simulation and stop conditions

time_sim        = t0:time_step:tf;                      % Simulation time [s]
voltage_limits  = [param.V_cutoff; param.V_cutover];    % Voltage window [V]

%% Initial conditions

y0 = spm_initial_conditions(SOC, initialState, input_density(0));

%% Spatial discretization

if param.EnableTemperature == 0
    [param.A, param.B] = spm_spatial_discretization(param.Dp, param.Dn);
end

%% Operating mode

switch param.OperatingMode
    case 1  % Constant current
        param.getCurrentDensity = input_density;
end

%% Model solver setup

stop_conditions = @(t,u) stop_conditions_fun(t, u, voltage_limits);
odefun          = @(t,u) spm_odefun_aging_new(t, u);

% Number of equations in the system
number_matrix = 2 * param.N + 2 + sum(param.aging_vector);

switch param.OperatingMode
    case {1, 2}
        if param.EnableSEI == 1 || param.EnableSEIoncracks == 1 || param.Enableplating == 1
            number_matrix = number_matrix + 1;  % jn
        end
        if param.EnableCEI == 1
            number_matrix = number_matrix + 1;  % jp
        end
end

if param.OperatingMode == 2
    number_matrix = number_matrix + 2;       % jn_tot and jp_tot
end

% Mass matrix initialization
M = diag(ones(number_matrix,1));

% Set algebraic rows to zero
switch param.OperatingMode
    case 1
        if param.EnableCEI == 1
            M(end,:) = 0; % jp
            if param.EnableSEI == 1 || param.EnableSEIoncracks == 1 || param.Enableplating == 1
                M(end-1,:) = 0; % jn
            end
        else
            if param.EnableSEI == 1 || param.EnableSEIoncracks == 1 || param.Enableplating == 1
                M(end,:) = 0; % jn
            end
        end

    case 2
        M(end-1,:) = 0;   % jn_tot
        M(end,:)   = 0;   % jp_tot

        if param.EnableCEI == 1
            M(end-2,:) = 0; % jp
            if param.EnableSEI == 1 || param.EnableSEIoncracks == 1 || param.Enableplating == 1
                M(end-3,:) = 0; % jn
            end
        else
            if param.EnableSEI == 1 || param.EnableSEIoncracks == 1 || param.Enableplating == 1
                M(end-2,:) = 0; % jn
            end
        end
end

% ODE options
opt = odeset( ...
    'Events', stop_conditions, ...
    'Mass', M, ...
    'RelTol', 1e-12, ...
    'AbsTol', 1e-12, ...
    'Refine', 5, ...
    'MassSingular', 'yes' ...
);

% Solve system
[time_sol, sol] = ode15s(odefun, time_sim, y0, opt);

%% Post-processing

t      = time_sol;  % Times for post-processing
result = spm_post_processing(sol, time_sol, t, linspace(0, 1, param.N + 1));

%% Simulation-end message 

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Simulation ended');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Reason:')

if time_sol(end) < time_sim(end)
    if ~isfield(result, 'current')
        disp('No saved results mode');
    else
        if result.current < 0
            switch param.OperatingMode
                case 1
                    disp('Battery voltage is above the cut-over Voltage');
                case 2
                    disp('Current is below the C/20 limit');
            end
        elseif result.current > 0
            disp('Battery voltage is below the cut-off Voltage');
        else
            disp('End of rest period');
        end
        disp('Battery SOC is');
        disp(result.SOC);
    end
end

end