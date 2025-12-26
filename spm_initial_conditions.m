function [initial_SPM] = spm_initial_conditions(SOC,initialState,current)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_INITIAL_CONDITIONS
% Determines the initial conditions for the SPM equations based on the 
% user-defined parameters and the initial state of charge (SOC).
%
%       initial_SPM = spm_initial_conditions(SOC, initialState, current)
%
% INPUT:
%   SOC              - Initial state of charge
%   initialState     - Vector of user-defined initial states (optional)
%   current          - Current at time t = 0
%
% OUTPUT:
%   initial_SPM      - Vector containing the initial conditions for the SPM
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

%% No user-defined initial state

if isempty(initialState)
    
    param.init_cell_soc = SOC / 100;  % SOC conversion into fraction [0,1]

    % Positive electrode initial concentration
    param.cs_p_init = (param.init_cell_soc * ...
        (param.theta_max_pos - param.theta_min_pos) + ...
        param.theta_min_pos) * param.cs_maxp;

    % Negative electrode initial concentration
    param.cs_n_init = (param.init_cell_soc * ...
        (param.theta_max_neg - param.theta_min_neg) + ...
        param.theta_min_neg) * param.cs_maxn;

    % Base state: concentrations (+ temperature if enabled)
    if param.EnableTemperature == 0
        initial_SPM = [ ...
            param.cs_p_init * ones(param.N+1,1); ...
            param.cs_n_init * ones(param.N+1,1)];
    else
        initial_SPM = [ ...
            param.cs_p_init * ones(param.N+1,1); ...
            param.cs_n_init * ones(param.N+1,1); ...
            param.Tref];
    end

    % Aging variables
    initial_state_aging = [];

    if param.EnableSEI == 1
        % SEI layer thickness
        initial_state_aging = [initial_state_aging; 0];

        % SEI on cracks
        if param.EnableSEIoncracks == 1
            switch param.CrackMechanism
                case {0,1}
                    A_crack = 2 * param.a_0 * param.N_cracks * param.w_crack;
                    initial_state_aging = [initial_state_aging; param.a_0];
                case 2
                    A_crack = param.A_crack;
            end
        else
            A_crack = 0;
        end
    else
        A_crack = 0;
    end

    if param.EnableCEI == 1
        initial_state_aging = [initial_state_aging; 0];
    end

    if param.Enableplating == 1
        initial_state_aging = [initial_state_aging; 0];
    end

    % Add fluxes depending on the operating mode
    if param.OperatingMode == 1
        if param.EnableSEI == 1 || param.Enableplating == 1
            jn0 = current / (param.F * (param.as_n * param.len_n + A_crack));
            initial_state_aging = [initial_state_aging; jn0];
        end
        if param.EnableCEI == 1
            jp0 = -current / (param.F * (param.as_p * param.len_p));
            initial_state_aging = [initial_state_aging; jp0];
        end
    else
        jn0 = current / (param.F * (param.as_n * param.len_n + A_crack));
        jp0 = -current / (param.F * (param.as_p * param.len_p));
        initial_state_aging = [initial_state_aging; jn0; jp0];
    end

    % Append aging state to base state
    initial_SPM = [initial_SPM; initial_state_aging];

else
%% User-defined initial state
    switch param.OperatingMode
        case 1
            initial_SPM = initialState;

        case 2
            if param.EnableSEI == 1 && param.EnableSEIoncracks == 1
                switch param.CrackMechanism
                    case 0
                        idx_a0 = 2*param.N + 2 + param.aging_vector(1) + param.aging_vector(2);
                        A_crack = 2 * initialState(idx_a0) * param.N_cracks * param.w_crack;
                    case 1
                        idx_a0 = 2*param.N + 2 + param.aging_vector(1) + param.aging_vector(2);
                        A_crack = 2 * initialState(idx_a0) * param.N_cracks * param.w_crack;
                    case 2
                        A_crack = param.A_crack;
                end
            else
                A_crack = 0;
            end

            jn0 = current / (param.F * (param.as_n * param.len_n + A_crack));
            jp0 = -current / (param.F * param.as_p * param.len_p);
            initial_SPM = [initialState; jn0; jp0];
    end
end

end
