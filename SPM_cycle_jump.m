function [initialState,cycle_count]=SPM_cycle_jump(delta_N,result,cycle_count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_CYCLE_JUMP
% Computes the optimal cycle jump size and updates the initial state for the
% next simulation segment, based on trends in degradation-related variables.
%
%       [initialState, cycle_count] = SPM_cycle_jump(delta_N, result, cycle_count)
%
% The cycle jump ΔN is estimated as:
%
%       ΔN = tol_jump_size * |Y_N - Y_{N-1}| / |Y_N - 2Y_{N-1} + Y_{N-2}|
%
% where Y is a degradation-related variable (e.g. L_sei, L_cei, a, etc.)
%
% The jump size is computed for each degradation variable. Then, the cycle
% jump is selected as the minimum among all ΔN computed. The new state
% vector is estimated by linear extrapolation. 
%
% INPUT:
%       delta_N         - Estimated cycle jump (can be overridden)
%       result          - Simulation result structure of the 3 previous cycles
%       cycle_count     - Current number of simulated cycles
%
% OUTPUT:
%       initialState    - Updated initial state for the next simulation
%       cycle_count     - Updated number of simulated cycles after jump
%
%   This file is part of POLIDEMO software
%
% 	Official GitHUB: 	https://github.com/Poli-ON/POLIDEMO
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
%   POLIDEMO is a free Matlab-based software distributed with BSD 3-Clause
%   license.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param

%% Extract relevant states from the last cycles

% Initial state
initial_state_N   = result(end,   3).initialState;
initial_state_N_1 = result(end-1, 3).initialState;

% L_sei
L_sei_N   = result(end,   3).L_sei(end);
L_sei_N_1 = result(end-1, 3).L_sei(end);
L_sei_N_2 = result(end-2, 3).L_sei(end);

% L_cei
L_cei_N   = result(end,   3).L_cei(end);
L_cei_N_1 = result(end-1, 3).L_cei(end);
L_cei_N_2 = result(end-2, 3).L_cei(end);

% L_plating
L_plating_N   = result(end,   3).L_plating(end);
L_plating_N_1 = result(end-1, 3).L_plating(end);
L_plating_N_2 = result(end-2, 3).L_plating(end);

% Stress
sigma_c_p_max_N   = result(end,   3).param.sigma_c_p_max;
sigma_c_p_max_N_1 = result(end-1, 3).param.sigma_c_p_max;
sigma_c_n_max_N   = result(end,   3).param.sigma_c_n_max;
sigma_c_n_max_N_1 = result(end-1, 3).param.sigma_c_n_max;
sigma_c_p_min_N   = result(end,   3).param.sigma_c_p_min;
sigma_c_p_min_N_1 = result(end-1, 3).param.sigma_c_p_min;
sigma_c_n_min_N   = result(end,   3).param.sigma_c_n_min;
sigma_c_n_min_N_1 = result(end-1, 3).param.sigma_c_n_min;

% Crack length
a_N   = result(end,   3).a(end);
a_N_1 = result(end-1, 3).a(end);
a_N_2 = result(end-2, 3).a(end);

% Volume fraction
volume_fraction_N   = result(end,   3).param.vol_fraction_solidphase;
volume_fraction_N_1 = result(end-1, 3).param.vol_fraction_solidphase;

%% Compute the jump size candidates

delta_N1 = param.tol_jump_size * abs(L_sei_N   - L_sei_N_1)   / abs(L_sei_N   - 2*L_sei_N_1   + L_sei_N_2);
delta_N2 = param.tol_jump_size * abs(L_cei_N   - L_cei_N_1)   / abs(L_cei_N   - 2*L_cei_N_1   + L_cei_N_2);
delta_N3 = param.tol_jump_size * abs(L_plating_N - L_plating_N_1) / abs(L_plating_N - 2*L_plating_N_1 + L_plating_N_2);
delta_N4 = 10*param.tol_jump_size * abs(a_N - a_N_1) / abs(a_N - 2*a_N_1 + a_N_2);
delta_N5 = param.tol_jump_size * abs(volume_fraction_N(1)) ./ abs(volume_fraction_N(1) - volume_fraction_N_1(1));
delta_N6 = param.tol_jump_size * abs(volume_fraction_N(2)) ./ abs(volume_fraction_N(2) - volume_fraction_N_1(2));

% Select the minimum jump after 50 cycles
if param.cycle_count > 50
    delta_N = min([delta_N1, delta_N2, delta_N3, delta_N4, delta_N5, delta_N6]);
end

% Apply constraints on minimum jump size
delta_N = ceil(delta_N);
if delta_N < 5
    delta_N = 5;
end

%% Update cycle count and initial state

cycle_count = cycle_count + delta_N;
initialState = initial_state_N + delta_N * (initial_state_N - initial_state_N_1);  % First-order extrapolation

%% Update mechanical degradation quantities if enabled

if param.EnableSEIoncracks == 1 || param.EnableLAM == 1
    param.sigma_c_p_max = sigma_c_p_max_N + delta_N * (sigma_c_p_max_N - sigma_c_p_max_N_1);
    param.sigma_c_n_max = sigma_c_n_max_N + delta_N * (sigma_c_n_max_N - sigma_c_n_max_N_1);
    param.sigma_c_p_min = sigma_c_p_min_N + delta_N * (sigma_c_p_min_N - sigma_c_p_min_N_1);
    param.sigma_c_n_min = sigma_c_n_min_N + delta_N * (sigma_c_n_min_N - sigma_c_n_min_N_1);

    delta_sigma_p = param.sigma_c_p_max - param.sigma_c_p_min;
    delta_sigma_n = param.sigma_c_n_max - param.sigma_c_n_min;
end

if param.EnableLAM == 1
    vol_fraction_solidphase_p = -1 / (param.n_LAM_p + 1) * ...
        (param.a_LAM_p * delta_sigma_p^param.m_LAM_p) * (cycle_count)^(param.n_LAM_p + 1) + ...
        param.vol_fraction_solidphase_0(1);

    vol_fraction_solidphase_n = -1 / (param.n_LAM_n + 1) * ...
        (param.a_LAM_n * delta_sigma_n^param.m_LAM_n) * (cycle_count)^(param.n_LAM_n + 1) + ...
        param.vol_fraction_solidphase_0(2);

    param.vol_fraction_solidphase = [vol_fraction_solidphase_p, vol_fraction_solidphase_n];
    param.as_p = 3 * param.vol_fraction_solidphase(1) / param.Rp;
    param.as_n = 3 * param.vol_fraction_solidphase(2) / param.Rn;
end

if param.EnableSEIoncracks == 1
    switch param.CrackMechanism
        case 1
            param.N_cracks = param.k_crack * ...
                (1 / (param.n_LAM_n + 1) * ...
                (param.a_LAM_n * delta_sigma_n^param.m_LAM_n) * ...
                (cycle_count)^(param.n_LAM_n + 1) + ...
                param.vol_fraction_solidphase_0(2));
        case 2
            param.A_crack = param.k_crack * ...
                ((-param.vol_fraction_solidphase(2) + param.vol_fraction_solidphase_0(2)) / ...
                 param.vol_fraction_solidphase_0(2));
    end
end

end








            



