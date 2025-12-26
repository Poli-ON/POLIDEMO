function dadt= crack_propagation(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRACK_PROPAGATION
% Computes the crack growth rate based on the Paris law using the stress
% intensity factor. Designed for use in the SPM aging model.
%
%       dadt = crack_propagation(a)
%
% INPUT:
%   a - current crack length [m] 
%
% OUTPUT:
%   dadt - crack growth rate [m/s] 
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

% Compute stress intensity factors at max and min loading
[K_n_max, K_n_min] = stress_intensity_factor_computation(a);

% Effective stress intensity factor range
delta_K = max(0, K_n_max - K_n_min);

% Crack growth rate using Paris law
dadt=param.C/param.tn*(delta_K).^param.m;

end