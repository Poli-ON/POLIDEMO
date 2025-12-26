function [K_max,K_min]=stress_intensity_factor_computation(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRESS_INTENSITY_FACTOR_COMPUTATION
% Computes the stress intensity factors for a particle with superficial 
% cracks. Designed for use in the SPM aging model and the Paris law. More
% details for the computation method can be found in this reference paper: 
% https://doi.org/10.1016/j.engfracmech.2023.109597
%
%       [K_max, K_min] = stress_intensity_factor_computation(a)
%
% INPUT:
%   a - current crack length [m] 
%
% OUTPUT:
%   K_max - maximum stress intensity factor [MPa*sqrt(m)]
%   K_min - minimum stress intensity factor [MPa*sqrt(m)] (currently set to 0)
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

% Geometric factors for particle with superficial cracks
p = [1.2231, 1.2231, -0.2373, -0.1111, -0.1440, -0.2040, -0.1500];
q = [0.1864, 0.4987,  0.5204, 0.3367, 0.3360, 0.3565, 0.3114];
r = [1.0210, 0.5692, 0.4305, 0.3833, 0.3266, 0.2828, 0.2567];

% Normalized crack length
a_frac_R = a / param.Rn;

% Geometric factor Y
Y = p .* (a_frac_R).^2 + q .* a_frac_R + r;

% Discretization points for fitting stress distribution
x = linspace(0, 1, param.N+1);

% Fit stress distribution
sigma_c_n_max_distribution_fitted=polyfit(x, param.sigma_c_n_max_distribution_crack,6);
sigma_c_n_max_distribution_fitted=fliplr(sigma_c_n_max_distribution_fitted);

% Prepare powers of a_frac_R for polynomial evaluation
a_frac_R_vector=[1,a_frac_R,a_frac_R^2,a_frac_R^3,a_frac_R^4,a_frac_R^5,a_frac_R^6];

% Compute stress intensity factor
K_max= sum(sigma_c_n_max_distribution_fitted/(10^6).*(1-a_frac_R_vector).*Y)*sqrt(a); %MPa * sqrt(m)
% K_min= sum(sigma_c_n_min_distribution_fitted/(10^6).*(1-a_fraq_R_vector).*Y)*sqrt(a);
K_min=0; % Ingnore compression 

end