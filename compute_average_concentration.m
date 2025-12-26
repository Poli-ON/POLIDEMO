function [cs_p_avg,cs_n_avg] = compute_average_concentration(cs_p,cs_n,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE_AVERAGE_CONCENTRATION
% This function computes the average lithium concentration in anode and
% cathode particles at a given time.
%
%       [cs_p_avg, cs_n_avg] = compute_average_concentration(cs_p, cs_n, x, param)
%
% INPUTS
% cs_p    : Lithium concentration in the cathode
% cs_n    : Lithium concentration in the anode
% x       : Spatial discretization
%
% OUTPUTS
% cs_p_avg : Average lithium concentration in the cathode
% cs_n_avg : Average lithium concentration in the anode
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

% Particle radius [m]
Rp = param.Rp;                                                   % Cathode
Rn = param.Rn;                                                   % Anode

% Average concentration
cs_p_avg = trapz(x*Rp,cs_p.*(x*Rp).^2)/trapz(x*Rp,(x*Rp).^2)/param.cs_maxp;   % Cathode
cs_n_avg =trapz(x*Rn,cs_n.*(x*Rn).^2)/trapz(x*Rn,(x*Rn).^2)/param.cs_maxn;    % Anode



end

