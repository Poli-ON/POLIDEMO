function SOC = compute_SOC(cs_p_avg,cs_n_avg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE_SOC
% This function computes the state of charge (SOC) of the battery.
%
%       SOC = compute_soc(cs_p_avg, cs_n_avg, param)
%
% INPUTS
%   cs_p_avg : Average lithium concentration in the cathode
%   cs_n_avg : Average lithium concentration in the anode
%   param    : Structure containing the parameters of the model
%
% OUTPUTS
%   SOC      : State of charge of the battery
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

% SOC
SOC=(cs_n_avg-param.theta_min_neg)/(param.theta_max_neg-param.theta_min_neg)*100;

end