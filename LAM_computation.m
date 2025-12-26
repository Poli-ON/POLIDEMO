function [LAM_p,LAM_n] = LAM_computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAM_COMPUTATION
% This function computes the loss of active material (LAM) in the electrodes
% based on the solutions from the Single Particle Model (SPM).
%
%       LAM = LAM_computation
%
% OUTPUTS
%   LAM : Loss of active material in the electrodes [%]
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


%% LAM Computation

global param

% Fresh state
Cp_0=param.vol_fraction_solidphase_0(1)*param.len_p*param.cs_maxp*param.F*param.pouch_area*param.no_of_layers/3600;
Cn_0=param.vol_fraction_solidphase_0(2)*param.len_n*param.cs_maxn*param.F*param.pouch_area*param.no_of_layers/3600;

% Aged state
Cp=param.vol_fraction_solidphase(1)*param.len_p*param.cs_maxp*param.F*param.pouch_area*param.no_of_layers/3600;
Cn=param.vol_fraction_solidphase(2)*param.len_n*param.cs_maxn*param.F*param.pouch_area*param.no_of_layers/3600;

LAM_p=(Cp_0-Cp)/Cp_0*100;
LAM_n=(Cn_0-Cn)/Cn_0*100;


end
