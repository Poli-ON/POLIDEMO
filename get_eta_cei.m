function eta_cei=get_eta_cei(ur,jp_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_ETA_CEI
% This function computes the cathode electrolyte interphase (CEI) overpotential
% based on the local electrochemical state variables.
%
%       eta_sei = get_eta_cei(ur, jp_int)
%
% INPUTS
% ur        : Vector of state variables at the current spatial location
% jp_int    : Intercalation current density at the positive electrode
%
% OUTPUTS
% eta_cei   : SEI overpotential
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

% Reaction rate constant [m^2.5/(mol^0.5 s)] - Anode
kp=param.kp;                    

% Surface lithium concentration [mol/m^3] - Cathode
cs_p = ur(1,:);                             

% Exchange current density [A/m^2]
i0_p = kp*param.F*sqrt(param.ce)*sqrt(cs_p(end))*sqrt(param.cs_maxp-cs_p(end)) ;

% Anode overpotential [V]
eta_p = 2*param.Rg*param.T/param.F*asinh(jp_int*param.F/(2*i0_p));

% OCV [V]
[U_p,U_n] =  openCircuitPotential(ur);

eta_cei = eta_p+U_p-param.U_cei;
end