function [Dp_eff,Dn_eff]= get_EffectiveDiffusion_coefficient(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_EFFECTIVE_DIFFUSION_COEFFICIENT
% This function computes the temperature-dependent effective diffusion 
% coefficients for the cathode and anode based on Arrhenius-type relations.
%
%       [Dp_eff, Dn_eff] = get_EffectiveDiffusion_coefficient(T)
%
% INPUTS
% T        : Temperature [K]
%
% OUTPUTS
% Dp_eff   : Effective diffusion coefficient in the cathode [m^2/s]
% Dn_eff   : Effective diffusion coefficient in the anode [m^2/s]
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

Dp_eff=param.Dp*exp(-param.Ea_Dp/param.Rg*(1/T-1/param.Tref));
Dn_eff=param.Dn*exp(-param.Ea_Dn/param.Rg*(1/T-1/param.Tref));

end