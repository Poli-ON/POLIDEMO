function [kp_eff,kn_eff]= get_Effectivek_coefficient(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_EFFECTIVEK_COEFFICIENT
% This function computes the temperature-dependent effective reaction rate 
% coefficients for the cathode and anode based on Arrhenius-type relations.
%
%       [kp_eff, kn_eff] = get_Effectivek_coefficient(T)
%
% INPUTS
% T        : Temperature [K]
%
% OUTPUTS
% kp_eff   : Effective reaction rate coefficient in the cathode
% kn_eff   : Effective reaction rate coefficient in the anode
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

kp_eff=param.kp*exp(-param.Ea_kp./param.Rg.*(1./T-1/param.Tref));
kn_eff=param.kn.*exp(-param.Ea_kn./param.Rg.*(1./T-1/param.Tref));

end