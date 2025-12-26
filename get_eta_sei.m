function eta_sei=get_eta_sei(ur,jn_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_ETA_SEI
% This function computes the overpotential associated with the SEI 
% (Solid Electrolyte Interphase) formation reaction at the anode.
%
%       eta_sei = get_eta_sei(ur, jn_int)
%
% INPUTS
% ur       : Unknown variable of the system
% jn_int   : Lithium intercalation flux at the anode
%
% OUTPUTS
% eta_sei  : SEI overpotential
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
kn=param.kn;                    

% Surface lithium concentration [mol/m^3] - Cathode
cs_n = ur(2,:);                             

% Exchange current density [A/m^2]
i0_n = kn*param.F*sqrt(param.ce)*sqrt(cs_n(end))*sqrt(param.cs_maxn-cs_n(end)) ;

% Anode overpotential [V]
eta_n = 2*param.Rg*param.T/param.F*asinh(jn_int*param.F/(2*i0_n));

% OCV [V]
[U_p,U_n] =  openCircuitPotential(ur);

eta_sei = eta_n+U_n-param.U_sei;
end