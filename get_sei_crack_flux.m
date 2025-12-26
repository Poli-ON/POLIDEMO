function j_crack_sei = get_sei_crack_flux(ur,L_sei_crack,jn_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_SEI_CRACK_FLUX
% This function computes the lithium flux associated with SEI growth through
% cracks in the electrode surface, based on the local SEI overpotential and
% transport properties.
%
%       j_crack_sei = get_sei_crack_flux(ur, L_sei_crack, jn_int)
%
% INPUTS
% ur            : Vector of state variables at the current spatial location
% L_sei_crack   : Effective SEI crack thickness
% jn_int        : Intercalation current density
%
% OUTPUTS
% j_crack_sei   : SEI growth flux through cracks
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

% Get SEI overpotential
eta_sei=get_eta_sei(ur,jn_int);

% Get SEI flux
j_crack_sei=-param.c_solv_electrolyte/(L_sei_crack/param.D_solv+1/param.k_sei*1/(exp(-0.5*param.F/(param.Rg*param.T)*eta_sei)));

end
