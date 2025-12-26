function j_sei=get_sei_flux(ur,L_sei,jn_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_SEI_FLUX
% This function calculates the flux of the Solid Electrolyte Interphase (SEI) 
% formation reaction at the anode.
%
%       j_sei=get_sei_flux(ur,L_sei,jn_int)
%
% INPUTS
% ur       : Unknown variable of the system
% jn_int   : Lithium intercalation flux at the anode
% L_sei    : Thickness of the SEI layer
%
% OUTPUTS
% j_sei    : SEI reaction flux
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
j_sei=-param.nSEI*param.c_solv_electrolyte./(L_sei./param.D_solv+1./param.k_sei*1./(exp(-0.5*param.nSEI*param.F./(param.Rg*param.T).*eta_sei)));
end