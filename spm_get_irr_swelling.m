%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_GET_IRR_SWELLING
% Computes the irreversible swelling thickness in a Li-ion pouch cell
% due to SEI growth, crack propagation, and gas generation effects.
%
%       irr_swelling = spm_get_irr_swelling(vol_fraction_solid_phase_n, L_sei, a)
%
% The irreversible swelling is calculated as the combination of:
%
%   - SEI-induced volume expansion
%   - Gas generation and accumulation inside the pouch
%   - Electrolyte consumption
%
% Depending on the selected crack mechanism (param.CrackMechanism),
% different degradation laws are applied to estimate SEI volume
% and crack growth contribution.
%
% INPUT:
%       vol_fraction_solid_phase_n - Volume fraction of the solid phase in
%                                    the negative electrode
%       L_sei                      - SEI layer thickness
%       a (optional)               - Crack geometry factor (only required
%                                    for mechanisms involving crack opening)
%
% OUTPUT:
%       irr_swelling               - Total irreversible swelling thickness [µm]
%
%   This file is part of POLIDEMO software
%
%   Official GitHUB:    https://github.com/Poli-ON/POLIDEMO
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
%   POLIDEMO is a free Matlab-based software distributed with BSD 3-Clause
%   license.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function irr_swelling = spm_get_irr_swelling(vol_fraction_solid_phase_n,L_sei,varargin)

global param

if nargin < 3 || isempty(varargin{1})
    a = []; 
else
    a = varargin{1};
end

switch param.CrackMechanism
    case 0
        epsv_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei/(param.pouch_area*param.len_n);
        n_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei*param.no_of_layers/param.V_sei;
    case 1
        epsv_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei/(param.pouch_area*param.len_n);
        n_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei*param.no_of_layers/param.V_sei;
    case 2
        epsv_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+param.A_crack)*L_sei/(param.pouch_area*param.len_n);
        n_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+param.A_crack)*L_sei*param.no_of_layers/param.V_sei;
end

param.kz=5284.89604*10^(-6);

switch param.CrackMechanism
    case 0
        volume_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei*param.no_of_layers;
    case 1
        volume_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+2*a*param.w_crack*param.N_cracks)*L_sei*param.no_of_layers;
    case 2
        volume_sei=(3*vol_fraction_solid_phase_n/param.Rn*param.len_n*param.pouch_area+param.A_crack)*L_sei*param.no_of_layers;
end

volume_electrolyte=n_sei*2*6.67e-5;
b_eqv=volume_electrolyte/(param.pouch_area)+param.kz*(101300)*param.pouch_thickness;
c_eqv=-param.kz*n_sei*param.Rg*param.T*param.pouch_thickness/(param.pouch_area)+...
    param.kz*(101300)*param.pouch_thickness*volume_electrolyte/param.pouch_area;
thk_gas=1000*(-b_eqv+sqrt(b_eqv^2-4*c_eqv))/2;
volume_gas=param.pouch_area*thk_gas+volume_electrolyte;

thk_sei=1000*epsv_sei*param.no_of_layers*param.len_n;
thk_electrolyte=1000*volume_electrolyte/param.pouch_area;

irr_swelling=thk_sei+thk_gas-thk_electrolyte;

