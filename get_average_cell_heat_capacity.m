function [mass_cell,cp_cell] = get_average_cell_heat_capacity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_AVERAGE_CELL_HEAT_CAPACITY - BETA VERSION
% This function computes the total mass and the average specific heat capacity
% of the pouch cell, to be used in the lumped thermal balance equation.
%
% OUTPUTS
% mass_cell    : Total mass of the cell
% cp_cell      : Average specific heat capacity of the cell
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

% The volume and the mass of the tab are ignored.

% Volume fraction occupied in electrodes and separator
vol_fraction_pos = param.vol_fraction_solidphase(1); 
vol_fraction_sep = (1 - param.eps_s); 
vol_fraction_neg = param.vol_fraction_solidphase(2); 

% Volume [m^3]
volume_al = param.pouch_area*ceil(0.5*param.no_of_layers)*param.len_al;
volume_pos = param.pouch_area*param.len_p*vol_fraction_pos*param.no_of_layers;
volume_sep = param.pouch_area*param.len_s*vol_fraction_sep*param.no_of_layers;
volume_neg = param.pouch_area*param.len_n*vol_fraction_neg*param.no_of_layers;
volume_cu = param.pouch_area*ceil(0.5*(param.no_of_layers+1))*param.len_cu;
volume_LiPF6 = param.pouch_area*(param.len_p*param.eps_p + param.len_s*param.eps_s + param.len_n*param.eps_n)*param.no_of_layers;
%volume_pouch = param.pouch_length*param.pouch_width*param.len_pouch;

% Mass [kg]
mass_al = param.rho_al*volume_al;
mass_pos = param.rho_p*volume_pos;
mass_sep = param.rho_s*volume_sep;
mass_neg = param.rho_n*volume_neg;
mass_cu = param.rho_cu*volume_cu;
mass_LiPF6 = param.rho_LiPF6*volume_LiPF6;
%mass_pouch = 2*param.rho_pouch*volume_pouch; % factor of two, for each of the top & bottom pouch covers
%mass_cell =  mass_al + mass_pos + mass_sep + mass_neg + mass_cu + mass_LiPF6 + mass_pouch;
mass_cell =  mass_al + mass_pos + mass_sep + mass_neg + mass_cu + mass_LiPF6;

% Cp of pouch is currently not available and not accounted for in
% calculation of Cp_avg. Also, note that tab material is not accounted for in Cp_avg calcs.
% Are these important?
% Commento fra: ho tolto cp_pouch

% Specific heat capacity, J/(kg K)
%cp_cell = (param.cp_al*mass_al + param.cp_p*mass_pos + param.cp_s*mass_sep + param.cp_n*mass_neg + param.cp_cu*mass_cu + param.cp_LiPF6*mass_LiPF6)/mass_cell;
mass_cell=360/1000;
%cp_cell = (param.cp_al*mass_al + param.cp_p*mass_pos + param.cp_s*mass_sep + param.cp_n*mass_neg + param.cp_cu*mass_cu + param.cp_LiPF6*mass_LiPF6)/mass_cell;  

cp_cell=(param.cp_al*param.len_al+param.cp_p*param.len_p+param.cp_s*param.len_s+param.cp_n*param.len_n+param.cp_cu*param.len_cu)/(param.len_al+param.len_p+param.len_s+param.len_n+param.len_cu);

end