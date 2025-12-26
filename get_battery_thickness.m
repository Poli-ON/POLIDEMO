function thk = get_battery_thickness(x,cs_p,cs_n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_BATTERY_THICKNESS
% This function computes the reversible thickness change of the battery at 
% a given time based on the lithium concentrations in the electrodes.
%
%       thk = get_battery_thickness(x, cs_p, cs_n)
%
% INPUTS
% cs_p    : Lithium concentration in the cathode at a given time
% cs_n    : Lithium concentration in the anode at a given time
% x       : Radial coordinate
%
% OUTPUTS
% thk     : Reversible thickness change of the battery at the given time
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

% For pouch cell
Ccase=1;

% Normalized concentration
cs_p_norm=cs_p/param.cs_maxp;
cs_n_norm=cs_n/param.cs_maxn;

% Get Omega*dC (Integration of partial molar volume)
OmdC_pf=param.OmdC_pf;
OmdC_nf=param.OmdC_nf;

% Displacement
u_p=param.Rp.*trapz(x,OmdC_pf(cs_p_norm).*x.^2);
u_n=param.Rn.*trapz(x,OmdC_nf(cs_n_norm).*x.^2);

% Volumetric deformation
epsv_p=(3*param.Rp^2.*u_p+3*param.Rp.*u_p.^2+u_p.^3)/param.Rp^3;
epsv_n=(3*param.Rn^2.*u_n+3*param.Rn.*u_n.^2+u_n.^3)/param.Rn^3;

% Reversible thickness change [mm]
thk=1000*Ccase*param.no_of_layers.*(param.vol_fraction_solidphase(1).*param.len_p.*epsv_p +...
            param.vol_fraction_solidphase(2).*param.len_n.*epsv_n);



end


