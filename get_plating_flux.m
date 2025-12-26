function j_plating=get_plating_flux(ur,cs_n_avg,jn_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_ETA_PLATING (Beta version)
% This function computes the flux associated with lithium plating
% at the anode.
%
%       eta_plating = get_eta_plating(ur, jn_int)
%
% INPUTS
%   ur         : Unknown variable of the system
%   cs_n_avg   : Average anode concentration
%   jn_int     : Lithium intercalation flux at the anode
%
% OUTPUTS
%   j_plating  : Flux of lithium plated
%
%   Beta version: This function is under development and may be updated in
%   future releases.
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
eta_plating=get_eta_plating(ur,jn_int);

cs_n_surf=ur(2);
cs_maxn=param.cs_maxn;
kdy=(cs_n_surf-cs_n_avg*param.cs_maxn)/cs_maxn;

if jn_int<0
    % Get SEI flux
    j_plating=-param.k_plating*param.ce*kdy*exp(-0.65*param.F/(param.Rg*param.T)*eta_plating);
else
    j_plating=0;
end

end