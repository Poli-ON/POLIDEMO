function [I] = get_Current_density(t,ur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_CURRENT_DENSITY
% Computes the applied current density based on the simulation operating mode
% (constant current or constant voltage).
%
%       I = get_Current_density(t, ur)
%
% INPUT:
%   t       - Current time [s]
%   ur      - State vector (contains unknown variables depending on the
%             operating mode)
%
% OUTPUT:
%   I       - Applied current density [A/m²]
%
%   Operating modes:
%     - Mode 1 (Constant Current): current density is provided directly by
%       param.getCurrentDensity(t)
%     - Mode 2 (Constant Voltage): current density is derived from the last
%       state variable (reaction current density at the positive electrode)
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

switch param.OperatingMode
    case 1  % Constant Current (CC)
        I = param.getCurrentDensity(t);

    case 2  % Constant Voltage (CV)
        I = -ur(end,:) * (param.F * (param.as_p * param.len_p));
end

end