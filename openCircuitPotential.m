function [U_p,U_n] = openCircuitPotential(ur,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPENCIRCUITPOTENTIAL
% Computes the Open Circuit Voltage (OCV) of the cathode and anode for a
% given battery chemistry and state of charge. It can include temperature
% dependence if enabled.
%
%       [U_p, U_n] = openCircuitPotential(ur, T)
%
% This function includes:
%   - Cathode OCV
%   - Anode OCV
%   - Temperature correction (optional)
%
% INPUTS:
%   ur  - Vector of surface lithium concentrations [mol/m^3], ur = [cs_p; cs_n]
%   T   - Temperature [K]
%
% OUTPUTS:
%   U_p - Cathode OCV [V]
%   U_n - Anode OCV [V]
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

% Spatial discretization
N = param.N;

% Lithium concentrations
cs_p = ur(1,:);  % Cathode
cs_n = ur(2,:);  % Anode

% Normalized surface concentrations (state of charge)
theta_p = cs_p./ param.cs_maxp;
theta_n = cs_n./ param.cs_maxn;

U_p = 4.04596+exp(-42.30027*theta_p+16.56714)-0.048.*atan(50.01833*theta_p-26.48897)-...
            0.05447.*atan(18.99678*theta_p-12.32362)-exp(78.240895*theta_p-78.68074);

if abs(param.Crate)>0.5 % High current
    U_n=0.6379+0.5416*exp(-305.5309*theta_n)+0.044*tanh(-(theta_n-0.1958)/0.1088)-...
        0.1978*tanh((theta_n-1.0571)/0.0854)-0.6875*tanh((theta_n+0.0117)/0.0529)-...
        0.0175*tanh((theta_n-0.5692)/0.0875);
else                    % Low current
    U_n=0.124+1.5*exp(-150.*theta_n)+0.0155*tanh((theta_n-0.205)./0.029)-...
        0.011*tanh((theta_n-0.124)./0.0226)-0.102*tanh((theta_n-0.194)./0.142)+...
        0.0347*tanh((theta_n-0.286)./0.083)-0.0147*tanh((theta_n-0.5)./0.034)-...
        0.0045*tanh((theta_n-0.9)./0.119)-0.022*tanh((theta_n-0.98)./0.0164)-...
        0.035*tanh((theta_n-0.99)./0.05);
end

% Temperature correction
if param.EnableTemperature == 1
    [dU_pdT,dU_ndT] = openCircuitPotentialDerivative(ur);
    U_p=U_p+(T-param.Tref).*dU_pdT;
    U_n=U_n+(T-param.Tref).*dU_ndT;

end


end