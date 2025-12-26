function [dU_pdT,dU_ndT]=openCircuitPotentialDerivative(ur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPENCIRCUITPOTENTIALDERIVATIVE
% Computes the entropic contribution, i.e., the derivative of the Open 
% Circuit Voltage (OCV) with respect to temperature for cathode and anode.
%
%       [dU_pdT, dU_ndT] = openCircuitPotentialDerivative(ur)
%
% INPUTS:
%   ur  - Surface lithium concentrations [mol/m^3], ur = [cs_p; cs_n]
%
% OUTPUTS:
%   dU_pdT - Derivative of cathode OCV with respect to temperature [V/K]
%   dU_ndT - Derivative of anode OCV with respect to temperature [V/K]
%   This file is part of the POLIDEMO software.
%
%   Official website:  [insert link]
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

% Lithium concentration [mol/m^3]
cs_p = ur(1,:);                      % Cathode
cs_n = ur(2,:);                      % Anode

%  Cathode
theta_p  = cs_p(end,:)./param.cs_maxp;
fit=[-6 -4 -2 2 5 6 -2 -2.5 -2.75 -3 -3.25 -3.5 -3.75 -4 -4.5 -5 -5.5 -6 -6 -6 -6].*10^-4; xfit=0.5:0.025:1;

dU_pdT=polyval(polyfit(xfit,fit,8),theta_p);

% Anode
theta_n  = cs_n(end,:)/param.cs_maxn;
p=[43284.7958739547,-311048.624060859,1023967.82305701,-2047611.40752125,...
    2779375.46621321,-2708548.72079468,1952437.13006254,-1055395.81834479,...
    428388.531061457,-129255.905213853,28333.5277851283,-4334.62976872029,...
    432.002906451329,-24.5563150153189,0.545059970471012,0.00498143840288525,...
    -8.04128007713281e-06];

dU_ndT=theta_n*0;

N=16;
for i=0:N
    dU_ndT=dU_ndT+p(i+1).*theta_n.^(N-i);
end

end