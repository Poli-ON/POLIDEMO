function [dTdt] = spm_lumped_thermal_model(t,ur,I,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_LUMPED_THERMAL_MODEL
% This function computes the lumped thermal balance of the cell and returns
% the time derivative of the temperature based on the applied current and
% internal heat generation.
%
%       dTdt = SPM_lumped_thermal_model(t, u, I, T)
%
% INPUTS
% t             Time
% u             State vector of the model
% I             Applied current density
% T             Cell temperature
%
% OUTPUTS
% dTdt          Time derivative of the cell temperature
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


V=get_Voltage(t,ur,I,T);

% OCV [V]
[U_p,U_n] = openCircuitPotential([ur(1,:);ur(2,:)],T);
[dU_pdT,dU_nT] = openCircuitPotentialDerivative([ur(1,:);ur(2,:)]);

Q_gen=abs(I)*param.pouch_area*param.no_of_layers*abs((V-(U_p-U_n)))-I*param.pouch_area*param.no_of_layers*T*(dU_pdT-dU_nT);
Q_loss=param.h_conv*param.A_conv*(T-param.Tref);

dTdt=(Q_gen-Q_loss)/(param.mass_cell*param.cp_cell);

end