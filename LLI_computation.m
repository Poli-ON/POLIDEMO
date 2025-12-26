function LLI = LLI_computation(cs_p_avg,cs_n_avg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LLI_COMPUTATION
% This function computes the loss of lithium inventory (LLI) in the battery
% based on the solutions from the Single Particle Model (SPM).
%
%       LLI = LLI_computation(cs_p_avg, cs_n_avg)
%
% INPUTS
%   cs_p_avg : Average lithium concentration in the cathode
%   cs_n_avg : Average lithium concentration in the anode
%
% OUTPUTS
%   LLI : Loss of lithium inventory [%]
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


%% LLI Computation

global param

% Computation of new Cp and Cn
Cp=param.vol_fraction_solidphase(1)*param.len_p*param.cs_maxp*param.F*param.pouch_area/3600;
Cn=param.vol_fraction_solidphase(2)*param.len_n*param.cs_maxn*param.F*param.pouch_area/3600;

theta_p=cs_p_avg;
theta_n=cs_n_avg;

param.theta_p=theta_p;
param.theta_n=theta_n;

n_Li=(theta_p*Cp+theta_n*Cn)*3600/(param.F); 

 if param.cycle_count ==1
     param.n_Li_initial=n_Li;
 end

%LLI=(param.n_Li_initial-n_Li)/param.n_Li_initial*100;

% n_loss1=2*(param.as_n*param.len_n*param.pouch_area+...
%     2*a*param.w_crack*param.N_cracks)*L_sei/param.V_sei;
% 
% n_loss1=2*(param.as_n*param.len_n*param.pouch_area+...
%     param.A_crack*param.pouch_area)*L_sei/param.V_sei;

n_loss2=(param.n_Li_initial-n_Li);
% 
LLI=n_loss2/param.n_Li_initial*100;
%LLI=n_loss1/param.n_Li_initial*100;


end
