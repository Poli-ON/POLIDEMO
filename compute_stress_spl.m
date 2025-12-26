function [stress_r_p,stress_c_p,stress_r_n,stress_c_n]= compute_stress_spl(x,cs_p,cs_n,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE_STRESS
%
% This function computes the stress distribution within the active material 
% particles based on the lithium concentration profile. 
% The implemented formulation follows the approach described in:
%   https://doi.org/10.3390/en13071717
%
%       compute_stress_spl(x,cs_p,cs_n,param)
%
% INPUTS:
%   x      - Radial spatial mesh (particle coordinates)
%   cs_p   - Lithium concentration in cathode particles
%   cs_n   - Lithium concentration in anode particles
%   param  - Structure of parameters
%
% OUTPUTS:
%   stress - Matrix containing:
%              stress(:,1): radial stress
%              stress(:,2): hoop (tangential) stress
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

%global param

% Fit coefficients
cs_p_fun_spl= spline(x,cs_p.*x.^2);
cs_n_fun_spl= spline(x,cs_n.*x.^2);

% Mechanical boundary conditions
mecBC_p=param.om_p/((1-param.nu_p))*2*(1-2*param.nu_p)*fnval(fnint(cs_p_fun_spl),x(end)); % Cathode
mecBC_n=param.om_n/((1-param.nu_n))*2*(1-2*param.nu_n)*fnval(fnint(cs_n_fun_spl),x(end)); % Anode

% Stress
const_tens_p=(param.om_p*param.Ep)/(3*(1-param.nu_p));           % Cathode
const_tensh_p=(2*param.om_p*param.Ep)/(9*(1-param.nu_p));        % Cathode

const_tens_n=(param.om_n*param.En)/(3*(1-param.nu_n));           % Anode
const_tensh_n=(2*param.om_n*param.En)/(9*(1-param.nu_n));        % Anode


% Cathode
stress_r_p(1)=((-const_tensh_p*cs_p(1))+param.Ep*mecBC_p/(3*(1-2*param.nu_p)));
stress_c_p(1)=(-const_tensh_p*cs_p(1)+((param.Ep*mecBC_p)/(3*(1-2*param.nu_p))));

stress_r_p=[stress_r_p,(param.Ep*mecBC_p/(3*(1-2*param.nu_p)))-(2*const_tens_p*x(2:end).^(-3)).*fnval(fnint(cs_p_fun_spl),x(2:end))];
stress_c_p=[stress_c_p,(const_tens_p*x(2:end).^-3).*fnval(fnint(cs_p_fun_spl),x(2:end))-const_tens_p*cs_p(2:end)+(param.Ep*mecBC_p/(3*(1-2*param.nu_p)))];


% Anode
stress_r_n(1)=((-const_tensh_n*cs_n(1))+param.En*mecBC_n/(3*(1-2*param.nu_n)));
stress_c_n(1)=(-const_tensh_n*cs_n(1)+((param.En*mecBC_n)/(3*(1-2*param.nu_n))));

stress_r_n=[stress_r_n,(param.En*mecBC_n/(3*(1-2*param.nu_n)))-(2*const_tens_n.*x(2:end).^(-3)).*fnval(fnint(cs_n_fun_spl),x(2:end))];
stress_c_n=[stress_c_n,(const_tens_n*x(2:end).^-3).*fnval(fnint(cs_n_fun_spl),x(2:end))-const_tens_n*cs_n(2:end)+(param.En*mecBC_n/(3*(1-2*param.nu_n)))];

end
