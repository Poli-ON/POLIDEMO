function [value,term,direction] = stop_conditions_fun(t,u,voltage_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOP_CONDITIONS_FUN
% Defines the event conditions to stop the integration of the electrochemical
% model when specific physical or operational limits are reached.
%
%       [value, term, direction] = stop_conditions_fun(t, u, voltage_limit)
%
% INPUT:
%   t                   - Current time [s]
%   u                   - State vector at the mesh points
%   voltage_limit       - Vector with the lower and upper voltage limits [V]
%
% OUTPUT:
%   value               - Event function values; integration stops when any
%                         element of 'value' is zero
%   term                - Vector indicating whether the integration stops when
%                         the corresponding event occurs (1 = stop, 0 = ignore)
%   direction           - Vector indicating the direction of zero-crossing:
%                           1  -> positive slope only
%                          -1  -> negative slope only
%                           0  -> any slope
%
%   Stopping conditions considered in this function:
%     - Voltage exceeding the limits imposed during constant current (CC)
%     - Current falling below the threshold imposed during constant voltage (CV)
%     - SOC outside its admissible range (0% in discharge, 100% in charge)
%     - Concentration values becoming negative (should always remain positive)
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

% Surface concentration in positive and negative electrodes
u1 = u(param.N+1);
u2 = u(2*param.N+2);

cs_p_value = u1;
cs_n_value = u2;

% Average concentration in the negative electrode (used for SOC calculation)
cs_n = u(param.N+2:2*param.N+2);
x = linspace(0,1,param.N+1);
cs_n_avg = 3/(param.Rn^3) * trapz(x*param.Rn, cs_n'.*(x*param.Rn).^2) / param.cs_maxn;

% State of charge (anode side)
SOC = (cs_n_avg - param.theta_min_neg) / (param.theta_max_neg - param.theta_min_neg) * 100;

% SOC stopping conditions
SOC_check_discharge = SOC - 0;    % Stop if SOC < 0%
SOC_check_charge    = 100 - SOC;  % Stop if SOC > 100%

if param.EnableTemperature == 0
    T=param.T;
else
    T=u(2*param.N+3);
end

if param.EnableTemperature == 0
    T=param.T;
else
    T=u(2*param.N+3);
end

% Current density
 I = get_Current_density(t,u(end));

% Volume fraction
vol_fraction_solidphase_p=param.vol_fraction_solidphase(1);
vol_fraction_solidphase_n=param.vol_fraction_solidphase(2);
as_p       = 3*vol_fraction_solidphase_p/param.Rp;
as_n       = 3*vol_fraction_solidphase_n/param.Rn;

% Aging variables (SEI, cracks, CEI, plating)
if param.EnableSEI == 1
    L_sei = u(2*param.N+2+param.aging_vector(1));
else
    L_sei = [];
end

if param.EnableSEIoncracks == 1 && param.CrackMechanism ~= 2
    a = u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2));
else
    a = [];
end

if param.EnableCEI == 1
    L_cei = u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+param.aging_vector(3));
else
    L_cei = [];
end

if param.Enableplating == 1
    L_plating = u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+ ...
                  param.aging_vector(3)+param.aging_vector(4));
else
    L_plating = [];
end

% Crack surface area
if param.EnableSEIoncracks == 1
    if param.CrackMechanism ~= 2
        a = u(2*param.N+2+ param.aging_vector(1)+ param.aging_vector(2));
        A_crack = param.N_cracks * a * param.w_crack * 2;
    else
        a = [];
        A_crack = param.A_crack;
    end
else
    a = [];
    A_crack = 0;
end

% Intercalation area
Ap_int=as_p*param.len_p;
An_int=param.as_n*param.len_n+A_crack;

% Voltage computation
if param.EnableSEI==1 || param.Enableplating==1
    if param.EnableCEI==1 % both jn e jp unknown
        switch param.OperatingMode
            case 1
                jn = u(end-1,:);
                jn_tot=I/(param.F*An_int);
                jp = u(end,:);
                jp_tot=-I/(param.F*Ap_int);
                V = get_Voltage(t,[u1(end);u2(end);L_sei;a;L_cei;L_plating;jn;jp],I);
            case 2
                jn = u(end-3,:);
                jn_tot=u(end-1,:);
                jp = u(end-2,:);
                jp_tot=u(end,:);
                V = get_Voltage(t,[u1(end);u2(end);L_sei;a;L_cei;L_plating;jn;jp;jn_tot;jp_tot],I);
        end
    else % only jn unknown
        switch param.OperatingMode
            case 1
                jn = u(end,:);
                jn_tot=I/(param.F*An_int);
                jp_tot = -I/(param.F*Ap_int);
                jp=jp_tot;
                 V = get_Voltage(t,[u1(end);u2(end);L_sei;a;L_plating;jn],I);
            case 2
                jn = u(end-2,:);
                jn_tot = u(end-1,:);
                jp_tot = u(end,:);
                jp=jp_tot;
                V = get_Voltage(t,[u1(end);u2(end);L_sei;a;L_plating;jn;jn_tot;jp_tot],I);
        end
    end
else % jn is known
    if param.EnableCEI==1 %  jp unknown
        switch param.OperatingMode
            case 1
                jn_tot =  I/(param.F*An_int);
                jn=jn_tot;
                jp = u(end,:);
                jp_tot=-I/(param.F*Ap_int);
                 V = get_Voltage(t,[u1(end);u2(end);L_cei;L_plating;jp],I);
            case 2
                jn = u(end-1,:);
                jn_tot=jn;
                jp = u(end-2,:);
                jp_tot=u(end,:);
                 V = get_Voltage(t,[u1(end);u2(end);L_cei;L_plating;jp;jn_tot;jp_tot],I);
        end
    else % no unknown
        switch param.OperatingMode
            case 1
                jn =  I/(param.F*An_int);
                jn_tot=jn;
                jp = -I/(param.F*Ap_int);
                jp_tot=jp;
                V = get_Voltage(t,[u1(end);u2(end);L_plating],I,T);
            case 2
                jn = u(end-1,:);
                jn_tot=jn;
                jp = u(end,:);
                jp_tot=jp;
                V = get_Voltage(t,[u1(end);u2(end);L_plating;jn_tot;jp_tot],I,T);
        end
    end
end

V = real(V);

% Event stopping conditions

% Voltage limits
V_check_discharge = V - voltage_limit(1);   % stop if below cutoff
V_check_charge    = voltage_limit(2) - V;   % stop if above cutoff

% Current stop for charge (stop when current < C/20)
I_check_charge = param.i_1C_density/20 - abs(I);

% Build event outputs
switch param.EnableSOCcontrol
    case 0
        % Events based only on voltage and concentrations
        switch param.OperatingMode
            case 1   % Constant current
                if param.Crate > 0 % Discharge
                    value = [cs_p_value, cs_n_value, V_check_discharge];
                    direction = [0,0,0];
                    term = [1,1,1];
                else   % Charge
                    value = [cs_p_value, cs_n_value, V_check_charge];
                    direction = [0,0,0];
                    term = [1,1,1];
                end
            case 2   % Constant voltage
                value = [cs_p_value, cs_n_value, V_check_discharge, I_check_charge];
                direction = [0,0,-1,1];
                term = [1,1,1,1];
        end

    case 1
        % Events based on voltage, concentrations and SOC
        switch param.OperatingMode
            case 1   % Constant current
                if param.Crate > 0 % Discharge
                    value = [cs_p_value, cs_n_value, V_check_discharge, SOC_check_discharge];
                    direction = [0,0,0,0];
                    term = [1,1,1,1];
                else   % Charge
                    value = [cs_p_value, cs_n_value, V_check_charge, SOC_check_charge];
                    direction = [0,0,0,0];
                    term = [1,1,1,1];
                end
            case 2   % Constant voltage
                value = [cs_p_value, cs_n_value, V_check_discharge, I_check_charge, SOC_check_charge];
                direction = [0,0,-1,1,1];
                term = [1,1,1,1,1];
        end
end

end














