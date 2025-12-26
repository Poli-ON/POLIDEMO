function V = get_Voltage(t,ur,I,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_VOLTAGE
% Computes the battery terminal voltage from the state vector of the SPM 
% model, accounting for overpotentials, interfacial resistances, and aging 
% mechanisms such as SEI, CEI, plating, and LAM.
%
%       V = get_Voltage(t, ur, I, T)
%
% INPUTS:
%   t   - Time [s]
%   ur  - State vector containing surface concentrations, aging variables, etc.
%   I   - Applied current [A]
%   T   - Temperature [K]
%
% OUTPUTS:
%   V   - Battery voltage [V]
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

% Surface lithium concentration [mol/m^3]
cs_p = ur(1,:);                       % Cathode
cs_n = ur(2,:);                       % Anode

% ================== INTERCALANT AREA ================== 
% Volume fraction
vol_fraction_solidphase_p=param.vol_fraction_solidphase(1);
vol_fraction_solidphase_n=param.vol_fraction_solidphase(2);
as_p       = 3*vol_fraction_solidphase_p/param.Rp;
as_n       = 3*vol_fraction_solidphase_n/param.Rn;

% Crack area
if param.EnableSEIoncracks == 1
    if param.CrackMechanism ~= 2
        if length(ur(:,1))<2*param.N
            a=ur(2+ param.aging_vector(1)+ param.aging_vector(2),:);
        else
            a=ur(2*param.N+2+ param.aging_vector(1)+ param.aging_vector(2),:);
        end
        A_crack=param.N_cracks*a*param.w_crack*2;
    else
        A_crack=param.A_crack;
    end
else
    a=[];
    A_crack=0;
end

% Intercalant area
Ap_int=as_p*param.len_p;
An_int=as_n*param.len_n+A_crack;

% ================== INTERCALATION FLUX ==================
if param.EnableSEI==1 || param.Enableplating==1
    if param.EnableCEI==1 % jn and jp unknown
        switch param.OperatingMode
            case 1
                jn = ur(end-1,:);
                jn_tot=I/(param.F*An_int);
                jp = ur(end,:);
                jp_tot=-I/(param.F*Ap_int);
            case 2
                jn = ur(end-3,:);
                jn_tot=ur(end-1,:);
                jp = ur(end-2,:);
                jp_tot=ur(end,:);
        end
    else % just jp unknown
        switch param.OperatingMode
            case 1
                jn = ur(end,:);
                jn_tot=I/(param.F*An_int);
                jp_tot = -I/(param.F*Ap_int);
                jp=jp_tot;
            case 2
                jn = ur(end-2,:);
                jn_tot = ur(end-1,:);
                jp_tot = ur(end,:);
                jp=jp_tot;
        end
    end
else % jn is known
    if param.EnableCEI==1 % jn and jp unknown
        switch param.OperatingMode
            case 1
                jn_tot =  I/(param.F*An_int);
                jn=jn_tot;
                jp = ur(end,:);
                jp_tot=-I/(param.F*Ap_int);
            case 2
                jn = ur(end-1,:);
                jn_tot=jn;
                jp = ur(end-2,:);
                jp_tot=ur(end,:);
        end
    else % just jp unknown
        switch param.OperatingMode
            case 1
                jn =  I/(param.F*An_int);
                jn_tot=jn;
                jp = -I/(param.F*Ap_int);
                jp_tot=jp;
            case 2
                jn = ur(end-1,:);
                jn_tot=jn;
                jp = ur(end,:);
                jp_tot=jp;
        end
    end
end

% ================== OCV AND OVERPOTENTIALS ==================
switch param.EnableTemperature
    case 0
        % Reaction rate constant [m^2.5/(mol^0.5 s)]
        kp=param.kp;                                     % Cathode
        kn=param.kn;                                     % Anode

        % Exchange current density [A/m^2]
        i0_p = kp*param.F*sqrt(param.ce)*sqrt(cs_p(:,:)).*sqrt(param.cs_maxp-cs_p(:,:));
        i0_n = kn*param.F*sqrt(param.ce)*sqrt(cs_n(:,:)).*sqrt(param.cs_maxn-cs_n(:,:));

        % Overpotential [V]
        eta_p=2*param.Rg*param.T/param.F*asinh(jp*param.F./(2*i0_p));
        eta_n=2*param.Rg*param.T/param.F*asinh(jn*param.F./(2*i0_n));

        % OCV [V]
        [U_p,U_n] =  openCircuitPotential([ur(1,:);ur(2,:)],param.T);
    case 1
        % Reaction rate constant [m^2.5/(mol^0.5 s)]
        [kp,kn]=get_Effectivek_coefficient(T);

        % Exchange current density [A/m^2]
        i0_p = kp'*param.F*sqrt(param.ce)*sqrt(cs_p(:,:)).*sqrt(param.cs_maxp-cs_p(:,:));
        i0_n = kn'*param.F*sqrt(param.ce)*sqrt(cs_n(:,:)).*sqrt(param.cs_maxn-cs_n(:,:));

        % Overpotential [V]
        eta_p = 2.*param.Rg.*T'./param.F.*asinh(jp.*param.F./(2.*i0_p));
        eta_n = 2.*param.Rg.*T'./param.F.*asinh(jn.*param.F./(2.*i0_n));

        % OCV [V]
        [U_p,U_n] =  openCircuitPotential([ur(1,:);ur(2,:)],T');

end

% ================== INTERNAL RESISTANCE ==================
R_int=param.Rf;

if param.EnableSEI == 1
    L_sei=ur(3,:);
    R_int=R_int+L_sei/param.sigma_sei;
end

if param.EnableCEI == 1
    L_cei=ur(2+ param.aging_vector(1)+ param.aging_vector(2)+param.aging_vector(3),:);
    R_int=R_int+L_cei/(param.sigma_cei);
end

if param.Enableplating == 1
    L_plating=ur(2+ param.aging_vector(1)+ param.aging_vector(2)+param.aging_vector(3)+param.aging_vector(4),:);
    R_int=R_int+L_plating/(param.sigma_plating);
end

if param.EnableLAM == 1
    R_int=R_int+param.k_LAM_p_resistance*(param.vol_fraction_solidphase_0(1)-param.vol_fraction_solidphase(1))/param.vol_fraction_solidphase_0(1)+...
        param.k_LAM_n_resistance*(param.vol_fraction_solidphase_0(2)-param.vol_fraction_solidphase(2))/param.vol_fraction_solidphase_0(2);
end

switch param.OperatingMode
    case 1
        V=U_p-U_n+eta_p-eta_n-R_int.*I;
    case 2
        V=U_p-U_n+eta_p-eta_n-R_int.*(jn_tot.*param.F.*An_int);
end

% ================== VOLTAGE COMPUTATION ==================
V=real(V);

end