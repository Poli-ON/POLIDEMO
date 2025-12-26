function result = spm_post_processing(sol,time_sol,t,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_POST_PROCESSING
% This function post-processes the results of the model and gives
% it back according to the variables to be stored.
%
%       result = spm_post_processing(sol,time_sol,t,x)
%
% INPUTS
% sol           Structure containing the solution computed with pdepe
% time_sol      Time at which the solution is computed
% t             Time at which the solution is post-processed
% x             Spatial mesh
%
% OUPUTS
% result        Structure with the results
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

% Preallocation
cs_p=zeros(length(time_sol),param.N+1);
cs_n=zeros(length(time_sol),param.N+1);
jp=zeros(length(time_sol),1);
jn=zeros(length(time_sol),1);
jn_tot=zeros(length(time_sol),1);
jp_tot=zeros(length(time_sol),1);
V=zeros(length(time_sol),1);
current=zeros(length(time_sol),1);
U_p=zeros(length(time_sol),1);
U_n=zeros(length(time_sol),1);
L_sei=zeros(length(time_sol),1);
j_sei=zeros(length(time_sol),1);
j_cei=zeros(length(time_sol),1);
j_plating=zeros(length(time_sol),1);
eta_p=zeros(length(time_sol),1);
eta_n=zeros(length(time_sol),1);
eta_sei=zeros(length(time_sol),1);
eta_cei=zeros(length(time_sol),1);
eta_plating=zeros(length(time_sol),1);
sigma_r_p=zeros(length(time_sol),param.N+1);
sigma_c_p=zeros(length(time_sol),param.N+1);
sigma_r_n=zeros(length(time_sol),param.N+1);
sigma_c_n=zeros(length(time_sol),param.N+1);
sigma_h_p=zeros(length(time_sol),param.N+1);
sigma_h_n=zeros(length(time_sol),param.N+1);
a=zeros(length(time_sol),1);
thk=zeros(length(time_sol),1);
cap=zeros(length(time_sol),1);
SOC=zeros(length(time_sol),1);
T=ones(length(time_sol),1)*param.Tref;
irr_swelling=0;

%--------------------------------------------------------------------------
% Temperature
%--------------------------------------------------------------------------
if param.EnableTemperature == 1
    T=sol(:,2*param.N+3);
end

%--------------------------------------------------------------------------
% Concentration at the specified time "t"
%--------------------------------------------------------------------------
[cs_p,cs_n]=compute_concentration(sol,time_sol,t);

%--------------------------------------------------------------------------
% Current density
%--------------------------------------------------------------------------
switch param.OperatingMode
    case 1
        current(:,1) = get_Current_density(t(:), [cs_p(:,end);cs_n(:,end)]);
    case 2
        current(:,1) = get_Current_density(t(:), sol(:,:)');
end

%--------------------------------------------------------------------------
% L_sei
%--------------------------------------------------------------------------
if param.EnableSEI==1
    L_sei=sol(:,2*param.N+2+param.aging_vector(1));
else
    L_sei=[];
end

%--------------------------------------------------------------------------
% L_CEI
%--------------------------------------------------------------------------
if param.EnableCEI==1
    L_cei=sol(:,2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+...
        param.aging_vector(3));
else
    L_cei=[];
end

%--------------------------------------------------------------------------
% L_plating
%--------------------------------------------------------------------------
if param.Enableplating==1
    L_plating=sol(:,2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+...
        param.aging_vector(3)+param.aging_vector(4));
else
    L_plating=[];
end

%--------------------------------------------------------------------------
% Crack length a
%--------------------------------------------------------------------------
if param.EnableSEIoncracks == 1
    if param.CrackMechanism ~= 2
        a=sol(:,2*param.N+2+ param.aging_vector(1)+ param.aging_vector(2));
        A_crack=param.N_cracks*a*param.w_crack*2;
    else
        a=[];
        A_crack=param.A_crack;
    end
else
    a=[];
    A_crack=0;
end
%--------------------------------------------------------------------------
% Intercalant area
%--------------------------------------------------------------------------
Ap_int=param.as_p*param.len_p;
An_int=param.as_n*param.len_n+A_crack;

%--------------------------------------------------------------------------
% Current density and Voltage
%--------------------------------------------------------------------------
if param.EnableSEI==1 || param.Enableplating==1
    if param.EnableCEI==1 % jn and jp uknown
        switch param.OperatingMode
            case 1
                jn = sol(:,end-1);
                jn_tot=current/(param.F*An_int);
                jp =sol(:,end);
                jp_tot=-current/(param.F*Ap_int);
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_sei';a';L_cei';L_plating';jn';jp'],current(:)',T);
            case 2
                jn =  sol(:,end-3);
                jn_tot= sol(:,end-1);
                jp =  sol(:,end-2);
                jp_tot= sol(:,end);
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_sei';a';L_cei';L_plating';jn';jp';jn_tot';jp_tot'],current(:)',T);
        end
    else % just jn unknown
        switch param.OperatingMode
            case 1
                jn = sol(:,end);
                jn_tot=current/(param.F*An_int);
                jp_tot = -current/(param.F*Ap_int);
                jp=jp_tot;
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_sei';a';L_plating';jn'],current(:)',T);

            case 2
                jn = sol(:,end-2);
                jn_tot = sol(:,end-1);
                jp_tot = sol(:,end);
                jp=jp_tot;
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_sei';a';L_plating';jn';jn_tot';jp_tot'],current(:)',T);
        end
    end
else % jn is known
    if param.EnableCEI==1 %  jp unkwown
        switch param.OperatingMode
            case 1
                jn_tot =  current/(param.F*An_int);
                jn=jn_tot;
                jp = sol(:,end);
                jp_tot=-current/(param.F*Ap_int);
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_cei';L_plating';jp'],current(:)',T);

            case 2
                jn = sol(:,end-1);
                jn_tot=jn;
                jp = sol(:,end-2);
                jp_tot=sol(:,end);
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_cei';L_plating';jp';jn_tot';jp_tot'],current(:)',T);
        end
    else
        switch param.OperatingMode
            case 1
                jn =  current/(param.F*An_int);
                jn_tot=jn;
                jp = -current/(param.F*Ap_int);
                jp_tot=jp;
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_plating'],current(:)',T);
            case 2
                jn = sol(:,end-1);
                jn_tot=jn;
                jp = sol(:,end);
                jp_tot=jp;
                V = get_Voltage(t(:),[cs_p(:,end)';cs_n(:,end)';L_plating';jn_tot';jp_tot'],current(:)',T);
        end
    end
end

V=real(V);

V=V';

%--------------------------------------------------------------------------
% OCV and Overpotential
%--------------------------------------------------------------------------
switch param.EnableTemperature
    case 0
        % Reaction rate constant [m^2.5/(mol^0.5 s)]
        kp=param.kp;                                     % Cathode
        kn=param.kn;                                     % Anode

        % Exchange current density [A/m^2]
        i0_p = kp*param.F*sqrt(param.ce)*sqrt(cs_p(:,end)).*sqrt(param.cs_maxp-cs_p(:,end));
        i0_n = kn*param.F*sqrt(param.ce)*sqrt(cs_n(:,end)).*sqrt(param.cs_maxn-cs_n(:,end));

        % Overpotential [V]
        eta_p=2*param.Rg*param.T/param.F*asinh(jp*param.F./(2*i0_p));
        eta_n=2*param.Rg*param.T/param.F*asinh(jn*param.F./(2*i0_n));

        % Open circuit potential
        for i=1:length(cs_p(:,end))
            [U_p(i,1),U_n(i,1)] =  openCircuitPotential([cs_p(i,end)';cs_n(i,end)'],T');
        end
    case 1
        % Reaction rate constant [m^2.5/(mol^0.5 s)]
        [kp,kn]=get_Effectivek_coefficient(T);

        % Exchange current density [A/m^2]
        i0_p = kp.*param.F*sqrt(param.ce).*sqrt(cs_p(:,end)).*sqrt(param.cs_maxp-cs_p(:,end)) ;
        i0_n = kn.*param.F*sqrt(param.ce).*sqrt(cs_n(:,end)).*sqrt(param.cs_maxn-cs_n(:,end)) ;

        % Overpotential [V]
        eta_p=2.*param.Rg.*T./param.F.*asinh(jp.*param.F./(2*i0_p));
        eta_n=2.*param.Rg.*T./param.F.*asinh(jn.*param.F./(2*i0_n));
end

%--------------------------------------------------------------------------
% SEI flux
%--------------------------------------------------------------------------
if param.EnableSEI == 1 || param.EnableSEIoncracks == 1
    eta_sei = eta_n+U_n-param.U_sei;
    j_sei=get_sei_flux([cs_p(:,end)';cs_n(:,end)'],L_sei',jn');
end

%--------------------------------------------------------------------------
% CEI flux
%--------------------------------------------------------------------------
if param.EnableCEI == 1
    eta_cei = eta_n+U_n-param.U_cei;
    j_cei=get_cei_flux([cs_p(:,end)';cs_n(:,end)'],L_cei',jp');
end

%--------------------------------------------------------------------------
% Plating flux
%--------------------------------------------------------------------------
% if param.Enableplating == 1
%     eta_plating = eta_n+U_n;
%     [cs_p_avg,cs_n_avg]=compute_average_concentration(u(1:param.N+1),u(param.N+2:2*param.N+2),linspace(0,1,param.N+1)');
%     j_plating=get_plating_flux([cs_p(:,end)';cs_n(:,enk_cd)'],cs_n_avg,jn_int);
% end

%--------------------------------------------------------------------------
% Internal resistance
%--------------------------------------------------------------------------
R_sei=(L_sei/param.sigma_sei);
R_cei=(L_cei/param.sigma_cei);
R_plating=L_plating/param.sigma_plating;
R_LAM=param.k_LAM_p_resistance*(param.vol_fraction_solidphase_0(1)-param.vol_fraction_solidphase(1))/param.vol_fraction_solidphase_0(1)+...
    param.k_LAM_n_resistance*(param.vol_fraction_solidphase_0(2)-param.vol_fraction_solidphase(2))/param.vol_fraction_solidphase_0(2);

%--------------------------------------------------------------------------
% Stress computation
%--------------------------------------------------------------------------
if param.EnableSEIoncracks ==1 || param.EnableLAM == 1
    for i=1:length(cs_p(:,1))
        [sigma_r_p(i,:),sigma_c_p(i,:),sigma_r_n(i,:),sigma_c_n(i,:)]= compute_stress_spl(x,cs_p(i,:),cs_n(i,:),param);
    end
end

%--------------------------------------------------------------------------
% Reversible battery thickness change
%--------------------------------------------------------------------------
%To speed up simulation, this computation can be avoided
 % for i=1:length(cs_p)
 % thk(i) = get_battery_thickness(x,cs_p(i,:),cs_n(i,:));
 % end
 % thk=thk-thk(1);
%--------------------------------------------------------------------------
% Capacity
%--------------------------------------------------------------------------
switch param.OperatingMode
    case 1
        cap= abs(param.Crate)*param.assumed_cell_capacity_Ah *time_sol/3600;
    case 2

        cap(1)=0;
        for i = 2:length(t)
            cap(i) =(param.pouch_area*param.no_of_layers)/3600*trapz(time_sol(1:i), abs(current(1:i)));
        end

end

%--------------------------------------------------------------------------
% Average concentration
%--------------------------------------------------------------------------
[cs_p_avg,cs_n_avg]=compute_average_concentration(cs_p(end,:),cs_n(end,:),x);

%--------------------------------------------------------------------------
% LAM and LLI
%--------------------------------------------------------------------------
LLI=LLI_computation(cs_p_avg,cs_n_avg);
[LAM_p,LAM_n]=LAM_computation;

%--------------------------------------------------------------------------
% Final SOC
%--------------------------------------------------------------------------
SOC=compute_SOC(cs_p_avg(:),cs_n_avg(:));

%--------------------------------------------------------------------------
% Compute irreversible swelling
%--------------------------------------------------------------------------

% Note: irreversibile swelling is computeg just at the end of the discharge

if all(current > 0) && (param.EnableSEI == 1)
    if ~isempty(a)
        irr_swelling=spm_get_irr_swelling(param.vol_fraction_solidphase(2),L_sei(end),a(end));
    else
        irr_swelling=spm_get_irr_swelling(param.vol_fraction_solidphase(2),L_sei(end));
    end

end

%--------------------------------------------------------------------------
% Build the structure of results
%--------------------------------------------------------------------------
result.x=x;
result.cs_p=cs_p;
result.cs_n=cs_n;
result.jp=jp;
result.jn=jn;
result.V=V;
result.current=current;
result.U_p=U_p;
result.U_n=U_n;
result.param=param;
result.time_sol=time_sol;
result.L_sei=L_sei;
result.L_cei=L_cei;
result.L_plating=L_plating;
result.j_sei=j_sei;
result.j_cei=j_cei;
result.j_plating=j_plating;
result.eta_plating=eta_plating;
result.eta_p=eta_p;
result.eta_n=eta_n;
result.eta_cei=eta_cei;
result.sigma_r_p=sigma_r_p;
result.sigma_c_p=sigma_c_p;
result.sigma_r_n=sigma_r_n;
result.sigma_c_n=sigma_c_n;
result.sigma_h_p=(sigma_r_p+2*sigma_c_p)/3;
result.sigma_h_n=(sigma_r_n+2*sigma_c_n)/3;
result.delta_sigma_p= param.sigma_c_p_max - param.sigma_c_p_min;
result.delta_sigma_n= param.sigma_c_n_max - param.sigma_c_n_min;
result.vol_fraction_solidphase_p=param.vol_fraction_solidphase(1);
result.vol_fraction_solidphase_n=param.vol_fraction_solidphase(2);
result.as_p = 3*result.vol_fraction_solidphase_p/param.Rp;
result.as_n = 3*result.vol_fraction_solidphase_n/param.Rn;
result.thk=thk;
result.cap=cap;
result.R_LAM=R_LAM;
result.T=T;
result.LLI=LLI;
result.LAM_n=LAM_n;
result.LAM_p=LAM_p;
result.jn_tot=jn_tot;
result.jp_tot=jp_tot;
result.SOC=SOC;
result.R_sei=R_sei;
result.R_cei=R_cei;
result.a=a;
result.irr_swelling=irr_swelling;
%result.R_plating=R_plating;

if isempty(L_sei) ==1
    result.L_sei=zeros(length(time_sol),1);
    result.R_sei=zeros(length(time_sol),1);
end

if isempty(L_cei) ==1
    result.L_cei=zeros(length(time_sol),1);
    result.R_cei=zeros(length(time_sol),1);
end

if isempty(L_plating) ==1
    result.L_plating=zeros(length(time_sol),1);
    result.R_plating=zeros(length(time_sol),1);
end

if isempty(a) ==1
    result.a=zeros(length(time_sol),1);
end

%--------------------------------------------------------------------------
% Update the initialState
%--------------------------------------------------------------------------

vector=[];

switch param.EnableTemperature
    case 0
        if isempty(L_sei) == 0
            vector=[vector;L_sei(end)];
        end

        if isempty(a) ==0
            vector=[vector;a(end)];
        end

        if isempty(L_cei) ==0
            vector=[vector;L_cei(end)];
        end

        if isempty(L_plating) ==0
            vector=[vector;L_plating(end)];
        end

        if isempty(L_sei)==0 || isempty(L_plating)==0
            vector=[vector;jn(end)];
        end

        if isempty(L_cei)==0
            vector=[vector;jp(end)];
        end

        result.initialState=[cs_p(end,:)';cs_n(end,:)';vector];

    case 1 % Thermal effects on - Beta versio
        switch param.EnableAging
            case 0
                result.initialState=[cs_p(end,:)';cs_n(end,:)';T(end)];
            case 1
                result.initialState=[cs_p(end,:)';cs_n(end,:)';L_sei(end);jn(end)];
            case 2
                result.initialState=[cs_p(end,:)';cs_n(end,:)';L_sei(end);a(end);jn(end)];
            case 3
                result.initialState=[cs_p(end,:)';cs_n(end,:)'];
            case 4
                result.initialState=[cs_p(end,:)';cs_n(end,:)';L_sei(end);a(end);jn(end)];
        end

end

%--------------------------------------------------------------------------
% Update cyclic stress
%--------------------------------------------------------------------------
% Maximum in discharge
if param.EnableSEIoncracks == 1 || param.EnableLAM == 1
    if current > 0
        param.sigma_c_p_max=max(sigma_c_p(2:end,end));
        [param.sigma_c_n_max,index]=max(sigma_c_n(2:end,end));
    end
    % Minimum in charge
    if current < 0
        param.sigma_c_p_min=min(sigma_c_p(2:end,end));
        [param.sigma_c_n_min,index]=min(sigma_c_n(2:end,end));
    end

    if param.OperatingMode == 2
        if min(sigma_c_n(2:end,end)) <  param.sigma_c_n_min
            param.sigma_c_p_min=min(sigma_c_p(2:end,end));
            [param.sigma_c_p_min,index]=min(sigma_c_n(2:end,end));
        end
        if min(sigma_c_p(2:end,end)) <  param.sigma_c_p_min
            param.sigma_c_p_min=min(sigma_c_p(2:end,end));
            [param.sigma_c_n_min,index]=min(sigma_c_n(2:end,end));
        end
    end

    if current >0 % at the end of discharge
        switch param.CrackMechanism
            case 1
                delta_sigma_p= param.sigma_c_p_max- param.sigma_c_p_min;
                delta_sigma_n= param.sigma_c_n_max- param.sigma_c_n_min;
                %param.N_cracks=param.k_crack*(1/(param.n_LAM_n+1).*(param.a_LAM_n*delta_sigma_n.^param.m_LAM_n)*(param.cycle_count).^(param.n_LAM_n+1)+param.vol_fraction_solidphase_0(2));

            case 2
                delta_sigma_p= param.sigma_c_p_max- param.sigma_c_p_min;
                delta_sigma_n= param.sigma_c_n_max- param.sigma_c_n_min;
                %param.A_crack=param.k_crack*(1/(param.n_LAM_n+1).*(param.a_LAM_n*delta_sigma_n.^param.m_LAM_n)*(param.cycle_count).^(param.n_LAM_n+1)+param.vol_fraction_solidphase_0(2));
        end
    end

end

%--------------------------------------------------------------------------
% Update cyclic stress distribution in the crack area
%--------------------------------------------------------------------------
if param.EnableSEIoncracks == 1
    if param.CrackMechanism ~= 2
        param.a=a(end);
        a_fraq_R=param.a/param.Rn;
        [~, pos] = min(abs(x - (1-a_fraq_R)));
        % Ora devo trovare l'istante di tempo in cui lo stress
        % all'apice della cricca è massimo

        if current <0
            %[param.sigma_c_n_max_crack_charge,index]=max(sigma_c_n(2:end,pos));
            if pos == 11
                [param.sigma_c_n_max_crack_charge,index]=max(sigma_c_n(2:end,pos));
                param.sigma_c_n_max_distribution_crack_charge=(sigma_c_n(index,:))';
            else
                [param.sigma_c_n_max_crack_charge,index]=max(trapz(x(pos:end),sigma_c_n(2:end,pos:end),2));
                param.sigma_c_n_max_distribution_crack_charge=(sigma_c_n(index,:))';
            end
        end

        if current >0
            if pos == 11
                [param.sigma_c_n_max_crack_discharge,index]=max(sigma_c_n(2:end,pos));
                param.sigma_c_n_max_distribution_crack_discharge=(sigma_c_n(index,:))';
            else
                [param.sigma_c_n_max_crack_discharge,index]=max(trapz(x(pos:end),sigma_c_n(2:end,pos:end),2));
                param.sigma_c_n_max_distribution_crack_discharge=(sigma_c_n(index,:))';
            end

            % Considero il valore massimo tra quello trovato in
            % carica e quello trovato in scarica

            if param.sigma_c_n_max_crack_discharge > param.sigma_c_n_max_crack_charge
                param.sigma_c_n_max_crack=param.sigma_c_n_max_crack_discharge;
                param.sigma_c_n_max_distribution_crack= param.sigma_c_n_max_distribution_crack_discharge;
            else
                param.sigma_c_n_max_crack=param.sigma_c_n_max_crack_charge;
                param.sigma_c_n_max_distribution_crack= param.sigma_c_n_max_distribution_crack_charge;

            end
        end
    end

end

%--------------------------------------------------------------------------
% Update volume fraction
%--------------------------------------------------------------------------
if param.EnableLAM == 1
    if current > 0
        delta_sigma_p= param.sigma_c_p_max- param.sigma_c_p_min;
        delta_sigma_n= param.sigma_c_n_max- param.sigma_c_n_min;
        vol_fraction_solidphase_p=-1/(param.n_LAM_p+1).*(param.a_LAM_p*delta_sigma_p.^param.m_LAM_p)*(param.cycle_count).^(param.n_LAM_p+1)+param.vol_fraction_solidphase_0(1);
        vol_fraction_solidphase_n=-1/(param.n_LAM_n+1).*(param.a_LAM_n*delta_sigma_n.^param.m_LAM_n)*(param.cycle_count).^(param.n_LAM_n+1)+param.vol_fraction_solidphase_0(2);
        param.as_p = 3*vol_fraction_solidphase_p/param.Rp;
        param.as_n = 3*vol_fraction_solidphase_n/param.Rn;
        param.vol_fraction_solidphase=[vol_fraction_solidphase_p;vol_fraction_solidphase_n];

        if param.EnableSEIoncracks == 1
            switch param.CrackMechanism
                case 1
                    param.N_cracks=param.k_crack*(1/(param.n_LAM_n+1).*(param.a_LAM_n*delta_sigma_n.^param.m_LAM_n)*(param.cycle_count).^(param.n_LAM_n+1)+param.vol_fraction_solidphase_0(2));
                case 2
                    param.A_crack=param.k_crack*((-vol_fraction_solidphase_n+param.vol_fraction_solidphase_0(2))/param.vol_fraction_solidphase_0(2));
            end
        end
    end

end


end
