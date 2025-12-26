function dudt = spm_odefun_aging_new(t,u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_PDEFUN_AGING_NEW
% Rewrites the SPM equations with aging mechanisms in the form expected by
% ODE solvers.
%
%       dudt = SPM_PDEFUN_AGING_NEW(t, u)
%
% This function includes:
%   - SEI growth
%   - SEI growth on crack surfaces
%   - Volume fraction variation due to fracture
%   - CEI growth (beta version)
%   - Lithium plating (beta version)
%
% INPUT:
%   t           - Independent time variable
%   u           - Vector of dependent variables
%
% OUTPUT:
%   dudt        - Time derivative of the state vector
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

    %% -----------------------------------------------------------------------
    %  MODE 1 - Aging during CC operation
    %% -----------------------------------------------------------------------

    case 1

        % Superficial concentrantion
        ur=[u(param.N+1);u(2*param.N+2)];

        % ================== SEI / CRACKS ==================
        switch param.EnableSEI
            case 0 % SEI off
                L_sei=[];
                A_crack=0;
                a=[];
            case 1 % SEI on
                L_sei=u(2*param.N+3);
                switch param.EnableSEIoncracks
                    case 0 % fracture off
                        a=[];
                        A_crack=0;
                    case 1
                        switch param.CrackMechanism
                            case 0
                                a=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2));
                                A_crack=param.N_cracks*a*param.w_crack*2;
                            case 1
                                a=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2));
                                A_crack=param.N_cracks*a*param.w_crack*2;
                            case 2
                                a=[];
                                A_crack=param.A_crack;
                        end
                end
        end

        % ================== CEI ================== (Beta version)
        switch param.EnableCEI
            case 0
                L_cei=[];
            case 1
                L_cei=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+param.aging_vector(3));
        end

        % ================== PLATING ================== (Beta version)
        switch param.Enableplating
            case 0
                L_plating=[];
            case 1
                L_plating=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+param.aging_vector(3)+param.aging_vector(4));
        end

        % Current density computation [A/m^2]
        I = get_Current_density(t,ur);

        % Volume fraction of solid phase
        vol_fraction_solidphase_p=param.vol_fraction_solidphase(1);
        vol_fraction_solidphase_n=param.vol_fraction_solidphase(2);

        % Update surface area to volume ratio
        as_p       = 3*vol_fraction_solidphase_p/param.Rp;
        as_n       = 3*vol_fraction_solidphase_n/param.Rn;

        % Intercalant area
        Ap_int=as_p*param.len_p;
        An_int=as_n*param.len_n+A_crack;

        % Molar flux [mol/(m^2*s)]
        jp_tot =  -I/(param.F*Ap_int);  % Cathode
        jn_tot = I/(param.F*An_int);    % Anode

        % ================== INTERCALATION FLUX ==================
        if isempty(L_sei)==0 || isempty(L_plating)==0
            if isempty(L_cei)==0
                jn_int=u(end-1);
            else
                jn_int=u(end);
            end
        else
            jn_int=jn_tot;
        end

        % ================== SEI FLUX ==================
        if isempty(L_sei)==0
            j_sei=get_sei_flux(ur,L_sei,jn_int);
        else
            j_sei=0;
        end

        % ================== PLATING FLUX ==================
        if isempty(L_plating)==0
            [cs_p_avg,cs_n_avg]=compute_average_concentration(u(1:param.N+1),u(param.N+2:2*param.N+2),linspace(0,1,param.N+1)');
            j_plating=get_plating_flux(ur,cs_n_avg,jn_int);
        else
            j_plating=0;
        end

        % ================== CEI FLUX ==================
        if isempty(L_cei)==0
            jp_int=u(end);
            j_cei=get_cei_flux(ur,L_cei,jp_int);
        else
            j_cei=0;
            jp_int=jp_tot;
        end

        % ================== CRACK PROPAGATION ==================
        if isempty(a)==0
            switch param.CrackMechanism
                case 0
                    dadt=crack_propagation(a);
                case 1
                    dadt=crack_propagation(a);
            end
        else
            dadt=[];
        end

        % ================== GOVERNING EQUATIONS ==================
        switch param.EnableTemperature
            case 0
                dudt= [param.A(1:param.N+1,:)*u(1:param.N+1)+param.B(1:param.N+1)*jp_int
                       param.A(param.N+2:2*param.N+2,:)*u(param.N+2:2*param.N+2)+param.B(param.N+2:2*param.N+2)*jn_int];
            case 1
                % Temperature
                T=u(end);

                % Temperature time derivative
                [Dp,Dn]=get_EffectiveDiffusion_coefficient(T);
                [param.A,param.B]=spm_spatial_discretization(Dp,Dn);
                dTdt= spm_lumped_thermal_model(t,ur,I,T);
                dudt=[param.A(1:param.N+1,:)*u(1:param.N+1)+param.B(1:param.N+1)*jp_int
                      param.A(param.N+2:2*param.N+2,:)*u(param.N+2:2*param.N+2)+param.B(param.N+2:2*param.N+2)*jn_int
                      dTdt];
        end

        % SEI and fracture
        if isempty(L_sei) == 0    % L_sei is not empty
            % Frattura
            if isempty(dadt) == 0 % dadt is not empty
                dudt= [dudt
                    (-1/2*param.V_sei*j_sei)*An_int/(An_int+2*param.N_cracks*dadt*param.w_crack)
                    dadt];
            else
                dudt= [dudt
                    (-1/2*param.V_sei*j_sei)];
            end
        end

        % CEI
        if isempty(L_cei)==0
            dudt=[dudt
                j_cei*param.V_cei];
        end

        % Plating
        if isempty(L_plating)==0
            dudt=[dudt
                -j_plating*param.V_plating];
        end

        % Add the flux
        if isempty(L_sei)==0 || isempty(L_plating)==0
            dudt=[dudt
                jn_int+j_sei+j_plating-jn_tot];
        end

        if isempty(L_cei)==0
            if jp_int >0
                dudt=[dudt
                    jp_int+j_cei-jp_tot];
            else
                dudt=[dudt
                    jp_int-j_cei-jp_tot];
            end
        end

        %% -----------------------------------------------------------------------
        %  MODE 2 - Aging during CV operation
        %% -----------------------------------------------------------------------

    case 2
        % The state vector u is [cp,cn,Lsei,a,lcei,lplating,jn_int,jp_int,jn_tot,jp_tot]

        % Superficial concentrantion
        ur=[u(param.N+1);u(2*param.N+2)];

        % ================== SEI / CRACKS ==================
        switch param.EnableSEI
            case 0 % SEI off
                L_sei=[];
                a=[];
                A_crack=0;
            case 1 % SEI on
                L_sei=u(2*param.N+3);
                switch param.EnableSEIoncracks
                    case 0 % fracture off
                        a=[];
                        A_crack=0;
                    case 1
                        switch param.CrackMechanism
                            case 0
                                a=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2));
                                A_crack=param.N_cracks*a*param.w_crack*2;
                            case 1
                                a=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2));
                                A_crack=param.N_cracks*a*param.w_crack*2;
                            case 2
                                a=[];
                                A_crack=param.A_crack;
                        end
                end
        end

        % ================== CEI ==================
        switch param.EnableCEI
            case 0
                L_cei=[];
            case 1
                L_cei=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+param.aging_vector(3));
        end

        % ================== PLATING ==================
        switch param.Enableplating
            case 0
                L_plating=[];
            case 1
                L_plating=u(2*param.N+2+param.aging_vector(1)+param.aging_vector(2)+param.aging_vector(3)+param.aging_vector(4));
        end

        % ===================== TOTAL FLUX =======================
        jp_tot = u(end);     % Cathode
        jn_tot = u(end-1);   % Anode

        % ================== INTERCALATION FLUX ==================
        if isempty(L_sei)==0 || isempty(L_plating)==0 % jn incongnito
            % Get lithium intercalation flux at boundary of anode
            if isempty(L_cei)==0 % jp incognito
                jn_int=u(end-3);
                jp_int=u(end-2);
                j_cei=get_cei_flux(ur,L_cei,jp_int);
            else
                jn_int=u(end-2);
                j_cei=0;
                jp_int=jp_tot;
            end
        else % jn is known
            jn_int=jn_tot;
            if isempty(L_cei)==0
                jp_int=u(end-2);
                j_cei=get_cei_flux(ur,L_cei,jp_int);
            else
                j_cei=0;
                jp_int=jp_tot;
            end
        end

        % ================== SEI FLUX ==================
        if isempty(L_sei)==0
            j_sei=get_sei_flux(ur,L_sei,jn_int);
        else
            j_sei=0;
        end

        % ================== PLATING FLUX ==================
        if isempty(L_plating)==0
            [cs_p_avg,cs_n_avg]=compute_average_concentration(u(1:param.N+1),u(param.N+2:2*param.N+2),linspace(0,1,param.N+1)');
            j_plating=get_plating_flux(ur,cs_n_avg,jn_int);
        else
            j_plating=0;
        end

        % ================== OCV AND OVERPOTENTIALS ==================
        [U_p,U_n] =  openCircuitPotential(ur(1:2));

        % Reaction rate constant [m^2.5/(mol^0.5 s)]
        kp=param.kp;                                     % Cathode
        kn=param.kn;                                     % Anode

        % Exchange current density [A/m^2]
        i0_p = kp*param.F*sqrt(param.ce)*sqrt(ur(1))*sqrt(param.cs_maxp-ur(1)) ;
        i0_n = kn*param.F*sqrt(param.ce)*sqrt(ur(2))*sqrt(param.cs_maxn-ur(2)) ;

        % Overpotential [V]
        eta_p = 2*param.Rg*param.T/param.F*asinh(jp_int*param.F/(2*i0_p));
        eta_n = 2*param.Rg*param.T/param.F*asinh(jn_int*param.F/(2*i0_n));

        % ================== INTERCALANT AREA ==================
        % Volume fraction of solid phase
        vol_fraction_solidphase_p=param.vol_fraction_solidphase(1);
        vol_fraction_solidphase_n=param.vol_fraction_solidphase(2);

        % Update surface area to volume ratio
        as_p       = 3*vol_fraction_solidphase_p/param.Rp;
        as_n       = 3*vol_fraction_solidphase_n/param.Rn;

        % Intercalant area
        Ap_int=as_p*param.len_p;
        An_int=as_n*param.len_n+A_crack;

        % ================== CRACK PROPAGATION ==================
        if isempty(a)==0
            switch param.CrackMechanism
                case 0
                    dadt=crack_propagation(a);
                case 1
                    dadt=crack_propagation(a);
            end
        else
            dadt=[];
        end

        % ================== GOVERNING EQUATIONS ==================
        dudt= [param.A(1:param.N+1,:)*u(1:param.N+1)+param.B(1:param.N+1)*jp_int
            param.A(param.N+2:2*param.N+2,:)*u(param.N+2:2*param.N+2)+param.B(param.N+2:2*param.N+2)*jn_int];

        % SEI and fracture
        if isempty(L_sei) == 0 % L_sei non è vuoto
            % Frattura
            if isempty(dadt) == 0 %dadt non è vuoto
                dudt= [dudt
                    (-1/2*param.V_sei*j_sei)*An_int/(An_int+2*param.N_cracks*dadt*param.w_crack)
                    dadt];
            else
                dudt= [dudt
                    (-1/2*param.V_sei*j_sei)];
            end
        end

        % CEI
        if isempty(L_cei)==0
            dudt=[dudt
                j_cei*param.V_cei];
        end

        % Plating
        if isempty(L_plating)==0
            dudt=[dudt
                -j_plating*param.V_plating];
        end

        % Add flux
        if isempty(L_sei)==0 || isempty(L_plating) == 0
            dudt=[dudt
                jn_int+j_sei+j_plating-jn_tot];
        end

        if isempty(L_cei)==0
            if jp_int >0
                dudt=[dudt
                    jp_int+j_cei-jp_tot];
            else
                dudt=[dudt
                    jp_int-j_cei-jp_tot];
            end
        end

        if isempty(L_sei)
            L_sei=0;
        end
        if isempty(L_cei)
            L_cei=0;
        end

        if isempty(L_plating)
            L_plating=0;
        end

        % Internal resistance
        R_int=param.Rf+(L_sei/param.sigma_sei)+ L_plating/param.sigma_plating+L_cei/param.sigma_cei+...
            param.k_LAM_p_resistance*(param.vol_fraction_solidphase_0(1)-param.vol_fraction_solidphase(1))/param.vol_fraction_solidphase_0(1)+...
            param.k_LAM_n_resistance*(param.vol_fraction_solidphase_0(2)-param.vol_fraction_solidphase(2))/param.vol_fraction_solidphase_0(2);

        dudt=[dudt
            U_p-U_n+eta_p-eta_n-R_int*(jn_tot*param.F*An_int)-...
            param.V_cutover
            jp_tot*(param.F*Ap_int)+jn_tot*(param.F*An_int)];

end

end
