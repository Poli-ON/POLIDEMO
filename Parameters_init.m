function param= Parameters_init(Crate,cycle_count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS_INIT
% Initializes and returns the model parameter structure.
%
%       param = parameters_init(Crate, cycle_count)
%
% INPUT:
%   Crate       - Charge/discharge rate (C-rate)
%   cycle_count - Current number of previously simulated cycles
%
% OUTPUT:
%   param       - Structure containing all model parameters
%
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

%% Simulation control

% Valid options are:
%                       0 - Disable thermal effects
%                       1 - Enable thermal effects

param.EnableTemperature=0;

% Valid options are:
%                       1 - Constant input current density (defined by the
%                       user)
%                       2 - Constant Voltage input

param.OperatingMode=1;


% Flags for degradation sub-models [Lsei, a, csei, plating]
% 0 = disabled, 1 = enabled
param.aging_vector=[0,0,0,0];

% Valid options are:
%                       0 - No SEI growth
%                       1 - Enable SEI growth (no fracture)
param.EnableSEI=1;

if param.EnableSEI==1
    param.aging_vector(1)=1;
end

% Valid options are:
%                       0 - No SEI on cracks
%                       1 - Enable SEI on cracks

param.EnableSEIoncracks=1;

% Add this control: if SEI is not enabled, also SEI on cracs in not
% enabled

if param.EnableSEI ==0
    param.EnableSEIoncracks=0;
end

% Valid options are:
%                       0 - Paris Law - constant number of cracks
%                       1 - Paris Law - number of cracks propotional to LAM
%                       2 - Cracked area proportional to LAM

param.CrackMechanism=2;

if param.EnableSEIoncracks==1 && param.CrackMechanism ~= 2
    param.aging_vector(2)=1;
end

% Valid options are:
%                       0 - No LAM
%                       1 - Enable LAM
param.EnableLAM=1;

% Valid options are:
%                       0 - No CEI growth
%                       1 - Enable CEI growth (no fracture)
param.EnableCEI=0;

if param.EnableCEI ==1
    param.aging_vector(3)=1;
end

% Valid options are:
%                       0 - No lithium plating
%                       1 - Enable lithium plating
param.Enableplating=0;

if param.Enableplating == 1
    param.aging_vector(4)=1;
end

% Valid options are:
%                       0 - No SOC control
%                       1 - SOC control

param.EnableSOCcontrol=0;

%% Spatial discretization

param.N=10;

%% Cycle jump

% Cycle count
param.cycle_count=cycle_count;

% Tolerance for cycle jump scheme
param.tol_jump_size=0.0008; % Lower tolerance means higher precision but higher computational cost

%% Universal constant

% Faraday Constant  [C/mol]
param.F         = 96487;
% Gas constant      [J / (mol K)]
param.Rg        = 8.314;

%% Voltage limits

% Cut-off voltage [V]
param.V_cutoff = 3;

% Cut-over voltage [V]
param.V_cutover = 4.2;

%% Current rate

% Initial SOC
param.initial_SOC=0;

% Crate
param.C_rate_charge=-1/2;
param.Crate_discharge=1;

% Time for a cycle
param.tn=3600/abs(param.C_rate_charge)+3600/abs(param.Crate_discharge);

%% Cell geometry

% Length of puch cell [m]
param.pouch_length 	= 174.5e-3;

%  Width of cell pouch [m]
param.pouch_width  	= 94.5e-3;

%  Thickness of cell pouch [m]
param.pouch_thickness=9.2e-3;

% Pouch Wrapper
param.len_pouch = 9.2e-3; %110e-6; % ref: "Li-Ion Pouch Cells for Vehicle Applications—Studies of Water Transmission and Packing Materials", Pontus Svens, Maria Hellqvist Kjell, Carl Tengstedt, Göran Flodberg and Göran Lindbergh, Energies 2013, 6, 400-410; doi:10.3390/en6010400

% Electrode area [m^2]
param.pouch_area    = param.pouch_length*param.pouch_width;

% Number of elmentary cells [-]
param.no_of_layers = 40;

% Stiffness of the case [1/Pa]
param.kz = 5284.89604*10^(-6);

%% Current density

% Input of Crate
param.Crate=Crate;

% Cell capacity [Ah]
param.assumed_cell_capacity_Ah = 22;

% Current to the single cell per unit of area [A/m^2]
param.i_1C_density=param.assumed_cell_capacity_Ah/(param.pouch_area*param.no_of_layers);

%param.i_1C_density=1.07*20; % CON QUESTO FUNZIONAVA

% Film resistance [Ohm]
param.Rf=40e-4;

%% Temperature Settings

% Environment (ambient) temperature [K]
param.Tref   = 21.45 + 273.15;

if param.EnableTemperature == 0
    param.T = param.Tref;
end

%% Density

% LiPF6 electrolyte [kg/m^3]
param.rho_LiPF6 = 1200; % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior Using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173, also supported by 'Thermal analysis of a cylindrical lithium-ion battery', Xiongwen Zhang, Electrochimica Acta,  56 (2011) 1246–1255, Table 3

% Pouch Material [kg/m^3]
param.rho_pouch = 2470; % 'Modeling for the scale-up of a lithium-ion polymer battery',Ui Seong Kim, Chee Burm Shin, Chi-Su Kim, Journal of Power Sources, 2008

% Filler/Binder [kg/m^3]
param.rho_pvdf = 1780;  % Table 1 from 'Characterization of Lithium-Ion Battery Thermal Abuse Behavior using Experimental and Computational Analysis',doi: 10.1149/2.0751510jes, J. Electrochem. Soc. 2015 volume 162, issue 10, A2163-A2173

% Aluminium current collector [kg/m^3]
param.rho_al = 2702;

% Positive electrode [kg/m^3]
param.rho_p  = 5180;

% Separator [kg/m^3]
param.rho_s  = 1200;

% Negative electrode [kg/m^3]
param.rho_n  = 2260;

% Copper current collector [kg/m^3]
param.rho_cu = 8933;

%% Sections thickness

% Positive Electrode [m]
param.len_p = 68e-6;

% Negative Electrode [m]
param.len_n = 88e-6;

% Separator [m]
param.len_s = 50e-6;

% Aluminium current collector [m]
param.len_al = 10e-6;

% Copper current collector [m]
param.len_cu = 10e-6;

%% Solid diffusion coefficients

% Positive electrode [m^2 / s]
param.Dp       = 6e-14;

% Negative electrode [m^2 / s]
switch abs(Crate)
    case 0
        param.Dn       = 3.e-14;
    case 1/20
        param.Dn       = 3e-14;
    case 1/5
        param.Dn       = 1.5e-14;
    case 1/2
        param.Dn       = 2e-14;
    case 1
        param.Dn       = 2.5e-14;
    case 2
        param.Dn       = 3.4e-14;
    case 3
        param.Dn       = 3.4e-14;
end

%% Stoichiometry limits [-]

param.theta_max_pos = 0.46;  % at 100% cell SOC
param.theta_max_neg = 0.7859;  % at 100% cell SOC
param.theta_min_pos = 0.91;  % at 0% cell SOC
param.theta_min_neg = 0.0075;  % at 0% cell SOC

%% Maximum concentration of Li-ions in the solid phase

% Positive electrode [mol/m^3]
param.cs_maxp   = 51597;

% Negative electrode [mol/m^3]
param.cs_maxn   = 28866;

%% Active material volume fraction

param.eps_fi = [0.028, 0.056];

% Positive electrode [-]
vol_fraction_solidphase_p = 0.76;

% Negative electrode [-]
vol_fraction_solidphase_n = 0.66;

if param.cycle_count == 1
    param.vol_fraction_solidphase = [vol_fraction_solidphase_p;vol_fraction_solidphase_n];
    param.vol_fraction_solidphase_0 = param.vol_fraction_solidphase;
end

% Volume fraction porosity
param.eps_p =  1 - param.vol_fraction_solidphase(1) - param.eps_fi(1);
param.eps_s = 0.55;
param.eps_n = 1 - param.vol_fraction_solidphase(2)- param.eps_fi(2);


%% Solid particle radius

% Positive electrode [m]
param.Rp     = 2e-6;

% Positive electrode [m]
param.Rn     = 10.e-6;

%% Particle surface area

% Positive electrode [m^2/m^3]
if param.cycle_count == 1
    param.as_p       = 3*param.vol_fraction_solidphase(1)/param.Rp;
    param.as_p0= param.as_p;
end

%Negative electrode [m^2/m^3]
if param.cycle_count == 1
    param.as_n       = 3*param.vol_fraction_solidphase(2)/param.Rn;
    param.as_n0= param.as_n;
end

%% Reaction rate constants

% Positive electrode [m^2.5/(mol^0.5s)]
param.kp       =  1e-11;

% Negative electrode [m^2.5/(mol^0.5s)]
param.kn       = 1e-11;

%%  Average concentration of Li-ions in the electrolyte

% Average electrolyte concentration [mol/m^3]
param.ce = 2000;

%% Thermal properties

% Battery specific heat capacity [J/(kg K)]
param.cp_al   = 903;        % Aluminium current collector
param.cp_p    = 700;        % Positive Electrode
param.cp_s    = 2050;       % Separator
param.cp_n    = 1437;       % Negative Electrode
param.cp_cu   = 385;        % Copper current collector
param.cp_LiPF6 =  2050;     % Electrolyte

% Assumption: Ignoring cp of binder/filler since they are negligible in content
% furthermore, the exterior pouch is also ignored in cp calculations
% (but it is accounted for in mass calculations)

if param.EnableTemperature == 1
    % Weighted calculation based on the constituents of the pouch material
    [param.mass_cell,param.cp_cell]=get_average_cell_heat_capacity;

    %% Convective properties

    % Area for convenction exchange flux [m^2]
    param.A_conv = (2*param.pouch_length*param.pouch_width+...
        2*param.pouch_width*param.pouch_thickness+...
        2*param.pouch_thickness*param.pouch_length);

    % Convective heat exchange coefficient [W/(m^2K)]
    param.h_conv = 0.001344*param.mass_cell*param.cp_cell/param.A_conv;

    % Commento Fra: questo coefficiente fa tornare le cose
    param.h_conv = 28;

end

%% Activation Energy

% Activation Energy for Temperature Dependent Solid Phase Diffusion [J/mol]

param.Ea_Dp  = 35e3;    % Positive electrode
param.Ea_Dn  = 35e3;    % Negative Electrode

% Activation Energy for Temperature Dependent Reaction Constant [J/mol]

param.Ea_kp  = 20e3;    % Positive electrode
param.Ea_kn  = 20e3;    % Negative electrode

%% Aging paramaters

% Concentration of EC in bulk electrolyte [mol/m^3]
param.c_solv_electrolyte=4541;

% SEI layer diffusivity [m^2/s]
param.D_solv=4e-16; % Rimpicciolisci se vuoi crescita esponenziale (prima era 4e-16)

% SEI kinetic rate constant [m/s]
param.k_sei=2.7e-18;

% SEI Molar Volume [m^3/mol]
param.V_sei=9.6e-5;

% Open circuit voltage for SEI side reaction [V]
param.U_sei          = 0.4;

% SEI conductivity [S/m]
param.sigma_sei= 1/30000;

% Stoichiometric coefficient for SEi reaction
param.nSEI=1;

% CEI kinetic rate constant [m/s]
param.k_cei=20e-16;

% Open circuit voltage for CEI side reaction [V]
param.U_cei=4.1;

% CEI conductivity [S/m]
param.sigma_cei=1/30000;

% CEI Molar Volume [m^3/mol]
param.V_cei=0.162/1.69e3;

% Plating
param.k_plating=0.4*2.54e-9;

% Open circuit voltage for plating side reaction [V]
param.U_plating=0;

% Plated lithium conductivity [S/m]
param.sigma_plating=1/30000;

% Plated lithium molar Volume [m^3/mol]
param.V_plating=1.30e-5;

%% Mechanical parameters

% Poisson ratio [-]
param.nu_p=0.24;
param.nu_n=0.3;

% Partial molar volume [m^3/mol]
param.om_p=-(0.0232/param.cs_maxp);
param.om_n=4.253e-6;

% Get Omega*dC (Integration of partial molar volume)
param.OmdC_pf=LatParLCO_Lion;
param.OmdC_nf=LatParGra_Lion;

% Young modulus [Pa]
param.Ep=191e9;
param.En=15e9;

%% Fracture parameters

% Paris coefficients
param.C=28.5e-7; %9e-6
param.m=2.15; % 3.1

% Initial crack length [m]
param.a_0=1.e-7;

% Crack width [m]
param.w_crack=1.5e-8;

% Crack density [m^(-2)]
param.rho_crack=4e20;

% Number of cracks (divido per pouch area perchè la corrente è
% normalizzata)
if param.cycle_count == 1
    switch param.CrackMechanism
        case 0
            param.N_cracks=4*pi*param.Rn^2*param.rho_crack;
        case 1
            param.k_crack=6.8*10^11;
            param.N_cracks=param.k_crack*(param.vol_fraction_solidphase_0(2));
        case 2
            param.k_crack=1.1;
            param.A_crack=param.k_crack*(param.vol_fraction_solidphase_0(2));
            param.A_crack=0; % no fracture at the fresh state
    end
end

if param.cycle_count == 1

    % Initial maximum and minimum circumferential stress [MPa]
    param.sigma_c_p_max=0;
    param.sigma_c_p_min=0;
    param.sigma_c_n_max=0;
    param.sigma_c_n_max_distribution=zeros(param.N+1,1);
    param.sigma_c_n_max_distribution_crack=zeros(param.N+1,1);
    param.sigma_c_n_min=0;
    param.sigma_c_n_min_distribution=zeros(param.N+1,1);
    param.delta_sigma_crack=0;

    % Initial maximum and minimum hoop stress [Pa]
    param.sigma_h_p_max=0;
    param.sigma_h_p_min=0;
    param.sigma_h_n_max=0;
    param.sigma_h_n_min=0;

end

% Initial cyclable lithium
if param.cycle_count ==1
    Cp=param.vol_fraction_solidphase(1)*param.len_p*param.cs_maxp*param.F*param.pouch_area/3600;
    Cn=param.vol_fraction_solidphase(2)*param.len_n*param.cs_maxn*param.F*param.pouch_area/3600;
    param.n_Li_initial=(param.theta_min_pos*Cp+param.theta_min_neg*Cn)*3600/(param.F);
end

%% LAM Parameters

param.a_LAM_p= 2.69e-6;    %2*90.6e-17*(1)
param.a_LAM_n= 12e-16;        %11.4e-16

param.n_LAM_p = 0;  %0.48
param.n_LAM_n = 5.55;  %  5.56

param.m_LAM_p = 0.19;    %  1.5*(1)
param.m_LAM_n =  -0.7;  % (-0.7)

param.k_LAM_p_resistance=0.05; %0.05
param.k_LAM_n_resistance=0.15;  %0.15


end