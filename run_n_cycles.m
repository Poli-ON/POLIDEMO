function result = run_n_cycles(n_cycles,initialState,cycle_count,Crate_charge,Crate_discharge,initial_SOC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN_N_CYCLES
% Runs a specified number of charge-discharge cycles using the SPM model.
%
%       result = run_n_cycles(n_cycles, EnableAgingsimplified, initialState, cycle_count)
%
% INPUT:
%   n_cycles              - Number of charge/discharge cycles to simulate
%   initialState          - Vector of initial states
%   cycle_count           - Current number of previously simulated cycles
%   Crate_charge          - Charge rate 
%   Crate_discharge       - Discharge rate
%   initial_SOC           - Initial SOC for cyclic simulations
%
% OUTPUT:
%   result                - Structure containing the simulation results
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

%% Initialization

SOC=initial_SOC;
param=Parameters_init(Crate_charge,cycle_count);

%% Cyclic simulation

for i=1:n_cycles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CC charge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CC charge
    param.OperatingMode = 1;   

    % Crate 
    param.Crate   = Crate_charge;

    switch abs(param.Crate)
        case 0
            param.Dn       = 3.e-14;
        case 1/20
            param.Dn       = 3e-14;
        case 1/5
            param.Dn       = 1.5e-14;  %1.5e-14 DAVE IO: 0.6e-14
        case 1/2
            param.Dn       = 2e-14;    %2e-14 DAVE IO 1e-14
        case 1
            param.Dn       = 2.5e-14;
        case 2
            param.Dn       = 3.4e-14;
        case 3
            param.Dn       = 3.4e-14;
    end

    % Current density [A/m^2]
    I       =  @(t) param.Crate*param.i_1C_density*ones(size(t));

    % Simulation time and stop conditions
    t0      = 0;                           % Initial time [s]
    tf      = abs(3600/param.Crate*2.2);   % Final time [s]
    time_step=10;                          % Time step [s]            

    % Simulation results
    result1{i} = spm_simulation(t0,tf,time_step,SOC,initialState,I);
    result1{i}.cycle_count = cycle_count;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CV charge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CV charge
    param.OperatingMode = 2; 

    % Initial conditions for the following step
    initialState=result1{i}.initialState;

    % Simulation results
    result2{i} =spm_simulation(t0,tf,time_step,SOC,initialState,I);

    % Initial conditions for the following step
    initialState=result2{i}.initialState;
    result2{i}.cycle_count = cycle_count;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CC discharge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CC discharge
    param.OperatingMode = 1; 

    % Crate
    param.Crate   = Crate_discharge;

    switch abs(param.Crate)
        case 0
            param.Dn       = 3.e-14;
        case 1/20
            param.Dn       = 3e-14;
        case 1/5
            param.Dn       = 1.5e-14; %1.5e-14 DAVE IO: 0.6e-14
        case 1/2
            param.Dn       = 2e-14;   %2e-14 DAVE IO 1e-14
        case 1
            param.Dn       = 2.5e-14;
        case 2
            param.Dn       = 3.4e-14;
        case 3
            param.Dn       = 3.4e-14;
    end

    % Current density [A/m^2]
    I       =  @(t) param.Crate*param.i_1C_density*ones(size(t));

    t0      = 0;                         % Initial time [s]
    tf      = abs(3600/param.Crate*2.2);       % Final time [s]
    time_step=10;                        % Time step [s]

    result3{i}=spm_simulation(t0,tf,time_step,SOC,initialState,I);
    result3{i}.cycle_count=cycle_count;

    % Initial conditions for the following step
    initialState=result3{i}.initialState;

    % Results of simulation
    result(i,:)= [result1{i},result2{i},result3{i}];

    cycle_count=cycle_count+1;
    param.cycle_count=cycle_count;
end

end

