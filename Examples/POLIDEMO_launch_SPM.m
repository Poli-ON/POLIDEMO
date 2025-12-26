%% POLIDEMO_launch_SPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POLIDEMO_launch_SPM.m
% This script reproduces the results presented in:
% https://doi.org/10.1016/j.apenergy.2025.126744
%
% The workflow demonstrates the use of POLIDEMO to simulate lithium-ion
% battery behavior during aging, including SPM solutions, and degradation
% indicators (LAM, LLI, irreversible swelling)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all

%% Add folder and subfolders to path

parentDir = fileparts(pwd);   % cartella madre della current folder
addpath(genpath(parentDir));

%% Define simulation parameters

global param

initialState=[];      % Define intialState
cycle_count=1;        % Define the initial number of cycles
Crate_discharge = -1/20;  % Define discharge rate
initial_SOC= 0;     % Define initial SOC

param=Parameters_init(Crate_discharge,cycle_count);

% Current density [A/m^2]
I       =  @(t) param.Crate*param.i_1C_density*ones(size(t));

% Simulation time and stop conditions
t0      = 0;                           % Initial time [s]
tf      = abs(3600/param.Crate*2.2);   % Final time [s]
time_step=10;                          % Time step [s]

% Simulation results
result = spm_simulation(t0,tf,time_step,initial_SOC,initialState,I);


%%
% CV charge
%param.OperatingMode = 2;

% Initial conditions for the following step
%initialState=result.initialState;

% Simulation results
%result2 =spm_simulation(t0,tf,time_step,initial_SOC,initialState,I);

%% Plot

figure()
plot(result.cap,result.V,'linewidth',1.5)
hold on
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Capaciy [Ah]', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('Voltage [V]','Fontsize',18.7,'Interpreter','Latex');

% Note: to get the results
figure()
plot(result.cap,result.thk,'linewidth',1.5)
hold on
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Capaciy [Ah]', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('Thickness change [mm]','Fontsize',18.7,'Interpreter','Latex');






