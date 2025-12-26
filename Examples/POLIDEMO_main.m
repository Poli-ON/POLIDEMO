%% POLIDEMO MAIN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POLIDEMO_main.m
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

%% Launch POLIDEMO with cycle jump

% Initialize simulations

n_cycles=1200;        % Define the number of cycles
initialState=[];      % Define intialState
cycle_count=1;        % Define the initial number of cycles
Crate_charge = -1/2;  % Define charge rate
Crate_discharge = 1;  % Define discharge rate
initial_SOC= 0;       % Define initial SOC
delta_N=50;           % Define the jump length

% Launch the first three cycles
result_cycle_jump(:,:)=run_n_cycles(3,initialState,...
                                             cycle_count,Crate_charge,...
                                             Crate_discharge,initial_SOC);
count=3;
cycle_count=result_cycle_jump(count,3).cycle_count+1; % Update the number of cycles

% Launch cyclic simulations
while cycle_count < 1200
    [initialState,cycle_count]=SPM_cycle_jump(delta_N,result_cycle_jump(count-2:count,:),cycle_count);
    count=count+3;
    result_cycle_jump(count-2:count,:)=run_n_cycles(3,initialState,...
                                             cycle_count,Crate_charge,...
                                             Crate_discharge,initial_SOC);
    if result_cycle_jump(count,3).cap(end) < 12
        disp('Capacity too low – stopping simulation');
        break
    end
    cycle_count=result_cycle_jump(count,3).cycle_count+1;
end

%% Get results

for i=1:length(result_cycle_jump)
    cycle_number(i)=result_cycle_jump(i,3).cycle_count;
    capacity_fade(i)=(result_cycle_jump(i,3).cap(end))/result_cycle_jump(1,3).cap(end);
    LLI(i)=result_cycle_jump(i,3).LLI(end); 
    LAM_n(i)=result_cycle_jump(i,3).LAM_n(end);
    LAM_p(i)=result_cycle_jump(i,3).LAM_p(end);
    irr_swelling(i)=result_cycle_jump(i,3).irr_swelling(end); % calcolato alla fine della caricae
end

%% Comparison with experimental results

% Load experimental data (complete aging dataset available on 10.5281/zenodo.14914172)
load('aging_data_set.mat')

figure()
plot(cycle_number,capacity_fade,'linewidth',1.5)
hold on
plot(aging_data_set.capacity_fade(1,:),aging_data_set.capacity_fade(2,:),'o','linewidth',1.2,'MarkerFaceColor',...
     'w','MarkerSize',7)
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Number of cycles', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('SOH [-]','Fontsize',18.7,'Interpreter','Latex');
grid on
legend('Model - RPT High current','Experimental - RPT High current',...
       'Interpreter','Latex','Location','Southwest')
xlim([0,1200])

figure()
plot(cycle_number,LLI,'linewidth',1.5)
hold on
plot(aging_data_set.LLI(1,:),aging_data_set.LLI(2,:),'o','linewidth',1.2,'MarkerFaceColor',...
     'w','MarkerSize',7)
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Number of cycles', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('LLI [$\%$]','Fontsize',18.7,'Interpreter','Latex');
grid on
legend('Model','Experimental','Interpreter','Latex','Location','Northwest')
xlim([0,1200])

figure()
plot(cycle_number,LAM_p,'linewidth',1.5)
hold on
plot(aging_data_set.LAM_p(1,:),aging_data_set.LAM_p(2,:),'o','linewidth',1.2,'MarkerFaceColor',...
     'w','MarkerSize',7)
legend('Model','Experimental','Interpreter','Latex','Location','Northwest')
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Number of cycles', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('LAM$_p$ [$\%$]','Fontsize',18.7,'Interpreter','Latex');
grid on
legend('Model','Experimental','Interpreter','Latex','Location','Northwest')
xlim([0,1200])

figure()
plot(cycle_number,LAM_n,'linewidth',1.5)
hold on
plot(aging_data_set.LAM_n(1,:),aging_data_set.LAM_n(2,:),'o','linewidth',1.2,'MarkerFaceColor',...
     'w','MarkerSize',7)
legend('Model','Experimental')
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Number of cycles', 'Fontsize',18.7,'Interpreter','Latex');
ylabel('LAM$_n$ [$\%$]','Fontsize',18.7,'Interpreter','Latex');
grid on
legend('Model','Experimental','Interpreter','Latex','Location','Northwest')
xlim([0,1200])

figure()
plot(cycle_number,irr_swelling,'-','linewidth',1.5,'MarkerFaceColor',...
     'w')
hold on
hold on
plot(aging_data_set.irreversible_swelling(1,:),aging_data_set.irreversible_swelling(2,:),'o','linewidth',1.2,'MarkerFaceColor',...
     'w','MarkerSize',7)
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',17);
xlabel('Number of cycles','Fontsize',18.7,'Interpreter','Latex');
ylabel('Irreversible swelling [mm]','Fontsize',18.7,'Interpreter','Latex');
grid on
legend('Model','Experimental','Interpreter','Latex','Location','Northwest')
xlim([0,1200])

 