%% POLIDEMO WO CYCLE JUMP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POLIDEMO_wo_cycle_jump.m
% This script lauches POLIDEMO simulations withotuh cycle jump:
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

clear
close all
clc

%% Add folder and subfolders to path

parentDir = fileparts(pwd);   % cartella madre della current folder
addpath(genpath(parentDir));

%% Launch POLIDEMO withouth cycle jump

% Initialize simulations

n_cycles=100;         % Define the number of cycles
initialState=[];      % Define intialState
cycle_count=1;        % Define the initial number of cycles
Crate_charge = -1/2;  % Define charge rate
Crate_discharge = 1;  % Define discharge rate
initial_SOC= 0;       % Define initial SOC

% Launch simulations
result_cycle_no_cycle_jump(:,:)=run_n_cycles(n_cycles,initialState,...
                                             cycle_count,Crate_charge,...
                                             Crate_discharge,initial_SOC);

% Get results
for i=1:length(result_cycle_no_cycle_jump)
    cycle_number(i)=result_cycle_no_cycle_jump(i,3).cycle_count;
    capacity_fade(i)=(result_cycle_no_cycle_jump(i,3).cap(end)-...
                      result_cycle_no_cycle_jump(1,3).cap(end))/result_cycle_no_cycle_jump(1,3).cap(end);
    LLI(i)=result_cycle_no_cycle_jump(i,3).LLI(end); 
    LAM_n(i)=result_cycle_no_cycle_jump(i,3).LAM_n(end);
    LAM_p(i)=result_cycle_no_cycle_jump(i,3).LAM_p(end);
    irr_swelling(i)=result_cycle_no_cycle_jump(i,3).irr_swelling(end); 
end