function [cs_p, cs_n] = compute_concentration(sol,time_sol,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE_CONCENTRATION
% This function computes the lithium concentration at each point of the spatial
% discretization at a specified time "t".
%
%       [cs_p, cs_n] = compute_concentration(sol, time_sol, t)
%
% INPUTS
% sol       : Structure containing the solution computed with pdepe
% time_sol  : Time points at which the solution is computed
% t         : Time at which the concentration has to be computed
%
% OUTPUTS
% cs_p      : Lithium concentration in the cathode
% cs_n      : Lithium concentration in the anode
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

index=find(time_sol==t);

cs_p=sol(index,1:param.N+1);
cs_n=sol(index,param.N+2:2*param.N+2);

end
    