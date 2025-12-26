function [A,B]= spm_spatial_discretization(Dp,Dn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM_SPATIAL_DISCRETIZATION
% Returns the matrices for the spatial discretization of the diffusion
% equations in the SPM using a finite difference scheme.
%
%       [A, B] = spm_spatial_discretization(Dp, Dn)
%
% The system of ODEs has the form:
%       c_dot = A * c + B * u,
% where:
%   - c is the concentration vector,
%   - u is the input (intercalation flux),
%   - A and B are the discretized system matrices.
%
% INPUTS:
%   Dp  - Diffusion coefficient in the positive electrode [m^2/s]
%   Dn  - Diffusion coefficient in the negative electrode [m^2/s]
%
% OUTPUTS:
%   A   - Block-diagonal matrix from finite difference discretization
%   B   - Input matrix (flux contribution)
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

% Spatial discretization increment
delta_rp = 1 / param.N;
delta_rn = 1 / param.N;

% Main diagonal 
main_diagonal       = -2 * ones(param.N + 1, 1);
main_diagonal(1)    = -6;

% Lower diagonal 
i                   = 1:param.N - 1;
bottom_diagonal     = (i - 1) ./ i;
bottom_diagonal     = [bottom_diagonal, 2];

% Upper diagonal 
i                   = 1:param.N - 1;
up_diagonal         = (i + 1) ./ i;
up_diagonal         = [6, up_diagonal];

% Cathode matrix (Ap)
Ap = Dp / (param.Rp^2 * delta_rp^2) * ...
     (diag(main_diagonal) + diag(up_diagonal, 1) + diag(bottom_diagonal, -1));

% Anode matrix (An)
An = Dn / (param.Rn^2 * delta_rn^2) * ...
     (diag(main_diagonal) + diag(up_diagonal, 1) + diag(bottom_diagonal, -1));

A=[Ap;An];

% Cathode input vector (Bp)
Bp         = zeros(param.N + 1, 1);
Bp(end)    = -2 * (param.N + 1) / param.N / (delta_rp * param.Rp);

% Anode input vector (Bn)
Bn         = zeros(param.N + 1, 1);
Bn(end)    = -2 * (param.N + 1) / param.N / (delta_rn * param.Rn);

B=[Bp;Bn];

end