function OmdC_f=LatParLCO_Lion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATPARLCO_LION
% This function computes the integral of the partial molar volume with respect
% to lithium concentration for the LCO (LiCoO2) cathode material:
%
%           ∫Ω(c) dc
%
% where Ω is the concentration-dependent partial molar volume.
%
% OUTPUTS
% OmdC_f      Value of the integrated partial molar volume as a function of
%             lithium concentration
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

%xf=[0.4,0.51,0.78,0.95,1];
xf=[0.4,0.53,0.76,0.95,1];
V=[0,0.18,-0.65,-1.67,-1.67]/100;


OmdC_f=@(x) ((xf(1)<=x) & (x<=xf(2))).*(((V(2)-V(1))/(xf(2)-xf(1))).*x-xf(1)*((V(2)-V(1))/(xf(2)-xf(1)))+V(1)) + ...
    ((xf(2)<=x) & (x<xf(3))).*(((V(3)-V(2))/(xf(3)-xf(2))).*x-xf(2)*((V(3)-V(2))/(xf(3)-xf(2)))+V(2)) + ...
    ((xf(3)<=x) & (x<=xf(4))).*(((V(4)-V(3))/(xf(4)-xf(3))).*x-xf(3)*((V(4)-V(3))/(xf(4)-xf(3)))+V(3))+...
    ((xf(4)<=x) & (x<=xf(5))).*(((V(5)-V(4))/(xf(5)-xf(4))).*x-xf(3)*((V(5)-V(4))/(xf(5)-xf(4)))+V(4));

end