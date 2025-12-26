function OmdC_f=LatParGra_Lion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATPARLCO_LION
% This function computes the integral of the partial molar volume with respect
% to lithium concentration for the graphite material:
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

%%%%%%%%     INPUT     %%%%%%%%
% cond:  0 scarica; 1 carica
% xs:    indice x
% Crate: Current rate (per fasi scarica)
if param.Crate>0 % Per me scarica è positiva
    cond=0;
    xf=[0 0.25 0.47 1]; %0.24 0.48
else
    cond=1;
    xf=[0 0.24 0.48 1]; %0.24 0.48
end
xs=linspace(0,1,1000);

switch cond
    case 0 %SCARICA
        switch param.Crate
            case 1/20
                % V=[0 0.058 0.059 0.130]; % SCHWEIDLER
                V=[0 0.053 0.0545 0.1240]; % DIDIER
            case 1/5
                % V=[0 0.057 0.061 0.130]; % SCHWEIDLER
                V=[0 0.05 0.0545 0.1240]; % DIDIER
            case 1/2
                % V=[0 0.053 0.061 0.130]; % SCHWEIDLER
                V=[0 0.0465 0.0545 0.1240]; % DIDIER
            case 1
                % V=[0 0.048 0.061 0.130]; % SCHWEIDLER
                V=[0 0.044 0.0545 0.1240]; % DIDIER
            case 2
                % V=[0 0.047 0.061 0.130]; % SCHWEIDLER
                V=[0 0.0430 0.0545 0.1240]; % DIDIER
            case 3
                % V=[0 0.045 0.061 0.130]; % SCHWEIDLER
                V=[0 0.0415 0.0545 0.1240]; % DIDIER

        end

    case 1 % CARICA
        V=[0 0.0415 0.051 0.124]; % DIDIER
end


OmdC_f=@(x) ((xf(1)<=x) & (x<=xf(2))).*((V(2)-V(1))/(xf(2)-xf(1))).*x + ((xf(2)<=x) & (x<xf(3))).*(((V(3)-V(2))/(xf(3)-xf(2))).*x-xf(2)*((V(3)-V(2))/(xf(3)-xf(2)))+V(2)) + ((xf(3)<=x) & (x<=xf(4))).*(((V(4)-V(3))/(xf(4)-xf(3))).*x-xf(3)*((V(4)-V(3))/(xf(4)-xf(3)))+V(3));

end