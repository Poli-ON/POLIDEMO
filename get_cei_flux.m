function j_cei=get_cei_flux(ur,L_cei,jp_int)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET_CEI_FLUX (Beta version)
% This function computes the flux of the CEI (Cathode Electrolyte Interphase)
% formation reaction at the cathode.
%
%       j_cei = get_cei_flux(ur, L_cei, jp_int)
%
% INPUTS
% ur       : Unknown variable of the system
% jp_int   : Lithium intercalation flux at the cathode
% L_cei    : Thickness of the CEI layer
%
% OUTPUTS
% j_cei    : CEI reaction flux
%
%   Beta version: This function is under development and may be updated in
%   future releases.
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

% Reaction rate constant [m^2.5/(mol^0.5 s)] - Anode
kp=param.kp;                    

% Surface lithium concentration [mol/m^3] - Cathode
cs_p = ur(1,:);                             

% Exchange current density [A/m^2]
i0_p = kp*param.F*sqrt(param.ce)*sqrt(cs_p(end))*sqrt(param.cs_maxp-cs_p(end)) ;

% Anode overpotential [V]
eta_p = 2*param.Rg*param.T/param.F*asinh(jp_int*param.F/(2*i0_p));

% OCV [V]
[U_p,U_n] =  openCircuitPotential(ur);

%Crea una maschera logica per gli elementi che soddisfano la condizione
mask = (U_p + eta_p) > 0;

% Prealloca i vettori di output
eta_cei = zeros(size(jp_int));  % Stesso size di jp_int per compatibilità
j_cei = zeros(size(jp_int));    % Stesso size di jp_int per compatibilità

valid_idx = find(mask); % calcola solo per i valori in cui mask >1 (condizione verificata)

if isempty(valid_idx) ==0
for i = valid_idx'
    eta_cei(i) = get_eta_cei(ur(:,i), jp_int(i));
    j_cei(i) = param.k_cei * param.ce * exp(-0.5 * param.F / (param.Rg * param.T) * eta_cei(i));
end

% eta_cei = get_eta_cei(ur, jp_int);
% j_cei = -param.k_cei * param.ce * exp(0.5 * param.F / (param.Rg * param.T) * eta_cei);
end

end