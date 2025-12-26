function [OmdC_f]=LatParLFP_Lion(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATPARLFP_LION
% This function computes the integration int_0^c\Omega dc for LFP cathode 
% material when the partial molar volume Omega is a function of the lithium 
% concentration c. 
%
%
% OUTPUTS
% OmdC_f          Function of int_0^c\Omega dc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global param

xs=linspace(0,1,1000);
ll=length(xs);
%om=(0.8*0.062/param.cs_maxp);
% om=(0.062/(param.cs_maxp*(0.062+1)));
om=(0.048/((0.048+1))); %Aggiustamento  23/11/24


OmdC=xs.*om;
OmdC_f = @(x) polyval(polyfit(xs,OmdC,1),x);


end





