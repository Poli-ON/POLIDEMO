function OmdC_f=LatParNMC_Lion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATPARLCO_LION
% This function computes the integration int_0^c\Omega dc for LCO cathode 
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

%om=param.om_p;

%OmdC=xs.*param.cs_maxp.*om;
OmdC=(-1.1e-2)*(1-xs);
OmdC_f = @(x) polyval(polyfit(xs,OmdC,1),x);

end