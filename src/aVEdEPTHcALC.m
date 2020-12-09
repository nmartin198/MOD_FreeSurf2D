function [BHTD,HTD,NSides] = aVEdEPTHcALC(BHTD,HTD,NSides,HTD1,HTD2,HTD3,HTD4)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% aVEdEPTHcALC calculates the average total depth for the node center location.
% In addition this function calculates the number of wet sides for each volume.
% This depth will be employed by plotting functions and by sediment transport
% calculations where the concentration is defined at the node center.
%
% Recieve:
% BHTD = BH [NUMNODES,1] enforces wet and dry faces.  Equals 1.0 for wet
%       0.0 for dry.
% HTD = H [NUMNODES,1] or total water depth defined at volume center.
% NSides = Sides [NUMNODES,1] or number of wet sides on each volume.
% HTD1 = h1 [NUMNODES,1] or the total water depth for side i+1/2,j.
% HTD1 = h2 [NUMNODES,1] or the total water depth for side i,j+1/2.
% HTD1 = h3 [NUMNODES,1] or the total water depth for side i-1/2,j.
% HTD1 = h4 [NUMNODES,1] or the total water depth for side i,j-1/2.
%
% Return
% BHTD = BH [NUMNODES,1] enforces wet and dry faces.  Equals 1.0 for wet
%       0.0 for dry.
% HTD = H [NUMNODES,1] or total water depth defined at volume center.
% NSides = Sides [NUMNODES,1] or number of wet sides on each volume.

global HCUTOFF HOld NUMNODES PREC

% local variables.
dph = double(zeros(NUMNODES,1));      % Sum of total depth for each side of each vol.
BTemp = zeros(NUMNODES,1);            % Boolean calculation variable.
Calc1 = double(zeros(NUMNODES,1));    % Calc. variable.
Calc2 = double(zeros(NUMNODES,1));    % Calc. variable.
Temp2 = double(zeros(NUMNODES,1));    % Calculation variable.
Temp3 = double(zeros(NUMNODES,1));    % Calculation variable.

% Save to old vector.
HOld = HTD;
% calculations.
NSides = ((HTD1 - double(0.0)) > HCUTOFF) + ((HTD2 - double(0.0)) > HCUTOFF) +...
   ((HTD3 - double(0.0)) > HCUTOFF) + ((HTD4 - double(0.0)) > HCUTOFF);
BTemp = ((NSides - double(0.0)) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
Temp2 = (Calc1.*NSides) + (Calc2.*PREC);
Temp3 = Calc1.*(1./Temp2);
dph = (HTD1 + HTD2 + HTD3 + HTD4);
HTD = dph.*Temp3;
BTemp = ((HTD - double(0.0)) > HCUTOFF);
BHTD = double(BTemp);
HTD = BHTD.*HTD;

clear dph BTemp Calc1 Calc2 Temp2 Temp3;
clear HTD1 HTD2 HTD3 HTD4;
return;
%EOF