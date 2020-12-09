function [Index1,Index2,Index3,Index4,newp] = yINDEXfIND(Col,Row,pp,...
   NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% yINDEXfIND finds the indices to be employed with v-velocity terms for
% bilinear interpolation.  For this function to work, the coordinate location
% for each particle needs to be known.  This function then determines the
% stencil to be employed in the bilinear calculation.
%       +
%  +3+++++++2++++++++        + = volume boundaries.
%   |   +   |   +            | = bilinear interpolation boundaries 
%   |   +   |   +                      
%   |   +   |   +
%  +4+++++++1++++++++      ->  positive directions.
%       +                \|/
%       +
% Col is the column location of each particle.
% Nodes is the node number of each particle.
% pp is the weight determined by the distance between n and X-particle
%    location coordinate.
%
% yIndex1, yIndex2, yIndex3, and yIndex4 are the indexes represented in
%     the schematic above.
% vp holds the horizontal weights to be employed in the v -
%     velocity bilinear interpolation.

global NUMCOLS PRECH

%local variables.
BTemp = zeros(NumInc,1);         % Boolean calculation variable.
Calc1 = double(zeros(NumInc,1)); % Calc variable.
Calc2 = double(zeros(NumInc,1)); % Calc variable.
CalcI1 = zeros(NumInc,1); % Int Calc variable.
CalcI2 = zeros(NumInc,1); % Int Calc variable.
Index1 = zeros(NumInc,1);
Index2 = zeros(NumInc,1);
Index3 = zeros(NumInc,1);
Index4 = zeros(NumInc,1);
TempCol = zeros(NumInc,1);
TempCol1 = zeros(NumInc,1);
newp = double(zeros(NumInc,1));

% Calculations.
BTemp = ((pp - 0.5) > PRECH);
CalcI1 = BTemp;
CalcI2 = (1 - BTemp);
TempCol = (CalcI1.*Col) + (CalcI2.*(Col + 1));
TempCol = cOLaDJ(TempCol,NumInc);
TempCol1 = (CalcI1.*(Col - 1)) + (CalcI2.*Col);
TempCol1 = cOLaDJ(TempCol1,NumInc);
Index1 = (((Row - 1).*NUMCOLS) + TempCol) + NUMCOLS;
Index4 = (((Row - 1).*NUMCOLS) + TempCol1) + NUMCOLS;
Index2 = (((Row - 1).*NUMCOLS) + TempCol);
Index3 = (((Row - 1).*NUMCOLS) + TempCol1);
% Now change the p weights to reflect the v velocity stencil.
Calc1 = double(CalcI1);
Calc2 = double(CalcI2);
newp = (Calc1.*(pp - 0.5)) + (Calc2.*(pp + 0.5));

clear BTemp Calc1 Calc2 CalcI1 CalcI2 Col NumInc pp Row TempCol TempCol1;
return;
%EOF