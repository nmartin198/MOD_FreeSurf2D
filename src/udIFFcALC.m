function DTerm = udIFFcALC(DTerm,PCol,PRow,pp,NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% udIFFcALC calculates the diffusion terms that contribute to the Fux term.
% A centered difference stencil is employed to calculate the diffusion 
% contribution in two dimensions.  Two different eddy vicosities, DiffH, can be
% employed.  If the parameter SMAGC has a value other than zero, a Smagorinsky
% diffusion model will be used to determine the horizontal eddy viscosity.
% If SMAGC is set to zero, then the kinematic viscosity will be used for the
% horizontal eddy viscosity.
%
% STENCIL:
%_|____|____|____|_
% |    |    |    |          0 = particle location.
% |    y2   |    | NRow     y2 = yIndex2
%_|____|____|____|_         y4 = yIndex4
% |    |    |    |          x1 = xIndex1
% x3  x20   x1   | PRow     x3 = xIndex3
%_|____|____|____|_         x2 = xIndex2
% |    |    |    |
% |   y4    |    | NRow1
%_|____|____|____|_
% |    |    |    |
% NCol1 NCol
%
% Recieved:
% DTerm [NumInc,1] is the calculated diffusion contribution.
% PCol = Col [NumInc,1] is the departure point column location.
% PRow = Row [NumInc,1] is the departure point row location.
% NumInc = scalar. Number of indexes.
%
% Returned:
% DTerm [NumInc,1] is the calculated diffusion contribution.

global fdt DX DY EVIS NUMCOLS NUMROWS PREC SMAGC u XINC

%local variables.
BTemp = zeros(NumInc,1);               % Boolean calcualtion variable.
CalcI1 = zeros(NumInc,1);              % Int Calc variable.
CalcI2 = zeros(NumInc,1);              % Int Calc variable.
NCol = zeros(NumInc,1);                % Column location.
NCol1 = zeros(NumInc,1);               % Column location.
NRow = zeros(NumInc,1);                % Row location.
NRow1 = zeros(NumInc,1);               % Row location.
xIndex1 = zeros(NumInc,1);             % Velocity index.
xIndex2 = zeros(NumInc,1);             % Velocity index.
xIndex3 = zeros(NumInc,1);             % Velocity index.
yIndex2 = zeros(NumInc,1);             % Velocity index.
yIndex4 = zeros(NumInc,1);             % Velocity index.

% First get Column and Row locations.
BTemp = ((pp - 0.5) > PREC);
CalcI1 = BTemp;
CalcI2 = (1 - BTemp);
NCol = (CalcI1.*PCol) + (CalcI2.*(PCol+1));
NCol = cOLaDJ(NCol,NumInc);
NCol1 = (CalcI1.*(PCol - 1)) + (CalcI2.*PCol);
NCol1 = cOLaDJ(NCol1,NumInc);
NRow = (PRow - 1);
NRow = rOWaDJ(NRow,NumInc);
NRow1 = (PRow + 1);
NRow1 = rOWaDJ(NRow1,NumInc);
% Now calculate Indexes.
xIndex2 = ((PRow - 1).*XINC) + NCol;
xIndex3 = ((PRow - 1).*XINC) + NCol1;
xIndex1 = ((PRow - 1).*XINC) + NCol + 1;
yIndex2 = ((NRow - 1).*XINC) + NCol;
yIndex4 = ((NRow1 - 1).*XINC) + NCol;
% Calculate diffusion part of the term.
DTerm = (EVIS*fdt).*(((u(xIndex1) - 2*u(xIndex2) + u(xIndex3))/(DX^2)) + ...
   ((u(yIndex4) - 2*u(xIndex2) + u(yIndex2))/(DY^2)));

clear BTemp CalcI1 CalcI2 NCol NCol1 NRow NRow1 NumInc PCol PRow pp xIndex1;
clear xIndex2 xIndex3 yIndex2 yIndex4;
return;

%EOF