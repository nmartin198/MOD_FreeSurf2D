function DTerm = vdIFFcALC(DTerm,PCol,PRow,pq,NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vdIFFcALC calculates the diffusion terms that contribute to the Fvy term.
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
% |    |    |    |          y2 = yIndex2
%_|____|_y2_|____|_NRow2    y4 = yIndex4
% |    |    |    |          y1 = yIndex1
% |    | 0  |    | PRow     x3 = xIndex3
%_|_x3_|_y1_|_x1_|_NRow     x1 = xIndex1
% |    |    |    |
% |    |    |    | 
%_|____|_y4_|____|_NRow1
% |    |    |    |
% NCol1 PCol NCol
%
% Recieved:
% DTerm [NumInc,1] is the calculated diffusion contribution.
% PCol = Col [NumInc,1] is the departure point column location.
% PRow = Row [NumInc,1] is the departure point row location.
% NumInc = scalar. Number of indexes.
%
% Returned:
% DTerm [NumInc,1] is the calculated diffusion contribution.

global fdt DX DY EVIS NUMCOLS NUMROWS PREC SMAGC v

% local variables.
BTemp = zeros(NumInc,1);               % Boolean calcualtion variable.
Calc1 = double(zeros(NumInc,1));       % Calc variable.
Calc2 = double(zeros(NumInc,1));       % Calc variable.
NCol = zeros(NumInc,1);
NCol1 = zeros(NumInc,1);
NRow = zeros(NumInc,1);
NRow1 = zeros(NumInc,1);
NRow2 = zeros(NumInc,1);
xIndex1 = zeros(NumInc,1);
xIndex3 = zeros(NumInc,1);
yIndex1 = zeros(NumInc,1);
yIndex2 = zeros(NumInc,1);
yIndex4 = zeros(NumInc,1);

% First get column and row locations.
BTemp = ((pq - 0.5) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
NRow = (Calc1.*(PRow - 1)) + (Calc2.*(PRow));
NRow = rOWaDJ(NRow,NumInc);
NRow1 = NRow + 1;
NRow1 = rOWaDJ(NRow1,NumInc);
NRow2 = NRow - 1;
NRow2 = rOWaDJ(NRow2,NumInc);
NCol = PCol + 1;
NCol = cOLaDJ(NCol,NumInc);
NCol1 = PCol - 1;
NCol1 = cOLaDJ(NCol1,NumInc);
% Now set indexes.
yIndex1 = (NRow.*NUMCOLS) + PCol;
yIndex2 = (NRow2.*NUMCOLS) + PCol;
yIndex4 = (NRow1.*NUMCOLS) + PCol;
xIndex3 = (NRow.*NUMCOLS) + NCol1;
xIndex1 = (NRow.*NUMCOLS) + NCol;
% Calculate diffusion part of the term.
DTerm = (EVIS*fdt).*(((v(xIndex1) - 2*v(yIndex1) + v(xIndex3))/(DX^2)) + ...
   ((v(yIndex4) - 2*v(yIndex1) + v(yIndex2))/(DY^2)));

clear BTemp Calc1 Calc2 NCol NCol1 NRow NRow1 NRow2 NumInc PCol pq PRow xIndex1;
clear xIndex3 yIndex1 yIndex2 yIndex4;
return;

%EOF