function [NCol,NRow,Nn,Nm,NNodes] = lOCATIONcALC(xen,yen,NumInc,NCol,...
    NRow,Nn,Nm,NNodes)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% lOCATIONcALC determines the location of a water particle in terms of 
% row, column, and node number.  The actual coordinate location is known.
%
% xen [NUMINC,1] = x-coordinate location for each particle.(Xe)
% yen [NUMINC,1] = y-coordinate location for each particle.(Ye)
% NumInc = number of index locations.
% NCol [NUMINC,1] = column index of particle location.
% NRow [NUMINC,1] = row index of particle location.
% Nn [NUMINC,1] = x-face boundary (greater) for particle location.
% Nm [NUMINC,1] = y-face boundary (greater) for particle location.
% NNodes [NUMINC,1] = volume index of particle location.

global DX DY NUMCOLS NUMROWS XINDEX YINDEX

% local variables.
BTemp = zeros(NumInc,1);        % Boolean.
CalcI1 = zeros(NumInc,1);       % Int Calc variable.
CalcI2 = zeros(NumInc,1);       % Int Calc variable.

% Calculate column locatin.
NCol = (floor(xen./DX)) + 1;
BTemp = (NCol > NUMCOLS);
CalcI1 = BTemp;
CalcI2 = (1 - BTemp);
NCol = (CalcI1.*NUMCOLS) + (CalcI2.*NCol);
% Calculate row location.
NRow = (floor(yen./DY)) + 1;
BTemp = (NRow > NUMROWS);
CalcI1 = BTemp;
CalcI2 = (1 - BTemp);
NRow = (CalcI1.*NUMROWS) + (CalcI2.*NRow);
% get value for max side of each volume.
Nn = (XINDEX(NCol + 1))';
Nm = (YINDEX(NRow + 1))';
% get node locations.
NNodes = (((NRow - 1).*NUMCOLS) + NCol);

clear BTemp CalcI1 CalcI2;
clear xen yen NumInc;
return;
%EOF