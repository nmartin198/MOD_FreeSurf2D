function [f1,f2,f3,f4] = aLLOCfACE(ValX,ValY,f1,f2,f3,f4)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% aLLOCfACE distributes the domain wide face indexed vectors for x
% and y faces to four vectors corresponding to particular faces for
% each volume.  The x-face vector is distributed to east (i+1/2,j) and
% west (i-1/2,j) faces for each volume.  The y-face vector is 
% allocated to north (i,j+1/2) and south (i,j-1/2) face vectors for
% each volume.
%
% ValX = x-face values. [NUMINCX,1].
% ValY = y-face values. [NUMINCY,1].
% f1 = east face values. [NUMNODES,1].
% f2 = south face values. [NUMNODES,1].
% f3 = west face values. [NUMNODES,1].
% f4 = north face values. [NUMNODES,1].

global X1 X3 Y2 Y4

% Calculations.
f1 = ValX(X1);
f2 = ValY(Y2);
f3 = ValX(X3);
f4 = ValY(Y4);

clear ValX ValY;
return;
%EOF