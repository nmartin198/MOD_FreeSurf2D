function [VolXP1,VolXM1,VolYP1,VolYM1] = aLLOCaDJvOL(MainVal,...
   VolXP1,VolXM1,VolYP1,VolYM1)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% aLLOCaDJvOL allocates volume indexed vectors to vectors that indexed 
% so that they correspond to adjacent volumes for velocity definition
% locations.
%
% MainVal is the original volume locations. [NUMNODES,1].
% VolXP1 is at (i+1,j) for the x-face at (i+1/2,j).
% VolXM1 is at (i,j) for the x-face at (i+1/2,j).
% VolYP1 is at (i,j+1) for the x-face at (i,j+1/2).
% VolYM1 is at (i,j) for the x-face at (i,j+1/2).

global UXM1 UXP1 VYP1 VYM1

VolXP1 = MainVal(UXP1);
VolXM1 = MainVal(UXM1);
VolYP1 = MainVal(VYP1);
VolYM1 = MainVal(VYM1);

clear MainVal;
return;
%EOF