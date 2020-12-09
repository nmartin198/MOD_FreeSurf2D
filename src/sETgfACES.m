function [GvX,GvY] = sETgfACES(GvX,GvY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETgfACES sets the g* vector values.  Each vector index represents the
% corresponding face for the corresponding volume.  g* is the same as 
% G* from Casulli and Cheng (1992).

global BHux BHvy DX DY EtaXM1 EtaXP1 EtaYP1 EtaYM1 fdt Fux Fvy G GAMMATX
global GAMMATY Hux Hvy NUMINCX NUMINCY THETA UA VA

% local variables.
MultX = double(0.0);
MultY = double(0.0);
PartX1 = double(zeros(NUMINCX,1));
PartX2 = double(0.0);
PartX3 = double(zeros(NUMINCX,1));
PartY1 = double(zeros(NUMINCY,1));
PartY2 = double(0.0);
PartY3 = double(zeros(NUMINCY,1));

% Calculations.
% gX.
% gX = Part1 + Part2 - Part3;
%     Part1 = Hux.*Fux;
%     Part2 = (dt*GAMMATX*UA);
%     Part3 = ((1 - THETA)*G*(dt/DX)).*Hux.*(EtaXP1 - EtaXM1);
MultX = ((1 - THETA)*G*(fdt/DX));
PartX1 = Hux.*Fux;
PartX2 = (fdt*GAMMATX*UA);
PartX3 = MultX.*Hux.*(EtaXP1 - EtaXM1);
GvX = PartX1 + PartX2 - PartX3;
% Make sure that only wet x-faces have G values.
GvX = BHux.*GvX;
% gY.
% gY = Part1 + Part2 - Part3;
%     Part1 = Hvy.*Fvy;
%     Part2 = (dt*GAMMATY*VA);
%     Part3 = ((1 - THETA)*G*(dt/DY)).*Hvy.*(EtaYP1 - EtaYM1);
MultY = ((1 - THETA)*G*(fdt/DY));
PartY1 = Hvy.*Fvy;
PartY2 = (fdt*GAMMATY*VA);
PartY3 = MultY.*Hvy.*(EtaYP1 - EtaYM1);
GvY = PartY1 + PartY2 - PartY3;
% Make sure that only wet x-faces have G values.
GvY = BHvy.*GvY;

clear PartX1 PartX2 PartX3 MultX;
clear PartY1 PartY2 PartY3 MultY;
return;
%EOF