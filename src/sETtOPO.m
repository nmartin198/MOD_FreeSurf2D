function [TopoX,TopoY] = sETtOPO(Topo,TopoX,TopoY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETtOPO interpolates topographic elevations from the volume centered
% topography file to the faces of each volume.  The interpolation 
% method is simple linear interpolation.
%
% Received:
%
% Topo [NUMNODES,1] = topographic elevations at volume centers.
% TopoX [NUMINCX,1] = topographic elevations at x-faces.
% TopoY [NUMINCY,1] = topographic elevations at y-faces.
%
% Returned:
% TopoX [NUMINCX,1] = topographic elevations at x-faces.
% TopoY [NUMINCY,1] = topographic elevations at y-faces.

global UXM1 UXP1 VYM1 VYP1

% calculations.
TopoX = 0.5.*(Topo(UXM1) + Topo(UXP1));
TopoY = 0.5.*(Topo(VYM1) + Topo(VYP1));

clear Topo;
return;
%EOF