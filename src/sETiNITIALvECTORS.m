function [FS,DCent,HeightX,HeightY] = sETiNITIALvECTORS(FS,DCent,HeightX,HeightY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETiNITIALvECTORS function sets the undisturbed water depth, HUX and HVY, and
% sets the free surface elevations for the domain.  At the end of this function,
% the Dirichlet source values are set for total water depth.
%
% FS [NUMNODES,1] = EtaNew or the free surface elevations for each volume in the
%                       domain.
% DCent [NUMNODES,1] = HNODE or the total water depth at the center of each volume.
% HeightX [NUMINCX,1] = HUX or the undisturbed water depth at each x-face.
% HeightY [NUMINCY,1] = HVY or the undistrubed water depth at each y-face.

global DATUM HDEPTH NUMINCX NUMINCY ZTopo ZtX ZtY

% Do the undisturbed water depths at the x-faces and y-faces first.
HeightX = DATUM - ZtX;
HeightY = DATUM - ZtY;
% Set the initial free surface.
FS = (HDEPTH + ZTopo) - DATUM;
DCent = DATUM - ZTopo;

return;
%EOF