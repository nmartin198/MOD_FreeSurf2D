function [ValueX,ValueY] = vOLtOfACE(MainVal,ValueX,ValueY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vOLtOfACE allocates volume values to face depending on velocity
% direction.  So upwinds volume center values for faces.
%
% Received:
% 
% MainVal [NUMNODES,1] = value for volume center
% ValueX [NUMINCX,1] = value upwinded to x-faces.
% ValueY [NUMINCY,1] = value upwinded to y-faces.
%
% Returned;
%
% ValueX [NUMINCX,1] = value upwinded to x-faces.
% ValueY [NUMINCY,1] = value upwinded to y-faces.

[ValueX] = vOLtOXfACE(MainVal,ValueX);
[ValueY] = vOLtOYfACE(MainVal,ValueY);

clear MainVal;
return;
%EOF