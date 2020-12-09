function [ValY] = vOLtOYfACE(MVal,ValY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vOLtOfACE allocates volume values to face depending on velocity
% direction.  So upwinds volume center values for faces.
%
% Received:
% 
% MVal [NUMNODES,1] = value for volume center
% ValY [NUMINCY,1] = value upwinded to y-faces.
%
% Returned;
%
% ValY [NUMINCY,1] = value upwinded to y-faces.

global NUMINCY VDirect VYM1 VYP1

% local variables.
BTempY = zeros(NUMINCY,1);          % boolean for y-faces.
CalcY1 = double(zeros(NUMINCY,1));  % Calc variable.
CalcY2 = double(zeros(NUMINCY,1));  % Calc variable.
ValYp1 = double(zeros(NUMINCY,1));  % Y-face value i,j+1/2.
ValYm1 = double(zeros(NUMINCY,1));  % Y-face value i,j-1/2.

% Distribute to vectors.
ValYp1 = MVal(VYP1);
ValYm1 = MVal(VYM1);
% Now allocate quanties to faces.
BTempY = (VDirect == 1);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
ValY = CalcY1.*ValYm1 + CalcY2.*ValYp1;

clear BTempY CalcY1 CalcY2 MVal ValYp1 ValYm1;
return;
%EOF