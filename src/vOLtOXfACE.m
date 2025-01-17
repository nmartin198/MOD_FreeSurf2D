function [ValX] = vOLtOXfACE(MVal,ValX)
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
% ValX [NUMINCX,1] = value upwinded to x-faces.
%
% Returned;
%
% ValX [NUMINCX,1] = value upwinded to x-faces.

global NUMINCX UDirect UXM1 UXP1

% local variables.
BTempX = zeros(NUMINCX,1);          % boolean for x-faces.
CalcX1 = double(zeros(NUMINCX,1));  % Calc variable.
CalcX2 = double(zeros(NUMINCX,1));  % Calc variable.
ValXp1 = double(zeros(NUMINCX,1));  % X-face value i+1/2,j
ValXm1 = double(zeros(NUMINCX,1));  % X-face value i-1/2,j

% Distribute to vectors.
ValXp1 = MVal(UXP1);
ValXm1 = MVal(UXM1);
% Now allocate quanties to faces.
BTempX = (UDirect == 1);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
ValX = CalcX1.*ValXm1 + CalcX2.*ValXp1;

clear BTempX CalcX1 CalcX2 MVal ValXp1 ValXm1;
return;
%EOF