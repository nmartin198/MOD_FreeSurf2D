function [XVel,YVel] = vELcALC(XVel,YVel)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vELcALC calculates the new velocities using the updated free surface values.
% The velocity calculation is explicit.
%
% Received and Returned:
%
%  XVel = u [NUMINCX,1] x-face velocities.
%  YVel = v [NUMINCY,1] y-face velocities.

global aXDen aYDen DX DY EtaXP1 EtaXM1 EtaYP1 EtaYM1 fdt G gX gY Hux
global Hvy PRECH NUMINCX NUMINCY RADVELBC THETA UOld VELDIRCBC VOld

%local variables.
BTempX = zeros(NUMINCX,1);          % Boolean to avoid divide by zero.
CalcX1 = double(zeros(NUMINCX,1));  % Calc variable.
BTempY = zeros(NUMINCY,1);          % Boolean to avoid divide by zero.
CalcY1 = double(zeros(NUMINCY,1));  % Calc variable.
PartX1 = double(zeros(NUMINCX,1));  % Calc variable - see below.
PartX2 = double(zeros(NUMINCX,1));  % Calc variable - see below.
PartX3 = double(zeros(NUMINCX,1));  % Calc variable - see below.
PartY1 = double(zeros(NUMINCY,1));  % Calc variable - see below.
PartY2 = double(zeros(NUMINCY,1));  % Calc variable - see below.
PartY3 = double(zeros(NUMINCY,1));  % Calc variable - see below.
Xmultp = double(0.0);               % Multiplier for x-face calculations.
Ymultp = double(0.0);               % Multipier for y faces.

% Save Old velocities.
UOld = XVel;
VOld = YVel;
% Calculate U.
% u = Part1 - Part2.*Part3
%     Part1 = gX./aX;
%     Part2 = Hux./aX;
%     Part3 = (THETA*G*(dt/DX)).*(EtaXP1 - EtaXM1);
Xmultp = (THETA*G*(fdt/DX));
PartX1 = gX.*aXDen;
PartX2 = Hux.*aXDen;
PartX3 = Xmultp.*(EtaXP1 - EtaXM1);
XVel = PartX1 - PartX2.*PartX3;
% Calculate V.
% v = Part1 - Part2.*Part3
%     Part1 = gY./aY;
%     Part2 = Hvy./aY;
%     Part3 = (THETA*G*(dt/DY)).*(EtaYP1 - EtaYM1);
Ymultp = (THETA*G*(fdt/DY));
PartY1 = gY.*aYDen;
PartY2 = Hvy.*aYDen;
PartY3 = Ymultp.*(EtaYP1 - EtaYM1);
YVel = PartY1 - PartY2.*PartY3;
% Set source values that have specified velocity components.
if (VELDIRCBC == 1);
   [XVel,YVel] = vELdIRCbc(XVel,YVel);
end
if (RADVELBC == 1)
   [XVel,YVel] = rADvELbc(XVel,YVel);
end
% check for precision.
BTempX = (abs(XVel - double(0.0)) > PRECH);
CalcX1 = double(BTempX);
XVel = CalcX1.*XVel;
BTempY = (abs(YVel - double(0.0)) > PRECH);
CalcY1 = double(BTempY);
YVel = CalcY1.*YVel;

clear PartX1 PartX2 PartX3 Xmultp;
clear PartY1 PartY2 PartY3 Ymultp;
clear BTempX CalcX1 BTempY CalcY1;
return;
%EOF