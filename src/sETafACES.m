function [AvX,AvY] = sETafACES(AvX,AvY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETafACES calculates the values for the A* vectors.  Each value
% represents the value of A from Casulli and Cheng (1992) for the
% respective face for every node in the domain.
%
% Received and Returned:
%
%  AvX = aXDen [NUMINCX,1] the denominator form of the A vector.
%  AvY = aYDen [NUMINCY,1] the denominator form of the A vector.

global BHux BHvy CzX CzY fdt G GAMMATX GAMMATY Hux Hvy NUMINCX
global NUMINCY PREC PRECH u UAve UDirect UXM1 UXP1 v VAve VDirect VYP1
global VYM1

% local variables.
BTempX = zeros(NUMINCX,1);          % Boolean.
BTempY = zeros(NUMINCY,1);          % Boolean.
CalcX1 = double(zeros(NUMINCX,1));  % Calc variable.
CalcX2 = double(zeros(NUMINCX,1));  % Calc variable.
CalcY1 = double(zeros(NUMINCY,1));  % Calc variable.
CalcY2 = double(zeros(NUMINCY,1));  % Calc variable.
Part2X = double(zeros(NUMINCX,1));  % Variable to store half of the A term.
Part3X = double(0.0);               % Variable to store part of the A term.
Part2Y = double(zeros(NUMINCY,1));  % Variable to store half of the A term.
Part3Y = double(0.0);               % Variable to store part of A term.
TempX = double(zeros(NUMINCX,1));   % Temporary calculation variable.
TempY = double(zeros(NUMINCY,1));   % Temporary calculation variable.
TCzX = double(zeros(NUMINCX,1));    % Temporary calculation variable to hold CzX.
TCzY = double(zeros(NUMINCY,1));    % Temporary calculation variable to hold CzX.
UforY = double(zeros(NUMINCY,1));   % Volume average u-velocity upwinded for each v-term.
VforX = double(zeros(NUMINCX,1));   % Volume average v-velocity upwinded for each u-term.

% calculations.
% aX
% aX = Hux + Part2 + Part3;
%     Part2 = (G*dt).*((sqrt(u.^2 + AveV.^2))./(Cz.^2));
%     Part3 = DT*GAMMATX
% First get Chezy coefficient in non-divide by zero format.
BTempX = ((CzX - double(0.0)) > PREC);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
TempX = (CalcX1.*(CzX.^2)) + (CalcX2.*PREC);
TCzX = CalcX1.*(1./TempX);
% Now get average v-velocity for each u location.
BTempX = (UDirect == 1);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
VforX = (CalcX1.*VAve(UXM1)) + (CalcX2.*VAve(UXP1));
% Now calculate velocity component of A term and store it in Temp.
TempX = sqrt(u.^2 + VforX.^2);
% Now calcualte the non-depth half of the A term.
Part2X = (fdt*G).*(TCzX.*TempX);
% Now calculate the wind contribution.
Part3X = fdt*GAMMATX;
% Finally calculate aX.
AvX = Hux + Part2X + Part3X;
% Now generate the denominator form.  A* should always be positive.
BTempX = ((AvX - double(0.0)) > PRECH);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
TempX = (CalcX1.*AvX) + (CalcX2.*PRECH);
AvX = 1./TempX;
% Make sure that only wet x-faces have AvX values.
AvX = BHux.*AvX;
% aY
% aY = Hvy + Part2;
%     Part2 = (G*dt).*((sqrt(AveU.^2 + v.^2))./(Cz.^2));
% calculations.
% First get Chezy coefficient in non-divide by zero format.
BTempY = ((CzY - double(0.0)) > PREC);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
TempY = (CalcY1.*(CzY.^2)) + (CalcY2.*PREC);
TCzY = CalcY1.*(1./TempY);
% Now get average v-velocity for each u location.
BTempY = (VDirect == 1);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
UforY = (CalcY1.*UAve(VYM1)) + (CalcY2.*UAve(VYP1));
% Now calculate velocity component of A term and store it in Temp.
TempY = sqrt(v.^2 + UforY.^2);
% Now calcualte the non-depth half of the A term.
Part2Y = (fdt*G).*(TCzY.*TempY);
% Now calculate the wind contribution.
Part3Y = fdt*GAMMATY;
% Finally calculate aY.
AvY = Hvy + Part2Y + Part3Y;
% Now generate the denominator form.  A* should always be positive.
BTempY = ((AvY - double(0.0)) > PRECH);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
TempY = (CalcY1.*AvY) + (CalcY2.*PRECH);
AvY = 1./TempY;
% Make sure that only wet y-faces have AvX values.
AvY = BHvy.*AvY;

clear BTempX CalcX1 CalcX2 Part2X Part3X TempX TCzX VforX;
clear BTempY CalcY1 CalcY2 Part2Y Part3Y TempY TCzY UforY;
return;
%EOF