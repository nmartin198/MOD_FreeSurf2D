function [UVel,VVel] = qiNfLUXbc(UVel,VVel)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% qiNfLUXbc sets up and enforces Dirichlet inflow flux boundaries.
% Each boundary flux value is set for the inflow face of the boundary volume.
% The same volume cannot be employed for both an x-face and y-face 
% Dirichlet inflow flux boundary.
%
% UVel [NUMINCX,1] = u or x-face depth-averaged velocity.
% VVel [NUMINCY,1] = v or y-face depth-averaged velocity.

global Hux Hvy NUMCOLS PRECH QINXVOL QINXFLUX QINYVOL QINYFLUX

% First do x-face sources.
if (QINXVOL(1) ~= 0)
% Local variables
   NumX = 0;
   NumX = size(QINXVOL,2);             % Number of x-face sources.
   BTemp = zeros(NumX,1);              % Boolean calculation variable.
   BTempC = zeros(NumX,1);             % Boolean calculation variable.
   Calc1 = double(zeros(NumX,1));      % Calc variable.
   Calc2 = double(zeros(NumX,1));      % Calc variable.
   CalcC1 = double(zeros(NumX,1));      % Calc variable.
   CalcC2 = double(zeros(NumX,1));      % Calc variable.
   Col = zeros(NumX,1);                % Col index for x-face sources.
   IndexUH = zeros(NumX,1);            % x-face index for vectors.
   NewNode = zeros(NumX,1);           % Row index for x-face sources.
   TempCalc1 = double(zeros(NumX,1));  % Temporary calculation variable.
   TempH = double(zeros(NumX,1));      % Total water depth values.
% Calculations
%  First get volume, row, column, and x-face vector indexes.
   NewNode = QINXVOL';
   Row = ceil(NewNode./NUMCOLS);
   TempCalc1 = mod(NewNode,NUMCOLS);
   BTempC = (TempCalc1 == 0);
   CalcC1 = double(BTempC);
   CalcC2 = double(~BTempC);
   Col = (CalcC1.*NUMCOLS) + (CalcC2.*TempCalc1);
   BTempC = (Col > 1);
   CalcC1 = double(BTempC);
   CalcC2 = double(~BTempC);
   IndexUH = (CalcC1.*((Row-1)+NewNode+1)) + ...
      (CalcC2.*((Row-1)+NewNode));
%  Now determine source volume face depth and put into denominator form.
   TempH = Hux(IndexUH);
   BTemp = ((TempH - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*TempH) + (Calc2.*PRECH);
   TempH = Calc1.*(1./TempCalc1);
%  Finally, determine inflow velocities given the specified flux.
   TempCalc1 = QINXFLUX'.*TempH;
   BTemp = ((TempCalc1 - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   TempCalc1 = Calc1.*TempCalc1;
   UVel(IndexUH) = (CalcC1.*(-1.*TempCalc1)) + (CalcC2.*TempCalc1);
   clear NumX BTemp BTempC Calc1 Calc2 CalcC1 CalcC2 Col IndexUH NewNode;
   clear Row TempCalc1 TempH;
end
% Then do y-face sources.
if (QINYVOL(1) ~= 0)
% Local variables
   NumY = 0;
   NumY = size(QINYVOL,2);             % Number of y-face sources.
   BTemp = zeros(NumY,1);              % Boolean calculation variable.
   BTempR = zeros(NumY,1);             % Boolean calculation variable.
   Calc1 = double(zeros(NumY,1));      % Calc variable.
   Calc2 = double(zeros(NumY,1));      % Calc variable.
   CalcR1 = double(zeros(NumY,1));      % Calc variable.
   CalcR2 = double(zeros(NumY,1));      % Calc variable.
   IndexVH = zeros(NumY,1);            % y-face index for vectors.
   NewNode = zeros(NumY,1);            % Volume index for y-faces sources.
   Row = zeros(NumY,1);                % Row index for y-face sources.
   TempCalc1 = double(zeros(NumY,1));  % Temporary calculation variable.
   TempH = double(zeros(NumY,1));      % Total water depth values.
% Calculations
%  First get volume, row, and y-face vector indexes.
   NewNode = QINYVOL';
   Row = ceil(NewNode./NUMCOLS);
   BTempR = (Row > 1);
   CalcR1 = double(BTempR);
   CalcR2 = double(~BTempR);
   IndexVH = (CalcR1.*(NewNode + NUMCOLS)) + (CalcR2.*NewNode);
%  Now determine source volume face depth and put into denominator form.
   TempH = Hvy(IndexVH);
   BTemp = ((TempH - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*TempH) + (Calc2.*PRECH);
   TempH = Calc1.*(1./TempCalc1);
%  Finally, determine inflow velocities given the specified flux.
   TempCalc1 = QINYFLUX'.*TempH;
   BTemp = ((TempCalc1 - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   TempCalc1 = Calc1.*TempCalc1;
   VVel(IndexVH) = (CalcR1.*(-1.*TempCalc1)) + (CalcR2.*TempCalc1);
   clear NumY BTemp BTempR Calc1 Calc2 CalcR1 CalcR2 IndexVH NewNode Row;
   clear TempCalc1 TempH;
end

return;
%EOF