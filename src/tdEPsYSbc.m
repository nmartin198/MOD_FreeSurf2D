function [FS1,FS2,FS3,FS4] = tdEPsYSbc(FS1,FS2,FS3,FS4)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% tdEPsYSbc sets Dirichlet boundaries for total water depth.
% The value of total water depth is set by adding the known
% free surface elevation at the boundary to the RHS of the
% free surface system of equations.
%
% FS1 [NUMROWS,1] = free surface value to be added to rhs.
% FS2 [NUMCOLS,1] = free surface value to be added to rhs.
% FS3 [NUMROWS,1] = free surface value to be added to rhs.
% FS4 [NUMCOLS,1] = free surface value to be added to rhs.

global DATUM NUMCOLS TDEPDXDEP TDEPDXVOL TDEPDYDEP TDEPDYVOL
global XINC ZtX ZtY

% First do x-face sources.
if (TDEPDXVOL(1) ~= 0)
% Local variables.
   NumX = 0;
   NumX = size(TDEPDXVOL,2);           % Number of x-face source volumes.
   BTempC = zeros(NumX,1);             % Boolean calculation variable.
   CalcI1 = zeros(NumX,1);             % Int Calc variable.
   CalcI2 = zeros(NumX,1);             % Int Calc variable.
   Col = zeros(NumX,1);                % Column indexes for sources.
   FSX = double(zeros(NumX,1));        % Specified free surface elevation.
   IndexX = zeros(NumX,1);             % Index of x-face boundary.
   NewNode = zeros(NumX,1);            % Volume indexes for sources.
   Row = zeros(NumX,1);                % Row index for sources.
   TempCalc1 = double(zeros(NumX,1));  % Temporary calculation variable.
% Get Indexes.
   NewNode = TDEPDXVOL';
   Row = ceil(NewNode./NUMCOLS);
   TempCalc1 = mod(NewNode,NUMCOLS);
   BTempC = (TempCalc1 == 0);
   CalcI1 = BTempC;
   CalcI2 = (1 - BTempC);
   Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
   BTempC = (Col > 1);
   CalcI1 = BTempC;
   CalcI2 = (1 - BTempC);
   IndexX = (CalcI1.*(((Row - 1).*XINC)+Col+1)) + ...
      (CalcI2.*(((Row - 1).*XINC)+Col));
   IndexF1 = find(BTempC == 1);
   IndexF3 = find(BTempC == 0);
   RowF1 = Row(IndexF1);
   RowF3 = Row(IndexF3);
% Calculate the specified free surface elevation in normalized form.
   FSX = (((TDEPDXDEP') + ZtX(IndexX,1)) - DATUM);
% Determine which face value to multiply the new free surface value by.
   FS1(RowF1) = FSX(IndexF1);
   FS3(RowF3) = FSX(IndexF3);
   clear NumX BTempC CalcI1 CalcI2 Col FSX IndexF1 IndexF3 IndexX NewNode
   clear Row RowF1 RowF3 TempCalc1;
end
% Then y-face sources.
if (TDEPDYVOL(1) ~= 0)
% Local variables.
   NumY = 0;
   NumY = size(TDEPDYVOL,2);           % Number of y-face source volumes.
   BTempR = zeros(NumY,1);             % Boolean calculation variable.
   CalcI1 = zeros(NumY,1);             % Int Calc variable.
   CalcI2 = zeros(NumY,1);             % Int Calc variable.
   FSY = double(zeros(NumY,1));        % Specified free surface elevation.
   IndexY = zeros(NumY,1);             % Index of y-face boundary.
   NewNode = zeros(NumY,1);            % Volume indexes for sources.
   Row = zeros(NumY,1);                % Row index for sources.
   TempCalc1 = double(zeros(NumY,1));  % Temporary calculation variable.
% Get Indexes.
   NewNode = TDEPDYVOL';
   TempCalc1 = mod(NewNode,NUMCOLS);
   BTempR = (TempCalc1 == 0);
   CalcI1 = BTempR;
   CalcI2 = (1 - BTempR);
   Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
   Row = ceil(NewNode./NUMCOLS);
   BTempR = (Row > 1);
   CalcI1 = BTempR;
   CalcI2 = (1 - BTempR);
   IndexY = (CalcI1.*(NewNode + NUMCOLS)) + (CalcI2.*NewNode);
   IndexF4 = find(BTempR == 1);
   IndexF2 = find(BTempR == 0);
   ColF4 = Col(IndexF4);
   ColF2 = Col(IndexF2);
% Calculate the specified free surface elevation in normalized form.
   FSY = (((TDEPDYDEP') + ZtY(IndexY,1)) - DATUM);
% Determine which face value to multiply the new free surface value by.
   FS4(ColF4) = FSY(IndexF4);
   FS2(ColF2) = FSY(IndexF2);
   clear NumY BTempR CalcI1 CalcI2 Col ColF2 ColF4 FSY IndexY IndexF4 IndexF2;
   clear NewNode Row TempCalc1;
end

return;
%EOF