function [DepX,DepY] = tdEPdIRCbc(DepX,DepY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% tdEPdIRCbc sets Dirichlet boundaries for total water depth at specified
% boundary volume faces.
%
% DepX [NUMINCX,1] = Hux or the total water depth at volume x-faces
%                       throughout the domain. 
% DepY [NUMINCY,1] = Hvy or the total water depth at volume y-faces
%                       throughout the domain.

global DATUM NUMCOLS TDEPDXDEP TDEPDXVOL TDEPDYDEP TDEPDYVOL
global XINC

% First do x-face sources.
if (TDEPDXVOL(1) ~= 0)
% local variables.
   NumX = 0;
   NumX = size(TDEPDXVOL,2);           % Number of x-face total depth sources.
   BTemp = zeros(NumX,1);              % Boolean calculation variable.
   CalcI1 = zeros(NumX,1);             % Int Calc variable.
   CalcI2 = zeros(NumX,1);             % Int Calc variable.
   Col = zeros(NumX,1);                % Column indexes for sources.
   IndexUH = zeros(NumX,1);            % x-face index for boundary.
   NewNode = zeros(NumX,1);            % Volume indexes
   Row = zeros(NumX,1);                % Row indexes.
   TempCalc1 = double(zeros(NumX,1));  % Calculation variable.
% Begin calculations.
   NewNode = TDEPDXVOL';   
   Row = ceil(NewNode./NUMCOLS);
   TempCalc1 = mod(NewNode,NUMCOLS);
   BTemp = (TempCalc1 == 0);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
   BTemp = (Col > 1);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   IndexUH = (CalcI1.*(((Row - 1).*XINC)+Col+1)) + ...
      (CalcI2.*(((Row - 1).*XINC)+Col));
   DepX(IndexUH,1) = TDEPDXDEP';
   clear NumX BTemp CalcI1 CalcI2 Col IndexUH NewNode Row TempCalc1;
end
% Next do y-face sources.
if (TDEPDYVOL(1) ~= 0)
% local variables.
   NumY = 0;
   NumY = size(TDEPDYVOL,2);        % Number of y-face total depth sources.
   BTemp = zeros(NumY,1);           % Boolean calculation variable.
   CalcI1 = zeros(NumY,1);          % Int Calc variable.
   CalcI2 = zeros(NumY,1);          % Int Calc variable.
   IndexVH = zeros(NumY,1);         % y-face index for boundary.
   NewNode = zeros(NumY,1);         % Volume indexes
   Row = zeros(NumY,1);             % Row indexes.
% Begin calculations.
   NewNode = TDEPDYVOL';   
   Row = ceil(NewNode./NUMCOLS);
   BTemp = (Row > 1);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   IndexVH = (CalcI1.*(NewNode + NUMCOLS)) + (CalcI2.*NewNode);
   DepY(IndexVH,1) = TDEPDYDEP';
   clear NumX BTemp CalcI1 CalcI2 IndexVH NewNode Row;
end

return;
%EOF