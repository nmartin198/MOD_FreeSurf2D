function rADvELsET
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% rADvELsET initializes the index vectors values for the radiation
% velocity boundary conditions.

global NUMCOLS RVELXVOL RVELYVOL u v
global VelIndexBX VelIndexBX1 VelROldXB VelROldXB1 
global VelRXMult VelSizeX
global VelIndexBY VelIndexBY1 VelROldYB VelROldYB1
global VelRYMult VelSizeY 

% First do x-faces.
if (RVELXVOL(1) ~= 0)
% Global variables.
   VelSizeX = 0;                             % Number of x-face source volumes.
   VelSizeX = size(RVELXVOL,2);
   VelIndexBX = zeros(VelSizeX,1);           % X-face boundary index.
   VelIndexBX1 = zeros(VelSizeX,1);          % X-face boundary - 1 index.
   VelROldXB = double(zeros(VelSizeX,1));    % Velocity value at boundary at N.
   VelROldXB1 = double(zeros(VelSizeX,1));   % Velocity value at boundary -1 at N.
   VelRXMult = zeros(VelSizeX,1);            % Outflow velocity direction.
% local variables.
   BTemp = zeros(VelSizeX,1);             % Boolean calc. variable.
   CalcI1 = zeros(VelSizeX,1);            % Int. calc. variable.
   CalcI2 = zeros(VelSizeX,1);            % Int. calc. variable.
   Col = zeros(VelSizeX,1);               % Column index.
   NewNode = zeros(VelSizeX,1);           % Volume index.
   Row = zeros(VelSizeX,1);               % Row index.
   TempCalc1 = double(zeros(VelSizeX,1)); % Calculation variable.
% Calculations
   NewNode = RVELXVOL';
   Row = ceil(NewNode./NUMCOLS);
   TempCalc1 = mod(NewNode,NUMCOLS);
   BTemp = (TempCalc1 == 0);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
   BTemp = (Col > 1);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   VelRXMult = CalcI1 + (CalcI2.*-1);
   VelIndexBX = (CalcI1.*((Row - 1)+NewNode+1)) + ...
      (CalcI2.*((Row - 1)+NewNode));
   VelIndexBX1 = (CalcI1.*((Row - 1)+NewNode)) + ...
      (CalcI2.*((Row - 1)+NewNode+1));
   VelROldXB = abs(u(VelIndexBX));
   VelROldXB1 = abs(u(VelIndexBX1));
   clear BTemp CalcI1 CalcI2 Col NewNode Row TempCalc1;
else
   VelSizeX = 0;
   VelIndexBX = 0;
   VelIndexBX1 = 0;
   VelROldXB = double(0.0);
   VelROldXB1 = double(0.0);
   VelRXMult = 0;   
end
% Next do y-faces.
if (RVELYVOL(1) ~= 0)
% Global variables.
   VelSizeY = 0;                             % Number of y-face source boundaries.
   VelSizeY = size(RVELYVOL,2);
   VelIndexBY = zeros(VelSizeY,1);           % Index for y-face at boundary.
   VelIndexBY1 = zeros(VelSizeY,1);          % Index for y-face at boundary-1.
   VelROldYB = double(zeros(VelSizeY,1));    % Velocity value at boundary at time N.
   VelROldYB1 = double(zeros(VelSizeY,1));   % Velocity value at boundary-1 at time N.
   VelRYMult = zeros(VelSizeY,1);            % Direction for outflow.
% local variables.
   BTemp = zeros(VelSizeY,1);                % Boolean calculation variable.
   CalcI1 = zeros(VelSizeY,1);            % Int. calc. variable.
   CalcI2 = zeros(VelSizeY,1);            % Int. calc. variable.
   NewNode = zeros(VelSizeY,1);              % Volume index.
   Row = zeros(VelSizeY,1);                  % Row index.
% Calculations
   NewNode = RVELYVOL';
   Row = ceil(NewNode./NUMCOLS);
   BTemp = (Row > 1);
   CalcI1 = BTemp;
   CalcI2 = (1 - BTemp);
   VelRYMult = CalcI1 + (CalcI2.*-1);
   VelIndexBY = (CalcI1.*(NewNode+NUMCOLS)) + (CalcI2.*NewNode);
   VelIndexBY1 = (CalcI1.*NewNode) + (CalcI2.*(NewNode+NUMCOLS));
   VelROldYB = abs(v(VelIndexBY));
   VelROldYB1 = abs(v(VelIndexBY1));
   clear BTemp CalcI1 CalcI2 NewNode Row TempCalc1;
else
   VelSizeY = 0;
   VelIndexBY = 0;
   VelIndexBY1 = 0;
   VelROldYB = double(0.0);
   VelROldYB1 = double(0.0);
   VelRYMult = 0;   
end

return;
%EOF