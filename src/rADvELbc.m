function [UVel,VVel] = rADvELbc(UVel,VVel)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% rADvELbc enforces radiation velocity boundaries.
% Each boundary value is set for the exit face of the boundary volume.
% The drift velocity employed is obtaine by simple velocity
% upwinding (Yu-Heng pers.comm.).
%
% UVel [NUMINCX,1] = XVel = u or updated x-face velocities.
% VVel [NUMINCY,1] = YVel = v or updated y-face velociites.

global DX DY fdt PREC RVELXVOL RVELYVOL
global VelIndexBX VelIndexBX1 VelROldXB VelROldXB1
global VelRXMult VelSizeX
global VelIndexBY VelIndexBY1 VelROldYB VelROldYB1
global VelRYMult VelSizeY

% First do x-face sources.
if (RVELXVOL(1) ~= 0)
% local variables.
   BTemp = zeros(VelSizeX,1);             % Boolean calc. variable.
   Calc1 = double(zeros(VelSizeX,1));     % Calc variable.
   Calc2 = double(zeros(VelSizeX,1));     % Calc variable.
   Ce = double(zeros(VelSizeX,1));        % Celerity term for radiation.
   VelCeX = double(zeros(VelSizeX,1));    % Upwinded velocity for celerity term.
% Calculations.
   VelCeX = abs(UVel(VelIndexBX1));
   Ce = VelCeX.*(fdt/DX);
   BTemp = ((Ce - 0.1) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   Ce = (Calc1.*Ce) + (Calc2.*1.0);
   BTemp = ((1.0 - Ce) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   Ce = (Calc1.*Ce) + (Calc2.*1.0);
   UVel(VelIndexBX) = VelRXMult.*((Ce.*VelROldXB1) + ...
      ((1.0 - Ce).*VelROldXB));
   VelROldXB = abs(UVel(VelIndexBX));
   VelROldXB1 = VelCeX;
   clear BTemp Calc1 Calc2 Ce VelCeX;
end
% Then do y-face sources.
if (RVELYVOL(1) ~= 0)
% local variables.
   BTemp = zeros(VelSizeY,1);             % Boolean calc. variable.
   Calc1 = double(zeros(VelSizeY,1));     % Calc variable.
   Calc2 = double(zeros(VelSizeY,1));     % Calc variable.
   Ce = double(zeros(VelSizeY,1));        % Celerity term.
   VelCeY = double(zeros(VelSizeY,1));    % Upwinded velocity for celerity term.
% Calculations.
   VelCeY = abs(VVel(VelIndexBY1));
   Ce = VelCeY.*(fdt/DY);
   BTemp = ((Ce - 0.1) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   Ce = (Calc1.*Ce) + (Calc2.*1.0);
   BTemp = ((1.0 - Ce) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   Ce = (Calc1.*Ce) + (Calc2.*1.0);
   VVel(VelIndexBY) = VelRYMult.*((Ce.*VelROldYB1) + ...
      ((1.0 - Ce).*VelROldYB));
   VelROldYB = abs(VVel(VelIndexBY));
   VelROldYB1 = VelCeY;
   clear BTemp Calc1 Calc2 Ce VelCeY;
end

return;
%EOF