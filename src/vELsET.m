function [UVel,VVel] = vELsET(UVel,VVel)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vELsET sets the initial velocity values.  If the files BU.txt exists, then
% UVel will be == BU.txt, else UVel = zeros.  For VVel, the file BV.txt needs
% to exist.  For both BU.txt and BV.txt, the file format should be ASCII column
% vector with one velocity value on each line.  The respective lengths should be
% NUMINCX for BU.txt and NUMINCY for BV.txt.
%
% Received and Returned:
%
% UVel [NUMINCX,1] = x-face velocity.(u)
% VVel [NUMINCY,1] = y-face velocity.(v)

% global variables.
global NUMINCX NUMINCY
% local variables.
Ux = 0;
Vy = 0;

% Try to open text files.
Ux = fopen('BU.txt');
Vy = fopen('BV.txt');

% If file exists, read.  Else close.
if (Ux ~= -1)
   UVel = fscanf(Ux,'%f',NUMINCX);
   fclose(Ux);
end
% If file exists, read.  Else close.
if (Vy ~= -1)
   VVel = fscanf(Vy,'%f',NUMINCY);
   fclose(Vy);
end

clear Ux Vy;
return;
%EOF