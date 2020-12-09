function fcORcALC
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fcORcALC calculates the Coriolis f-parameter using an f-plane Coriolis model.
% The f-parameter enters into the calculations in the inertial term --- see the
% Fv or Fu term from Casulli and Cheng (1992).

global FCOR COROMEGA CENTRALLATITUDE

FCOR = 2*COROMEGA*(sin((pi/180)*CENTRALLATITUDE));

return;
%EOF