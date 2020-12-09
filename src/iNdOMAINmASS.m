function IMass = iNdOMAINmASS(IMass)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% iNdOMAINmASS calculates the mass of water in the simulation domain.
% Mass of water is calculated by taking the free surface elevation, 
% calculating the representative water depth, and deterimining the mass
% in the cell from the this height.
%
% Received and Returned:
%
% IMass [1] = mass of water in kg in simulation domain. (BMass).

global DX DY EtaNew HNODE NUMNODES RHOW 

%local variables.
TempHeight = double(zeros(NUMNODES,1));

% Calculations.
TempHeight = HNODE + EtaNew;
IMass = (ones(1,NUMNODES)*TempHeight)*(DX*DY)*RHOW;

clear TempHeight;
return;
%EOF