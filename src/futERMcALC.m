function Vel = futERMcALC(Vel,NCol,NRow,pq,pp,NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fvtERMcALC calculates the Coriolis term that will be subtracted from 
% the particle velocity at the previous time step.  For the Coriolis term
% an average v-velocity for the particle location will be computed with
% bilinear interpolation.
% Received:
% 
% Vel [NumInc,1] = Coriolis term. (UVel)
% NCol [NumInc,1] = column index of particle location.
% NRow [NumInc,1] = row index of particle location.
% pq [NumInc,1] = x-direction dimensionless distance to greater boundary.
% pp [NumInc,1] = y-direction dimensionless distance to greater boundary.
% NumInc [1] = number of indexes.
%
% Returned:
%
% Vel [NumInc,1] = Coriolis term. (VVel)

global FCOR fdt u

% local variables.
xIndex1 = zeros(NumInc,1);
xIndex2 = zeros(NumInc,1);
xIndex3 = zeros(NumInc,1);
xIndex4 = zeros(NumInc,1);
uq = double(zeros(NumInc,1));

[xIndex1,xIndex2,xIndex3,xIndex4,uq] = xINDEXfIND(NCol,NRow,pq,NumInc);
% Calculate the average u-velocity at the new locations using bilinear
% bilinear interpolation to use to generate the Coriolis term.  The Coriolis
% term is now stored in UVel.
Vel = (((1 - pp).*((1-uq).*u(xIndex1) + uq.*u(xIndex2))) + ...
   pp.*((1-uq).*u(xIndex4) + uq.*(u(xIndex3))));
Vel = (FCOR*fdt).*Vel;

clear NCol NRow pp pq NumInc xIndex1 xIndex2 xIndex3 xIndex4 uq;
return;
%EOF