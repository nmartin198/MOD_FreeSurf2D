function Vel = fvtERMcALC(Vel,NCol,NRow,pq,pp,NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fvtERMcALC calculates the Coriolis term that will be subtracted from 
% the particle velocity at the previous time step.  For the Coriolis term
% bilinear interpolation will be employed to calculate the v - velocity
% at the new location.
%
% Received:
% 
% Vel [NumInc,1] = Coriolis term. (VVel)
% NCol [NumInc,1] = column index of particle location.
% NRow [NumInc,1] = row index of particle location.
% pq [NumInc,1] = x-direction dimensionless distance to greater boundary.
% pp [NumInc,1] = y-direction dimensionless distance to greater boundary.
% NumInc [1] = number of indexes.
%
% Returned:
%
% Vel [NumInc,1] = Coriolis term. (VVel)


global FCOR fdt v

%local variables.
yIndex1 = zeros(NumInc,1);
yIndex2 = zeros(NumInc,1);
yIndex3 = zeros(NumInc,1);
yIndex4 = zeros(NumInc,1);
vp = double(zeros(NumInc,1));

% calculations.
[yIndex1,yIndex2,yIndex3,yIndex4,vp] = yINDEXfIND(NCol,NRow,pp,NumInc);
% Calculate the average v-velocity at the new locations using
% bilinear interpolation to generate the Coriolis term.  The Coriolis
% term is now stored in UVel.
Vel = (((1 - vp).*((1-pq).*v(yIndex1) + pq.*v(yIndex2))) + ...
   vp.*((1-pq).*v(yIndex4) + pq.*(v(yIndex3))));
Vel = (FCOR*fdt).*Vel;

clear NCol NRow pp pq NumInc yIndex1 yIndex2 yIndex3 yIndex4 vp;
return;
%EOF