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
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copyright and License
    %
    % Copyright 2021 Nick Martin
    %
    % This file is part of MOD_FreeSurf2D.
    %
    % MOD_FreeSurf2d is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % MOD_FreeSurf2D is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU Affero General Public License for more details.
    %
    % You should have received a copy of the GNU Affero General Public License
    % along with MOD_FreeSurf2D.  If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % not needed
    %clear NCol NRow pp pq NumInc yIndex1 yIndex2 yIndex3 yIndex4 vp;
end
%EOF