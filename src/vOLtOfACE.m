function [ValueX,ValueY] = vOLtOfACE(MainVal,ValueX,ValueY)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vOLtOfACE allocates volume values to face depending on velocity
    % direction.  So upwinds volume center values for faces.
    %
    % Received:
    % 
    % MainVal [NUMNODES,1] = value for volume center
    % ValueX [NUMINCX,1] = value upwinded to x-faces.
    % ValueY [NUMINCY,1] = value upwinded to y-faces.
    %
    % Returned;
    %
    % ValueX [NUMINCX,1] = value upwinded to x-faces.
    % ValueY [NUMINCY,1] = value upwinded to y-faces.
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
    % along with  MOD_FreeSurf2D.  If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ValueX] = vOLtOXfACE(MainVal,ValueX);
    [ValueY] = vOLtOYfACE(MainVal,ValueY);
    % not needed
    %clear MainVal;
end
%EOF