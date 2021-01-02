function [FS,DCent,HeightX,HeightY] = sETiNITIALvECTORS(FS,DCent,HeightX,HeightY)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % sETiNITIALvECTORS function sets the undisturbed water depth, HUX and HVY, and
    % sets the free surface elevations for the domain.  At the end of this function,
    % the Dirichlet source values are set for total water depth.
    %
    % FS [NUMNODES,1] = EtaNew or the free surface elevations for each volume in the
    %                       domain.
    % DCent [NUMNODES,1] = HNODE or the total water depth at the center of each volume.
    % HeightX [NUMINCX,1] = HUX or the undisturbed water depth at each x-face.
    % HeightY [NUMINCY,1] = HVY or the undistrubed water depth at each y-face.
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

    global DATUM HDEPTH ZTopo ZtX ZtY

    % Do the undisturbed water depths at the x-faces and y-faces first.
    HeightX = DATUM - ZtX;
    HeightY = DATUM - ZtY;
    % Set the initial free surface.
    FS = (HDEPTH + ZTopo) - DATUM;
    DCent = DATUM - ZTopo;

end
%EOF