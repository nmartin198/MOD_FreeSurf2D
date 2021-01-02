function [TopoX,TopoY] = sETtOPO(Topo,TopoX,TopoY)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % sETtOPO interpolates topographic elevations from the volume centered
    % topography file to the faces of each volume.  The interpolation 
    % method is simple linear interpolation.
    %
    % Received:
    %
    % Topo [NUMNODES,1] = topographic elevations at volume centers.
    % TopoX [NUMINCX,1] = topographic elevations at x-faces.
    % TopoY [NUMINCY,1] = topographic elevations at y-faces.
    %
    % Returned:
    % TopoX [NUMINCX,1] = topographic elevations at x-faces.
    % TopoY [NUMINCY,1] = topographic elevations at y-faces.
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

    global UXM1 UXP1 VYM1 VYP1

    % calculations.
    TopoX = 0.5.*(Topo(UXM1) + Topo(UXP1));
    TopoY = 0.5.*(Topo(VYM1) + Topo(VYP1));
    % no longer needed
    %clear Topo;
end
%EOF