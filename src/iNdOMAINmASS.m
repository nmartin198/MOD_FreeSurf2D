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

    global DX DY EtaNew HNODE NUMNODES RHOW 

    %local variables.
    TempHeight = double(zeros(NUMNODES,1));

    % Calculations.
    TempHeight = HNODE + EtaNew;
    IMass = (ones(1,NUMNODES)*TempHeight)*(DX*DY)*RHOW;
    % no longer needed
    %clear TempHeight;
end
%EOF