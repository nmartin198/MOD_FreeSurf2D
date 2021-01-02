function [VolXP1,VolXM1,VolYP1,VolYM1] = aLLOCaDJvOL(MainVal,...
   VolXP1,VolXM1,VolYP1,VolYM1)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % aLLOCaDJvOL allocates volume indexed vectors to vectors that indexed 
    % so that they correspond to adjacent volumes for velocity definition
    % locations.
    %
    % MainVal is the original volume locations. [NUMNODES,1].
    % VolXP1 is at (i+1,j) for the x-face at (i+1/2,j).
    % VolXM1 is at (i,j) for the x-face at (i+1/2,j).
    % VolYP1 is at (i,j+1) for the x-face at (i,j+1/2).
    % VolYM1 is at (i,j) for the x-face at (i,j+1/2).
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

    global UXM1 UXP1 VYP1 VYM1

    VolXP1 = MainVal(UXP1);
    VolXM1 = MainVal(UXM1);
    VolYP1 = MainVal(VYP1);
    VolYM1 = MainVal(VYM1);
    % not needed anymore
    %clear MainVal;
end
%EOF