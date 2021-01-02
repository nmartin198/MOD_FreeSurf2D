function [f1,f2,f3,f4] = aLLOCfACE(ValX,ValY,f1,f2,f3,f4)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % aLLOCfACE distributes the domain wide face indexed vectors for x
    % and y faces to four vectors corresponding to particular faces for
    % each volume.  The x-face vector is distributed to east (i+1/2,j) and
    % west (i-1/2,j) faces for each volume.  The y-face vector is 
    % allocated to north (i,j+1/2) and south (i,j-1/2) face vectors for
    % each volume.
    %
    % ValX = x-face values. [NUMINCX,1].
    % ValY = y-face values. [NUMINCY,1].
    % f1 = east face values. [NUMNODES,1].
    % f2 = south face values. [NUMNODES,1].
    % f3 = west face values. [NUMNODES,1].
    % f4 = north face values. [NUMNODES,1].
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

    global X1 X3 Y2 Y4

    % Calculations.
    f1 = ValX(X1);
    f2 = ValY(Y2);
    f3 = ValX(X3);
    f4 = ValY(Y4);
    % no longer needed
    %clear ValX ValY;
end
%EOF