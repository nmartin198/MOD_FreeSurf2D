function [ValY] = vOLtOYfACE(MVal,ValY)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vOLtOfACE allocates volume values to face depending on velocity
    % direction.  So upwinds volume center values for faces.
    %
    % Received:
    % 
    % MVal [NUMNODES,1] = value for volume center
    % ValY [NUMINCY,1] = value upwinded to y-faces.
    %
    % Returned;
    %
    % ValY [NUMINCY,1] = value upwinded to y-faces.
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

    global NUMINCY VDirect VYM1 VYP1

    % local variables.
    BTempY = zeros(NUMINCY,1);          % boolean for y-faces.
    CalcY1 = double(zeros(NUMINCY,1));  % Calc variable.
    CalcY2 = double(zeros(NUMINCY,1));  % Calc variable.
    ValYp1 = double(zeros(NUMINCY,1));  % Y-face value i,j+1/2.
    ValYm1 = double(zeros(NUMINCY,1));  % Y-face value i,j-1/2.

    % Distribute to vectors.
    ValYp1 = MVal(VYP1);
    ValYm1 = MVal(VYM1);
    % Now allocate quanties to faces.
    BTempY = (VDirect == 1);
    CalcY1 = double(BTempY);
    CalcY2 = double(~BTempY);
    ValY = CalcY1.*ValYm1 + CalcY2.*ValYp1;
    %  not needed
    %clear BTempY CalcY1 CalcY2 MVal ValYp1 ValYm1;
end
%EOF