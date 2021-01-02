function rhs = rhscALC(rhs,Del)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rhscALC calculates the right hand side (RHS) of the system of equations
    % for free surface elevation.  q is the vector employed to store the terms
    % in the RHS that are not included in Delta.  The notation and numeric
    % representation is very similar to that of Casulli (1999) and Casulli and
    % Cheng (1992).
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

    global aXDen aYDen DX DY fdt gX gY Hux Hvy NUMNODES NUMINCX NUMINCY
    global THETA

    %local variables.
    qX = double(zeros(NUMINCX,1));         % q term for x-faces.
    qY = double(zeros(NUMINCY,1));         % q term for y-faces.
    q1 = double(zeros(NUMNODES,1));        % q term for east faces.
    q2 = double(zeros(NUMNODES,1));        % q term for south faces.
    q3 = double(zeros(NUMNODES,1));        % q term for west faces.
    q4 = double(zeros(NUMNODES,1));        % q term for north faces.
    TempX = double(zeros(NUMINCX,1));      % X-face temp. calc. variable.
    TempXq = double(zeros(NUMINCX,1));     % X-face temp. calc. variable.
    TempYq = double(zeros(NUMINCY,1));     % Y-face temp. calc. variable.
    TempY = double(zeros(NUMINCY,1));      % Y-face temp. calc. variable.
    XMult = double(0.0);                   % Multiplier for x-faces.
    YMult = double(0.0);                   % Multiplier for y-faces.

    % q X-face calculations.
    XMult = THETA*(fdt/DX);
    TempXq = gX.*Hux;
    TempX = TempXq.*aXDen;
    qX = XMult.*TempX;
    % q Y-face cacluations.
    YMult = THETA*(fdt/DY);
    TempYq = gY.*Hvy;
    TempY = TempYq.*aYDen;
    qY = YMult.*TempY;
    % alloc qX and qY to respective faces.
    [q1,q2,q3,q4] = aLLOCfACE(qX,qY,q1,q2,q3,q4);
    % calculate the RHS.
    rhs = Del - q1 + q3 - q4 + q2;
    % not needed
    %clear q1 q2 q3 q4 qX qY;
    %clear Del;
    %clear TempX TempXq TempY TempYq XMult YMult;
end
%EOF