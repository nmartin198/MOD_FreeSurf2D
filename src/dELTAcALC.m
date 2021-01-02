function Del = dELTAcALC(Del)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % dELTAcALC calculates the part of the RHS of the system of equations for
    % free surface elevations that corresponds to the delta(i,j) term from 
    % Casulli and Cheng (1992) and Casulli (1999).
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

    global DX DY Eta fdt Hux Hvy NUMINCX NUMINCY NUMNODES THETA
    global u v

    % local variables.
    DeltaX = double(zeros(NUMINCX,1));% Delta_{i+-1/2,j} term.
    DeltaY = double(zeros(NUMINCY,1));% Delta_{i,j+-1/2} term.
    TempX = double(zeros(NUMINCX,1)); % Temporary calculation variable for x-faces.
    TempY = double(zeros(NUMINCY,1)); % Temporary calculation variable for y-faces.
    Delta1 = double(zeros(NUMNODES,1));% Delta_{i+1/2,j}
    Delta2 = double(zeros(NUMNODES,1));% Delta_{i,j-1/2}
    Delta3 = double(zeros(NUMNODES,1));% Delta_{i-1/2,j}
    Delta4 = double(zeros(NUMNODES,1));% Delta_{i,j+1/2}
    XMult = double(0.0);              % Factor to multiply all x-faces by.
    YMult = double(0.0);              % Factor to multiply all y-faces by.

    % calculations.
    XMult = (1 - THETA)*(fdt/DX);
    YMult = (1 - THETA)*(fdt/DY);
    % X -face calculations.
    % DeltaX = ((1 - THETA)*(dt/DX)).*Hux.*u;
    TempX = Hux.*u;
    DeltaX = XMult.*TempX;
    % Y -face calculations.
    % DeltaY = ((1 - THETA)*(dt/DY)).*Hvy.*v;
    TempY = Hvy.*v;
    DeltaY = YMult.*TempY;
    % allocate DeltaX and DeltaY to the faces of each node.
    [Delta1,Delta2,Delta3,Delta4] = aLLOCfACE(DeltaX,DeltaY,Delta1,...
        Delta2,Delta3,Delta4);
    % calculate full delta values.
    % Delta = Eta - [Delta1 - Delta3] - [Delta4 - Delta2];
    Del = Eta - Delta1 + Delta3 - Delta4 + Delta2;
    % not needed
    %clear Delta1 Delta2 Delta3 Delta4 DeltaX DeltaY XMult YMult;
    %clear TempX TempY;
end
%EOF