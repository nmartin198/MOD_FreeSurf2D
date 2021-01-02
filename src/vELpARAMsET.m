function [AveU,DirectU,AveV,DirectV] = vELpARAMsET(AveU,DirectU,AveV,DirectV,UXVel,VYVel)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vELpARAMsET sets the average component velocities for each volume and determines the
    % direction (positive or negative) for each volume face velocity.
    %
    % Received:
    %
    %  AveU = UAve [NUMNODES,1] volume average x-direction velocity component.
    %  AveV = VAve [NUMNODES,1] volume average y-direction velocity component.
    %  DirectU = UDirect [NUMINCX,1] x-face velocity direction (+,-).
    %  DirectV = VDirect [NUMINCX,1] y-face velocity direction (+,-).
    %  XVel = u [NUMINCX,1] x-face velocities.
    %  YVel = v [NUMINCY,1] y-face velocities.
    %
    % Returned:
    %
    %  AveU = UAve [NUMNODES,1] volume average x-direction velocity component.
    %  AveV = VAve [NUMNODES,1] volume average y-direction velocity component.
    %  DirectU = UDirect [NUMINCX,1] x-face velocity direction (+,-).
    %  DirectV = VDirect [NUMINCX,1] y-face velocity direction (+,-).
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

    global NUMINCX NUMINCY NUMNODES PREC

    % local variables.
    BTempX = zeros(NUMINCX,1);          % Boolean for if positive.
    BTempY = zeros(NUMINCY,1);          % Boolean for if positive.
    CalcXI1 = double(zeros(NUMINCX,1));  % Calc variable.
    CalcXI2 = double(zeros(NUMINCX,1));  % Calc variable.
    CalcYI1 = double(zeros(NUMINCY,1));  % Calc variable.
    CalcYI2 = double(zeros(NUMINCY,1));  % Calc variable.
    u1 = double(zeros(NUMNODES,1));     % velocity at (i+1/2,j) for each volume.
    v2 = double(zeros(NUMNODES,1));     % velocity at (i,j-1/2) for each volume.
    u3 = double(zeros(NUMNODES,1));     % velocity at (i-1/2,j) for each volume.
    v4 = double(zeros(NUMNODES,1));     % velocity at (i,j+1/2) for each volume.

    % get face values.
    [u1,v2,u3,v4] = aLLOCfACE(UXVel,VYVel,u1,v2,u3,v4);
    % Get averages.
    AveU = 0.5.*(u1 + u3);
    AveV = 0.5.*(v2 + v4);
    % Determine direction.
    BTempX = ((UXVel - double(0.0)) > -PREC);
    CalcXI1 = (BTempX);
    CalcXI2 = (1 - BTempX);
    BTempY = ((VYVel - double(0.0)) > -PREC);
    CalcYI1 = (BTempY);
    CalcYI2 = (1 - BTempY);
    DirectU = (CalcXI1.*1.0) + (CalcXI2.*-1.0);
    DirectV = (CalcYI1.*1.0) + (CalcYI2.*-1.0);
    % not needed
    %clear BTempX BTempY CalcXI1 CalcXI2 CalcYI1 CalcYI2 UXVel VYVel;
    %clear u1 v2 u3 v4;
end
%EOF