function [XBnd,YBnd,UMaj,UMin,VMaj,VMin,GDx,GDy,LDx,LDy] = vOLUMEfIND(Xp,Yp,UVel,...
   VVel,VelX,VelY,pp,pq,NRow,NCol,NumInc)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vOLUMEfIND determines the current volume, that a particle will backwards
    % across as part of the semi-analytical pathline tracing method.  This
    % function is needed because this method of pathline tracing is intended
    % to end up at volume boundaries where extra decision information is needed.
    %
    % Xp [NumInc,1] is the particle current x-coordinate.  == XLoc or Xe.
    % Yp [NumInc,1] is the particle current y-coordinate.  == YLoc or Ye.
    % UVel [NumInc,1] is bilinearly interpolated u-velocity at new location.
    %     == Vxp
    % VVel [NumInc,1] is bilinearly interpolated v-velocity at new location.
    %     == Vyp
    % VelX [NUMINCX,1] is x-face velocity.
    % VelY [NUMINCY,1] is y-face velocity.
    % pp [NumInc,1] is the normalized distance in x-direction to greatest magnitude
    %     x-face. == p
    % pq [NumInc,1] is the normalized distance in y-direction to greatest magnitude
    %     y-face. == q
    % NRow [NumInc,1] is the row index of the current location. == PRow
    % NCol [NumInc,1] is the column index of the current location.  == PCol.
    % NumInc is the number of indexes. == NUMINC
    %
    %
    % _|____|VMin|____|_          X is the particle location.  Velocity
    %  |    |    X    |              directions are -> and \|/
    %  |  UMin   UMaj |                                         
    % _|____|____|____|_ Y2
    %  |    |VMaj|    |
    %           X2
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

    global DX DY PRECH NUMCOLS PREC XINC XINDEX YINDEX
    % local variables.
    BTempU = zeros(NumInc,1);           % Boolean variable for x-quantites.
    BTempV = zeros(NumInc,1);           % Boolean variable for y-quantites.
    BTempX = zeros(NumInc,1);           % Boolean variable for x-quantites.
    BTempY = zeros(NumInc,1);           % Boolean variable for y-quantites.
    CalcU1 = double(zeros(NumInc,1));   % Calc variable.
    CalcU2 = double(zeros(NumInc,1));   % Calc variable.
    CalcV1 = double(zeros(NumInc,1));   % Calc variable.
    CalcV2 = double(zeros(NumInc,1));   % Calc variable.
    CalcX1 = double(zeros(NumInc,1));   % Calc variable.
    CalcX2 = double(zeros(NumInc,1));   % Calc variable.
    CalcY1 = double(zeros(NumInc,1));   % Calc variable.
    CalcY2 = double(zeros(NumInc,1));   % Calc variable.
    CalcIU1 = zeros(NumInc,1);          % Int Calc variable.
    CalcIU2 = zeros(NumInc,1);          % Int Calc variable.
    CalcIV1 = zeros(NumInc,1);          % Int Calc variable.
    CalcIV2 = zeros(NumInc,1);          % Int Calc variable.
    CalcIX1 = zeros(NumInc,1);          % Int Calc variable.
    CalcIX2 = zeros(NumInc,1);          % Int Calc variable.
    CalcIY1 = zeros(NumInc,1);          % Int Calc variable.
    CalcIY2 = zeros(NumInc,1);          % Int Calc variable.
    GDx = double(zeros(NumInc,1));      % Gradient delta X for Vxp calculations.
    GDy = double(zeros(NumInc,1));      % Gradient delta Y for Vyp calculations.
    LDx = double(zeros(NumInc,1));      % Distance delta X for time of flight calcs.
    LDy = double(zeros(NumInc,1));      % Distance delta Y for time fo flight calcs.
    NewCol = zeros(NumInc,1);           % New column location.
    NewRow = zeros(NumInc,1);           % New row location.
    NewNode = zeros(NumInc,1);          % New node index.
    NewUBase = zeros(NumInc,1);         % base for u-indexes.
    NPrec = double(0.0);                % Precision for normalized variables.
    UMaj = double(zeros(NumInc,1));     % U-value at downwind definition location.
    UMin = double(zeros(NumInc,1));     % U-value for upwind definition location.
    VMaj = double(zeros(NumInc,1));     % V-value at downwind definition location.
    VMin = double(zeros(NumInc,1));     % V-value at upwind definition location.
    XBnd = double(zeros(NumInc,1));     % X-coordinate of downwind boundary.
    YBnd = double(zeros(NumInc,1));     % Y-coordinate of downwing boundary.

    % calculations.
    NPrec = PREC^(1/5);
    % Set-up initial boolean quantities.
    BTempU = ((UVel - double(0.0)) > -PREC);
    CalcIU1 = BTempU;
    CalcIU2 = (1 - BTempU);
    BTempV = ((VVel - double(0.0)) > -PREC);
    CalcIV1 = BTempV;
    CalcIV2 = (1 - BTempV);
    % First determine the column index of the new volume.
    BTempX = (abs(pp - 1.0) <= NPrec);
    CalcIX1 = BTempX;
    CalcIX2 = (1 - BTempX);
    NewCol = (CalcIX1.*((CalcIU1.*(NCol - 1)) + (CalcIU2.*NCol))) + ...
       (CalcIX2.*NCol);
    BTempX = (abs(pp - 0.0) <= NPrec);
    CalcIX1 = BTempX;
    CalcIX2 = (1 - BTempX);
    NewCol = (CalcIX1.*((CalcIU1.*NCol) + (CalcIU2.*(NCol + 1)))) + ...
       (CalcIX2.*NewCol);
    NewCol = cOLaDJ(NewCol,NumInc);
    % Then determine the row index of the new volume.
    BTempY = (abs(pq - 1.0) <= NPrec);
    CalcIY1 = BTempY;
    CalcIY2 = (1 - BTempY);
    NewRow = (CalcIY1.*((CalcIV1.*(NRow - 1)) + (CalcIV2.*NRow))) + ...
       (CalcIY2.*NRow);
    BTempY = (abs(pq - 0.0) <= NPrec);
    CalcIY1 = BTempY;
    CalcIY2 = (1 - BTempY);
    NewRow = (CalcIY1.*((CalcIV1.*NRow) + (CalcIV2.*(NRow + 1)))) + ...
       (CalcIY2.*NewRow);
    NewRow = rOWaDJ(NewRow,NumInc);
    % Now get the volume index of the new volume.
    NewNode = ((NewRow - 1).*NUMCOLS) + NewCol;
    NewUBase = ((NewRow - 1).*XINC) + NewCol;
    % Now that have the location - generate the quantities of interest.
    CalcU1 = double(CalcIU1);
    CalcU2 = double(CalcIU2);
    CalcV1 = double(CalcIV1);
    CalcV2 = double(CalcIV2);
    XBnd = (CalcU1.*(XINDEX(NewCol + 1))') + (CalcU2.*(XINDEX(NewCol))');
    YBnd = (CalcV1.*(YINDEX(NewRow + 1))') + (CalcV2.*(YINDEX(NewRow))');
    UMaj = (CalcU1.*VelX(NewUBase + 1)) + (CalcU2.*VelX(NewUBase));
    UMin = (CalcU1.*VelX(NewUBase)) + (CalcU2.*VelX(NewUBase + 1));
    VMaj = (CalcV1.*VelY(NewNode + NUMCOLS)) + (CalcV2.*VelY(NewNode));
    VMin = (CalcV1.*VelY(NewNode)) + (CalcV2.*VelY(NewNode + NUMCOLS));
    GDx = XBnd - Xp;
    BTempX = ((abs(GDx) - 0.0) > PRECH);
    CalcX1 = double(BTempX);
    CalcX2 = double(~BTempX);
    GDx = (CalcX1.*GDx) + (CalcX2.*0.0);
    BTempX = ((DX - abs(GDx)) > PRECH);
    CalcX1 = double(BTempX);
    CalcX2 = double(~BTempX);
    GDx = (CalcX1.*GDx) + (CalcX2.*DX); 
    GDy = YBnd - Yp;
    BTempY = ((abs(GDy) - 0.0) > PRECH);
    CalcY1 = double(BTempY);
    CalcY2 = double(~BTempY);
    GDy = (CalcY1.*GDy) + (CalcY2.*0.0);
    BTempY = ((DY - abs(GDy)) > PRECH);
    CalcY1 = double(BTempY);
    CalcY2 = double(~BTempY);
    GDy = (CalcY1.*GDy) + (CalcY2.*DY); 
    LDx = DX - abs(GDx);
    LDy = DY - abs(GDy);
    BTempX = ((LDx - 0.0) > PRECH);
    CalcX1 = double(BTempX);
    CalcX2 = double(~BTempX);
    LDx = (CalcX1.*LDx) + (CalcX2.*0.0);
    BTempX = ((DX - LDx) > PRECH);
    CalcX1 = double(BTempX);
    CalcX2 = double(~BTempX);
    LDx = (CalcX1.*LDx) + (CalcX2.*DX); 
    BTempY = ((LDy - 0.0) > PRECH);
    CalcY1 = double(BTempY);
    CalcY2 = double(~BTempY);
    LDy = (CalcY1.*LDy) + (CalcY2.*0.0);
    BTempY = ((DY - LDy) > PRECH);
    CalcY1 = double(BTempY);
    CalcY2 = double(~BTempY);
    LDy = (CalcY1.*LDy) + (CalcY2.*DY); 
    % not needed
    %clear BTempU BTempV BTempX BTempY CalcU1 CalcU2 CalcV1 CalcV2;
    %clear CalcX1 CalcX2 CalcY1 CalcY2 CalcIU1 CalcIU2 CalcIV1 CalcIV2;
    %clear CalcIX1 CalcIX2 CalcIY1 CalcIY2 NewCol NewRow NewNode NewUBase;
end
%EOF