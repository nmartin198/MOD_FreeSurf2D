function [XLoc,YLoc,XPos,YPos,PCol,PRow,PNode] = eULiNT(XLoc,YLoc,tUvel,tVvel,UMain,VMain,...
   DTMain,MaxN,MinN,NumInc)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % eULiNT uses Eulerian integration to follow a particle trajectory backwards.
    %
    % Received:
    %
    % XLoc [NumInc,1] = x-coordinate location or Xe.
    % YLoc [NumInc,1] = y-coordinate location or Ye.
    % tUvel [NumInc,1] = interpolated u velocity.
    % tVvel [NumInc,1] = interpolated v velocity.
    % UMain [NUMINCX,1] = u velocities in simulation domain.
    % VMain [NUMINCY,1] = v velocities in the simulation domain.
    % DTMain = time duration for path line tracing.
    % MaxN = max number of partial steps backwards.
    % MinN = minimum number of partial steps backwards.
    % NumInc = number of location indexes.
    %
    % Returned:
    %
    % XLoc [NumInc,1] = departure point x-location.
    % YLoc [NumInc,1] = departure point y-location.
    % XPos [NumInc,1] = x-coordinate of max. x boundary.
    % YPos [NumInc,1] = y-coordinate of max. y boundary.
    % PCol [NumInc,1] = departure point column index.
    % PRow [NumInc,1] = departure point row index.
    % PNode [NumInc,1] = departure point node index.
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

    global DX DY NUMCOLS NUMROWS PREC PRECH XINC XINDEX YINDEX

    % local variables.
    ax = double(zeros(NumInc,1)); % Weight for distance from defined location.
    bx = double(zeros(NumInc,1)); % Weight for distance from defined location.
    ay = double(zeros(NumInc,1)); % Weight for distance from defined location.
    by = double(zeros(NumInc,1)); % Weight for distance from defined location.
    Bound = zeros(NumInc,1);      % Integer to determine if on a boundary.
    BoundX = zeros(NumInc,1);     % Boolean to determine if on a boundary.
    Bound1 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    Bound3 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    BoundY = zeros(NumInc,1);     % Boolean to determine if on a boundary.
    Bound2 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    Bound4 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    BTemp = zeros(NumInc,1);      % Boolean calc variable.
    CalcB1 = double(zeros(NumInc,1)); % Calc variable.
    CalcB2 = double(zeros(NumInc,1)); % Calc variable.
    DeltaT = double(0.0);         % Partial time step.
    maxTau = double(0.0);         % Maximum calculated partial time step.
    NewCol = zeros(NumInc,1);     % Column location tracking variables.
    NewCol1 = zeros(NumInc,1);    % Column location tracking variable.
    NewRow = zeros(NumInc,1);     % Row location tracking variable.
    NewRow1 = zeros(NumInc,1);    % Row location tracking variable.
    Number = zeros(NumInc,1);     % Number of partial time steps.
    XPos = double(zeros(NumInc,1));% X direction volume boundary locations.
    YPos = double(zeros(NumInc,1));% Y direction volume boudnary locations.
    PCol = zeros(NumInc,1);       % Column locations.
    PRow = zeros(NumInc,1);       % Row locations.
    PNode = zeros(NumInc,1);      % Node Locations.
    uIndex1 = zeros(NumInc,1);    % Index of u velocity location 1.
    uIndex2 = zeros(NumInc,1);    % Index of u velocity location 2.
    uIndex3 = zeros(NumInc,1);    % Index of u velocity location 3.
    uIndex4 = zeros(NumInc,1);    % Index of u velocity location 4.
    vIndex1 = zeros(NumInc,1);    % Index of v velocity location 1.
    vIndex2 = zeros(NumInc,1);    % Index of v velocity location 2.
    vIndex3 = zeros(NumInc,1);    % Index of v velocity location 3.
    vIndex4 = zeros(NumInc,1);    % Index of v velocity location 4.
    XPart = double(0.0);          % Min time to cross volume in x-direction.
    YPart = double(0.0);          % Min time to cross volume in y-direction.

    % Calculations.
    % Determine the partial time step DeltaT to employ in the integration.
    XPart = max(abs(UMain));
    BTemp = ((XPart - double(0.0)) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    XPart = (CalcB1.*XPart) + (CalcB2.*PREC);
    XPart = DX/XPart;
    YPart = max(abs(VMain));
    BTemp = ((YPart - double(0.0)) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    YPart = (CalcB1.*YPart) + (CalcB2.*PREC);
    YPart = (0.5*DY)/YPart;
    maxTau = min(XPart,YPart);
    maxTau = max(maxTau,(DTMain/MaxN));
    % Now set dTau so that divides into dt evenly.
    Number = ceil(DTMain/maxTau);
    Number = max(Number,MinN);
    DeltaT = DTMain/Number;
    % Calculations.
    % Adjust for boundary locations.
    Bound4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
       ((tVvel - double(0.0)) <= PRECH));
    Bound2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
       ((tVvel - double(0.0)) >= -PRECH));
    BoundY = Bound4 + Bound2;
    Bound1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
       ((tUvel - double(0.0)) <= PRECH));
    Bound3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
       ((tUvel - double(0.0)) >= -PRECH));
    BoundX = Bound1 + Bound3;
    Bound = BoundX + BoundY;
    BTemp = (Bound <= 0);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    % Now call calculate first location.
    XLoc = (CalcB1.*(XLoc - (DeltaT.*tUvel))) + (CalcB2.*XLoc);
    YLoc = (CalcB1.*(YLoc - (DeltaT.*tVvel))) + (CalcB2.*YLoc);
    [XLoc,YLoc] = iNTbOUND(XLoc,YLoc,NumInc);
    % Call routine to trace back the particle path through the complete time
    %     step.
    % First calculate new location.
    [PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NumInc,...
        PCol,PRow,XPos,YPos,PNode);
    % Now time loop.
    for i=2:Number
       % First calculate new u velocity.
       pw = (XPos - XLoc)./DX;
       qw = (YPos - YLoc)./DY;
       BTemp = ((qw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewRow = (CalcB1.*PRow) + (CalcB2.*(PRow+1));
       NewRow = rOWaDJ(NewRow,NumInc);
       NewRow1 = (CalcB1.*(PRow - 1)) + (CalcB2.*(PRow));
       NewRow1 = rOWaDJ(NewRow1,NumInc);
       uIndex1 = ((NewRow - 1).*XINC) + (PCol + 1);
       uIndex4 = ((NewRow - 1).*XINC) + (PCol);
       uIndex2 = ((NewRow1 - 1).*XINC) + (PCol + 1);
       uIndex3 = ((NewRow1 - 1).*XINC) + (PCol);
       ax = pw;
       bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
       tUvel = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
          (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
       % Next calculate new v velocity.
       BTemp = ((pw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewCol = (CalcB1.*PCol) + (CalcB2.*(PCol + 1));
       NewCol = cOLaDJ(NewCol,NumInc);
       NewCol1 = (CalcB1.*(PCol - 1)) + (CalcB2.*(PCol));
       NewCol1 = cOLaDJ(NewCol1,NumInc);
       by = qw;
       ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
       vIndex1 = ((PRow - 1).*NUMCOLS) + (NewCol + NUMCOLS);
       vIndex2 = ((PRow - 1).*NUMCOLS) + (NewCol);
       vIndex4 = ((PRow - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
       vIndex3 = ((PRow - 1).*NUMCOLS) + (NewCol1);
       tVvel = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
          (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
       % Adjust for boundary locations.
       Bound4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
           ((tVvel - double(0.0)) <= PRECH));
       Bound2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
           ((tVvel - double(0.0)) >= -PRECH));
       BoundY = Bound4 + Bound2;
       Bound1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
           ((tUvel - double(0.0)) <= PRECH));
       Bound3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
           ((tUvel - double(0.0)) >= -PRECH));
       BoundX = Bound1 + Bound3;
       Bound = BoundX + BoundY;
       BTemp = (Bound <= 0);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       % Finally calculate the new locations.
       XLoc = (CalcB1.*(XLoc - (DeltaT.*tUvel))) + (CalcB2.*XLoc);
       YLoc = (CalcB1.*(YLoc - (DeltaT.*tVvel))) + (CalcB2.*YLoc);
       [XLoc,YLoc] = iNTbOUND(XLoc,YLoc,NumInc);
       % Determine the indexes for the new location.
       [PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NumInc,...
           PCol,PRow,XPos,YPos,PNode);
    % End of loop.
    end
    % not needed
    %clear ax bx ay by Bound BoundX Bound1 Bound3 BoundY Bound2 Bound4;
    %clear BTemp CalcB1 CalcB2 DeltaT maxTau NewCol NewCol1 NewRow NewRow1;
    %clear Number uIndex1 uIndex2 uIndex3 uIndex4 vIndex1 vIndex2 vIndex3;
    %clear vIndex4 XPart YPart;
    %clear tUvel tVvel UMain VMain DTMain MaxN MinN NumInc;
end
%EOF
