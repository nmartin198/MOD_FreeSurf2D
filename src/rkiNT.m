function [XLoc,YLoc,XPos,YPos,PCol,PRow,PNode] = rkiNT(XLoc,YLoc,tU1,tV1,...
   UMain,VMain,DTMain,MaxN,MinN,NumInc)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rkiNT uses Runge-Kutta integration to follow a particle trajectory
    % backwards along a pathline.  Specifically, this is classical, four-step,
    % explicit Runge-Kutta integration.
    %
    % Received:
    %
    % XLoc [NumInc,1] = x-coordinate location or Xe.
    % YLoc [NumInc,1] = y-coordinate location or Ye.
    % tU1 [NumInc,1] = u velocity at initial point.
    % tV1 [NumInc,1] = v velocity at initial point.
    % UMain [NUMINCX,1] = u velocity for Eulerian Grid.
    % VMain [NUMINCY,1] = v velocity for Eulerian Grid.
    % DTMain = total time step == DT.
    % MaxN = maximum number of partial steps.
    % MinN = minimum number of partial steps.
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

    global DX DY PRECH NUMCOLS NUMROWS PREC XINC XINDEX YINDEX

    % local variables.
    ax = double(zeros(NumInc,1));       % Weight for distance from defined location.
    bx = double(zeros(NumInc,1));       % Weight for distance from defined location.
    ay = double(zeros(NumInc,1));       % Weight for distance from defined location.
    by = double(zeros(NumInc,1));       % Weight for distance from defined location.
    Bound = zeros(NumInc,1);            % Boolean boundary location variable.
    BoundX = zeros(NumInc,1);           % Boolean boundary location variable.
    Bound1 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    Bound3 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    BoundY = zeros(NumInc,1);           % Boolean boundary location variable.
    Bound2 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    Bound4 = zeros(NumInc,1);    % Boolean to determine if on a boundary.
    BTemp = zeros(NumInc,1);      % Boolean calc variable.
    CalcB1 = double(zeros(NumInc,1)); % Calc variable.
    CalcB2 = double(zeros(NumInc,1)); % Calc variable.
    DeltaT = double(0.0);               % Partial time step.
    maxTau = double(0.0);               % Maximum calculated partial time step.
    NewCol = zeros(NumInc,1);           % Column location tracking variables.
    NewCol1 = zeros(NumInc,1);          % Column location tracking variable.
    NewRow = zeros(NumInc,1);           % Row location tracking variable.
    NewRow1 = zeros(NumInc,1);          % Row location tracking variable.
    Number = 0;                         % Number of partial steps.
    XPos = double(zeros(NumInc,1));     % X direction volume boundary locations.
    YPos = double(zeros(NumInc,1));     % Y direction volume boudnary locations.
    PCol = zeros(NumInc,1);             % Column locations.
    PRow = zeros(NumInc,1);             % Row locations.
    PNode = zeros(NumInc,1);            % Node Locations.
    X1Pos = double(zeros(NumInc,1));    % X direction volume boundary locations.
    Y1Pos = double(zeros(NumInc,1));    % Y direction volume boudnary locations.
    P1Col = zeros(NumInc,1);            % Column locations for location 1.
    P1Row = zeros(NumInc,1);            % Row locations for location 1.
    P1Node = zeros(NumInc,1);           % Node Locations for location 1.
    X2Pos = double(zeros(NumInc,1));    % X direction volume boundary locations.
    Y2Pos = double(zeros(NumInc,1));    % Y direction volume boudnary locations.
    P2Col = zeros(NumInc,1);            % Column locations for location 2.
    P2Row = zeros(NumInc,1);            % Row locations for location 2.
    P2Node = zeros(NumInc,1);           % Node Locations for location 2.
    X3Pos = double(zeros(NumInc,1));    % X direction volume boundary locations.
    Y3Pos = double(zeros(NumInc,1));    % Y direction volume boudnary locations.
    P3Col = zeros(NumInc,1);            % Column locations for location 3
    P3Row = zeros(NumInc,1);            % Row locations for location 3.
    P3Node = zeros(NumInc,1);           % Node Locations for location 3.
    tU2 = double(zeros(NumInc,1));      % Interpolated u velocity value 2.
    tU3 = double(zeros(NumInc,1));      % Interpolated u velocity value 3.
    tU4 = double(zeros(NumInc,1));      % Interpolated u velocity value 4.
    tV2 = double(zeros(NumInc,1));      % Interpolated v velocity value 2.
    tV3 = double(zeros(NumInc,1));      % Interpolated v velocity value 3.
    tV4 = double(zeros(NumInc,1));      % Interpolated v velocity value 4.
    uIndex1 = zeros(NumInc,1);          % Index of u velocity location 1.
    uIndex2 = zeros(NumInc,1);          % Index of u velocity location 2.
    uIndex3 = zeros(NumInc,1);          % Index of u velocity location 3.
    uIndex4 = zeros(NumInc,1);          % Index of u velocity location 4.
    vIndex1 = zeros(NumInc,1);          % Index of v velocity location 1.
    vIndex2 = zeros(NumInc,1);          % Index of v velocity location 2.
    vIndex3 = zeros(NumInc,1);          % Index of v velocity location 3.
    vIndex4 = zeros(NumInc,1);          % Index of v velocity location 4.
    X1 = double(zeros(NumInc,1));       % X-location for velocity 2.
    X2 = double(zeros(NumInc,1));       % X-location for velocity 3.
    X3 = double(zeros(NumInc,1));       % X-location for velocity 4.
    XPart = double(0.0);                % Min time to cross volume in x-direction.
    Y1 = double(zeros(NumInc,1));       % Y-location for velocity 2.
    Y2 = double(zeros(NumInc,1));       % Y-location for velocity 3.
    Y3 = double(zeros(NumInc,1));       % Y-location for velocity 4.
    YPart = double(0.0);                % Min time to cross volume in y-direction.

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
    % Now set dTau so that divides into dt evenly.
    Number = ceil(DTMain/maxTau);
    Number = max(Number,MinN);
    DeltaT = DTMain/Number;
    % Have tU1 and tV1 so calculate location 1 and get Velocity 2.
    X1 = XLoc - ((0.5*DeltaT).*tU1);
    Y1 = YLoc - ((0.5*DeltaT).*tV1);
    [X1,Y1] = iNTbOUND(X1,Y1,NumInc);
    [P1Col,P1Row,X1Pos,Y1Pos,P1Node] = lOCATIONcALC(X1,Y1,NumInc,...
        P1Col,P1Row,X1Pos,Y1Pos,P1Node);
    % Now start calculating tU2 and tV2
    pw = (X1Pos - X1)./DX;
    qw = (Y1Pos - Y1)./DY;
    BTemp = ((qw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewRow = (CalcB1.*P1Row) + (CalcB2.*(P1Row+1));
    NewRow = rOWaDJ(NewRow,NumInc);
    NewRow1 = (CalcB1.*(P1Row - 1)) + (CalcB2.*(P1Row));
    NewRow1 = rOWaDJ(NewRow1,NumInc);
    uIndex1 = ((NewRow - 1).*XINC) + (P1Col + 1);
    uIndex4 = ((NewRow - 1).*XINC) + (P1Col);
    uIndex2 = ((NewRow1 - 1).*XINC) + (P1Col + 1);
    uIndex3 = ((NewRow1 - 1).*XINC) + (P1Col);
    ax = pw;
    bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
    tU2 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
       (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
    % Next calculate new v velocity.
    BTemp = ((pw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewCol = (CalcB1.*P1Col) + (CalcB2.*(P1Col + 1));
    NewCol = cOLaDJ(NewCol,NumInc);
    NewCol1 = (CalcB1.*(P1Col - 1)) + (CalcB2.*(P1Col));
    NewCol1 = cOLaDJ(NewCol1,NumInc);
    by = qw;
    ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
    vIndex1 = ((P1Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
    vIndex2 = ((P1Row - 1).*NUMCOLS) + (NewCol);
    vIndex4 = ((P1Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
    vIndex3 = ((P1Row - 1).*NUMCOLS) + (NewCol1);
    tV2 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
       (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
    % Next calculate location 2 and get Velocity 3.
    X2 = XLoc - ((0.5*DeltaT).*tU2);
    Y2 = YLoc - ((0.5*DeltaT).*tV2);
    [X2,Y2] = iNTbOUND(X2,Y2,NumInc);
    [P2Col,P2Row,X2Pos,Y2Pos,P2Node] = lOCATIONcALC(X2,Y2,NumInc,...
        P2Col,P2Row,X2Pos,Y2Pos,P2Node);
    % Now start calculating tU3 and tV3
    pw = (X2Pos - X2)./DX;
    qw = (Y2Pos - Y2)./DY;
    BTemp = ((qw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewRow = (CalcB1.*P2Row) + (CalcB2.*(P2Row+1));
    NewRow = rOWaDJ(NewRow,NumInc);
    NewRow1 = (CalcB1.*(P2Row - 1)) + (CalcB2.*(P2Row));
    NewRow1 = rOWaDJ(NewRow1,NumInc);
    uIndex1 = ((NewRow - 1).*XINC) + (P2Col + 1);
    uIndex4 = ((NewRow - 1).*XINC) + (P2Col);
    uIndex2 = ((NewRow1 - 1).*XINC) + (P2Col + 1);
    uIndex3 = ((NewRow1 - 1).*XINC) + (P2Col);
    ax = pw;
    bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
    tU3 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
       (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
    % Next calculate new v velocity.
    BTemp = ((pw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewCol = (CalcB1.*P2Col) + (CalcB2.*(P2Col + 1));
    NewCol = cOLaDJ(NewCol,NumInc);
    NewCol1 = (CalcB1.*(P2Col - 1)) + (CalcB2.*(P2Col));
    NewCol1 = cOLaDJ(NewCol1,NumInc);
    by = qw;
    ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
    vIndex1 = ((P2Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
    vIndex2 = ((P2Row - 1).*NUMCOLS) + (NewCol);
    vIndex4 = ((P2Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
    vIndex3 = ((P2Row - 1).*NUMCOLS) + (NewCol1);
    tV3 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
       (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
    % Next calculate location 3 and get Velocity 4.
    X3 = XLoc - (DeltaT.*tU3);
    Y3 = YLoc - (DeltaT.*tV3);
    [X3,Y3] = iNTbOUND(X3,Y3,NumInc);
    [P3Col,P3Row,X3Pos,Y3Pos,P3Node] = lOCATIONcALC(X3,Y3,NumInc,...
        P3Col,P3Row,X3Pos,Y3Pos,P3Node);
    % Finally calculate tU4 and tV4
    pw = (X3Pos - X3)./DX;
    qw = (Y3Pos - Y3)./DY;
    BTemp = ((qw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewRow = (CalcB1.*P3Row) + (CalcB2.*(P3Row+1));
    NewRow = rOWaDJ(NewRow,NumInc);
    NewRow1 = (CalcB1.*(P3Row - 1)) + (CalcB2.*(P3Row));
    NewRow1 = rOWaDJ(NewRow1,NumInc);
    uIndex1 = ((NewRow - 1).*XINC) + (P3Col + 1);
    uIndex4 = ((NewRow - 1).*XINC) + (P3Col);
    uIndex2 = ((NewRow1 - 1).*XINC) + (P3Col + 1);
    uIndex3 = ((NewRow1 - 1).*XINC) + (P3Col);
    ax = pw;
    bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
    tU4 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
       (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
    % Next calculate new v velocity.
    BTemp = ((pw - 0.5) > PREC);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    NewCol = (CalcB1.*P3Col) + (CalcB2.*(P3Col + 1));
    NewCol = cOLaDJ(NewCol,NumInc);
    NewCol1 = (CalcB1.*(P3Col - 1)) + (CalcB2.*(P3Col));
    NewCol1 = cOLaDJ(NewCol1,NumInc);
    by = qw;
    ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
    vIndex1 = ((P3Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
    vIndex2 = ((P3Row - 1).*NUMCOLS) + (NewCol);
    vIndex4 = ((P3Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
    vIndex3 = ((P3Row - 1).*NUMCOLS) + (NewCol1);
    tV4 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
       (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
    % Adjust for boundary locations.
    Bound4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
       ((tV1 - double(0.0)) <= PRECH));
    Bound2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
       ((tV1 - double(0.0)) >= -PRECH));
    BoundY = Bound4 + Bound2;
    Bound1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
       ((tU1 - double(0.0)) <= PRECH));
    Bound3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
       ((tU1 - double(0.0)) >= -PRECH));
    BoundX = Bound1 + Bound3;
    Bound = BoundX + BoundY;
    BTemp = (Bound <= 0);
    CalcB1 = double(BTemp);
    CalcB2 = double(~BTemp);
    % Finally calculate the new locations using all four velocities.
    XLoc = (CalcB1.*(XLoc - ((DeltaT/6).*(tU1 + (2.*tU2) + (2.*tU3) + tU4)))) + ...
       (CalcB2.*XLoc);
    YLoc = (CalcB1.*(YLoc - ((DeltaT/6).*(tV1 + (2.*tV2) + (2.*tV3) + tV4)))) + ...
       (CalcB2.*YLoc);
    [XLoc,YLoc] = iNTbOUND(XLoc,YLoc,NumInc);
    % First calculate new location one step back.
    [PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NumInc,...
        PCol,PRow,XPos,YPos,PNode);
    % Now loop.
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
       tU1 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
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
       tV1 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
          (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
       % Next calculate location 1 and get Velocity 2.
       X1 = XLoc - ((0.5*DeltaT).*tU1);
       Y1 = YLoc - ((0.5*DeltaT).*tV1);
       [X1,Y1] = iNTbOUND(X1,Y1,NumInc);
       [P1Col,P1Row,X1Pos,Y1Pos,P1Node] = lOCATIONcALC(X1,Y1,NumInc,...
           P1Col,P1Row,X1Pos,Y1Pos,P1Node);
       % Now start calculating tU2 and tV2
       pw = (X1Pos - X1)./DX;
       qw = (Y1Pos - Y1)./DY;
       BTemp = ((qw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewRow = (CalcB1.*P1Row) + (CalcB2.*(P1Row+1));
       NewRow = rOWaDJ(NewRow,NumInc);
       NewRow1 = (CalcB1.*(P1Row - 1)) + (CalcB2.*(P1Row));
       NewRow1 = rOWaDJ(NewRow1,NumInc);
       uIndex1 = ((NewRow - 1).*XINC) + (P1Col + 1);
       uIndex4 = ((NewRow - 1).*XINC) + (P1Col);
       uIndex2 = ((NewRow1 - 1).*XINC) + (P1Col + 1);
       uIndex3 = ((NewRow1 - 1).*XINC) + (P1Col);
       ax = pw;
       bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
       tU2 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
          (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
       % Next calculate new v velocity.
       BTemp = ((pw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewCol = (CalcB1.*P1Col) + (CalcB2.*(P1Col + 1));
       NewCol = cOLaDJ(NewCol,NumInc);
       NewCol1 = (CalcB1.*(P1Col - 1)) + (CalcB2.*(P1Col));
       NewCol1 = cOLaDJ(NewCol1,NumInc);
       by = qw;
       ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
       vIndex1 = ((P1Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
       vIndex2 = ((P1Row - 1).*NUMCOLS) + (NewCol);
       vIndex4 = ((P1Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
       vIndex3 = ((P1Row - 1).*NUMCOLS) + (NewCol1);
       tV2 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
          (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
       % Next calculate location 2 and get Velocity 3.
       X2 = XLoc - ((0.5*DeltaT).*tU2);
       Y2 = YLoc - ((0.5*DeltaT).*tV2);
       [X2,Y2] = iNTbOUND(X2,Y2,NumInc);
       [P2Col,P2Row,X2Pos,Y2Pos,P2Node] = lOCATIONcALC(X2,Y2,NumInc,...
           P2Col,P2Row,X2Pos,Y2Pos,P2Node);
       % Now start calculating tU3 and tV3
       pw = (X2Pos - X2)./DX;
       qw = (Y2Pos - Y2)./DY;
       BTemp = ((qw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewRow = (CalcB1.*P2Row) + (CalcB2.*(P2Row+1));
       NewRow = rOWaDJ(NewRow,NumInc);
       NewRow1 = (CalcB1.*(P2Row - 1)) + (CalcB2.*(P2Row));
       NewRow1 = rOWaDJ(NewRow1,NumInc);
       uIndex1 = ((NewRow - 1).*XINC) + (P2Col + 1);
       uIndex4 = ((NewRow - 1).*XINC) + (P2Col);
       uIndex2 = ((NewRow1 - 1).*XINC) + (P2Col + 1);
       uIndex3 = ((NewRow1 - 1).*XINC) + (P2Col);
       ax = pw;
       bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
       tU3 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
          (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
       % Next calculate new v velocity.
       BTemp = ((pw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewCol = (CalcB1.*P2Col) + (CalcB2.*(P2Col + 1));
       NewCol = cOLaDJ(NewCol,NumInc);
       NewCol1 = (CalcB1.*(P2Col - 1)) + (CalcB2.*(P2Col));
       NewCol1 = cOLaDJ(NewCol1,NumInc);
       by = qw;
       ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
       vIndex1 = ((P2Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
       vIndex2 = ((P2Row - 1).*NUMCOLS) + (NewCol);
       vIndex4 = ((P2Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
       vIndex3 = ((P2Row - 1).*NUMCOLS) + (NewCol1);
       tV3 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
          (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
       % Next calculate location 3 and get Velocity 4.
       X3 = XLoc - (DeltaT.*tU3);
       Y3 = YLoc - (DeltaT.*tV3);
       [X3,Y3] = iNTbOUND(X3,Y3,NumInc);
       [P3Col,P3Row,X3Pos,Y3Pos,P3Node] = lOCATIONcALC(X3,Y3,NumInc,...
           P3Col,P3Row,X3Pos,Y3Pos,P3Node);
       % Finally calculate tU4 and tV4
       pw = (X3Pos - X3)./DX;
       qw = (Y3Pos - Y3)./DY;
       BTemp = ((qw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewRow = (CalcB1.*P3Row) + (CalcB2.*(P3Row+1));
       NewRow = rOWaDJ(NewRow,NumInc);
       NewRow1 = (CalcB1.*(P3Row - 1)) + (CalcB2.*(P3Row));
       NewRow1 = rOWaDJ(NewRow1,NumInc);
       uIndex1 = ((NewRow - 1).*XINC) + (P3Col + 1);
       uIndex4 = ((NewRow - 1).*XINC) + (P3Col);
       uIndex2 = ((NewRow1 - 1).*XINC) + (P3Col + 1);
       uIndex3 = ((NewRow1 - 1).*XINC) + (P3Col);
       ax = pw;
       bx = (CalcB1.*(qw - 0.5)) + (CalcB2.*(qw + 0.5));
       tU4 = ((1 - ax).*((1 - bx).*UMain(uIndex1) + bx.*UMain(uIndex2))) + ...
          (ax.*((1 - bx).*UMain(uIndex4) + bx.*UMain(uIndex3)));
       % Next calculate new v velocity.
       BTemp = ((pw - 0.5) > PREC);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       NewCol = (CalcB1.*P3Col) + (CalcB2.*(P3Col + 1));
       NewCol = cOLaDJ(NewCol,NumInc);
       NewCol1 = (CalcB1.*(P3Col - 1)) + (CalcB2.*(P3Col));
       NewCol1 = cOLaDJ(NewCol1,NumInc);
       by = qw;
       ay = (CalcB1.*(pw - 0.5)) + (CalcB2.*(pw + 0.5));
       vIndex1 = ((P3Row - 1).*NUMCOLS) + (NewCol + NUMCOLS);
       vIndex2 = ((P3Row - 1).*NUMCOLS) + (NewCol);
       vIndex4 = ((P3Row - 1).*NUMCOLS) + (NewCol1 + NUMCOLS);
       vIndex3 = ((P3Row - 1).*NUMCOLS) + (NewCol1);
       tV4 = ((1 - ay).*((1 - by).*VMain(vIndex1) + by.*VMain(vIndex2))) + ...
          (ay.*((1 - by).*VMain(vIndex4) + by.*VMain(vIndex3)));
       % Adjust for boundary locations.
       Bound4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
           ((tV1 - double(0.0)) <= PRECH));
       Bound2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
           ((tV1 - double(0.0)) >= -PRECH));
       BoundY = Bound4 + Bound2;
       Bound1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
           ((tU1 - double(0.0)) <= PRECH));
       Bound3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
           ((tU1 - double(0.0)) >= -PRECH));
       BoundX = Bound1 + Bound3;
       Bound = BoundX + BoundY;
       BTemp = (Bound <= 0);
       CalcB1 = double(BTemp);
       CalcB2 = double(~BTemp);
       % Finally calculate the new locations using all four velocities.
       XLoc = (CalcB1.*(XLoc - ((DeltaT/6).*(tU1 + (2.*tU2) + (2.*tU3) + tU4)))) + ...
          (CalcB2.*XLoc);
       YLoc = (CalcB1.*(YLoc - ((DeltaT/6).*(tV1 + (2.*tV2) + (2.*tV3) + tV4)))) + ...
          (CalcB2.*YLoc);
       [XLoc,YLoc] = iNTbOUND(XLoc,YLoc,NumInc);
       % Determine the indexes for the new location.
       [PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NumInc,...
           PCol,PRow,XPos,YPos,PNode);
    % End of loop.
    end
    % not needed
    %clear ax bx ay by Bound BoundX Bound1 Bound3 BoundY Bound2 Bound4;
    %clear BTemp CalcB1 CalcB2 DeltaT maxTau NewCol NewCol1 NewRow NewRow1;
    %clear Number X1Pos Y1Pos P1Col P1Row P1Node;
    %clear X2Pos Y2Pos P2Col P2Row P2Node X3Pos Y3Pos P3Col P3Row P3Node;
    %clear tU2 tU3 tU4 tV2 tV3 tV4 uIndex1 uIndex2 uIndex3 uIndex4 vIndex1;
    %clear vIndex2 vIndex3 vIndex4 X1 X2 X3 XPart Y1 Y2 Y3 YPart;
    %clear tU1 tV1 UMain VMain DTMain MaxN MinN NumInc;
end
%EOF