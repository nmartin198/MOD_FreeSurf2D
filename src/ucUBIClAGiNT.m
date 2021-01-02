function Vel = ucUBIClAGiNT(Vel,XLoc,YLoc,XPos,YPos,PCol,PNode,PRow,pq,NumInc)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % ucUBIClAGiNT uses cubic Lagrangian interpolation to calculate the 
    % u velocity terms at the "new" location.  The "new" location is obtained
    % by tracing back individual path lines to find the particle location
    % at the beginning of the time step.  This calculation algorithim is 
    % set up directionally.  So, 1.1 is always downwind of the actual
    % particle location.
    %                                                              |
    % Stencil (assuming positive u and v at new location): -> and \|/
    %                                                              -
    %_|____|____|____|____|____|____|____|_
    % |    |    |    |    |    |    |    |   The stencil is shown here
    % |    |    |    |    |    |    |    |    using a.b where a is the
    %_|____|____|____|____|____|____|____|_   row number and b is the
    % |    |    |    |    |    |    |    |    index number of each row.
    % |   4.4  4.3 D4.2  4.1   |    |    |    b will correspond to xIndex*.
    %_|____|____|____|____|____|____|____|_   
    % |    |    |    |    |    |    |    |   After A, B,C, and D are calculated
    % |   3.4  3.3 C3.2  3.1   |    |    |    the interpolated value will be
    %_|____|____|____|____|__ _|____|____|_   obtained from one last interpolation
    % |    |    |  X |    |    |    |    |    between these four values.
    % |   2.4  2.3 B2.2  2.1   |    |    |
    %_|____|____|____|____|____|____|____|_  The final interpolated value
    % |    |    |    |    |    |    |    |    represents the velocity at location X
    % |   1.4  1.3 A1.2  1.1   |    |    |    which is the new location obtained
    %_|____|____|____|____|____|____|____|_   by tracking back the streamline.
    % |    |    |    |    |    |    |    |
    %
    % Received:
    %
    % Vel = x-face velocity == u.
    % XLoc = x-coordinate particle location (Xe).
    % YLoc = y-coordinate particle location (Ye).
    % XPos = x-coordinate of greater volume boundary (n).
    % YPos = y-coordinate of greater volume boundary (m).
    % PCol = column index of particle location.
    % PRow = row index of volume location.
    % pq = normalized distance from YPos to YLoc.
    % NumInc = number of face indexes.
    %
    % Returned:
    %
    % Vel = interpolated velocity.
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
    %
    % globals.
    global BHux DX DY PRECH PREC UAve VAve XINC

    % First allocate new variables.
    A = double(zeros(NumInc,1));        % Interpolated velocity at location A.
    B = double(zeros(NumInc,1));        % Interpolated velocity at location B.
    BTempN = zeros(NumInc,1);           % Boolean temporary variable.
    BTempU = zeros(NumInc,1);           % Boolean temporary variable.
    BTempV = zeros(NumInc,1);           % Boolean temporary variable.
    CalcN1 = double(zeros(NumInc,1));   % Calc variable.
    CalcN2 = double(zeros(NumInc,1));   % Calc variable.
    CalcU1 = double(zeros(NumInc,1));   % Calc variable.
    CalcU2 = double(zeros(NumInc,1));   % Calc variable.
    CalcV1 = double(zeros(NumInc,1));   % Calc variable.
    CalcV2 = double(zeros(NumInc,1));   % Calc variable.
    C = double(zeros(NumInc,1));        % Interpolated velocity at location C.
    Col1 = zeros(NumInc,1);             % Column index for column location 1.
    Col2 = zeros(NumInc,1);             % Column index for column location 1.
    Col3 = zeros(NumInc,1);             % Column index for column location 1.
    Col4 = zeros(NumInc,1);             % Column index for column location 1.
    D = double(zeros(NumInc,1));        % Interpolated velocity at location D.
    DirectU = zeros(NumInc,1);          % u velocity direction.
    DirectV = zeros(NumInc,1);          % v velocity direction.
    LWeight1 = double(zeros(NumInc,1)); % Weight for Index1.
    LWeight2 = double(zeros(NumInc,1)); % Weight for Index2.
    LWeight3 = double(zeros(NumInc,1)); % Weight for Index3.
    LWeight4 = double(zeros(NumInc,1)); % Weight for Index4.
    Row1 = zeros(NumInc,1);             % Row index for row 1.
    Row2 = zeros(NumInc,1);             % Row index for row 2.
    Row3 = zeros(NumInc,1);             % Row index for row 3.
    Row4 = zeros(NumInc,1);             % Row index for row 4.
    RowLoc = zeros(NumInc,1);           % Volume index of the last column in the previous row.
    X1 = double(zeros(NumInc,1));       % X - coordinate of column 1.
    X2 = double(zeros(NumInc,1));       % X - coordinate of column 2.
    X3 = double(zeros(NumInc,1));       % X - coordinate of column 3.
    X4 = double(zeros(NumInc,1));       % X - coordinate of column 4.
    Y1 = double(zeros(NumInc,1));       % Y - coordinate of row 1.
    Y2 = double(zeros(NumInc,1));       % Y - coordinate of row 2.
    Y3 = double(zeros(NumInc,1));       % Y - coordinate of row 3.
    Y4 = double(zeros(NumInc,1));       % Y - coordinate of row 4.
    xIndex1 = zeros(NumInc,1);          % u Index for *.1.
    xIndex2 = zeros(NumInc,1);          % u Index for *.2.
    xIndex3 = zeros(NumInc,1);          % u Index for *.3.
    xIndex4 = zeros(NumInc,1);          % u Index for *.4.

    % Get directions for stencil.
    BTempN = ((UAve(PNode) - double(0.0)) >= -PREC);
    CalcN1 = double(BTempN);
    CalcN2 = double(~BTempN);
    DirectU = (CalcN1.*1.0) + (CalcN2.*-1.0);
    BTempN = ((VAve(PNode) - double(0.0)) >= -PREC);
    CalcN1 = double(BTempN);
    CalcN2 = double(~BTempN);
    DirectV = (CalcN1.*1.0) + (CalcN2.*-1.0);
    % Start with row 2 location 2 -> 2.2;
    % First get the row and column locations for Row 2 and Column 2.
    BTempU = (DirectU == 1);
    CalcU1 = double(BTempU);
    CalcU2 = double(~BTempU);
    Col2 = (CalcU1.*(PCol+1)) + (CalcU2.*PCol);
    X2 = (CalcU1.*XPos) + (CalcU2.*(XPos - DX));
    BTempN = ((pq - 0.5) > PREC);
    CalcN1 = double(BTempN);
    CalcN2 = double(~BTempN);
    BTempV = (DirectV == 1);
    CalcV1 = double(BTempV);
    CalcV2 = double(~BTempV);
    Row2 = (CalcV1.*((CalcN1.*PRow) + (CalcN2.*(PRow+1)))) + ...
       (CalcV2.*((CalcN1.*(PRow-1)) + (CalcN2.*(PRow))));
    Y2 = (CalcV1.*((CalcN1.*(YPos - (0.5*DY))) + (CalcN2.*(YPos + (0.5*DY))))) +...
       (CalcV2.*((CalcN1.*(YPos - (1.5*DY))) + (CalcN2.*(YPos - (0.5*DY)))));
    Row2 = rOWaDJ(Row2,NumInc);
    % Go ahead and determine all row and column locations.
    Col1 = (CalcU1.*(Col2 + 1)) + (CalcU2.*(Col2 - 1));
    Col1 = cOLbaDJ(Col1,NumInc);
    Col3 = (CalcU1.*(Col2 - 1)) + (CalcU2.*(Col2 + 1));
    Col3 = cOLbaDJ(Col3,NumInc);
    Col4 = (CalcU1.*(Col3 - 1)) + (CalcU2.*(Col3 + 1));
    Col4 = cOLbaDJ(Col4,NumInc);
    Row1 = (CalcV1.*(Row2 + 1)) + (CalcV2.*(Row2 - 1));
    Row1 = rOWaDJ(Row1,NumInc);
    Row3 = (CalcV1.*(Row2 - 1)) + (CalcV2.*(Row2 + 1));
    Row3 = rOWaDJ(Row3,NumInc);
    Row4 = (CalcV1.*(Row3 - 1)) + (CalcV2.*(Row3 + 1));
    Row4 = rOWaDJ(Row4,NumInc);
    X1 = (CalcU1.*(X2 + DX)) + (CalcU2.*(X2 - DX));
    X3 = (CalcU1.*(X2 - DX)) + (CalcU2.*(X2 + DX));
    X4 = (CalcU1.*(X3 - DX)) + (CalcU2.*(X3 + DX));
    Y1 = (CalcV1.*(Y2 + DY)) + (CalcV2.*(Y2 - DY));
    Y3 = (CalcV1.*(Y2 - DY)) + (CalcV2.*(Y2 + DY));
    Y4 = (CalcV1.*(Y3 - DY)) + (CalcV2.*(Y3 + DY));
    % Now calculate weight values for the cubic Lagrange interpolation along
    % all rows.  Can do this because have the same column location for each
    % row.
    LWeight1 = ((XLoc - X2).*(XLoc - X3).*(XLoc - X4))./...
       ((X1 - X2).*(X1 - X3).*(X1 - X4));
    LWeight2 = ((XLoc - X1).*(XLoc - X3).*(XLoc - X4))./...
       ((X2 - X1).*(X2 - X3).*(X2 - X4));
    LWeight3 = ((XLoc - X1).*(XLoc - X2).*(XLoc - X4))./...
       ((X3 - X1).*(X3 - X2).*(X3 - X4));
    LWeight4 = ((XLoc - X1).*(XLoc - X2).*(XLoc - X3))./...
       ((X4 - X1).*(X4 - X2).*(X4 - X3));
    % Now get the index locations to calculate B.
    RowLoc = ((Row2 - 1).*XINC);
    xIndex1 = RowLoc + Col1;
    xIndex2 = RowLoc + Col2;
    xIndex3 = RowLoc + Col3;
    xIndex4 = RowLoc + Col4;
    % Now interpolate the value at B.
    B = (LWeight1.*Vel(xIndex1)) + (LWeight2.*Vel(xIndex2)) + ...
       (LWeight3.*Vel(xIndex3)) + (LWeight4.*Vel(xIndex4));
    % Next do A.
    % Now get the index locations to calculate B.
    RowLoc = ((Row1 - 1).*XINC);
    xIndex1 = RowLoc + Col1;
    xIndex2 = RowLoc + Col2;
    xIndex3 = RowLoc + Col3;
    xIndex4 = RowLoc + Col4;
    % Now interpolate the value at B.
    A = (LWeight1.*Vel(xIndex1)) + (LWeight2.*Vel(xIndex2)) + ...
       (LWeight3.*Vel(xIndex3)) + (LWeight4.*Vel(xIndex4));
    % Next do C.
    % Now get the index locations to calculate B.
    RowLoc = ((Row3 - 1).*XINC);
    xIndex1 = RowLoc + Col1;
    xIndex2 = RowLoc + Col2;
    xIndex3 = RowLoc + Col3;
    xIndex4 = RowLoc + Col4;
    % Now interpolate the value at B.
    C = (LWeight1.*Vel(xIndex1)) + (LWeight2.*Vel(xIndex2)) + ...
       (LWeight3.*Vel(xIndex3)) + (LWeight4.*Vel(xIndex4));
    % Finally do D.
    % Now get the index locations to calculate B.
    RowLoc = ((Row4 - 1).*XINC);
    xIndex1 = RowLoc + Col1;
    xIndex2 = RowLoc + Col2;
    xIndex3 = RowLoc + Col3;
    xIndex4 = RowLoc + Col4;
    % Now interpolate the value at B.
    D = (LWeight1.*Vel(xIndex1)) + (LWeight2.*Vel(xIndex2)) + ...
       (LWeight3.*Vel(xIndex3)) + (LWeight4.*Vel(xIndex4));
    % Now calculate the weights for the interpolation along the column of interpolated
    % values.
    LWeight1 = ((YLoc - Y2).*(YLoc - Y3).*(YLoc - Y4))./...
       ((Y1 - Y2).*(Y1 - Y3).*(Y1 - Y4));
    LWeight2 = ((YLoc - Y1).*(YLoc - Y3).*(YLoc - Y4))./...
       ((Y2 - Y1).*(Y2 - Y3).*(Y2 - Y4));
    LWeight3 = ((YLoc - Y1).*(YLoc - Y2).*(YLoc - Y4))./...
       ((Y3 - Y1).*(Y3 - Y2).*(Y3 - Y4));
    LWeight4 = ((YLoc - Y1).*(YLoc - Y2).*(YLoc - Y3))./...
       ((Y4 - Y1).*(Y4 - Y2).*(Y4 - Y3));
    % Finally calculate the velocity.  Then adjust velocity for areas near
    % the domain boundary.
    Vel = (LWeight1.*A) + (LWeight2.*B) + (LWeight3.*C) + (LWeight4.*D);
    % Ensure positive depth.
    Vel = BHux.*Vel;
    % Check for precision.
    BTempN = (abs(Vel - double(0.0)) > PRECH);
    CalcN1 = double(BTempN);
    Vel = CalcN1.*Vel;
    % not needed
    %clear A B BTempN BTempU BTempV C CalcN1 CalcN2 CalcU1 CalcU2 CalcV1;
    %clear CalcV2Col1 Col2 Col3 Col4 D LWeight1 LWeight2;
    %clear LWeight3 LWeight4 Row1 Row2 Row3 Row4 RowLoc DirectU DirectV;
    %clear X1 X2 X3 X4 Y1 Y2 Y3 Y4 xIndex1 xIndex2 xIndex3 xIndex4;
    %clear XLoc YLoc XPos YPos PCol PNode PRow pq NumInc;
end
%EOF