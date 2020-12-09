function Vel = vcUBIClAGiNT(Vel,XLoc,YLoc,XPos,YPos,PCol,PNode,PRow,pp,NumInc)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% vcUBIClAGiNT uses cubic Lagrangian interpolation to calculate the 
% v velocity terms at the "new" location.  The "new" location is obtained
% by tracing back individual path lines to find the particle location
% at the beginning of the time step.  The stencil below is always 
% oriented according to velocity direction.  So, 1.1 will always be
% downwind of the current particle location.
%                                                              |
% Stencil (assuming positive u and v at new location): -> and \|/
%                                                              -
%_|____|____|____|____|____|____|____|_  
% |    |    |    |    |    |    |    |  
% |    |    |    |    |    |    |    |  
%_|____|____|_4.4|_4.3|D4.2|_4.1|____|_
% |    |    |    |    |    |    |    |   The stencil is shown here
% |    |    |    |    |    |    |    |    using a.b where a is the
%_|____|____|_3.4|_3.3|C3.2|_3.1|____|_   column number and b is the
% |    |    |    |    |    |    |    |    index number of each row.
% |    |    |    |    |X   |    |    |    b will correspond to yIndex*.
%_|____|____|_2.4|_2.3|B2.2|_2.1|____|_   
% |    |    |    |    |    |    |    |   After A, B,C, and D are calculated
% |    |    |    |    |    |    |    |    the interpolated value will be
%_|____|____|_1.4|_1.3|A1.2|_1.1|____|_   obtained from one last interpolation
% |    |    |    |    |    |    |    |    between these four values.
% |    |    |    |    |    |    |    |
%_|____|____|____|____|____|____|____|_  The final interpolated value
% |    |    |    |    |    |    |    |    represents the velocity at location X
% |    |    |    |    |    |    |    |    which is the new location obtained
%_|____|____|____|____|____|____|____|_   by tracking back the streamline.
% |    |    |    |    |    |    |    |
%
% Received:
%
% Vel = y-face velocity == v.
% XLoc = x-coordinate particle location (Xe).
% YLoc = y-coordinate particle location (Ye).
% XPos = x-coordinate of greater volume boundary (n).
% YPos = y-coordinate of greater volume boundary (m).
% PCol = column index of particle location.
% PNode = volume index of particle location.
% PRow = row index of volume location.
% pp = normalized distance from XPos to XLoc.
% NumInc = number of face indexes.
%
% Returned:
%
% Vel = interpolated velocity.
%
% globals
global BHvy DX DY PRECH NUMCOLS PREC UAve VAve

% First allocate new variables.
A = double(zeros(NumInc,1));        % Interpolated velocity at location A.
B = double(zeros(NumInc,1));        % Interpolated velocity at location B.
BTempN = zeros(NumInc,1);           % Boolean calculation variable.
BTempU = zeros(NumInc,1);           % Boolean calculation variable.
BTempV = zeros(NumInc,1);           % Boolean calculation variable.
C = double(zeros(NumInc,1));        % Interpolated velocity at location C.
CalcN1 = double(zeros(NumInc,1));   % Calc variable.
CalcN2 = double(zeros(NumInc,1));   % Calc variable.
CalcU1 = double(zeros(NumInc,1));   % Calc variable.
CalcU2 = double(zeros(NumInc,1));   % Calc variable.
CalcV1 = double(zeros(NumInc,1));   % Calc variable.
CalcV2 = double(zeros(NumInc,1));   % Calc variable.
Col1 = zeros(NumInc,1);             % Column index for column location 1.
Col2 = zeros(NumInc,1);             % Column index for column location 2.
Col3 = zeros(NumInc,1);             % Column index for column location 3.
Col4 = zeros(NumInc,1);             % Column index for column location 4.
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
RowLoc = zeros(NumInc,1);           % Index at start of current row.
X1 = double(zeros(NumInc,1));       % X - coordinate of column 1.
X2 = double(zeros(NumInc,1));       % X - coordinate of column 2.
X3 = double(zeros(NumInc,1));       % X - coordinate of column 3.
X4 = double(zeros(NumInc,1));       % X - coordinate of column 4.
Y1 = double(zeros(NumInc,1));       % Y - coordinate of row 1.
Y2 = double(zeros(NumInc,1));       % Y - coordinate of row 2.
Y3 = double(zeros(NumInc,1));       % Y - coordinate of row 3.
Y4 = double(zeros(NumInc,1));       % Y - coordinate of row 4.
yIndex1 = zeros(NumInc,1);          % v Index for *.1.
yIndex2 = zeros(NumInc,1);          % v Index for *.2.
yIndex3 = zeros(NumInc,1);          % v Index for *.3.
yIndex4 = zeros(NumInc,1);          % v Index for *.4.

% Calculations.
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
BTempV = (DirectV == 1);
CalcV1 = double(BTempV);
CalcV2 = double(~BTempV);
BTempU = (DirectU == 1);
CalcU1 = double(BTempU);
CalcU2 = double(~BTempU);
Row2 = (CalcV1.*(PRow+1)) + (CalcV2.*PRow);
Y2 = (CalcV1.*YPos) + (CalcV2.*(YPos - DY));
BTempN = ((pp - 0.5) > PREC);
CalcN1 = double(BTempN);
CalcN2 = double(~BTempN);
Col2 = (CalcV1.*((CalcN1.*PCol) + (CalcN2.*(PCol+1)))) + ...
   (CalcV2.*((CalcN1.*(PCol - 1)) + (CalcN2.*PCol)));
X2 = (CalcV1.*((CalcN1.*(XPos - (0.5*DX))) + ...
   (CalcN2.*(XPos + (0.5*DX))))) + (CalcV2.*...
   ((CalcN1.*(XPos - (1.5*DX))) + (CalcN2.*(XPos - (0.5*DX)))));
Col2 = cOLaDJ(Col2,NumInc);
% Go ahead and determine all row and column locations.
Row1 = (CalcV1.*(Row2 + 1)) + (CalcV2.*(Row2 - 1));
Row1 = rOWbaDJ(Row1,NumInc);
Row3 = (CalcV1.*(Row2 - 1)) + (CalcV2.*(Row2 + 1));
Row3 = rOWbaDJ(Row3,NumInc);
Row4 = (CalcV1.*(Row3 - 1)) + (CalcV2.*(Row3 + 1));
Row4 = rOWbaDJ(Row4,NumInc);
Col1 = (CalcU1.*(Col2 + 1)) + (CalcU2.*(Col2 - 1));
Col1 = cOLaDJ(Col1,NumInc);
Col3 = (CalcU1.*(Col2 - 1)) + (CalcU2.*(Col2 + 1));
Col3 = cOLaDJ(Col3,NumInc);
Col4 = (CalcU1.*(Col3 - 1)) + (CalcU2.*(Col3 + 1));
Col4 = cOLaDJ(Col4,NumInc);
Y1 = (CalcV1.*(Y2 + DY)) + (CalcV2.*(Y2 - DY));
Y3 = (CalcV1.*(Y2 - DY)) + (CalcV2.*(Y2 + DY));
Y4 = (CalcV1.*(Y3 - DY)) + (CalcV2.*(Y3 + DY));
X1 = (CalcU1.*(X2 + DX)) + (CalcU2.*(X2 - DX));
X3 = (CalcU1.*(X2 - DX)) + (CalcU2.*(X2 + DX));
X4 = (CalcU1.*(X3 - DX)) + (CalcU2.*(X3 + DX));
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
RowLoc = ((Row2 - 1).*NUMCOLS);
yIndex1 = RowLoc + Col1;
yIndex2 = RowLoc + Col2;
yIndex3 = RowLoc + Col3;
yIndex4 = RowLoc + Col4;
% Now interpolate the value at B.
B = (LWeight1.*Vel(yIndex1)) + (LWeight2.*Vel(yIndex2)) + ...
   (LWeight3.*Vel(yIndex3)) + (LWeight4.*Vel(yIndex4));
% Next do A.
RowLoc = ((Row1 - 1).*NUMCOLS);
yIndex1 = RowLoc + Col1;
yIndex2 = RowLoc + Col2;
yIndex3 = RowLoc + Col3;
yIndex4 = RowLoc + Col4;
% Now interpolate the value at A.
A = (LWeight1.*Vel(yIndex1)) + (LWeight2.*Vel(yIndex2)) + ...
   (LWeight3.*Vel(yIndex3)) + (LWeight4.*Vel(yIndex4));
% Next do C.
RowLoc = ((Row3 - 1).*NUMCOLS);
yIndex1 = RowLoc + Col1;
yIndex2 = RowLoc + Col2;
yIndex3 = RowLoc + Col3;
yIndex4 = RowLoc + Col4;
% Now interpolate the value at C.
C = (LWeight1.*Vel(yIndex1)) + (LWeight2.*Vel(yIndex2)) + ...
   (LWeight3.*Vel(yIndex3)) + (LWeight4.*Vel(yIndex4));
% Finally do D.
RowLoc = ((Row4 - 1).*NUMCOLS);
yIndex1 = RowLoc + Col1;
yIndex2 = RowLoc + Col2;
yIndex3 = RowLoc + Col3;
yIndex4 = RowLoc + Col4;
% Now interpolate the value at D.
D = (LWeight1.*Vel(yIndex1)) + (LWeight2.*Vel(yIndex2)) + ...
   (LWeight3.*Vel(yIndex3)) + (LWeight4.*Vel(yIndex4));
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
Vel = BHvy.*Vel;
% Check for precision.
BTempN = (abs(Vel - double(0.0)) > PRECH);
CalcN1 = double(BTempN);
Vel = CalcN1.*Vel;

clear A B BTempN BTempU BTempV C CalcN1 CalcN2 CalcU1 CalcU2 CalcV1 CalcV2;
clear Col1 Col2 Col3 Col4 D LWeight1 LWeight2;
clear LWeight3 LWeight4 Row1 Row2 Row3 Row4 RowLoc DirectU DirectV;
clear X1 X2 X3 X4 Y1 Y2 Y3 Y4 yIndex1 yIndex2 yIndex3 yIndex4;
clear XLoc YLoc XPos YPos PCol PNode PRow pp NumInc;
return;
%EOF