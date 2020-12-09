function [XLoc,YLoc,XPos,YPos,PCol,PRow,PNode] = sEMIapATH(XLoc,YLoc,Vxp,Vyp,Vx2,...
   Vy2,Vx1,Vy1,UVel,VVel,Ax,Ay,X2,Y2,Ddx,Ddy,DT,Bh,NUMINC)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sEMIapATH employs the semianalytical path line tracing method of Pollock (1988)
% to trace back the characteristic from a particular volume location to the 
% equivalent particle location at the beginning of the previous time step.
% Can be used for either volume center or volume face locations.  The semi-
% analytical pathline tracing method allows the movement back across one
% complete volume at a time.  This entails fewer calculations for large
% Courant numbers.
%
% In addition to finding that cubic interpolation can be considered a 
% good compromise between accuracy and efficiency for the quantity of 
% interest at the foot (or previous time step) location, Staniforth and 
% Cote (1991) also determined that bilinear interpolation (which is
% employed in the method of Pollock (1988)) is sufficient for interpolation
% when tracing back the trajectory or characteristic.
%
% XLoc = Xe or the x-coordinate of the initial location.
% YLoc = Ye or the y-coordinate of the initial location.
% Vxp = x-direction particle velccity.
% Vyp = y-direction particle velocity.
% Vx2 = u velocity at downwind x-face.
% Vy2 = v velocity at downwind y-face.
% Vx1 = u velocity at upwind x-face.
% Vy1 = v velocity at upwind y-face.
% UVel = u velocity at x-faces (u) [NUMINCX,1].
% VVel = v velocity at y-faces (v) [NUMINCY,1].
% Ax = x-direction velocity gradient in current volume.
% Ay = y-direction velocity gradient in current volume.
% X2 = x-coordinate of downstream x-face.
% Y2 = y-coordinate of downstream y-face.
% Ddx = x direction distance from the particle to volume boundary.
% Ddy = y direction distance from the particle to volume boundary.
% DT = dt or the fluid time step.
% Bh = boolean for wet face. 1 == wet. 0 == dry.
% NUMINC = the number of indexes of the face that are dealing with.

% global variables.
global DX DY PRECH MAXCR NUMCOLS NUMROWS PREC XINC XINDEX YINDEX

% Initialize local variables.
Adx = double(zeros(NUMINC,1));         % X-direction distance from particle to downstream boundary.
Ady = double(zeros(NUMINC,1));         % Y-direction distance from particle to downstream boundary.
BTemp = zeros(NUMINC,1);               % Boolean temporary variable.
Calc1 = double(zeros(NUMINC,1));       % Calc variable.
Calc2 = double(zeros(NUMINC,1));       % Calc variable.
CalcAx1 = double(zeros(NUMINC,1));       % Calc variable.
CalcAx2 = double(zeros(NUMINC,1));       % Calc variable.
CalcAy1 = double(zeros(NUMINC,1));       % Calc variable.
CalcAy2 = double(zeros(NUMINC,1));       % Calc variable.
CalcX1 = double(zeros(NUMINC,1));       % Calc variable.
CalcX2 = double(zeros(NUMINC,1));       % Calc variable.
CalcY1 = double(zeros(NUMINC,1));       % Calc variable.
CalcY2 = double(zeros(NUMINC,1));       % Calc variable.
BoundX = zeros(NUMINC,1);              % Boolean for if on domain boundary.
BoundX1 = zeros(NUMINC,1);              % Boolean for if on domain boundary.
BoundX3 = zeros(NUMINC,1);              % Boolean for if on domain boundary.
BoundY = zeros(NUMINC,1);              % Boolean for if on domain boundary.
BoundY2 = zeros(NUMINC,1);              % Boolean for if on domain boundary.
BoundY4 = zeros(NUMINC,1);              % Boolean for if on domain boundary.
DeltaTE = double(zeros(NUMINC,1));     % Smallest tracking vector.
DeltaTX = double(zeros(NUMINC,1));     % Time increment to cross the current node in x.
DeltaTY = double(zeros(NUMINC,1));     % Time increment to cross the current node in y.
LOOPINC = 0;                           % Boolean looping criteria.
PartTX = double(zeros(NUMINC,1));      % Part of the time increment calculation.
PartTY = double(zeros(NUMINC,1));      % Part of the time increment calculation.
PCol = zeros(NUMINC,1);                % particle column index.
PNode = zeros(NUMINC,1);               % particle node index.
PRow = zeros(NUMINC,1);                % particle row index.
TempCalc1 = double(zeros(NUMINC,1));   % Temporary calculation variable.
TempCalc2 = double(zeros(NUMINC,1));   % Temporary calculation variable.
TempCalc3 = double(zeros(NUMINC,1));   % Temporary calculation variable.
TempTE = double(zeros(NUMINC,1));      % Temporary time step.
uq = double(zeros(NUMINC,1));          % q weights for u terms.
vp = double(zeros(NUMINC,1));          % p weights for v terms.
xIndex1 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
xIndex2 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
xIndex3 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
xIndex4 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
XPos = double(zeros(NUMINC,1));        % X-coordinate of larger volume boundary.
yIndex1 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
yIndex2 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
yIndex3 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
yIndex4 = zeros(NUMINC,1);             % Index vector for bilinear interpolations.
YPos = double(zeros(NUMINC,1));        % Y-coordinate of larger volume boundary.

% Now calculate DeltaTE.  Use DeltaTE to determine whether will need to
% go back through more than one volume.  If DeltaTE is less than DT than 
% will need to keep tracking backwards until DeltaTE = DT.

% PartTX and PartTY are the parts of DeltaTE that are inside the log().
% Need to ensure that the divisor is not equal to zero, and to set the 
% value to zero if the divisor is equal to zero.  Because are going 
% backwards, flip the values from in relation to Pollock (1988).
BTemp = (abs(Vx1 - 0.0) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*Vx1) + (Calc2.*PREC);
TempCalc2 = 1./TempCalc1;
PartTX = TempCalc2.*Vxp;
BTemp = (abs(Vy1 - 0.0) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*Vy1) + (Calc2.*PREC);
TempCalc2 = 1./TempCalc1;
PartTY = TempCalc2.*Vyp;
% Calculate the time needed to reach a volume boundary in the x-direction.
% Make sure to avoid divide by zero and log(0) errors.  Do the same for y.
BTemp = (abs(Ax - double(0.0)) > PREC);
CalcAx1 = double(BTemp);
CalcAx2 = double(~BTemp);
TempCalc1 = (CalcAx1.*Ax) + (CalcAx2.*PREC);
TempCalc2 = 1./TempCalc1;
BTemp = ((PartTX - double(0.0)) > PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*PartTX) + (Calc2.*PRECH);
DeltaTX = abs(TempCalc2.*log(TempCalc1));
BTemp = (abs(Ay - double(0.0)) > PREC);
CalcAy1 = double(BTemp);
CalcAy2 = double(~BTemp);
TempCalc1 = (CalcAy1.*Ay) + (CalcAy2.*PREC);
TempCalc2 = 1./TempCalc1;
BTemp = ((PartTY - double(0.0)) > PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*PartTY) + (Calc2.*PRECH);
DeltaTY = abs(TempCalc2.*log(TempCalc1));
% Now adjust for the situations where GradX == 0, PartTX = 0, and PartTX == 1.
%  When PartTX == 0, then the exit velocity for the volume is zero.  In this case,
%  the particle will not leave the volume in that direction, so DeltaT = DT.
%  When GradX == 0, then the velocity gradient is zero, and the particle will
%  travel in a straight line.  When PartTX == 1, then ln(1) == 0 so the calculation
%  will be screwed up.  Generally this occurs close to a boundary.
%  Do the same for y quantities.
BTemp = (abs(Vxp - double(0.0)) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*abs(Vxp)) + (Calc2.*PREC);
TempCalc2 = 1./TempCalc1;
DeltaTX = (CalcAx2.*(TempCalc2.*Ddx)) + (CalcAx1.*DeltaTX);
BTemp = (abs(PartTX - double(0.0)) <= PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTX = (Calc1.*DT) + (Calc2.*DeltaTX);
BTemp = (abs(PartTX - 1.0) <= PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTX = (Calc1.*(TempCalc2.*Ddx)) + (Calc2.*DeltaTX);
BTemp = (abs(Vyp - double(0.0)) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempCalc1 = (Calc1.*abs(Vyp)) + (Calc2.*PREC);
TempCalc2 = 1./TempCalc1;
DeltaTY = (CalcAy2.*(TempCalc2.*Ddy)) + (CalcAy1.*DeltaTY);
BTemp = (abs(PartTY - double(0.0)) <= PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTY = (Calc1.*DT) + (Calc2.*DeltaTY);
BTemp = (abs(PartTY - 1.0) <= PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTY = (Calc1.*(TempCalc2.*Ddy)) + (Calc2.*DeltaTY);
% Now calculate the time traveling backwards along the pathline for each particle
% originally located at the x-face.  Use the smaller time of DeltaTX and DeltaTY.
% Now create the new time step to move back through the next cell.
BTemp = ((DeltaTX - DeltaTY) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempTE = (Calc2.*DeltaTX) + (Calc1.*DeltaTY);
BTemp = (((TempTE + DeltaTE) - DT) < -PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
TempTE = (Calc1.*TempTE) + (Calc2.*(TempTE - ((TempTE + DeltaTE)- DT)));
% Adjust TempTE so that if have no depth than do not step backwards
%   - ignore dry cells.
TempTE = (Bh.*TempTE);
% Now calculate new value of DeltaTE to determine if will loop again.
BTemp = ((DeltaTE - DT) < -PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTE = (Calc1.*(DeltaTE + TempTE)) + (Calc2.*DT);
% Adjust TempTe so that if have no depth than do not step backwards
%   - ignore dry cells.
DeltaTE = ((1 - Bh).*DT) + (Bh.*DeltaTE);
% Adjust so that DeltaTE is not equal to zero.
BTemp = ((DeltaTE - double(0.0)) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
DeltaTE = (Calc1.*DeltaTE) + (Calc2.*PREC);
% Now adjust the time step if at domain boundaries.  Want the particle to 
%  stop when it hits a boundary.
BoundY2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
   ((Vyp - double(0.0)) >= -PRECH));
BoundY4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
   ((Vyp - double(0.0)) <= PRECH));
BoundY = BoundY4 + BoundY2;
CalcY1 = double(BoundY);
CalcY2 = double(1 - BoundY);
TempTE = (CalcY2.*TempTE);
DeltaTE = (CalcY1.*DT) + (CalcY2.*DeltaTE);
BoundX3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
   ((Vxp - double(0.0)) >= -PRECH));
BoundX1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
   ((Vxp - double(0.0)) <= PRECH));
BoundX = BoundX1 + BoundX3;
CalcX1 = double(BoundX);
CalcX2 = double(1 - BoundX);
TempTE = (CalcX2.*TempTE);
DeltaTE = (CalcX1.*DT) + (CalcX2.*DeltaTE);
% Now have the velocities, time steps, and gradients so calculate
% the new position.  First calculate position and avoid divide by
% zero erros.  If will have a zero gradient correct the value of
% XLoc to reflect this.  The divide by zero correction will give an
% incorrect value for XLoc.
TempCalc1 = (CalcAx1.*Ax) + (CalcAx2.*PREC);
% TempCalc2 contains a correction for zero gradient.
TempCalc2 = 1./TempCalc1;
% Set-up TempCalc1 to represent solution with a velocity gradient and
%  TempCalc2 to provide the solution in the absence of a gradient.
TempCalc1 = (CalcAx1.*(X2 - (TempCalc2.*(Vx2 - (Vxp./exp(Ax.*TempTE))))));
TempCalc2 = (CalcAx2.*(XLoc - (Vxp.*TempTE)));
XLoc = (CalcX2.*(TempCalc2 + TempCalc1)) + (CalcX1.*XLoc);
TempCalc1 = (CalcAy1.*Ay) + (CalcAy2.*PREC);
TempCalc2 = 1./TempCalc1;
TempCalc1 = (CalcAy1.*(Y2 - (TempCalc2.*(Vy2 - (Vyp./exp(Ay.*TempTE))))));
TempCalc2 = (CalcAy2.*(YLoc - (Vyp.*TempTE)));
YLoc = (CalcY2.*(TempCalc2 + TempCalc1)) + (CalcY1.*YLoc);
% Next adjust the location to be within the x-velocity face framework.  The location
% calculations may permit the Ye value to fall outside of the x-face framework.
% The Xe correction is included, but should not be needed.
BTemp = ((YINDEX(1) - YLoc) < -PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
YLoc = (Calc2.*YINDEX(1)) + (Calc1.*YLoc);
BTemp = ((YINDEX(NUMROWS+1) - YLoc) > PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
YLoc = (Calc2.*YINDEX(NUMROWS+1)) + (Calc1.*YLoc);
BTemp = ((XINDEX(1) - XLoc) < -PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
XLoc = (Calc2.*XINDEX(1)) + (Calc1.*XLoc);
BTemp = ((XINDEX(NUMCOLS+1) - XLoc) > PRECH);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
XLoc = (Calc2.*XINDEX(NUMCOLS+1)) + (Calc1.*XLoc);
% First need to get the present locations from the simulation grid.
[PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NUMINC,PCol,PRow,...
    XPos,YPos,PNode);
% The particle location is now defined for the initial step backwards along the
% pathline.  Particles will either have moved back a distance equivalent to the
% entire time step, or will have moved back to a volume boundary.  The semi-
% analytic pathline solution only applies within a volume.  If need to traverse
% more than one volume, then need to track the pathline in a piecewise manner.
sacounter = 1;
% For particles that have moved backwards only to a volume boundary and not
% completely integrated for the time step.  Move back to the next volume boundary.
% Loop back piecewise across volumes until integrate over the entire time-step.
BTemp = (abs(DT - DeltaTE) > PRECH);
LOOPINC = sum(BTemp);
while ((LOOPINC > 0) & (sacounter < (MAXCR+1)))
% Could use a DT to DeltaTE check to limit calculations.  However, in Matlab
% vector format - this check doubles the calculations.  In essence calculate
% two vectors rather than one.
% Set up time increment vectors for moving back across multiple volumes.
   TempTE = double(zeros(NUMINC,1));
% Now calculate the weights to be used.  The weight location is u-velocity
% in the positive direction and v-velocity in the positive direction.  
%  X-velocity in the positive direction is located
% at the top corner of the bilinear interpolation grid.  Determines percentage
% of the volume separating actual location from defined velocity location.
   p = (abs((XPos - XLoc)))./DX;
   q = (abs((YPos - YLoc)))./DY;
% Determine y indices for v- velocity bilinear interpolation.
   [yIndex1,yIndex2,yIndex3,yIndex4,vp] = yINDEXfIND(PCol,PRow,p,NUMINC);
% Use bilinear interpolation to get the v-velocity at the new location.  Will
%  be employed to determine which volume to calculate next.
   Vyp = ((1 - vp).*((1-q).*VVel(yIndex1) + q.*VVel(yIndex2))) + ...
      vp.*((1-q).*VVel(yIndex4) + q.*(VVel(yIndex3)));
% Determine the x indices for u - velocity bilinear interpolation.
   [xIndex1,xIndex2,xIndex3,xIndex4,uq] = xINDEXfIND(PCol,PRow,q,NUMINC);
% Calculate the average v-velocity at the new locations for particles that still
% need to be traced backwards.  Bilinear interpolation will be employed to 
% determine the next volume of note.
   Vxp = (((1 - p).*((1-uq).*UVel(xIndex1) + uq.*UVel(xIndex2))) + ...
      p.*((1-uq).*UVel(xIndex4) + uq.*(UVel(xIndex3))));
% Now have the bilinear values for velocity and the initial volume, row, column
%  locations, so adjust these locations to reflect the next volume traversed.
%  Then, get from this volume the velocity, location, and distance values needed
%  to complete the calculations for the next step.
   [X2,Y2,Vx2,Vx1,Vy2,Vy1,Adx,Ady,Ddx,Ddy] = vOLUMEfIND(XLoc,YLoc,Vxp,Vyp,UVel,VVel,...
      p,q,PRow,PCol,NUMINC);
% Now calculate the new x-direction gradient values with consideration for
% divide by zero.
   BTemp = ((Vxp - double(0.0)) > -PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*(Vx2 - Vx1)) + (Calc2.*(Vx1 - Vx2));
   Ax = TempCalc1./DX;
% Now calculate new y-direction gradient with consideration for divide by zero.
   BTemp = ((Vyp - double(0.0)) > -PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*(Vy2 - Vy1)) + (Calc2.*(Vy1 - Vy2));
   Ay = TempCalc1./DY;
% Now determine the u-direction velocity value to apply across the current 
%     volume.  Linear interpolation.
   Vxp = Vx2 - Ax.*Adx;
% Now determine the u-direction velocity value to apply across the current 
%     volume.  Linear interpolation.
   Vyp = Vy2 - Ay.*Ady;
% Now are ready for the time step calculations.  PartTX = Vxp/Vx1.
%     PartTY = Vyp/Vy1.
   BTemp = (abs(Vx1 - 0.0) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*Vx1) + (Calc2.*PREC);
   TempCalc2 = 1./TempCalc1;
   PartTX = TempCalc2.*Vxp;
   BTemp = (abs(Vy1 - 0.0) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*Vy1) + (Calc2.*PREC);
   TempCalc2 = 1./TempCalc1;
   PartTY = TempCalc2.*Vyp;
% Next calculate the time step in the x-direction.  Do the same for y.
   BTemp = (abs(Ax - double(0.0)) > PREC);
   CalcAx1 = double(BTemp);
   CalcAx2 = double(~BTemp);
   TempCalc1 = (CalcAx1.*Ax) + (CalcAx2.*PREC);
   TempCalc2 = 1./TempCalc1;
   BTemp = ((PartTX - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*PartTX) + (Calc2.*PRECH);
   DeltaTX = abs(TempCalc2.*log(TempCalc1));
   BTemp = (abs(Ay - double(0.0)) > PREC);
   CalcAy1 = double(BTemp);
   CalcAy2 = double(~BTemp);
   TempCalc1 = (CalcAy1.*Ay) + (CalcAy2.*PREC);
   TempCalc2 = 1./TempCalc1;
   BTemp = ((PartTY - double(0.0)) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*PartTY) + (Calc2.*PRECH);
   DeltaTY = abs(TempCalc2.*log(TempCalc1));
% Now adjust for the situations where GradX == 0, PartTX = 0, and PartTX == 1.
%  When PartTX == 0, then the exit velocity for the volume is zero.  In this case,
%  the particle will not leave the volume in that direction, so DeltaT = DT.
%  When GradX == 0, then the velocity gradient is zero, and the particle will
%  travel in a straight line.  When PartTX == 1, then ln(1) == 0 so the calculation
%  will be screwed up.  Generally this occurs close to a boundary.
%  Do the same for y quantities.
   BTemp = (abs(Vxp - double(0.0)) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*abs(Vxp)) + (Calc2.*PREC);
   TempCalc2 = 1./TempCalc1;
   DeltaTX = (CalcAx2.*(TempCalc2.*Ddx)) + (CalcAx1.*DeltaTX);
   BTemp = (abs(PartTX - double(0.0)) <= PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTX = (Calc1.*DT) + (Calc2.*DeltaTX);
   BTemp = (abs(PartTX - 1.0) <= PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTX = (Calc1.*(TempCalc2.*Ddx)) + (Calc2.*DeltaTX);
   BTemp = (abs(Vyp - double(0.0)) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1.*abs(Vyp)) + (Calc2.*PREC);
   TempCalc2 = 1./TempCalc1;
   BTemp = (abs(Ay - double(0.0)) <= PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTY = (Calc1.*(TempCalc2.*Ddy)) + (Calc2.*DeltaTY);
   BTemp = (abs(PartTY - double(0.0)) <= PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTY = (Calc1.*DT) + (Calc2.*DeltaTY);
   BTemp = (abs(PartTY - 1.0) <= PREC);
   BTemp = (abs(PartTY - double(0.0)) <= PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTY = (Calc1.*(TempCalc2.*Ddy)) + (Calc2.*DeltaTY);
% Now create the new time step to move back through the next cell.
   BTemp = ((DeltaTX - DeltaTY) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);  
   TempTE = (Calc2.*DeltaTX) + (Calc1.*DeltaTY);
   BTemp = (((TempTE + DeltaTE) - DT) < PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);  
   TempTE = (Calc1.*TempTE) + (Calc2.*(TempTE - ((TempTE + DeltaTE)- DT)));
   BTemp = ((DeltaTE - DT) < PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);     
   TempTE = (Calc1.*TempTE) + (Calc2.*0.0);
% The new time increment for calculations is now stored in TempTE.
% Now modify the applicable time step if the particle has reached 
%  a domain boundary.  Want the particle to stop when it hits a
%  boundary.
   BoundY2 = (((YINDEX(1) - YLoc) >= -PRECH) & ...
      ((Vyp - double(0.0)) >= -PRECH));
   BoundY4 = (((YLoc - YINDEX(NUMROWS+1)) >= -PRECH) & ...
      ((Vyp - double(0.0)) <= PRECH));
   BoundY = BoundY2 + BoundY4;
   CalcY1 = double(BoundY);
   CalcY2 = double(1 - BoundY);
   TempTE = (CalcY2.*TempTE);
   DeltaTE = (CalcY1.*DT) + (CalcY2.*DeltaTE);
   BoundX3 = (((XINDEX(1) - XLoc) >= -PRECH) & ...
      ((Vxp - double(0.0)) >= -PRECH));
   BoundX1 = (((XLoc - XINDEX(NUMCOLS+1)) >= -PRECH) & ...
      ((Vxp - double(0.0)) <= PRECH));
   BoundX = BoundX1 + BoundX3;
   CalcX1 = double(BoundX);
   CalcX2 = double(1 - BoundX);
   TempTE = (CalcX2.*TempTE);
   DeltaTE = (CalcX1.*DT) + (CalcX2.*DeltaTE);
% Now do new location calculations.
   TempCalc1 = (CalcAx1.*Ax) + (CalcAx2.*PREC);
   TempCalc2 = 1./TempCalc1;
% Set-up TempCalc1 to represent solution with a velocity gradient and
%  TempCalc2 to provide the solution in the absence of a gradient.
   TempCalc1 = (CalcAx1.*(X2 - (TempCalc2.*(Vx2 - (Vxp./exp(Ax.*TempTE))))));
   TempCalc2 =  (CalcAx2.*(XLoc - (Vxp.*TempTE)));
   XLoc = (CalcX2.*(TempCalc2 + TempCalc1)) + (CalcX1.*XLoc);
   BTemp = (abs(Ay - double(0.0)) > PREC);
   TempCalc1 = (CalcAy1.*Ay) + (CalcAy2.*PREC);
   TempCalc2 = 1./TempCalc1;
   TempCalc1 = (CalcAy1.*(Y2 - (TempCalc2.*(Vy2 - (Vyp./exp(Ay.*TempTE))))));
   TempCalc2 = (CalcAy2.*(YLoc - (Vyp.*TempTE)));
   YLoc = (CalcY2.*(TempCalc2 + TempCalc1)) + (CalcY1.*YLoc);
% Now adjust location to be within the y-face and x-face stencils.
   BTemp = ((YINDEX(1) - YLoc) < -PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   YLoc = (Calc2.*YINDEX(1)) + (Calc1.*YLoc);
   BTemp = ((YINDEX(NUMROWS+1) - YLoc) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   YLoc = (Calc2.*YINDEX(NUMROWS+1)) + (Calc1.*YLoc);
   BTemp = ((XINDEX(1) - XLoc) < -PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   XLoc = (Calc2.*XINDEX(1)) + (Calc1.*XLoc);
   BTemp = ((XINDEX(NUMCOLS+1) - XLoc) > PRECH);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   XLoc = (Calc2.*XINDEX(NUMCOLS+1)) + (Calc1.*XLoc);
% Now determine the new locations calculated in XLoc and YLoc.
   [PCol,PRow,XPos,YPos,PNode] = lOCATIONcALC(XLoc,YLoc,NUMINC,...
       PCol,PRow,XPos,YPos,PNode);
% Now calculate new value of DeltaTE to determine if will loop again.
   BTemp = ((DeltaTE - DT) < -PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   DeltaTE = (Calc1.*(DeltaTE + TempTE)) + (Calc2.*DT);
% Now adjust the time step if at domain boundaries.
   DeltaTE = (CalcY1.*DT) + (CalcY2.*DeltaTE);
   DeltaTE = (CalcX1.*DT) + (CalcX2.*DeltaTE);
% Now have the particle locations determined so loop/ check if have integrated
% all particles over the current time step.
   sacounter = sacounter + 1;
   BTemp = (abs(DT - DeltaTE) > PRECH);
   LOOPINC = sum(BTemp);
end

clear Adx Ady Vxp Vyp Vx1 Vx2 Vy1 Vy2 Ax Ay Ddx Ddy DT NUMINC h;
clear BTemp BoundX BoundY DeltaTE DeltaTX DeltaTY PartTX PartTY;
clear TempCalc1 TempCalc2 TempTE uq vp xIndex1 xIndex2 xIndex3 xIndex4;
clear yIndex1 yIndex2 yIndex3 yIndex4;
return;
%EOF