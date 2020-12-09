function cOMPsTATgEN
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cOMPsTATgEN generates statistics for model accuracy evaluation.  This function 
% is designed to work with the Kootenai River Reach 1 data.  This function 
% is also designed to be called by cpLOT.m.  The format of all input files
% is X   Y  Data Value.
%  
global DX DY NUMROWS NUMCOLS NUMNODES PREC

% Parameters
SIZE = 1591;                           % Number of measurement locations.
XPlotMin = 543380;
XPlotMax = 543890;
YPlotMin = 5403660;
YPlotMax = 5403990;
XMax = XPlotMin + (NUMCOLS*DX) + (0.5*DX);
XMin = XPlotMin - (0.5*DX);              % X-coordinate of lower left-hand corner.
YMax = YPlotMin + (NUMROWS*DY) + (0.5*DY);
YMin = YPlotMin - (0.5*DY);             % Y-coordinate of lower left-hand corner.
FDVel = 'Rch1_Comp_AVel.dat';   % Measured velocity file.
FDDep = 'Rch1_Comp_Depth.dat';    % Measured depth file.
FCVel = 'AVelxyz.dat'; % Calculated velocity file.
FCDep = 'Depxyz.dat';  % Calculated depth file.
FOVel = 'VErrxyz.dat';  % Velocity residuals
FODep = 'DErrxyz.dat';  % Depth residuals
FOVel1 = 'V1Errxyz.dat'; % Velocity residuals normalized by max.
FODep1 = 'D1Errxyz.dat'; % Depth residuals normalized by max.
FOVel2 = 'V2Errxyz.dat'; % Velocity residuals normalized by median.
FODep2 = 'D2Errxyz.dat'; % Depth residuals normalized by median.

% Variables.
AErVel = zeros(SIZE,1);                % Absolute velocity error.
AErDep = zeros(SIZE,1);                % Absolute error depth.
AErVel1 = zeros(SIZE,1);               % Absolute velocity error normalized by median of data.
AErDep1 = zeros(SIZE,1);               % Absolute depth error normalized by median of data.
AErVel2 = zeros(SIZE,1);               % Absolute velocity error normalized by data.
AErDep2 = zeros(SIZE,1);               % Absolute depth error normalized by data.
AErVstd = zeros(1,1);                  % STD of absolute velocity error.
AErVstd1 = zeros(1,1);                 % STD of absolute velocity error normalized by median.
AErVstd2 = zeros(1,1);                 % STD of absolute velocity error normalized by data.
AErDstd = zeros(1,1);                  % STD of absolute depth error.
AErDstd1 = zeros(1,1);                 % STD of absolute depth error normalized by median.
AErDstd2 = zeros(1,1);                 % STD of absolute depth error normalized by data.
CDep = zeros(SIZE,1);                  % Interpolated depth values.
CVel = zeros(SIZE,1);                  % Interplated velocity values.
CalDep = zeros(NUMNODES,3);            % Matrix holding calculated depth values.
CalVel = zeros(NUMNODES,3);            % Matrix holding calculated velocity values.
Col = zeros(SIZE,1);                   % Column location of each measurement.
Col1 = zeros(SIZE,1);                  % Column location of interpolation point 1.
Col2 = zeros(SIZE,1);                  % Column location of interpolation point 2.
Col3 = zeros(SIZE,1);                  % Column location of interpolation point 3.
Col4 = zeros(SIZE,1);                  % Column location of interpolation point 4.
DataDep = zeros(SIZE,3);               % Matrix holding depth data.
DataVel = zeros(SIZE,3);               % Matrix holding velocity data.
DDep = zeros(SIZE,1);                  % Measured depth values.
DVel = zeros(SIZE,1);                  % Measured velocity values.
IntP = zeros(SIZE,1);                  % Interpolation x-distance weight.
IntQ = zeros(SIZE,1);                  % Interpolation y-distance weight.
m = zeros(SIZE,1);                     % Distance from greater y-volume boundary to
                                       %     measurement location.
n = zeros(SIZE,1);                     % Distance from greater x-volume boundary to
                                       %     measurement location.
Node = zeros(SIZE,1);                  % Node location of each measurement.
Node1 = zeros(SIZE,1);                 % Node location of interpolation point 1.
Node2 = zeros(SIZE,1);                 % Node location of interpolation point 2.
Node3 = zeros(SIZE,1);                 % Node location of interpolation point 3.
Node4 = zeros(SIZE,1);                 % Node location of interpolation point 4.
p = zeros(SIZE,1);                     % Normalized distance from location to n.
q = zeros(SIZE,1);                     % Normalized distance from location to m.
Quad1 = zeros(SIZE,1);                 % Boolean for quadrant location.
Quad2 = zeros(SIZE,1);                 % Boolean for quadrant location.
Quad3 = zeros(SIZE,1);                 % Boolean for quadrant location.
Quad4 = zeros(SIZE,1);                 % Boolean for quadrant location.
Row = zeros(SIZE,1);                   % Row location of each measurement.
Row1 = zeros(SIZE,1);                  % Row location of interpolation point 1.
Row2 = zeros(SIZE,1);                  % Row location of interpolation point 2.
Row3 = zeros(SIZE,1);                  % Row location of interpolation point 3.
Row4 = zeros(SIZE,1);                  % Row location of interpolation point 4.
XLoc = zeros(SIZE,1);                  % Adjusted X-coordinate of measurement locations.
YLoc = zeros(SIZE,1);                  % Adjusted Y-coordinate of measurement locations.
XPlot = zeros(SIZE,1);                 % X-coordinate of measurement locations.
YPlot = zeros(SIZE,1);                 % Y-coordinate of measurement locations.


%  Calculations.
%     Load data or measured values and calculated or model values.
DataDep = load(FDDep,'-ASCII');
DataVel = load(FDVel,'-ASCII');
CalDep = load(FCDep,'-ASCII');
CalVel = load(FCVel,'-ASCII');
%     Find the volume location of each data location.  Depth and velocity
%        measurements were taken at same locations.
XPlot = DataDep(:,1);
XLoc = XPlot - XMin;
YPlot = DataDep(:,2);
YLoc = YPlot - YMin;
[Col,Row,n,m,Nodes] = lOCATIONcALC(XLoc,YLoc,SIZE);
%     Next determine the quadrant of the measurement location.
p = (n - XLoc)/DX;
q = (m - YLoc)/DY;
Quad1 = (((p - 0.5) <= -PREC) & ((q - 0.5) > -PREC));
Quad2 = (((p - 0.5) > -PREC) & ((q - 0.5) > -PREC));
Quad3 = (((p - 0.5) > -PREC) & ((q - 0.5) <= -PREC));
Quad4 = (((p - 0.5) <= -PREC) & ((q - 0.5) <= -PREC));
%     Employ the quadrant to determine the IntP, IntQ weights for
%        bilinear interpolation.
IntP = (Quad1.*(p + 0.5)) + ((1 - Quad1).*IntP);
IntP = (Quad4.*(p + 0.5)) + ((1 - Quad4).*IntP);
IntP = (Quad2.*(p - 0.5)) + ((1 - Quad2).*IntP);
IntP = (Quad3.*(p - 0.5)) + ((1 - Quad3).*IntP);
IntQ = (Quad4.*(q + 0.5)) + ((1 - Quad4).*IntQ);
IntQ = (Quad3.*(q + 0.5)) + ((1 - Quad3).*IntQ);
IntQ = (Quad1.*(q - 0.5)) + ((1 - Quad1).*IntQ);
IntQ = (Quad2.*(q - 0.5)) + ((1 - Quad2).*IntQ);
%     Employ the quadrant to determine the Row location for each bilinear
%        interpolation point.
Row1 = (Quad1.*(Row - 1)) + ((1 - Quad1).*Row1);
Row1 = (Quad2.*(Row - 1)) + ((1 - Quad2).*Row1);
Row1 = (Quad3.*(Row)) + ((1 - Quad3).*Row1);
Row1 = (Quad4.*(Row)) + ((1 - Quad4).*Row1);
Row1 = ((Row1 >= 1).*Row1) + ((Row1 < 1).*1);
Row1 = ((Row1 <= NUMROWS).*Row1) + ((Row1 > NUMROWS).*NUMROWS);
Row2 = (Quad1.*(Row - 1)) + ((1 - Quad1).*Row2);
Row2 = (Quad2.*(Row - 1)) + ((1 - Quad2).*Row2);
Row2 = (Quad3.*(Row)) + ((1 - Quad3).*Row2);
Row2 = (Quad4.*(Row)) + ((1 - Quad4).*Row2);
Row2 = ((Row2 >= 1).*Row2) + ((Row2 < 1).*1);
Row2 = ((Row2 <= NUMROWS).*Row2) + ((Row2 > NUMROWS).*NUMROWS);
Row3 = (Quad1.*(Row)) + ((1 - Quad1).*Row3);
Row3 = (Quad2.*(Row)) + ((1 - Quad2).*Row3);
Row3 = (Quad3.*(Row+1)) + ((1 - Quad3).*Row3);
Row3 = (Quad4.*(Row+1)) + ((1 - Quad4).*Row3);
Row3 = ((Row3 >= 1).*Row3) + ((Row3 < 1).*1);
Row3 = ((Row3 <= NUMROWS).*Row3) + ((Row3 > NUMROWS).*NUMROWS);
Row4 = (Quad1.*(Row)) + ((1 - Quad1).*Row4);
Row4 = (Quad2.*(Row)) + ((1 - Quad2).*Row4);
Row4 = (Quad3.*(Row+1)) + ((1 - Quad3).*Row4);
Row4 = (Quad4.*(Row+1)) + ((1 - Quad4).*Row4);
Row4 = ((Row4 >= 1).*Row4) + ((Row4 < 1).*1);
Row4 = ((Row4 <= NUMROWS).*Row4) + ((Row4 > NUMROWS).*NUMROWS);
%     Employ the quadrant location to determine the column location for
%        each bilinear interpolation point.
Col1 = (Quad1.*(Col + 1)) + ((1 - Quad1).*Col1);
Col1 = (Quad2.*(Col)) + ((1 - Quad2).*Col1);
Col1 = (Quad3.*(Col)) + ((1 - Quad3).*Col1);
Col1 = (Quad4.*(Col + 1)) + ((1 - Quad4).*Col1);
Col1 = ((Col1 >= 1).*Col1) + ((Col1 < 1).*1);
Col1 = ((Col1 <= NUMCOLS).*Col1) + ((Col1 > NUMCOLS).*NUMCOLS);
Col2 = (Quad1.*(Col)) + ((1 - Quad1).*Col2);
Col2 = (Quad2.*(Col - 1)) + ((1 - Quad2).*Col2);
Col2 = (Quad3.*(Col - 1)) + ((1 - Quad3).*Col2);
Col2 = (Quad4.*(Col)) + ((1 - Quad4).*Col2);
Col2 = ((Col2 >= 1).*Col2) + ((Col2 < 1).*1);
Col2 = ((Col2 <= NUMCOLS).*Col2) + ((Col2 > NUMCOLS).*NUMCOLS);
Col3 = (Quad1.*(Col)) + ((1 - Quad1).*Col3);
Col3 = (Quad2.*(Col - 1)) + ((1 - Quad2).*Col3);
Col3 = (Quad3.*(Col - 1)) + ((1 - Quad3).*Col3);
Col3 = (Quad4.*(Col)) + ((1 - Quad4).*Col3);
Col3 = ((Col3 >= 1).*Col3) + ((Col3 < 1).*1);
Col3 = ((Col3 <= NUMCOLS).*Col3) + ((Col3 > NUMCOLS).*NUMCOLS);
Col4 = (Quad1.*(Col + 1)) + ((1 - Quad1).*Col4);
Col4 = (Quad2.*(Col)) + ((1 - Quad2).*Col4);
Col4 = (Quad3.*(Col)) + ((1 - Quad3).*Col4);
Col4 = (Quad4.*(Col + 1)) + ((1 - Quad4).*Col4);
Col4 = ((Col4 >= 1).*Col4) + ((Col4 < 1).*1);
Col4 = ((Col4 <= NUMCOLS).*Col4) + ((Col4 > NUMCOLS).*NUMCOLS);
%     Now get the volume index for each of the four bilinear interpolation
%        locations.
Node1 = ((Row1 - 1).*NUMCOLS) + Col1;
Node2 = ((Row2 - 1).*NUMCOLS) + Col2;
Node3 = ((Row3 - 1).*NUMCOLS) + Col3;
Node4 = ((Row4 - 1).*NUMCOLS) + Col4;
%     Finally bilinearly interpolate calculated values to compare to the
%        measured values.
CVel = ((1 - IntP).*((IntQ.*CalVel(Node1,3)) + ((1 - IntQ).*CalVel(Node4,3)))) + ...
   (IntP.*((IntQ.*CalVel(Node2,3)) + ((1 - IntQ).*CalVel(Node3,3))));
CDep = ((1 - IntP).*((IntQ.*CalDep(Node1,3)) + ((1 - IntQ).*CalDep(Node4,3)))) + ...
   (IntP.*((IntQ.*CalDep(Node2,3)) + ((1 - IntQ).*CalDep(Node3,3))));
%     Make data values into vectors.
DVel = DataVel(:,3);
DDep = DataDep(:,3);
%     Pass values to mODsTATgEN.m for statistics generation and plotting.
[AErVel,AErVel1,AErVel2,AErDep,AErDep1,AErDep2,AErVstd,AErDstd,...
    AErVstd1,AErDstd1,AErVstd2,AErDstd2] = mODsTATgEN(CDep,CVel,DDep,...
    DVel,SIZE);
%     Finally output the data to an ascii file.
DataVel(:,3) = AErVel;
DataDep(:,3) = AErDep;
save(FOVel,'DataVel','-ASCII');
save(FODep,'DataDep','-ASCII');
DataVel(:,3) = AErVel1;
DataDep(:,3) = AErDep1;
save(FOVel1,'DataVel','-ASCII');
save(FODep1,'DataDep','-ASCII');
DataVel(:,3) = AErVel2;
DataDep(:,3) = AErDep2;
save(FOVel2,'DataVel','-ASCII');
save(FODep2,'DataDep','-ASCII');


clear SIZE XMin YMin FDVel FDDep FCVel FCDep CDep CVel CalDep CalVel Col;
clear Col1 Col2 Col3 Col4 DataDep DataVel DDep DVel IntP IntQ m;
clear NewXindex NewYindex n Node Node1 Node2 Node3 Node4 p q Quad1;
clear Quad2 Quad3 Quad4 Row Row1 Row2 Row3 Row4 XLoc YLoc XPlot YPlot;
clear AErVel AErDep AErVel1 AErDep1 AErVel2 AErDep2 FOVel FODep;
clear AErVstd AErVstd1 AErVstd2 AErDstd AErDstd1 AErDstd2;
clear FOVel1 FOVel2 FODep1 FODep2;
return;
%EOF
