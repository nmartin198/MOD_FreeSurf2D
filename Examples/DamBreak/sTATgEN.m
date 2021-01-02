function sTATgEN
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sTATgEN generates statistics for model accuracy evaluation.  This function 
% is designed to work with the dambreak data.  This format of all input files
% is data values in [# of time increments,1]. This script is for use 
% when simulating the dambreak style flume experiment of Bellos et al.
% (1992).
%  
global ENDTIME INFO PREC

% First determine the locations of the values to be employed in the 
% comparison.

% Open Latex format output file.

% Get calculation time increments.
TimeInc = load('TC.txt','-ASCII');
TotTimes = size(TimeInc,1);
% Get calculated/simulated depths.
CalVal1 = load('DL1.txt','-ASCII');
CalVal2 = load('DL2.txt','-ASCII');
CalVal3 = load('DL3.txt','-ASCII');
CalVal4 = load('DL4.txt','-ASCII');
CalVal5 = load('DL5.txt','-ASCII');
CalVal6 = load('DL6.txt','-ASCII');
CalVal8 = load('DL8.txt','-ASCII');
% Get measured values.
DLoc1 = load('loc1.txt','-ASCII');
Len1 = size(DLoc1,1);
DLoc2 = load('loc2.txt','-ASCII');
Len2 = size(DLoc2,1);
DLoc3 = load('loc3.txt','-ASCII');
Len3 = size(DLoc3,1);
DLoc4 = load('loc4.txt','-ASCII');
Len4 = size(DLoc4,1);
DLoc5 = load('loc5.txt','-ASCII');
Len5 = size(DLoc5,1);
DLoc6 = load('loc6.txt','-ASCII');
Len6 = size(DLoc6,1);
DLoc8 = load('loc8.txt','-ASCII');
Len8 = size(DLoc8,1);
% Set-up Time-Averaged depth files for output.
TAWrite1 = zeros(Len1,2);
TAWrite2 = zeros(Len2,2);
TAWrite3 = zeros(Len3,2);
TAWrite4 = zeros(Len4,2);
TAWrite5 = zeros(Len5,2);
TAWrite6 = zeros(Len6,2);
TAWrite8 = zeros(Len8,2);
% pad the bottom of the measured values matrixes to enable for loop
%     to complete.
Ender = ceil(ENDTIME*(24*3600)) + 150.0;
DLoc1 = [DLoc1;[Ender 0.0]];
DLoc2 = [DLoc2;[Ender 0.0]];
DLoc3 = [DLoc3;[Ender 0.0]];
DLoc4 = [DLoc4;[Ender 0.0]];
DLoc5 = [DLoc5;[Ender 0.0]];
DLoc6 = [DLoc6;[Ender 0.0]];
DLoc8 = [DLoc8;[Ender 0.0]];
% Initialize Times* indexing vectors and initialize vectors to
%     hold the time averaged depth values.
Times1 = zeros(Len1,1);
TADep1 = zeros(Len1,1);
Times2 = zeros(Len2,1);
TADep2 = zeros(Len2,1);
Times3 = zeros(Len3,1);
TADep3 = zeros(Len3,1);
Times4 = zeros(Len4,1);
TADep4 = zeros(Len4,1);
Times5 = zeros(Len5,1);
TADep5 = zeros(Len5,1);
Times6 = zeros(Len6,1);
TADep6 = zeros(Len6,1);
Times8 = zeros(Len8,1);
TADep8 = zeros(Len8,1);
% Fist do time = 0.0.
Times1(1) = 1;
Times2(1) = 1;
Times3(1) = 1;
Times4(1) = 1;
Times5(1) = 1;
Times6(1) = 1;
Times8(1) = 1;
% Now do the remaining times.
LocCounter1 = 2;
SumD1 = 0.0;
TpCounter1 = 0;
LocCounter2 = 2;
SumD2 = 0.0;
TpCounter2 = 0;
LocCounter3 = 2;
SumD3 = 0.0;
TpCounter3 = 0;
LocCounter4 = 2;
SumD4 = 0.0;
TpCounter4 = 0;
LocCounter5 = 2;
SumD5 = 0.0;
TpCounter5 = 0;
LocCounter6 = 2;
SumD6 = 0.0;
TpCounter6 = 0;
LocCounter8 = 2;
SumD8 = 0.0;
TpCounter8 = 0;
% Loop through the simulation time increments to find the indexes
%  to use to compare the simulated values to measured values.
%  Also, generate time-averaged depth values.
for k=1:TotTimes
   if ((TimeInc(k) - DLoc1(LocCounter1,1)) > -PREC)
      Times1(LocCounter1) = k;
      SumD1 = SumD1 + CalVal1(k);
      TpCounter1 = TpCounter1 + 1;
      TADep1(LocCounter1) = SumD1/TpCounter1;
      LocCounter1 = LocCounter1 + 1;
      SumD1 = 0.0;
      TpCounter1 = 0;
   else
      SumD1 = SumD1 + CalVal1(k);
      TpCounter1 = TpCounter1 + 1;
   end
   if ((TimeInc(k) - DLoc2(LocCounter2,1)) > -PREC)
      Times2(LocCounter2) = k;
      SumD2 = SumD2 + CalVal2(k);
      TpCounter2 = TpCounter2 + 1;
      TADep2(LocCounter2) = SumD2/TpCounter2;
      LocCounter2 = LocCounter2 + 1;
      SumD2 = 0.0;
      TpCounter2 = 0;
   else
      SumD2 = SumD2 + CalVal2(k);
      TpCounter2 = TpCounter2 + 1;
   end
   if ((TimeInc(k) - DLoc3(LocCounter3,1)) > -PREC)
      Times3(LocCounter3) = k;
      SumD3 = SumD3 + CalVal3(k);
      TpCounter3 = TpCounter3 + 1;
      TADep3(LocCounter3) = SumD3/TpCounter3;
      LocCounter3 = LocCounter3 + 1;
      SumD3 = 0.0;
      TpCounter3 = 0;
   else
      SumD3 = SumD3 + CalVal3(k);
      TpCounter3 = TpCounter3 + 1;
   end
   if ((TimeInc(k) - DLoc4(LocCounter4,1)) > -PREC)
      Times4(LocCounter4) = k;
      SumD4 = SumD4 + CalVal4(k);
      TpCounter4 = TpCounter4 + 1;
      TADep4(LocCounter4) = SumD4/TpCounter4;
      LocCounter4 = LocCounter4 + 1;
      SumD4 = 0.0;
      TpCounter4 = 0;
   else
      SumD4 = SumD4 + CalVal4(k);
      TpCounter4 = TpCounter4 + 1;
   end
   if ((TimeInc(k) - DLoc5(LocCounter5,1)) > -PREC)
      Times5(LocCounter5) = k;
      SumD5 = SumD5 + CalVal5(k);
      TpCounter5 = TpCounter5 + 1;
      TADep5(LocCounter5) = SumD5/TpCounter5;
      LocCounter5 = LocCounter5 + 1;
      SumD5 = 0.0;
      TpCounter5 = 0;
   else
      SumD5 = SumD5 + CalVal5(k);
      TpCounter5 = TpCounter5 + 1;
   end
   if ((TimeInc(k) - DLoc6(LocCounter6,1)) > -PREC)
      Times6(LocCounter6) = k;
      SumD6 = SumD6 + CalVal6(k);
      TpCounter6 = TpCounter6 + 1;
      TADep6(LocCounter6) = SumD6/TpCounter6;
      LocCounter6 = LocCounter6 + 1;
      SumD6 = 0.0;
      TpCounter6 = 0;
   else
      SumD6 = SumD6 + CalVal6(k);
      TpCounter6 = TpCounter6 + 1;
   end
   if ((TimeInc(k) - DLoc8(LocCounter8,1)) > -PREC)
      Times8(LocCounter8) = k;
      SumD8 = SumD8 + CalVal8(k);
      TpCounter8 = TpCounter8 + 1;
      TADep8(LocCounter8) = SumD8/TpCounter8;
      LocCounter8 = LocCounter8 + 1;
      SumD8 = 0.0;
      TpCounter8 = 0;
   else
      SumD8 = SumD8 + CalVal8(k);
      TpCounter8 = TpCounter8 + 1;
   end
end
% Write time averaged depth values to files.
TAWrite1 = [TimeInc(Times1) TADep1];
TAWrite2 = [TimeInc(Times2) TADep2];
TAWrite3 = [TimeInc(Times3) TADep3];
TAWrite4 = [TimeInc(Times4) TADep4];
TAWrite5 = [TimeInc(Times5) TADep5];
TAWrite6 = [TimeInc(Times6) TADep6];
TAWrite8 = [TimeInc(Times8) TADep8];
f1 = 'TADL1.txt';
f2 = 'TADL2.txt';
f3 = 'TADL3.txt';
f4 = 'TADL4.txt';
f5 = 'TADL5.txt';
f6 = 'TADL6.txt';
f8 = 'TADL8.txt';
save(f1,'TAWrite1','-ASCII');
save(f2,'TAWrite2','-ASCII');
save(f3,'TAWrite3','-ASCII');
save(f4,'TAWrite4','-ASCII');
save(f5,'TAWrite5','-ASCII');
save(f6,'TAWrite6','-ASCII');
save(f8,'TAWrite8','-ASCII');
% Calculations.
% Generate statistics for Location 1.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len1,1);               % Calculated values for comparison times.
Diff = zeros(Len1,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len1,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len1,1);          % Calculation Variable.
TempCalc2 = zeros(Len1,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc1(1:Len1,2);
Maxer = max(MVal);
CVal = CalVal1(Times1);
Diff = CVal - MVal;
N = 1/Len1;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 1 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error:  %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep1;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 2.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len2,1);               % Calculated values for comparison times.
Diff = zeros(Len2,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len2,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len2,1);          % Calculation Variable.
TempCalc2 = zeros(Len2,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc2(1:Len2,2);
Maxer = max(MVal);
CVal = CalVal2(Times2);
Diff = CVal - MVal;
N = 1/Len2;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 2 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep2;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 3.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len3,1);               % Calculated values for comparison times.
Diff = zeros(Len3,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len3,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len3,1);          % Calculation Variable.
TempCalc2 = zeros(Len3,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc3(1:Len3,2);
Maxer = max(MVal);
CVal = CalVal3(Times3);
Diff = CVal - MVal;
N = 1/Len3;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 3 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep3;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 4.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len4,1);               % Calculated values for comparison times.
Diff = zeros(Len4,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len4,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len4,1);          % Calculation Variable.
TempCalc2 = zeros(Len4,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc4(1:Len4,2);
Maxer = max(MVal);
CVal = CalVal4(Times4);
Diff = CVal - MVal;
N = 1/Len4;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 4 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep4;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 5.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len5,1);               % Calculated values for comparison times.
Diff = zeros(Len5,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len5,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len5,1);          % Calculation Variable.
TempCalc2 = zeros(Len5,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc5(1:Len5,2);
Maxer = max(MVal);
CVal = CalVal5(Times5);
Diff = CVal - MVal;
N = 1/Len5;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 5 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep5;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 6.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len6,1);               % Calculated values for comparison times.
Diff = zeros(Len6,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len6,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len6,1);          % Calculation Variable.
TempCalc2 = zeros(Len6,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc6(1:Len6,2);
Maxer = max(MVal);
CVal = CalVal6(Times6);
Diff = CVal - MVal;
N = 1/Len6;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 6 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep6;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

% Generate statistics for Location 8.
CMean = 0.0;                        % Mean of calculated values.
Cstd = 0.0;                         % Standard deviation of calculated values.
CVal = zeros(Len8,1);               % Calculated values for comparison times.
Diff = zeros(Len8,1);               % Difference between calculated and measured.
EMean = 0.0;                        % Error mean.
Estd = 0.0;                         % Error standard deviation.
N = 0;                              % Number of comparison values.
NMEBias = 0.0;                      % Normalized bias.
NRMSE = 0.0;                        % Normalized RMSE.
MMean = 0.0;                        % Mean of measured values.
Mstd = 0.0;                         % Standard deviation of measured values.
MVal = zeros(Len8,1);               % Measured values for comparison times.
MEBias = 0.0;                       % Mean Error or Bias.
MinAE = 0.0;                        % Minimum error.
MaxAE = 0.0;                        % Maximum error.
RMSE = 0.0;                         % Root mean square error.
TempCalc1 = zeros(Len8,1);          % Calculation Variable.
TempCalc2 = zeros(Len8,1);          % Calculation Variable.
% First do actual simulated values.
MVal = DLoc8(1:Len8,2);
Maxer = max(MVal);
CVal = CalVal8(Times8);
Diff = CVal - MVal;
N = 1/Len8;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'WATER DEPTH COMPARISON for LOC. 8 \n\n');
fprintf(INFO,'SIMULATED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n',RMSE,NRMSE);
% Next do the time-averaged values.
CVal = TADep8;
Diff = CVal - MVal;
MEBias = N*(sum(Diff));
NMEBias = MEBias*(1/Maxer);
%  Next Root Mean Square Error.
TempCalc1 = Diff.^2;
TempCalc2 = sum(TempCalc1);
TempCalc1 = N*TempCalc2;
RMSE = sqrt(TempCalc1);
NRMSE = RMSE*(1/Maxer);
MinAE = min(Diff);
MaxAE = max(Diff);
% Standard deviation and mean for calculated values, measured values, and
%  absolute error.
Cstd = std(CVal);
Mstd = std(MVal);
CMean = mean(CVal);
MMean = mean(MVal);
Estd = std(Diff);
EMean = mean(Diff);
% Output the calculated quantities.
fprintf(INFO,'SIMULATED, TIME-AVERAGED and MEASURED STATS\n');
fprintf(INFO,'Max Sim:   %6.4f\t Min Sim:   %6.4f\t Mean Sim:   %6.4f\n',...
   max(CVal),min(CVal),mean(CVal));
fprintf(INFO,'Max Meas.: %6.4f\t Min Meas.: %6.4f\t Mean Meas.: %6.4f\n\n',...
   max(MVal),min(MVal),mean(MVal));
fprintf(INFO,'DATA and SIMULATION, TIME-AVERAGED MEANS and STANDARD DEVS.\n');
fprintf(INFO,'Measured mean:  %9.5f\t Measured STD:  %9.5f\n',MMean,Mstd);
fprintf(INFO,'Simulated mean: %9.5f\t Simulated STD: %9.5f\n\n',CMean,Cstd);
fprintf(INFO,'TIME-AVERAGED ERROR ANALYSIS \n');
fprintf(INFO,'Max error: %6.4f\t\t\t Min error: %6.4f\n',MaxAE,MinAE);
fprintf(INFO,'Mean error: %6.4f\t\t\t Error Standard Dev.: %6.4f\n',EMean,Estd);
fprintf(INFO,'Mean Error (BIAS): %6.4f\t\t Normalized M.E.: %6.4f\n',MEBias,NMEBias);
fprintf(INFO,'Root Mean Square Error (RMSE): %6.4f\t Normalized RMSE: %6.4f\n\n\n',RMSE,NRMSE);

clear CMean Cstd CVal Diff EMean Estd N NMEBias NRMSE MMean Mstd MVal MEBias;
clear MinAE MaxAE RMSE TempCalc1 TempCalc2;

clear TimeInc TotTimes CalVal1 CalVal2 CalVal3 CalVal4 CalVal5 CalVal6;
clear CalVal8 Len1 Len2 Len3 Len4 Len5 Len6 Len8 DLoc1 DLoc2 DLoc3 DLoc4;
clear DLoc5 DLoc6 DLoc8 Times1 TADep1 Times2 TADep2 Times3 TADep3 Times4;
clear TADep4 Times5 TADep5 Times6 TADep6 Times8 TADep8 LocCounter1 SumD1;
clear TpCounter1 LocCounter2 SumD2 TpCounter2 LocCounter3 SumD3 TpCounter3;
clear LocCounter4 SumD4 TpCounter4 LocCounter5 SumD5 TpCounter5 LocCounter6;
clear SumD6 TpCounter6 LocCounter8 SumD8 TpCounter8;
clear f1 f2 f3 f4 f5 f6 f8 TAWrite1 TAWrite2 TAWrite3 TAWrite4 TAWrite5;
clear TAWrite6 TAWrite8 Ender;

return;
%EOF
