function [AbsErrorVel,AbsErrorVel1,AbsErrorVel2,AbsErrorDep,AbsErrorDep1,...
    AbsErrorDep2,AEVstd,AEDstd,N1AEVstd,N1AEDstd,N2AEVstd,N2AEDstd] = ...
    mODsTATgEN(CDep,CVel,DDep,DVel,N)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% mODsTATgEN generates statistics for model accuracy evaluation.  This function 
% is designed to work with the Kootenai River Reach 1 data.  This function 
% is also designed to be called by cpLOT.m.  The format of all input files
% is X   Y  Data Value.
%  
global DX DY HPREC Hux Hvy INFO NUMROWS NUMCOLS PREC u v

% local variables.
AbsErrorVel = zeros(N,1); % Absolute error.
AbsErrorDep = zeros(N,1); % Absolute error.
AbsErrorVel1 = zeros(N,1); % Absolute error normalized by data median.
AbsErrorDep1 = zeros(N,1); % Absolute error normalized by data median.
AbsErrorVel2 = zeros(N,1); % Absolute error normalized by data.
AbsErrorDep2 = zeros(N,1); % Absolute error normalized by data.
AEVmax = 0.0;       % Velocity abs. error max.
AEVmin = 0.0;       % Velocity abs. error min.
AEVmean = 0.0;      % Velocity abs. error mean.
AEVmed = 0.0;       % Velocity abs. error median.
AEVstd = 0.0;       % Velocity abs. error std.
AEDmax = 0.0;       % Depth abs. error max.
AEDmin = 0.0;       % Depth abs. error min.
AEDmean = 0.0;      % Depth abs. error mean.
AEDmed = 0.0;       % Depth abs. error median.
AEDstd = 0.0;       % Depth abs. error std.
AEVmax = 0.0;       % Velocity abs. error max.
AEVmin = 0.0;       % Velocity abs. error min.
AEVmean = 0.0;      % Velocity abs. error mean.
AEVmed = 0.0;       % Velocity abs. error median.
AEVstd = 0.0;       % Velocity abs. error std.
N1AEDmax = 0.0;       % Depth abs. error max.
N1AEDmin = 0.0;       % Depth abs. error min.
N1AEDmean = 0.0;      % Depth abs. error mean.
N1AEDmed = 0.0;       % Depth abs. error median.
N1AEDstd = 0.0;       % Depth abs. error std.
N2AEVmax = 0.0;       % Velocity abs. error max.
N2AEVmin = 0.0;       % Velocity abs. error min.
N2AEVmean = 0.0;      % Velocity abs. error mean.
N2AEVmed = 0.0;       % Velocity abs. error median.
N2AEVstd = 0.0;       % Velocity abs. error std.
N2AEDmax = 0.0;       % Depth abs. error max.
N2AEDmin = 0.0;       % Depth abs. error min.
N2AEDmean = 0.0;      % Depth abs. error mean.
N2AEDmed = 0.0;       % Depth abs. error median.
N2AEDstd = 0.0;       % Depth abs. error std.
BTemp = 0;          % Boolean calculation variable.
BTemp1 = zeros(N,1);% Boolean calculation variable.
CDepmax = 0.0;      % calculated depth std.
CDepmin = 0.0;      % calculated depth min.
CDepmean = 0.0;     %  calculated depth mean.
CDepmed = 0.0;      % calculated depth median.
CDepstd = 0.0;      % calculated depth std.
DDepmax = 0.0;      % data depth std.
DDepmin = 0.0;      % data depth min.
DDepmean = 0.0;     %  data depth mean.
DDepmed = 0.0;      % data depth median.
DDepstd = 0.0;      % data depth std.
CVelmax = 0.0;      % calculated velocity std.
CVelmin = 0.0;      % calculated velocity min.
CVelmean = 0.0;     %  calculated velocity mean.
CVelmed = 0.0;      % calculated velocity median.
CVelstd = 0.0;      % calculated velocity std.
DVelmax = 0.0;      % data velocity std.
DVelmin = 0.0;      % data velocity min.
DVelmean = 0.0;     %  data velocity mean.
DVelmed = 0.0;      % data velocity median.
DVelstd = 0.0;      % data velocity std.
Denom = 0.0;
Denom1 = zeros(N,1);
Huxmax = 0.0;       % Depth abs. error max.
Huxmin = 0.0;       % Depth abs. error min.
Huxmean = 0.0;      % Depth abs. error mean.
Huxmed = 0.0;       % Depth abs. error median.
Huxstd = 0.0;       % Depth abs. error std.
Hvymax = 0.0;       % Depth abs. error max.
Hvymin = 0.0;       % Depth abs. error min.
Hvymean = 0.0;      % Depth abs. error mean.
Hvymed = 0.0;       % Depth abs. error median.
Hvystd = 0.0;       % Depth abs. error std.
MEDep = 0.0;                % Mean error depth.
MEVel = 0.0;                % Mean error velocity.
NMEDep = 0.0;               % Normalized Mean error depth.
NMEVel = 0.0;               % Normalized Mean error velocity.
NRMSEDep = 0.0;             % Normalized Root mean square error depth.
NRMSEVel = 0.0;             % Normalized Root mean square error velocity.
N1MEDep = 0.0;               % Normalized Mean error depth.
N1MEVel = 0.0;               % Normalized Mean error velocity.
N1RMSEDep = 0.0;             % Normalized Root mean square error depth.
N1RMSEVel = 0.0;             % Normalized Root mean square error velocity.
N2MEDep = 0.0;               % Normalized Mean error depth.
N2MEVel = 0.0;               % Normalized Mean error velocity.
N2RMSEDep = 0.0;             % Normalized Root mean square error depth.
N2RMSEVel = 0.0;             % Normalized Root mean square error velocity.
RMSEDep = 0.0;              % Root mean square error depth.
RMSEVel = 0.0;              % Root mean square error velocity.
SqErrorVel = zeros(N,1);  % Square of absolute error.
SqErrorDep = zeros(N,1);  % Square of absolute error.
SqErrorVel1 = zeros(N,1);  % Square of absolute error.
SqErrorDep1 = zeros(N,1);  % Square of absolute error.
SqErrorVel2 = zeros(N,1);  % Square of absolute error.
SqErrorDep2 = zeros(N,1);  % Square of absolute error.
Temp = 0.0;
Temp1 = zeros(N,1);
Umax = 0.0;       % Velocity abs. error max.
Umin = 0.0;       % Velocity abs. error min.
Umean = 0.0;      % Velocity abs. error mean.
Umed = 0.0;       % Velocity abs. error median.
Ustd = 0.0;       % Velocity abs. error std.
Vmax = 0.0;       % Depth abs. error max.
Vmin = 0.0;       % Depth abs. error min.
Vmean = 0.0;      % Depth abs. error mean.
Vmed = 0.0;       % Depth abs. error median.
Vstd = 0.0;       % Depth abs. error std.

% Calculations.
% Calculate errors.
AbsErrorVel = CVel - DVel;
AbsErrorDep = CDep - DDep;
SqErrorVel = AbsErrorVel.^2;
SqErrorDep = AbsErrorDep.^2;
PlotXVel = AbsErrorVel;
PlotXDep = AbsErrorDep;
% Now get max,min,mean,median, and STD for all values of interest.
AEVmax = max(AbsErrorVel);
AEVmin = min(AbsErrorVel);
AEVmean = mean(AbsErrorVel);
AEVmed = median(AbsErrorVel);
AEVstd = std(AbsErrorVel);
AEDmax = max(AbsErrorDep);
AEDmin = min(AbsErrorDep);
AEDmean = mean(AbsErrorDep);
AEDmed = median(AbsErrorDep);
AEDstd = std(AbsErrorDep);
Umax = max(u);
Umin = min(u);
Umean = mean(u);
Umed = median(u);
Ustd = std(u);
Vmax = max(v);
Vmin = min(v);
Vmean = mean(v);
Vmed = median(v);
Vstd = std(v);
Huxmax = max(Hux);
Huxmin = min(Hux);
Huxmean = mean(Hux);
Huxmed = median(Hux);
Huxstd = std(Hux);
Hvymax = max(Hvy);
Hvymin = min(Hvy);
Hvymean = mean(Hvy);
Hvymed = median(Hvy);
Hvystd = std(Hvy);
CVelmax = max(CVel);
CVelmin = min(CVel);
CVelmean = mean(CVel);
CVelmed = median(CVel);
CVelstd = std(CVel);
DVelmax = max(DVel);
DVelmin = min(DVel);
DVelmean = mean(DVel);
DVelmed = median(DVel);
DVelstd = std(DVel);
CDepmax = max(CDep);
CDepmin = min(CDep);
CDepmean = mean(CDep);
CDepmed = median(CDep);
CDepstd = std(CDep);
DDepmax = max(DDep);
DDepmin = min(DDep);
DDepmean = mean(DDep);
DDepmed = median(DDep);
DDepstd = std(DDep);
% Finally calculate RMSE and Mean Error.  Also calculate Normalized RMSE and
%   Normalized Mean Error.  Where the median data value of either depth or
%   velocity is employed for normalization.
RMSEVel = sqrt((1/N)*sum(SqErrorVel));
RMSEDep = sqrt((1/N)*sum(SqErrorDep));
MEVel = (1/N)*sum(AbsErrorVel);
MEDep = (1/N)*sum(AbsErrorDep);
% Lastly get different normalized versions of absolute error.
%  1st do normalized by data median.
BTemp = (abs(DVelmed - double(0.0)) > PREC);
Temp = (BTemp*DVelmed) + (1 - BTemp);
Denom = (BTemp*(1/Temp)) + (1 - BTemp);
AbsErrorVel1 = Denom.*AbsErrorVel;
BTemp = (abs(DDepmed - double(0.0)) > PREC);
Temp = (BTemp*DDepmed) + (1 - BTemp);
Denom = (BTemp*(1/Temp)) + (1 - BTemp);
AbsErrorDep1 = Denom.*AbsErrorDep;
SqErrorVel1 = AbsErrorVel1.^2;
SqErrorDep1 = AbsErrorDep1.^2;
N1RMSEVel = sqrt((1/N)*sum(SqErrorVel1));
N1RMSEDep = sqrt((1/N)*sum(SqErrorDep1));
N1MEVel = (1/N)*sum(AbsErrorVel1);
N1MEDep = (1/N)*sum(AbsErrorDep1);
N1AEVmax = max(AbsErrorVel1);
N1AEVmin = min(AbsErrorVel1);
N1AEVmean = mean(AbsErrorVel1);
N1AEVmed = median(AbsErrorVel1);
N1AEVstd = std(AbsErrorVel1);
N1AEDmax = max(AbsErrorDep1);
N1AEDmin = min(AbsErrorDep1);
N1AEDmean = mean(AbsErrorDep1);
N1AEDmed = median(AbsErrorDep1);
N1AEDstd = std(AbsErrorDep1);
% Next do normalized by data value.
BTemp1 = (abs(DVel - double(0.0)) > PREC);
Temp1 = (BTemp1.*DVel) + (1 - BTemp1);
Denom1 = (BTemp1.*(1./Temp1)) + (1 - BTemp1);
AbsErrorVel2 = Denom1.*AbsErrorVel;
BTemp1 = (abs(DDep - double(0.0)) > PREC);
Temp1 = (BTemp1.*DDep) + (1 - BTemp1);
Denom1 = (BTemp1.*(1./Temp1)) + (1 - BTemp1);
AbsErrorDep2 = Denom1.*AbsErrorDep;
SqErrorVel2 = AbsErrorVel2.^2;
SqErrorDep2 = AbsErrorDep2.^2;
N2RMSEVel = sqrt((1/N)*sum(SqErrorVel2));
N2RMSEDep = sqrt((1/N)*sum(SqErrorDep2));
N2MEVel = (1/N)*sum(AbsErrorVel2);
N2MEDep = (1/N)*sum(AbsErrorDep2);
N2AEVmax = max(AbsErrorVel2);
N2AEVmin = min(AbsErrorVel2);
N2AEVmean = mean(AbsErrorVel2);
N2AEVmed = median(AbsErrorVel2);
N2AEVstd = std(AbsErrorVel2);
N2AEDmax = max(AbsErrorDep2);
N2AEDmin = min(AbsErrorDep2);
N2AEDmean = mean(AbsErrorDep2);
N2AEDmed = median(AbsErrorDep2);
N2AEDstd = std(AbsErrorDep2);

% Output the calculated quantities.
fprintf(INFO,'N: %9.1f \n',N);
fprintf(INFO,'Velocity Comparisons ====================== \n\n');
fprintf(INFO,'Statistics for U and V values \n\n');
fprintf(INFO,'Max U: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
   Umax,Umin,Umean,Umed,Ustd);
fprintf(INFO,'Max V: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n\n',...
   Vmax,Vmin,Vmean,Vmed,Vstd);
fprintf(INFO,'Means and standard deviations of simulated and measured velocity \n');
fprintf(INFO,'magnitude values \n\n');
fprintf(INFO,'Measured max:  %6.3f min: %6.3f mean: %6.3f med: %6.3f STD: %6.3f \n',...
    DVelmax,DVelmin,DVelmean,DVelmed,DVelstd);
fprintf(INFO,'Simulated max: %6.3f min: %6.3f mean: %6.3f med: %6.3f STD: %6.3f \n\n',...
    CVelmax,CVelmin,CVelmean,CVelmed,CVelstd);
fprintf(INFO,'Error statistics \n\n');
fprintf(INFO,'Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    AEVmax,AEVmin,AEVmean,AEVmed,AEVstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    RMSEVel,MEVel);
fprintf(INFO,'Normalized by median data value Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    N1AEVmax,N1AEVmin,N1AEVmean,N1AEVmed,N1AEVstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    N1RMSEVel,N1MEVel);
fprintf(INFO,'Normalized by data values Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    N2AEVmax,N2AEVmin,N2AEVmean,N2AEVmed,N2AEVstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    N2RMSEVel,N2MEVel);
fprintf(INFO,'Depth Comparisons ======================== \n\n');
fprintf(INFO,'Statistics for Hux and Hvy values \n\n');
fprintf(INFO,'Hux Max: %5.2f Min: %5.2f Mean: %5.2f Med: %5.2f STD: %5.2f \n',...
   Huxmax,Huxmin,Huxmean,Huxmed,Huxstd);
fprintf(INFO,'Hvy Max: %5.2f Min: %5.2f Mean: %5.2f Med: %5.2f STD: %5.2f \n\n',...
   Hvymax,Hvymin,Hvymean,Hvymed,Hvystd);
fprintf(INFO,'Means and standard deviations of simulated and measured volume \n');
fprintf(INFO,'\t center depth values \n\n');
fprintf(INFO,'Measured max: %6.3f min: %6.3f mean: %6.3f med: %6.3f STD: %6.3f \n',...
    DDepmax,DDepmin,DDepmean,DDepmed,DDepstd);
fprintf(INFO,'Simulated max: %6.3f min: %6.3f mean: %6.3f med.: %6.3f STD: %6.3f \n\n',...
    CDepmax,CDepmin,CDepmean,CDepmed,CDepstd);
fprintf(INFO,'Error statistics \n\n');
fprintf(INFO,'Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    AEDmax,AEDmin,AEDmean,AEDmed,AEDstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    RMSEDep,MEDep);
fprintf(INFO,'Normalized by median data value Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    N1AEDmax,N1AEDmin,N1AEDmean,N1AEDmed,N1AEDstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    N1RMSEDep,N1MEDep);
fprintf(INFO,'Normalized by data values Absolute Error Stats \n');
fprintf(INFO,'Max: %6.3f Min: %6.3f Mean: %6.3f Med: %6.3f STD: %6.3f \n',...
    N2AEDmax,N2AEDmin,N2AEDmean,N2AEDmed,N2AEDstd);
fprintf(INFO,'RMSE:       %9.5f \t Mean Error:       %9.5f \n\n',...
    N2RMSEDep,N2MEDep);
fprintf(INFO,'=========================================== \n\n');


clear AEVmax AEVmin AEVmean AEVmed AEDmax;
clear AEDmin AEDmean AEDmed BTemp CDep;
clear CDepmax CDepmin CDepmean CDepmed CDepstd DDep DDepmax DDepmin DDepmean;
clear DDepmed DDepstd CVel CVelmax CVelmin CVelmean CVelmed CVelstd DVel;
clear DVelmax DVelmin DVelmean DVelmed DVelstd Huxmax Huxmin Huxmean Huxmed;
clear Huxstd Hvymax Hvymin Hvymean Hvymed Hvystd;
clear MEDep MEVel N NMEDep NMEVel NRMSEDep NRMSEVel RMSEDep RMSEVel;
clear SqErrorVel SqErrorDep Umax Umin Umean Umed Ustd Vmax Vmin
clear Vmean Vmed Vstd;
clear N1AEVmax N1AEVmin N1AEVmean N1AEVmed N1AEDmax N1AEDmin;
clear N1AEDmean N1AEDmed SqErrorVel1 SqErrorDep1 N1RMSEDep;
clear N1RMSEVel N1MEVel N1MEDep;
clear N2AEVmax N2AEVmin N2AEVmean N2AEVmed N2AEDmax N2AEDmin;
clear N2AEDmean N2AEDmed SqErrorVel2 SqErrorDep2 N2RMSEDep;
clear N2RMSEVel N2MEVel N2MEDep;

return;
%EOF