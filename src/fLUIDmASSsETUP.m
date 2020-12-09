function fLUIDmASSsETUP
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fLUIDmASSsETUP sets up variables and files for fluid mass balance
% calculations.  This function must be called befor sEDmASSsETUP.m to
% provide FILENAME.

global BMass FLUID_DT DX DY EMass1 ENDTIME FIDHISTORY INFO 
global MFlux1 OUTP STARTTIME THETA SIMNAME TotalFlux TotalFlux1
global TotalMass TotalMass1

%global variable setup.
EMass1 = double(0.0);           % End of time step mass in domain.
FIDHISTORY = 0;                 % File identifier.
INFO = 0;                       % File identifier.
MFlux1 = double(0.0);           % End of time step net mass flux into the domain.
OUTP = 0;                       % File identifier.
TotalFlux = double(0.0);        % Total flux into domain including current time step.
TotalFlux1 = double(0.0);       % Total flux into domain as of previous time step.
TotalMass = double(0.0);        % Total mass of water in domain including current time step.
TotalMass1 = double(0.0);       % Total mass of water in domain as of previous time step.

% local variable setup.
TotTime = double(0.0); 

% calculations.
TotTime = ENDTIME - STARTTIME;
% set file info.
FIDHISTORY = fopen('Mass.txt','w+');

% write file headers.
fprintf(FIDHISTORY, 'Fluid v. 2.7 \n\n');
fprintf(FIDHISTORY, 'Topo [m]= %15s\t Start [hrs]= %5.3f\t End [hrs]= %5.3f\n',...
   SIMNAME,STARTTIME,ENDTIME);
fprintf(FIDHISTORY,...
   'DT [s]= %5.3f\t DX [m]= %8.2f\t DY [m]= %8.2f\t',...
   FLUID_DT,DX,DY);
fprintf(FIDHISTORY,'Beginning Mass [kg]= %15.5f\t\n', BMass);
fprintf(FIDHISTORY,'\n\n');
fprintf(FIDHISTORY,'  T [hr]   MFlux [kg]   EMass [kg]  DBalance [kg]    TotalMBS [kg]   TotalMFaB [kg]    Q3 [m3/s]    Q1 [m3/s]    Q2 [m3/s]    Q4 [m3/s]\n');
fprintf(FIDHISTORY,'\n');


% Open the global information file.
INFO = fopen('Info.txt','w+');
fprintf(INFO,'Fluid v. 2.7 \n\n');
% Open the global output file.
OUTP = fopen('Output.txt','w+');
fprintf(OUTP,'Fluid v. 2.7 \n\n');

clear TotTime;
return;
%EOF