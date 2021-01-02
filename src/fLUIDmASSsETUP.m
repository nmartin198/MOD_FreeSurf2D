function fider = fLUIDmASSsETUP
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fLUIDmASSsETUP sets up variables and files for fluid mass balance
    % calculations.  This function must be called befor sEDmASSsETUP.m to
    % provide FILENAME.
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

    global BMass FLUID_DT DX DY EMass1 ENDTIME FIDHISTORY 
    global MFlux1 OUTP STARTTIME SIMNAME TotalFlux TotalFlux1
    global TotalMass TotalMass1
    global INFO

    %global variable setup.
    EMass1 = double(0.0);           % End of time step mass in domain.
    FIDHISTORY = 0;                 % File identifier.
    fider = 0;                      % File identifier, returned
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
    fprintf(FIDHISTORY, 'Fluid v. 2.0 \n\n');
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
    fider = fopen('Info.txt','w+');
    INFO = fider;
    fprintf(INFO,'Fluid v. 2.0 \n\n');
    % Open the global output file.
    OUTP = fopen('Output.txt','w+');
    fprintf(OUTP,'Fluid v. 2.0 \n\n');
    % not needed anymore
    %clear TotTime;
end
%EOF