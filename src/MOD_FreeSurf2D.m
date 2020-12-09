function MOD_FreeSurf2D
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MOD_FreeSurf2D                                              %
%                                                               %
%       Nick Martin                                             %
%       Start: May 2, 2002                                      %
%       Most Recent: April 13, 2004                             %
%                                                               %
%    MOD_FreeSurf2D.m provides the "main" program functionality %
% of the collection of Matlab scripts (or *.m files) composing  %
% the program MOD_FreeSurf2D.  The "MOD" part of the name comes %
% from the intention to incorporate this fluid flow model into a%
% larger program to simulate sediment transport over short      %
% geologic time scales.  The purpose of the larger modsed pro-  %
% gram will be to simulate depositional systems creation on the %
% spatial scale of tens of meters.  This program,               %
% MOD_FreeSurf2D, will simulate surface fluid flow to provide   %
% the velocity field to a coupled sediment transport model.     %
%                                                               %
%    MOD_FreeSurf2D simulates free surface flows in rivers      %
% and shallow estuaries. The governing equations are the depth- %
% averaged free surface equations.  The numerical representation%
% of the depth-averaged equations is based on the TRIM method   %
% for 2-D from Casulli and Cheng (1992).  Th basic TRIM 2-D     %
% model is augmented by the inclusion of the THETA method of    %
% time integration which allows the user to specify explicit,   %
% semi-implicit, or fully implicit time treatment.  The         %
% integration of the THETA method with the TRIM method for 3-D  %
% is presented in Casulli (1999).                               %
%                                                               %
%    Several departures from the TRIM method for 2-D outlined   %
% in Casulli and Cheng (1992), in addition to the THETA method, %
% have been incorporated.  One of the departures, or more       %
% correctly one thing that is not mentioned by Casulli and Cheng%
% (1992), is the treatment of domain edge boundaries.  Another  %
% slight difference is numerical representation of the inertial %
% terms.                                                        %
%                                                               %
%    Two different types of domain edge boundaries are a        %
% vailable in this program.  One type is Dirichlet or constant  %
% value.  Constant values of total water depth, velocity, and   %
% flux (into the domain) can be set in the program.  The other  %
% type of boundary condition available is radiation outflow.    %
% Three different radiation outflow boundary conditions can be  %
% implemented.  One is a radiation conditions for outflow       %
% velocity where the drift velocity term is simply the upwind   %
% velocity (Yu-Heng pers.comm.).  Another radiation boundary    %
% condition is an absorbing radiation boundary condition for    %
% free surface elevation that is taken from Orlanski (1976).    %
% The final radiation boundary condition is a radiation         %
% condition for flux.  The radiation condition for flux employs %
% the Orlanski boundary for free surface elevation which allows %
% the calculation of total depth at the outflow volume of       %
% interest.  The total depth is used to adjust the outflow      %
% velocities so that the total flux equals a specified value.   %
% The radiation flux boundary condition is technically over-    %
% specified.                                                    %
%                                                               %
%    The inertial terms from the depth-averaged shallow water   %
% equations are modeled with a semi-Lagrangian (sL) method.  In %
% the sL methods, the current location of interest (in this     %
% case the center of each volume face) is traced back the       %
% characteristic to the location at the previous time step.     %
% The quantity of interest (in this case velocity) is then      %
% interpolated for the calculated position.  So, the sL methods %
% consist of a Lagrangian treatment on an Eulerian grid.  The sL%
% treatment of inertial term is part of the original TRIM       %
% formulation from Casulli and Cheng (1992), but the sL         %
% formulation in the original TRIM follows the bilinear         %
% interpolation method of Cheng, Casulli, and Milford (1984).   %
% Analysis of sL methods applied to atmospheric models by       %
% Staniforth and Cote (1991) demonstrates that a cubic          %
% interpolation method will provide better accuracy with a sL   %
% method.  So, cubic Lagrangian polynomial interpolation in 2-D %
% is employed to interpolate velocity at the calculated location%
% along the characteristic.  To move back along the             %
% characteristic, two methods are availabe for path line        %
% tracing.  One is the classical, four-step, explicit Runge-    %
% Kutta.  The other method is the semi-analytical path line     %
% tracing method of Pollock (1988).                             %
%                                                               %
% References:                                                   %
%                                                               %
% Casulli, V. and Cheng, R.T. (1992) "Semi-implicit finite      %
%     difference methods for three-dimensional shallow water    %
%     flow". International Journal for Numerical Methods in     %
%     Fluids, v. 15, pp. 629-648.                               %
%                                                               %
% Casulli, V. (1999) "A semi-implicit finite difference method  %
%     for non-hydrostatic, free-surface flows". International   %
%     Journal for Numerical Methods in Fluids, v. 30, pp. 425-  %
%     440.                                                      %
%                                                               %
% Cheng, R.T., Casulli, V., and Milford, S.N. (1984). "Eulerian-%
%     Lagrangian solution of the convection-dispersion equation %
%     in natural coordinates".  Water Resources Research, v. 20,%
%     n. 7, pp. 944-952.                                        %
%                                                               %
% Orlanski, I. (1976). "A simple boundary condition for         %
%     unbounded hyperbolic flows". Journal of Computational     %
%     Physics, v. 21, pp. 251-269.                              %
%                                                               %
% Pollock, D.W. (1988) "Semi-analytical computation of path     %
%     lines for finite difference models". Ground Water, v. 26, %
%     n. 6, pp. 743-750.                                        %
%                                                               %
% Staniforth, A. and Cote, J. (1991). "Semi-Lagrangian          %
%     integration schemes for atmospheric models - a review".   %
%     Monthly Weather Review, v. 119. pp. 2206-2223.            %
%                                                               %
% Yu-Heng, T. (2003). On the Development of a Ghost-Cell        %
%     Immersed Boundary Method and Its Application to Large     %
%     Eddy Simulation and Geophysical Fluid Dynamics. [PhD      %
%     Thesis], Stanford University.                             %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start the timer.
tic;
% global variables.
global PREC PRECH HCUTOFF
% local variables.
InfFile = 0;
TotTime = double(0.0);
% Declare global variables - should save memory.
rEADiNPUT;
% set numeric precision values.
[PREC,PRECH,HCUTOFF] = sETpRECISION(PREC,PRECH,HCUTOFF);
% DOMAIN setup.
sETdOMAIN;
% Initialize the remainder of the global variables.
fLUIDiNIT;
% set Coriolis - location specific.
fcORcALC;
% would start reading overall time and going through main time loop here.
InfFile = fluid;
% end the time counter.
TotalTime = toc;
% Output the total time to the information file.
fprintf(InfFile,'Total Elapsed time in min. for the simulation is: %20.4f\n',...
   (TotalTime/60));
fclose(InfFile);
% output profile information.
%profile report preporti
clear InfFile TotalTime;
clear global;
clear all;
return;
%EOF