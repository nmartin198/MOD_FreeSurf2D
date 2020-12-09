function [cntr,initcntr] = fLUIDmASSwRITE(cntr,initcntr)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fLUIDmASSwRITE writes mass balance information to a file for analysis.

global BMass EMass1 FIDHISTORY fluid_t MFlux1 TotalFlux
global TotalFlux1 TotalMass TotalMass1
global Hux Hvy XINC ULASTROW DY DX u v NUMINCX NUMINCY NUMCOLS NUMROWS
global ROWBEGIN ROWEND

% local variables.
DeltaFlux = double(0.0);     % Change in flux.
DeltaMass = double(0.0);     % Change in mass of water.
EMass = double(0.0);         % Current total mass.
MFlux = double(0.0);         % Current mass flux.
MassBalance = double(0.0);   % mass balance.
TempMass = double(0.0);

% mass in domain calculations.
EMass = iNdOMAINmASS(EMass);
DeltaMass = EMass - EMass1;
TotalMass = TotalMass + DeltaMass;
EMass1 = EMass;
TotalMass1 = TotalMass;
% flux calculations.
if (initcntr == 0)
   TempMass = mASSfLUXcALC(TempMass);
   MFlux = TempMass + BMass;
   DeltaFlux = MFlux - MFlux1;
   TotalFlux = TotalFlux + MFlux;
   MFlux1 = MFlux;
   TotalFlux1 = TotalFlux;
   initcntr = 1;
else
   MFlux = mASSfLUXcALC(MFlux);
   DeltaFlux = MFlux - MFlux1;
   TotalFlux = TotalFlux + MFlux;
   MFlux1 = MFlux;
   TotalFlux1 = TotalFlux;
   cntr = 1;
end
MassBalance = MFlux - DeltaMass;
% special Q in and Q out calculations.
Q3 = (Hux(1:XINC:ULASTROW+1)'*u(1:XINC:ULASTROW+1))*DY;
Q1 = (Hux(NUMCOLS+1:XINC:NUMINCX)'*u(NUMCOLS+1:XINC:NUMINCX))*DY;
Q2 = (Hvy(ROWBEGIN(1):1:ROWEND(1))'*v(ROWBEGIN(1):1:ROWEND(1)))*DX;
Q4 = (Hvy(ROWEND(NUMROWS)+1:1:NUMINCY)'*v(ROWEND(NUMROWS)+1:1:NUMINCY))*DX;
% print values.
fprintf(FIDHISTORY,'%8.4f %12.5E %12.5E %14.5E %15.10E %15.10E %12.5E %12.5E %12.5E %12.5E\n',...
   (fluid_t/3600),MFlux,DeltaMass,MassBalance,TotalMass,TotalFlux,Q3,Q1,Q2,Q4);

clear DeltaFlux DeltaMass MassBalance Q3 Q1 Q2 Q4 TempMass;
return;
%EOF