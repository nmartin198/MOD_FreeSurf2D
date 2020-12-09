function FMass = mASSfLUXcALC(FMass)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% mASSfLUXcALC calculates the mass leaving the domain during a time step.

global DX DY fdt Hux Hvy NUMCOLS NUMINCX NUMINCY NUMNODES NUMROWS RHOW
global RINC ROWBEGIN ROWEND u ULASTROW v XINC

% local variables.
Temp1 = double(zeros(NUMROWS,1));
Temp2 = double(zeros(NUMCOLS,1));
Temp3 = double(zeros(NUMROWS,1));
Temp4 = double(zeros(NUMCOLS,1));
Tempu = double(zeros(NUMROWS,1));
Tempuh = double(zeros(1,NUMROWS));
Tempv = double(zeros(NUMCOLS,1));
Tempvh = double(zeros(1,NUMCOLS));

% calcs.
% u3 side.
Temp3 = u(1:XINC:ULASTROW+1);
Tempu = fdt.*Temp3;
Temp3 = Hux(1:XINC:ULASTROW+1)';
Tempuh = DY.*Temp3;
FMass = FMass + ((Tempuh*Tempu)*RHOW);

%u1 side.
Temp1 = u(NUMCOLS+1:XINC:NUMINCX);
Tempu = fdt.*Temp1;
Temp1 = Hux(NUMCOLS+1:XINC:NUMINCX)';
Tempuh = DY.*Temp1;
FMass = FMass - ((Tempuh*Tempu)*RHOW);

% v2 side.
Temp2 = v(ROWBEGIN(1):1:ROWEND(1));
Tempv = fdt.*Temp2;
Temp2 = Hvy(ROWBEGIN(1):1:ROWEND(1))';
Tempvh = DX.*Temp2;
FMass = FMass + ((Tempvh*Tempv)*RHOW);

% v4 side.
Temp4 = v(ROWEND(NUMROWS)+1:1:NUMINCY);
Tempv = fdt.*Temp4;
Temp4 = Hvy(ROWEND(NUMROWS)+1:1:NUMINCY)';
Tempvh = DX.*Temp4;
FMass = FMass - ((Tempvh*Tempv)*RHOW);

clear Tempuh Tempu Tempvh Tempv;
clear Temp1 Temp2 Temp3 Temp4;
return;
%EOF