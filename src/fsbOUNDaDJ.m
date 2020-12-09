function [FS1,FS3,FS4,FS2] = fsbOUNDaDJ(FS,FS1,FS3,FS4,FS2)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fsbOUNDaDJ sets the free surface gradient across open boundaries.
% The free surface gradient is originally set to enforce a constant
% slope.  All closed boundaries have already been set to a constant
% slope = 0 with the function aLLOCaDJvOL.m.
%
% FS  [NUMNODES,1] = EtaNew or free surface elevations.
% FS1 [NUMINCX,1] = EtaXP1. Adjacent to the east.
% FS2 [NUMINCY,1] = EtaYM1. Adjacent to the south.
% FS3 [NUMINCX,1] = EtaXM1. Adjacent to the west.
% FS4 [NUMINCY,1] = EtaYP1. Adjacent to the north.


global Hux Hvy PRECH LASTROW NUMCOLS NUMINCX NUMINCY NUMNODES NUMROWS
global PREC ROWBEGIN ROWEND TDEPDIRCBC ULASTROW XINC

% local variables.
Calc1A = double(zeros(NUMROWS,1)); % Calc form of Side 1.
Calc1B = double(zeros(NUMROWS,1)); % Calc form of Side 1.
Calc2A = double(zeros(NUMCOLS,1)); % Calc form of Side 2.
Calc2B = double(zeros(NUMCOLS,1)); % Calc form of Side 2.
Calc3A = double(zeros(NUMROWS,1)); % Calc form of Side 3.
Calc3B = double(zeros(NUMROWS,1)); % Calc form of Side 3.
Calc4A = double(zeros(NUMCOLS,1)); % Calc form of Side 4.
Calc4B = double(zeros(NUMCOLS,1)); % Calc form of Side 4.
CFS1 = double(zeros(NUMROWS,1));      % Constant slope eelvation east side.
CFS2 = double(zeros(NUMCOLS,1));      % Constant slope eelvation south side.
CFS3 = double(zeros(NUMROWS,1));      % Constant slope eelvation west side.
CFS4 = double(zeros(NUMCOLS,1));      % Constant slope eelvation north side.
Side1 = zeros(NUMROWS,1);             % Boolean for wet east side.
Side2 = zeros(NUMCOLS,1);             % Boolean for wet south side.
Side3 = zeros(NUMROWS,1);             % Boolean for wet west side.
Side4 = zeros(NUMCOLS,1);             % Boolean for wet north side.

% First determine which domain faces are open on each boundary.
Side1 = ((Hux(XINC:XINC:NUMINCX) - double(0.0)) > PRECH);
Calc1A = double(Side1);
Calc1B = double(~Side1);
Side2 = ((Hvy(ROWBEGIN(1):1:ROWEND(1)) - double(0.0)) > PRECH);
Calc2A = double(Side2);
Calc2B = double(~Side2);
Side3 = ((Hux(1:XINC:ULASTROW+1) - double(0.0)) > PRECH);
Calc3A = double(Side3);
Calc3B = double(~Side3);
Side4 = ((Hvy(ROWEND(NUMROWS)+1:1:NUMINCY) - double(0.0)) > PRECH);
Calc4A = double(Side4);
Calc4B = double(~Side4);
% Calculate constant slope values.
CFS1 = (2.*FS(NUMCOLS:NUMCOLS:NUMNODES)) - ...
   FS((NUMCOLS-1):NUMCOLS:(NUMNODES-1));
CFS2 = (2.*FS(ROWBEGIN(1):1:ROWEND(1))) - FS(ROWBEGIN(2):1:ROWEND(2));
CFS3 = (2.*FS(1:NUMCOLS:LASTROW+1)) - FS(2:NUMCOLS:LASTROW+2);
CFS4 = (2.*FS(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS))) - ...
   FS(ROWBEGIN(NUMROWS-1):1:ROWEND(NUMROWS-1));
if (TDEPDIRCBC == 1)
   [CFS1,CFS2,CFS3,CFS4] = tdEPsYSbc(CFS1,CFS2,CFS3,CFS4);
end
% Allocate constant slope values to adjacent volume vectors if boundary face
%  is wet.
FS1(XINC:XINC:NUMINCX) = (Calc1A.*CFS1) + ...
   (Calc1B.*FS1(XINC:XINC:NUMINCX));
FS2(ROWBEGIN(1):1:ROWEND(1)) = (Calc2A.*CFS2) + ...
   (Calc2B.*FS2(ROWBEGIN(1):1:ROWEND(1)));
FS3(1:XINC:ULASTROW+1) = (Calc3A.*CFS3) + ...
   (Calc3B.*FS3(1:XINC:ULASTROW+1));
FS4(ROWEND(NUMROWS)+1:1:NUMINCY) = (Calc4A.*CFS4) + ...
   (Calc4B.*FS4(ROWEND(NUMROWS)+1:1:NUMINCY));
%
clear Calc1A Calc1B Calc2A Calc2B Calc3A Calc3B Calc4A Calc4B;
clear CFS1 CFS2 CFS3 CFS4 FS Side1 Side2 Side3 Side4;
return;
%EOF