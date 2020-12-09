function [Side1,Side2,Side3,Side4,Xplus1,Yminus1,Xminus1,Yplus1,UXplus1,...
      UXminus1,UVplus1,UVminus1,VYplus1,VYminus1,VVplus1,VVminus1] = gENfACEiNDEX(Side1,...
      Side2,Side3,Side4,Xplus1,Yminus1,Xminus1,Yplus1,UXplus1,UXminus1,UVplus1,...
      UVminus1,VYplus1,VYminus1,VVplus1,VVminus1);
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% gENfACEiNDEX generates vectors of indexes that correspond to specific
% stencil location.  For a description of the stencil locations see
% sETdOMAIN.m
%

global LASTROW NUMCOLS NUMINCX NUMINCY NUMNODES NUMROWS RINC ROWBEGIN ROWEND XINC

% Calculations. Loop through domain by row.
for i = 1:NUMROWS
   Side1(ROWBEGIN(i):1:ROWEND(i),1) = ((((i-1)*XINC)+2):1:(i*XINC))';
   Side2(ROWBEGIN(i):1:ROWEND(i)) = (ROWBEGIN(i):1:ROWEND(i))';
   Side3(ROWBEGIN(i):1:ROWEND(i)) = ((((i-1)*XINC)+1):1:(((i-1)*XINC)+NUMCOLS))';
   Side4(ROWBEGIN(i):1:ROWEND(i)) = (ROWBEGIN(i)+NUMCOLS:1:ROWEND(i)+NUMCOLS)';
   Xplus1(ROWBEGIN(i):1:ROWEND(i)-1,1) = (ROWBEGIN(i)+1:1:ROWEND(i))';
   if (i > 1)
      Yminus1(ROWBEGIN(i):1:ROWEND(i),1) = (ROWBEGIN(i-1):1:ROWEND(i-1))';
   end
   Xminus1(ROWBEGIN(i)+1:1:ROWEND(i),1) = (ROWBEGIN(i):1:ROWEND(i)-1)';
   if (i < NUMROWS)
      Yplus1(ROWBEGIN(i):1:ROWEND(i),1) = (ROWBEGIN(i+1):1:ROWEND(i+1))';
   end
   UXplus1((((i-1)*XINC)+1):1:((i*XINC)-1),1) = (ROWBEGIN(i):1:ROWEND(i))';
   UXplus1((i*XINC),1) = ROWEND(i);
   UXminus1((((i-1)*XINC)+1),1) = ROWBEGIN(i);
   UXminus1((((i-1)*XINC)+2):1:(i*XINC),1) = (ROWBEGIN(i):1:ROWEND(i))';
   UVplus1((((i-1)*XINC)+1):1:((i*XINC)-1),1) = ((((i-1)*XINC)+2):1:((i*XINC)))';
   UVplus1((i*XINC),1) = (i*XINC);
   UVminus1((((i-1)*XINC)+1),1) = (((i-1)*XINC)+1);
   UVminus1((((i-1)*XINC)+2):1:(i*XINC),1) = ((((i-1)*XINC)+1):1:((i*XINC)-1))';
end
Xplus1(NUMCOLS:RINC:NUMNODES,1) = (NUMCOLS:RINC:NUMNODES)';
Yminus1(ROWBEGIN(1):1:ROWEND(1),1) = (ROWBEGIN(1):1:ROWEND(1))';
Xminus1(1:RINC:LASTROW+1,1) = (1:RINC:LASTROW+1)';
Yplus1(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS),1) = (ROWBEGIN(NUMROWS):1:...
   ROWEND(NUMROWS))';
% finally do the v velocity indexes.
VYminus1 = [Yminus1;Yplus1(ROWBEGIN(NUMROWS-1):1:ROWEND(NUMROWS-1))];
VYplus1 = [Yminus1(ROWBEGIN(2):1:ROWEND(2));Yplus1];
VVplus1 = [Yplus1;(ROWEND(NUMROWS)+1:1:NUMINCY)'];
VVplus1(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS),1) = (ROWEND(NUMROWS)+1:1:NUMINCY)';
VVminus1 = [Yminus1;(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS))'];
return;
%EOF