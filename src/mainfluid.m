function mAINfLUID(cntr)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% The function fluid provides the main body of the caluclations for the 
% depth-averaged surface flow program.

global aXDen aYDen BH BHux BHvy Cz CzX CzY EPSILON Eta EtaNew EtaXM1 EtaXP1
global EtaYP1 EtaYM1 EToEta fluid_t Fux Fvy gX gY H Hux Hvy h1 h2 h3 h4
global MN MnX MnY NUMINCX NUMINCY NUMNODES QINBC RADFLUXBC RADORLFSBC
global Sides u UAve UDirect v VAve VDirect

% Calculations.
% Calculate the Chezy coefficients.
[MnX,MnY] = vOLtOfACE(MN,MnX,MnY);
[Cz,CzX,CzY] = cHEZYcALC(MN,MnX,MnY,Cz,CzX,CzY);
%Set the values for the a* vectors.
[aXDen,aYDen] = sETafACES(aXDen,aYDen);
%Calculate the convective and viscous terms using a semi-Lagrange
% method.
Fux = fUXcALC(Fux,NUMINCX);
Fvy = fVYcALC(Fvy,NUMINCY);
%Set the values for the g* vectors.
[gX,gY] = sETgfACES(gX,gY);
% Get the new free surface elevations.
EtaNew = fREEsURF(EtaNew,fluid_t,cntr);
% allocate free surface values to vectors for adjacent nodes.
[EtaXP1,EtaXM1,EtaYP1,EtaYM1] = aLLOCaDJvOL(EtaNew,EtaXP1,EtaXM1,...
   EtaYP1,EtaYM1);
% adjust boundary values for open boundaries.
[EtaXP1,EtaXM1,EtaYP1,EtaYM1] = fsbOUNDaDJ(EtaNew,EtaXP1,EtaXM1,...
   EtaYP1,EtaYM1);
% Set specified free surface values for radiation boundaries.
% Also, specify the free surface be equal to adjacent volume for
% flux and velocity boundaries.
if (RADORLFSBC == 1)
   [EtaXP1,EtaXM1,EtaYP1,EtaYM1] = rADoRLfSbc(EtaXP1,EtaXM1,EtaYP1,EtaYM1);
end
% now calculate velocities with the new free surface values and
% adjust velocities at "point source" locations.
[u,v] = vELcALC(u,v);
% set new total depths.
[Hux,Hvy,h1,h2,h3,h4,BHux,BHvy] = sETtOTALdEPTH(Hux,Hvy,h1,h2,h3,h4,BHux,BHvy);
[BH,H,Sides] = aVEdEPTHcALC(BH,H,Sides,h1,h2,h3,h4);
% set specified flux values now that have the new depth.
if (QINBC == 1)
   [u,v] = qiNfLUXbc(u,v);
end
if ((RADFLUXBC == 1) & (RADORLFSBC == 1))
   [u,v] = rADfLUXbc(u,v);
end
% get velocity parameters.
[UAve,UDirect,VAve,VDirect] = vELpARAMsET(UAve,UDirect,VAve,VDirect,u,v);
% return to fluid.m
clear cntr;
return;
%EOF
