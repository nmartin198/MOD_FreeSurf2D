function INFO = fluid
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% fluid.m is a depth-averaged surface flow model.  The numerics are based
% on the depth-averaged TRIM method presented in Casulli and Cheng (1992).

global BH BHux BHvy BMass Eta EtaNew EtaXM1 EtaXP1 EtaYM1 EtaYP1 fdt FLUID_DT
global ENDTIME FIDHISTORY FILENAME fluid_t H HNODE HUX Hux HVY Hvy
global h1 h2 h3 h4 INFO MN NUMNODES OUTINT PREC QINBC RADFLUXBC RADORLFSBC
global RADVELBC STARTTIME Sides TIMESTEP u UDirect UAve v VAve VDirect
global VELDIRCBC

% Set time increment to start of zero and initialize local.
fluid_t = double(0.0);
TotTime = double(0.0);
NumSteps = 0;
% Set the undisturbed depths, free surface heights, and initial mass
% balance values.  Also set the values for source locations.  Finally
% create the vectors corresponding to free surface depths at adjacent
% nodes.
[EtaNew,HNODE,HUX,HVY] = sETiNITIALvECTORS(EtaNew,HNODE,HUX,HVY);
% Now allocate the free surface values to the vectors representing adjacent cell
% free surface values.
[EtaXP1,EtaXM1,EtaYP1,EtaYM1] = aLLOCaDJvOL(EtaNew,EtaXP1,EtaXM1,...
   EtaYP1,EtaYM1);
% Now set the total depth vectors, HUX and HVY, and the corresponding
% total depth vectors for each face of each volume, h1,h2,h3,h4.
[Hux,Hvy,h1,h2,h3,h4,BHux,BHvy] = sETtOTALdEPTH(Hux,Hvy,h1,h2,h3,h4,...
    BHux,BHvy);
% Now set initial velocities.
[u,v] = vELsET(u,v);
% Now calculate the initial mass in the domain for mass balance checks.
BMass = iNdOMAINmASS(BMass);
% Set the mass calculation variable equal to the beginning mass.  Also,
% set-up the mass balance output file.
fLUIDmASSsETUP;
% calculate initial average depth for each volume.
[BH,H,Sides] = aVEdEPTHcALC(BH,H,Sides,h1,h2,h3,h4);
% Boundary conditions which need an initial set-up.
if (VELDIRCBC == 1);
   [u,v] = vELdIRCbc(u,v);
end
if (QINBC == 1)
   [u,v] = qiNfLUXbc(u,v);
end
if (RADORLFSBC == 1)
   rADoRLsET;
end
if (RADVELBC == 1)
   rADvELsET;
end
[UAve,UDirect,VAve,VDirect] = vELpARAMsET(UAve,UDirect,VAve,VDirect,u,v);
% get value for Manning's n.
MN = tOTmANNINGcALC(MN);
% Setup and initialization is done.
% Start Main Loop = Time Loop.
fdt = FLUID_DT;
counter = 1;
ocounter = 1;
inocntr = 0;
fluid_t = (STARTTIME*(60*60));
TotTime = ((ENDTIME - STARTTIME)*(60*60));
NumSteps = ceil(TotTime/fdt);
% special dambreak lines %%%%%%%%%%%%%%%%%%%%%%%%%
dbinit;
% end special lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:1:NumSteps
   fluid_t = fluid_t + fdt;   
   Eta = EtaNew;
   mainfluid(ocounter);
   if ((ocounter == OUTINT) | (inocntr == 0))
      [ocounter,inocntr] = fLUIDmASSwRITE(ocounter,inocntr);
   else
      ocounter = ocounter + 1;
   end
   counter = counter+1;
% special dambreak lines %%%%%%%%%%%%%%%%%%%%%%%%%
   dbset(counter,fluid_t);
% end special lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
[ocounter,inocntr] = fLUIDmASSwRITE(ocounter,inocntr);
% Close open files.
fclose(FIDHISTORY);
% Function to save velocities (U,V,W) and total water depth (Hux,Hvy) and 
% topography (ZTOPO) at the end - potentially once "steady state" has been
% achieved.
rIToUTPUT;
% special dambreak lines %%%%%%%%%%%%%%%%%%%%%%%%%%
dbwrite;
sTATgEN;
aFTERpLOT;
% end special lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear counter inocntr ocounter;
clear NumSteps TotTime t;
return;
%EOF
