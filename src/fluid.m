function fider = fluid(H5FileP)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fluid.m is a depth-averaged surface flow model.  The numerics are based
    % on the depth-averaged TRIM method presented in Casulli and Cheng (1992).
    % 
    % H5FileP : HDF5 file name with path
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

    global BH BHux BHvy BMass Eta EtaNew EtaXM1 EtaXP1 EtaYM1 EtaYP1 fdt FLUID_DT
    global ENDTIME FIDHISTORY fluid_t H HNODE HUX Hux HVY Hvy
    global h1 h2 h3 h4 MN OUTINT QINBC RADORLFSBC
    global RADVELBC STARTTIME Sides u UDirect UAve v VAve VDirect
    global VELDIRCBC NUMILOCS

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
        BHux,BHvy,0.0+STARTTIME);
    % Now set initial velocities - not needed always start at 0
    %[u,v] = vELsET(u,v);
    % Now calculate the initial mass in the domain for mass balance checks.
    BMass = iNdOMAINmASS(BMass);
    % Set the mass calculation variable equal to the beginning mass.  Also,
    % set-up the mass balance output file.
    fider = fLUIDmASSsETUP;
    % calculate initial average depth for each volume.
    [BH,H,Sides] = aVEdEPTHcALC(BH,H,Sides,h1,h2,h3,h4);
    % Boundary conditions which need an initial set-up.
    if (VELDIRCBC == 1)
       [u,v] = vELdIRCbc(u,v,0.0+STARTTIME);
    end
    if (QINBC == 1)
       [u,v] = qiNfLUXbc(u,v,0.0+STARTTIME);
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
    outIntCnt = 0;
    fluid_t = (STARTTIME*(24.0 * 60.0 * 60.0));
    TotTime = ((ENDTIME - STARTTIME)*(24.0 * 60.0 * 60.0 ));
    NumSteps = ceil(TotTime/fdt);
    % set up the outputs
    sETh5Outputs(H5FileP,NumSteps);
    for t=1:1:NumSteps
       fluid_t = fluid_t + fdt;
       flt_days = fluid_t / ( 24.0 * 60.0 * 60.0 );
       Eta = EtaNew;
       mainfluid(ocounter, flt_days);
       if ((ocounter == OUTINT) || (inocntr == 0))
          [ocounter,inocntr] = fLUIDmASSwRITE(ocounter,inocntr);
          outIntCnt = outIntCnt + 1;
          rITh5oUTPUT(H5FileP,outIntCnt);
       else
          ocounter = ocounter + 1;
       end
       if (NUMILOCS > 0)
           % need to output for this time step and specified locations
           rITh5iLOCS(H5FileP,counter);
       end
       counter = counter+1;
    end
    [ocounter,inocntr] = fLUIDmASSwRITE(ocounter,inocntr);
    outIntCnt = outIntCnt + 1;
    rITh5oUTPUT(H5FileP,outIntCnt);
    % Close open files.
    fclose(FIDHISTORY);
    % probably not needed anymore
    %clear counter inocntr ocounter;
    %clear NumSteps TotTime t;
end
%EOF
