function mainfluid(cntr,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % The function fluid provides the main body of the caluclations for the 
    % depth-averaged surface flow program.
    %
    % cntr is the output counter to see if need to output
    % cETime is current elapsed time in decimal days
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

    global aXDen aYDen BH BHux BHvy Cz CzX CzY EtaNew EtaXM1 EtaXP1
    global EtaYP1 EtaYM1 fluid_t Fux Fvy gX gY H Hux Hvy h1 h2 h3 h4
    global MN MnX MnY NUMINCX NUMINCY QINBC RADFLUXBC RADORLFSBC
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
    EtaNew = fREEsURF(EtaNew,fluid_t,cntr,cETime);
    % allocate free surface values to vectors for adjacent nodes.
    [EtaXP1,EtaXM1,EtaYP1,EtaYM1] = aLLOCaDJvOL(EtaNew,EtaXP1,EtaXM1,...
       EtaYP1,EtaYM1);
    % adjust boundary values for open boundaries.
    [EtaXP1,EtaXM1,EtaYP1,EtaYM1] = fsbOUNDaDJ(EtaNew,EtaXP1,EtaXM1,...
       EtaYP1,EtaYM1,cETime);
    % Set specified free surface values for radiation boundaries.
    % Also, specify the free surface be equal to adjacent volume for
    % flux and velocity boundaries.
    if (RADORLFSBC == 1)
       [EtaXP1,EtaXM1,EtaYP1,EtaYM1] = rADoRLfSbc(EtaXP1,EtaXM1,EtaYP1,EtaYM1);
    end
    % now calculate velocities with the new free surface values and
    % adjust velocities at "point source" locations.
    [u,v] = vELcALC(u,v,cETime);
    % set new total depths.
    [Hux,Hvy,h1,h2,h3,h4,BHux,BHvy] = sETtOTALdEPTH(Hux,Hvy,h1,h2,h3,h4,...
                                      BHux,BHvy,cETime);
    [BH,H,Sides] = aVEdEPTHcALC(BH,H,Sides,h1,h2,h3,h4);
    % set specified flux values now that have the new depth.
    if (QINBC == 1)
       [u,v] = qiNfLUXbc(u,v,cETime);
    end
    if ((RADFLUXBC == 1) && (RADORLFSBC == 1))
       [u,v] = rADfLUXbc(u,v);
    end
    % get velocity parameters.
    [UAve,UDirect,VAve,VDirect] = vELpARAMsET(UAve,UDirect,VAve,VDirect,u,v);
    % return to fluid.m
    % not needed
    %clear cntr;
end
%EOF
