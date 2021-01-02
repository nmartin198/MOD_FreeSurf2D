function fLUIDiNIT
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fLUIDiNIT initializes the global variables needed for fluid flow
    % calculations.
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

    global aXDen aYDen BH BHux BHvy Cz CzX CzY Eta EtaNew EtaXP1 EtaXM1 EtaYM1
    global EtaYP1 Fux Fvy gX gY H HOld h1 h2 h3 h4 HUX HVY HNODE Hux HuxOld Hvy
    global HvyOld MnX MnY NUMINCX NUMINCY NUMNODES qRHS Sides
    global sLHS u UAve UDirect UOld v VAve VDirect VOld

    % Free surface variables for the domain.  Eta is the free surface elevation for
    % the old time step (N).  EtaNew is the free surface elevation for the new time
    % step (N+1).  Eta1 is the free surface elevation for the node adjacent to 
    % the current node but at (i,j+1).  Eta2 is for (i-1,j).  Eta3 is for (i,j-1).
    % and Eta4 is for (i+1,j);  EtaXP1 and EtaXM1 are for use when determining
    % x-face values.  EtaXP1 is the free surface on the i+1 side of the face and
    % EtaMP1 is the free surface on the i side of the face.
    % All values are in meters [m].
    Eta = double(zeros(NUMNODES,1));
    EtaNew = double(zeros(NUMNODES,1));
    EtaXP1 = double(zeros(NUMINCX,1));
    EtaXM1 = double(zeros(NUMINCX,1));
    EtaYM1 = double(zeros(NUMINCY,1));
    EtaYP1 = double(zeros(NUMINCY,1));

    % Allocate memory for the velocity vectors.  The u vector is for
    % velocities in the x-direction, v vector is for velocities in the Y-direction.
    % Velocity values are in meters/second [m/s].

    u = double(zeros(NUMINCX,1));    % x-direction velocity defined at x-faces.
    UAve = double(zeros(NUMNODES,1));% volume average x-component velocity.
    UDirect = zeros(NUMINCX,1);      % u velocity direction (+1,-1);
    UOld = double(zeros(NUMINCX,1)); % x-direction velocity for previous time step.
    v = double(zeros(NUMINCY,1));    % y-direction velocity defined at y-faces.
    VAve = double(zeros(NUMNODES,1));% volume average y-component velocity.
    VDirect = zeros(NUMINCY,1);      % u velocity direction (+1,-1);
    VOld = double(zeros(NUMINCY,1)); % y-direction velocity at previous time step.

    % Convective and viscous term vectors.  Will calculate these values by integrating
    % back the particle path line.
    Fux = double(zeros(NUMINCX,1));    % Fu term from Casulli and Cheng (1992).
    Fvy = double(zeros(NUMINCY,1));    % Fv term from Casulli and Cheng (1992).

    % Now vectors for the undisturbed water height at every place where velocity is 
    % defined in 2-D planview. HUX correspondes to U locations and HVY corresponds
    % to V locations.

    HUX = double(zeros(NUMINCX,1));    % undisturbed water depth at x-faces.
    HVY = double(zeros(NUMINCY,1));    % undisturbed water depth at y-faces.
    HNODE = double(zeros(NUMNODES,1)); % undisturbed depth at volume center.

    % Now vectors for total water depth (HUX or HVY + Eta) at every place where
    % velocity is defined in 2-D planview.  Hux corresponds to U locations and 
    % Hvy corresponds to V locations.

    BH = double(zeros(NUMNODES,1));   % wet volume.
    BHux = double(zeros(NUMINCX,1));   % wet face.
    BHvy = double(zeros(NUMINCY,1));   % wet face.
    Hux = double(zeros(NUMINCX,1));    % total water depth at x-faces.
    Hvy = double(zeros(NUMINCY,1));    % total water depth at y-faces.
    HuxOld = double(zeros(NUMINCX,1)); % total water depth at x-faces previous time.
    HvyOld = double(zeros(NUMINCY,1)); % total water depth at y-faces previous time.
    H = double(zeros(NUMNODES,1));     % average total water depth - at volume center.
    HOld = double(zeros(NUMNODES,1));  % average total water depth previous time.
    h1 = double(zeros(NUMNODES,1));    % water depth at i+1/2,j.
    h2 = double(zeros(NUMNODES,1));    % water depth at i,j-1/2.
    h3 = double(zeros(NUMNODES,1));    % water depth at i-1/2,j.
    h4 = double(zeros(NUMNODES,1));    % water depth at i,j+1/2.
    Sides = zeros(NUMNODES,1);         % Number of wet sides for each volume.

    % Now allocate for calculation vectors.
    % a* corresponds to A from Casulli and Cheng (1992).
    aXDen = double(zeros(NUMINCX,1));  % Denominator form.
    aYDen = double(zeros(NUMINCY,1));  % Denominator form.

    % Chezy friction factors.
    Cz = double(zeros(NUMNODES,1));    % Chezy friction factor defined at volume center.
    CzX = double(zeros(NUMINCX,1));    % Chezy friciton factors defined at x-faces
    CzY = double(zeros(NUMINCY,1));    % Chezy friction factors defined at y-faces.
    MnX = double(zeros(NUMINCX,1));    % Manning's n upwinded for x-faces.
    MnY = double(zeros(NUMINCY,1));    % Manning's n upwinded for y-faces.

    % g* corresponds to G from Casulli and Cheng (1992).
    gX = double(zeros(NUMINCX,1));
    gY = double(zeros(NUMINCY,1));

    % Free surface system of equation variables.
    sLHS = sparse(1:NUMNODES,1:NUMNODES,1,NUMNODES,NUMNODES,5*NUMNODES);
    qRHS = double(zeros(NUMNODES,1));

end
%EOF