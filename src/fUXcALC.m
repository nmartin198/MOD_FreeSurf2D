function FVel = fUXcALC(FVel,NUMINC)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fUXcALC calculates the inertial terms for the x-faces of volumes.  The inertial
    % terms include the advective, diffusive, and Coriolis terms.  This component,
    % FVel, is equivalent to the Fu term from Casulli and Cheng (1992).   The Fu or
    % FVel term is obtained with a semi-Lagrangian algorithm.  In this method, 
    % a particle, originally located at the x-face of interest, is traced back
    % along the particle pathline (or streamline) to find the particles position
    % at time step N (current time step would be N+1) which can also be called
    % the particle's departure point.  The departure point provides the location for
    % the generation of the advective, diffusive, and Coriolis terms.  As a result,
    % the Fu terms are generated with an explicit semi-Lagrangian method.
    %
    % Two different methods are available to trace the particle pathline.  One
    % method is the classical, four-step, explicit Runge-Kutta method.  The
    % other method is the semi-analytical path line tracing method of Pollock
    % (1988).
    %
    % At the departure point, several different interpolation schemes are employed
    % to obtain the values of interest.  For the advective contribution, cubic 
    % Lagrange polynomials in 2-D are used to determine the velocity from the
    % surrounding Eulerian grid.  For the diffusive and Coriolis terms, bilinear
    % interpolation is employed.  Also for the diffusive term, either a Smagorinsky
    % model eddy diffusivity or a specified constant eddy diffusivity may be
    % employed.
    %
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

    % global variables.
    global BHux DX DY fdt OUTP MAXSTEPS MINSTEPS NUMNODES
    global PATHTRAC PREC u UDirect UVM1 UVP1 UXM1
    global UXP1 v VAve XXeI XYeI

    % Initialize local variables.

    BTemp = zeros(NUMINC,1);                  % Boolean temporary calculation variable.
    CalcU1 = double(zeros(NUMINC,1));          % Calc variable.
    CalcU2 = double(zeros(NUMINC,1));          % Calc variable.
    CalcV1 = double(zeros(NUMINC,1));          % Calc variable.
    CalcV2 = double(zeros(NUMINC,1));          % Calc variable.
    Col = zeros(NUMINC,1);                    % Column locations for each particle.
    DiffTerm = double(zeros(NUMINC,1));       % Diffusion terms.
    DirectV = zeros(NUMINC,1);                % Average v direction.
    GradX = double(zeros(NUMINC,1));          % vel. gradient, x, that fluid particle will see.
    GradY = double(zeros(NUMINC,1));          % vel. gradient, y, that fluid particle will see.
    m = double(zeros(NUMINC,1));              % YINDEX locations for each particle.
    n = double(zeros(NUMINC,1));              % XINDEX locations for each particle.
    Nodes = zeros(NUMINC,1);                  % Node Number for each particle.
    p = double(zeros(NUMINC,1));              % X weights for each particle - distance from n.
    q = double(zeros(NUMINC,1));              % Y weights for each particle - distance from m.
    Row = zeros(NUMINC,1);                    % Row locations for each particle.
    TempCalc1 = double(zeros(NUMINC,1));      % Temporary calculation variable.
    TempCalc2 = double(zeros(NUMINC,1));      % Temporary calculation variable.
    TempDX = DX.*double(ones(NUMINC,1));      % Adjust DX for x-face location.
    TempDY = (0.5*DY).*double(ones(NUMINC,1));% Adjust DY for x-face location.
    u1 = double(zeros(NUMNODES,1));           % velocity at (i+1/2,j).
    u3 = double(zeros(NUMNODES,1));           % velocity at (i-1/2,j).
    UVel = double(u);                         % averaged v - velocity value for each x-face.
    Uxit = double(zeros(NUMINC,1));           % potential x-exit velocity.
    v2 = double(zeros(NUMNODES,1));           % velocity at (i,j-1/2).
    v4 = double(zeros(NUMNODES,1));           % velocity at (i,j+1/2).
    VVel = double(zeros(NUMINC,1));           % v velocity at initial location.
    VBnd = double(zeros(NUMINC,1));           % y-coordinate of downstream y-face.
    Vent = double(zeros(NUMINC,1));           % v velocity at y-direction entrance.
    Vxit = double(zeros(NUMINC,1));           % potential y-exit velocity.
    Xe = double(zeros(NUMINC,1));             % X location of the particle org. at face.
    Ye = double(zeros(NUMINC,1));             % Y location of particle org. at face.

    Xe = XXeI;
    Ye = XYeI;

    %*****************************************************************************
    % Use a pathline tracing method to follow the particle, currently located at the
    % u - velocity location back to the particle location at the beginning of the
    % time step.
    %
    % The initial values change depending on which face, u or v, is being 
    % evaluated.  Calculate the values needed by the program for the first
    % step back.  To do this calculation - loop through the domain by column 
    % of u -velocity locations.  There are NUMCOLS + 1 rows of u-velocity
    % locations.
    % set-up.
    [u1,v2,u3,v4] = aLLOCfACE(u,v,u1,v2,u3,v4);
    BTemp = (UDirect == 1);
    CalcU1 = double(BTemp);
    CalcU2 = double(~BTemp);
    Uxit = (CalcU1.*u(UVM1)) + (CalcU2.*u(UVP1));
    VVel = (CalcU1.*VAve(UXM1)) + (CalcU2.*VAve(UXP1));
    DirectV = (((VVel - double(0.0)) >= -PREC).*1) + ...
       (((VVel - double(0.0)) < -PREC).*-1);
    BTemp = (DirectV == 1);
    CalcV1 = double(BTemp);
    CalcV2 = double(~BTemp);
    VBnd = (CalcV1.*(Ye + (0.5*DY))) + (CalcV2.*(Ye - (0.5*DY)));
    Vxit = (CalcU1.*((CalcV1.*v2(UXM1)) + (CalcV2.*v4(UXM1)))) + ...
       (CalcU2.*((CalcV1.*v2(UXP1)) + (CalcV2.*v4(UXP1))));
    Vent = (CalcU1.*((CalcV1.*v4(UXM1)) + (CalcV2.*v2(UXM1)))) + ...
       (CalcU2.*((CalcV1.*v4(UXP1)) + (CalcV2.*v2(UXP1))));
    % adjust VVel for when depth is zero.
    VVel = BHux.*VVel;
    % Now Calculate the gradients.
    BTemp = ((UVel - double(0.0)) > -PREC);
    CalcU1 = double(BTemp);
    CalcU2 = double(~BTemp);
    TempCalc1 = (Uxit - UVel)./DX;
    TempCalc2 = (UVel - Uxit)./DX;
    GradX = (CalcU1.*TempCalc2) + (CalcU2.*TempCalc1);
    BTemp = ((VVel - double(0.0)) > -PREC);
    CalcV1 = double(BTemp);
    CalcV2 = double(~BTemp);
    TempCalc1 = (Vxit - Vent)./DY;
    TempCalc2 = (Vent - Vxit)./DY;
    GradY = (CalcV1.*TempCalc2) + (CalcV2.*TempCalc1);
    % Now trace back the path lines.
    if (PATHTRAC == 1)
       [Xe,Ye,n,m,Col,Row,Nodes] = sEMIapATH(Xe,Ye,UVel,VVel,UVel,Vent,Uxit,...
          Vxit,u,v,GradX,GradY,Xe,VBnd,TempDX,TempDY,fdt,BHux,NUMINC);
    elseif (PATHTRAC == 2)
       [Xe,Ye,n,m,Col,Row,Nodes] = rkiNT(Xe,Ye,UVel,VVel,u,v,fdt,MAXSTEPS,...
          MINSTEPS,NUMINC);
    elseif (PATHTRAC == 3)
       [Xe,Ye,n,m,Col,Row,Nodes] = eULiNT(Xe,Ye,UVel,VVel,u,v,fdt,MAXSTEPS,...
          MINSTEPS,NUMINC);
    else
       disp('Invalid value for PATHTRAC in fUXcALC.m');
       fprintf(OUTP,'Invalid value for PATHTRAC in fUXcALC.m \n');
    end
    %%Calculate velocity, diffusion terms and Coriolis terms for integrated locations.
    %
    % Calculate volume weights to be employed in setting up the stencil
    % for cubic Lagrange interpolation, and for calculating the velocity term
    % for Coriolis calculations.
    p = (abs((n - Xe)))./DX;
    q = (abs((m - Ye)))./DY;
    % Coriolis velocity calculations
    VVel = fvtERMcALC(VVel,Col,Row,q,p,NUMINC);
    % Velocity interpolation for advective terms with Lagrange cubic polynomials.
    %  UVel that is passed is u.  UVel returned is interpolated velocity at the 
    %  foot of the characteristic.
    UVel = ucUBIClAGiNT(UVel,Xe,Ye,n,m,Col,Nodes,Row,q,NUMINC);
    %Diffusion term calculation.
    DiffTerm = udIFFcALC(DiffTerm,Col,Row,p,NUMINC);
    %
    %Calculate Fv term from TRIM method.
    FVel = UVel + DiffTerm + VVel;
    % not needed
    %clear BTemp NUMINC Col DiffTerm GradX GradY m n Nodes p q Row TempCalc1;
    %clear TempCalc2 TempDX TempDY UVel Uxit VVel Vxit Vent Xe Ye VBnd DirectV;
    %clear u1 v2 u3 v4;
    %clear CalcU1 CalcU2 CalcV1 CalcV2;
end
%EOF