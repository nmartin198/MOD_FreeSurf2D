function FVel = fVYcALC(FVel,NUMINC)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % fVYcALC calculates the inertial terms for the y-faces of volumes.  The inertial
    % terms include the advective, diffusive, and Coriolis terms.  This component,
    % FVel, is equivalent to the Fv term from Casulli and Cheng (1992).   The Fv or
    % FVel term is obtained with a semi-Lagrangian algorithm.  In this method, 
    % a particle, originally located at the y-face of interest, is traced back
    % along the particle pathline (or streamline) to find the particles position
    % at time step N (current time step would be N+1) which can also be called
    % the particle's departure point.  The departure point provides the location for
    % the generation of the advective, diffusive, and Coriolis terms.  As a result,
    % the Fv terms are generated with an explicit semi-Lagrangian method.
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
    global BHvy DX DY fdt MAXSTEPS MINSTEPS NUMNODES
    global OUTP PATHTRAC PREC u UAve v VDirect
    global VYM1 VYP1 VVM1 VVP1 YXeI YYeI

    % Initialize variables.

    AveU = double(zeros(NUMNODES,1));   % Average u - velocity.
    BTemp = zeros(NUMINC,1);            % Boolean temporary calculation variable.
    CalcU1 = double(zeros(NUMINC,1));   % Calc variable.
    CalcU2 = double(zeros(NUMINC,1));   % Calc variable.
    CalcV1 = double(zeros(NUMINC,1));   % Calc variable.
    CalcV2 = double(zeros(NUMINC,1));   % Calc variable.
    Col = zeros(NUMINC,1);              % Column locations for each particle.
    DiffTerm = double(zeros(NUMINC,1)); % Diffusion terms.
    DirectU = zeros(NUMINC,1);          % u velocity direction.
    GradX = double(zeros(NUMINC,1));    % vel. gradient, x, that fluid particle will see.
    GradY = double(zeros(NUMINC,1));    % vel. gradient, y, that fluid particle will see.
    m = double(zeros(NUMINC,1));        % YINDEX locations for each particle.
    n = double(zeros(NUMINC,1));        % XINDEX locations for each particle.
    Nodes = zeros(NUMINC,1);            % Node Number for each particle.
    p = double(zeros(NUMINC,1));        % X weights for each particle - distance from n.
    q = double(zeros(NUMINC,1));        % Y weights for each particle - distance from m.
    Row = zeros(NUMINC,1);              % Row locations for each particle.
    TempCalc1 = double(zeros(NUMINC,1));% Temporary calculation variable.
    TempCalc2 = double(zeros(NUMINC,1));% Temporary calculation variable.
    TempDX = 0.5*DX;                    % Adjust DX for y-face location.
    TempDY = DY;                        % Adjust DY for y-face location.
    u1 = double(zeros(NUMNODES,1));     % velocity at (i+1/2,j).
    u3 = double(zeros(NUMNODES,1));     % velocity at (i-1/2,j).
    Uent = double(zeros(NUMINC,1));     % downstream u velocity.
    UVel = double(zeros(NUMINC,1));     % averaged v - velocity value for each x-face.
    Uxit = double(zeros(NUMINC,1));     % upwind u - velocity.
    v2 = double(zeros(NUMNODES,1));     % velocity at (i,j-1/2).
    v4 = double(zeros(NUMNODES,1));     % velocity at (i,j+1/2).
    VVel = double(v);                   % v velocity at initial location.
    Vxit = double(zeros(NUMINC,1));     % upwind v velocity.
    XBnd = double(zeros(NUMINC,1));     % downwind x-face coordinates.
    Xe = double(zeros(NUMINC,1));       % X location of the particle org. at face.
    Ye = double(zeros(NUMINC,1));       % Y location of particle org. at face.

    Xe = YXeI;
    Ye = YYeI;
    % Use a pathline tracing method to follow the particle, currently located at the
    % v - velocity location back to the particle location at the beginning of the
    % time step.
    %
    % The initial values change depending on which face, u or v, is being 
    % evaluated.  Calculate the values needed by the program for the first
    % step back.  To do this calculation - loop through the domain by row 
    % of v-velocity locations.  There are NUMROWS + 1 rows of v-velocity
    % locations.
    [u1,v2,u3,v4] = aLLOCfACE(u,v,u1,v2,u3,v4);
    BTemp = (VDirect == 1);
    CalcV1 = double(BTemp);
    CalcV2 = double(~BTemp);
    Vxit = (CalcV1.*v(VVM1)) + (CalcV2.*v(VVP1));
    UVel = (CalcV1.*UAve(VYM1)) + (CalcV2.*UAve(VYP1));
    DirectU = (((UVel - double(0.0)) >= -PREC).*1) + ...
       (((UVel - double(0.0)) < -PREC).*-1);
    BTemp = (DirectU == 1);
    CalcU1 = double(BTemp);
    CalcU2 = double(~BTemp);
    XBnd = (CalcU1.*(Xe + (0.5*DX))) + (CalcU2.*(Xe - (0.5*DX)));
    Uxit = (CalcV1.*((CalcU1.*u3(VYM1)) + (CalcU2.*u1(VYM1)))) + ...
       (CalcV2.*((CalcU1.*u3(VYP1)) + (CalcU2.*u1(VYP1))));
    Uent = (CalcV1.*((CalcU1.*u1(VYM1)) + (CalcU2.*u3(VYM1)))) + ...
       (CalcV2.*((CalcU1.*u1(VYP1)) + (CalcU2.*u3(VYP1))));
    % Correct UVel for when depth at V-location is zero.
    UVel = BHvy.*UVel;
    % Now Calculate the gradients at the initial locations.
    BTemp = ((UVel - double(0.0)) > -PREC);
    CalcU1 = double(BTemp);
    CalcU2 = double(~BTemp);
    TempCalc1 = (Uxit - Uent)./DX;
    TempCalc2 = (Uent - Uxit)./DX;
    GradX = (CalcU1.*TempCalc2) + (CalcU2.*TempCalc1);
    BTemp = ((VVel - double(0.0)) > -PREC);
    CalcV1 = double(BTemp);
    CalcV2 = double(~BTemp);
    TempCalc1 = (Vxit - VVel)./DY;
    TempCalc2 = (VVel - Vxit)./DY;
    GradY = (CalcV1.*TempCalc2) + (CalcV2.*TempCalc1);
    % Now trace back the path lines.
    if (PATHTRAC == 1)
       [Xe,Ye,n,m,Col,Row,Nodes] = sEMIapATH(Xe,Ye,UVel,VVel,Uent,VVel,Uxit,...
          Vxit,u,v,GradX,GradY,XBnd,Ye,TempDX,TempDY,fdt,BHvy,NUMINC);
    elseif (PATHTRAC == 2)
       [Xe,Ye,n,m,Col,Row,Nodes] = rkiNT(Xe,Ye,UVel,VVel,u,v,fdt,MAXSTEPS,...
          MINSTEPS,NUMINC);
    elseif (PATHTRAC == 3)
       [Xe,Ye,n,m,Col,Row,Nodes] = eULiNT(Xe,Ye,UVel,VVel,u,v,fdt,MAXSTEPS,...
          MINSTEPS,NUMINC);
    else
       disp('Invalid value for PATHTRAC in fVYcALC.m');
       fprintf(OUTP,'Invalid value for PATHTRAC in fVYcALC.m \n');
    end
    %Calculate velocity, diffusion terms and Coriolis terms for integrated locations.
    % Calculate volume weights to be employed in setting up the stencil
    % for cubic Lagrange interpolation, and for calculating the velocity term
    % for Coriolis calculations.
    p = (abs((n - Xe)))./DX;
    q = (abs((m - Ye)))./DY;
    % Coriolis velocity calculations
    UVel = futERMcALC(UVel,Col,Row,q,p,NUMINC);
    % interpolate convective velocity values with cubic Lagrange polynomials.
    % VVel is passed as v, but returned as the interpolated velocity at the
    % foot of the characteristic.
    VVel = vcUBIClAGiNT(VVel,Xe,Ye,n,m,Col,Nodes,Row,p,NUMINC);
    %Diffusion term calculation.
    DiffTerm = vdIFFcALC(DiffTerm,Col,Row,p,NUMINC);
    %Calculate Fv term from TRIM method.
    FVel = VVel + DiffTerm - UVel;
    % not needed
    %clear BTemp CalcU1 CalcU2 CalcV1 CalcV2 Col DiffTerm GradX GradY m n Nodes;
    %clear NUMINC p q Row TempCalc1 TempCalc2 TempDX TempDY UVel VVel Xe Ye Uxit;
    %clear Uent Vxit XBnd DirectU u1 v2 u3 v4;
end
%EOF