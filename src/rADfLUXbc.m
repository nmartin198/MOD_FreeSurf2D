function [UVel,VVel] = rADfLUXbc(UVel,VVel)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% rADfLUXbc sets up and enforces radiation outflow flux boundaries.
% This script works in conjunction with the rADoRLfSbc.m script to
% create a specified flux, radiation outflow boundary condition.
% This is an over-specified boundary condition.  First, the free surface
% elevation immediately adjacent to the simulation domain will be
% calculated using the radiation boundary condition of Orlanski.
% This free surface elevation will be employed to calculate total
% water depth at the outflow boundary and to generate initial 
% outflow velocities.  The outflow velocities are then updated in
% this script so that the total flux across the specified free
% surface boundary volumes is equal to a given value.
%
% UVel [NUMINCX,1] = u or x-face depth-averaged velocity.
% VVel [NUMINCY,1] = v or y-face depth-averaged velocity.

global DX DY Hux Hvy PREC RFLUXXFLUX RFLUXYFLUX RORLFSXVOL RORLFSYVOL
global OrlSizeX FSXCols FSXIndexBX
global OrlSizeY FSYRows FSYIndexBY

% First do x-face sources.
if (RORLFSXVOL(1) ~= 0)
% Local variables.
   BTemp = 0;                             % Boolean calculation variable.
   BTempV = zeros(OrlSizeX,1);            % Boolean calculation variable.
   BTempV1 = zeros(OrlSizeX,1);           % Boolean calculation variable.
   Calc1 = double(0.0);                   % Calc. variable.
   Calc2 = double(0.0);                   % Calc. variable.
   CalcV1 = double(zeros(OrlSizeX,1));    % Calc. variable.
   CalcV2 = double(zeros(OrlSizeX,1));    % Calc. variable.
   CalcVv1 = double(zeros(OrlSizeX,1));    % Calc. variable.
   CalcVv2 = double(zeros(OrlSizeX,1));    % Calc. variable.
   CFlux = double(0.0);                   % Calculated flux.
   CDenom = double(0.0);                  % Calculated flux in denominator form.
   Depth = double(zeros(1,OrlSizeX));     % Depth a x-face sources.
   Multi = zeros(OrlSizeX,1);             % Direction for outflow.
   Vel = double(zeros(OrlSizeX,1));       % Calculated velocity.
   TempCalc1 = double(0.0);               % Temporary calculation variable.
% Calculations.
%   Get depth at specified boundary locations.
   Depth = Hux(FSXIndexBX)';
%   Get velocity at specified boundary locations.
   Vel = UVel(FSXIndexBX);
   BTempV = (FSXCols > 1);
   CalcV1 = double(BTempV);
   CalcV2 = double(~BTempV);
   BTempV1 = ((Vel - double(0.0)) > PREC);
   CalcVv1 = double(BTempV1);
   CalcVv2 = double(~BTempV1);
   Vel = ((CalcV1.*CalcVv1).*Vel) + ((CalcV2.*CalcVv2).*Vel);
   Vel = abs(Vel);
%   Determine flux at specified boundary locations and put in denominator form.
   CFlux = DY*(Depth*Vel);
   BTemp = ((CFlux - double(0.0)) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1*CFlux) + (Calc2*PREC);
   CDenom = Calc1*(1/TempCalc1);
%   Get flux ratio.
   Ratio = RFLUXXFLUX*CDenom;
% Update velocities to reflect calculations.
   BTempV = (FSXCols > 1);
   CalcV1 = double(BTempV);
   CalcV2 = double(~BTempV);
   Multi = (CalcV1.*1.0) + (CalcV2.*-1.0);
   UVel(FSXIndexBX) = Multi.*(Vel.*Ratio);
% clear variables.
   clear BTemp BTempV BTempV1 Calc1 Calc2 CalcV1 CalcV2 CalcVv1 CalcVv2;
   clear CFlux CDenom Depth Multi Vel TempCalc1;
end
% Next y-face sources.
if (RORLFSYVOL(1) ~= 0)
% Local variables.
   BTemp = 0;                             % Boolean calculation variable.
   BTempV = zeros(OrlSizeY,1);            % Boolean calculation variable.
   BTempV1 = zeros(OrlSizeY,1);           % Boolean calculation variable.
   Calc1 = double(0.0);                   % Calc. variable.
   Calc2 = double(0.0);                   % Calc. variable.
   CalcV1 = double(zeros(OrlSizeY,1));    % Calc. variable.
   CalcV2 = double(zeros(OrlSizeY,1));    % Calc. variable.
   CalcVv1 = double(zeros(OrlSizeY,1));    % Calc. variable.
   CalcVv2 = double(zeros(OrlSizeY,1));    % Calc. variable.
   CFlux = double(0.0);                   % Calculated flux.
   CDenom = double(0.0);                  % Calculated flux in denominator form.
   Depth = double(zeros(1,OrlSizeY));     % Depth a x-face sources.
   Multi = zeros(OrlSizeY,1);             % Direction for outflow.
   Vel = double(zeros(OrlSizeY,1));       % Calculated velocity.
   TempCalc1 = double(0.0);               % Temporary calculation variable.
% Calculations.
%   Get depth at specified boundary locations.
   Depth = Hvy(FSYIndexBY)';
%   Get velocity at specified boundary locations.
   Vel = VVel(FSYIndexBY);
   BTempV = (FSYRows > 1);
   CalcV1 = double(BTempV);
   CalcV2 = double(~BTempV);
   BTempV1 = ((Vel - double(0.0)) > PREC);
   CalcVv1 = double(BTempV1);
   CalcVv2 = double(~BTempV1);
   Vel = ((CalcV1.*CalcVv1).*Vel) + ((CalcV2.*CalcVv2).*Vel);
   Vel = abs(Vel);
%   Determine flux at specified boundary locations and put in denominator form.
   CFlux = DX*(Depth*Vel);
   BTemp = ((CFlux - double(0.0)) > PREC);
   Calc1 = double(BTemp);
   Calc2 = double(~BTemp);
   TempCalc1 = (Calc1*CFlux) + (Calc2*PREC);
   CDenom = Calc1*(1/TempCalc1);
%   Get flux ratio.
   Ratio = RFLUXYFLUX*CDenom;
% Update velocities to reflect calculations.
   BTempV = (FSYRows > 1);
   CalcV1 = double(BTempV);
   CalcV2 = double(~BTempV);
   Multi = (CalcV1.*1.0) + (CalcV2.*-1.0);
   VVel(FSYIndexBY) = Multi.*(Vel.*Ratio);
% clear variables.
   clear BTemp BTempV BTempV1 Calc1 Calc2 CalcV1 CalcV2 CalcVv1 CalcVv2 CFlux;
   clear CDenom Depth Multi Vel TempCalc1;
end  

return;
%EOF
