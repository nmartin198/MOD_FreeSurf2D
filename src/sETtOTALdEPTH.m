function [HTDepX,HTDepY,HTD1,HTD2,HTD3,HTD4,BHTx,BHTy] = sETtOTALdEPTH(HTDepX,...
    HTDepY,HTD1,HTD2,HTD3,HTD4,BHTx,BHTy)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETtOTALdEPTH sets the total depth (undisturbed height plus free surface)
%  for the x-faces (Hux) and the y-faces (Hvy).  First the free surface or
%  water elevation across the face is determined.  Then this water elevation is
%  employed with the face undisturbed water depth (HUX or HVY) to get the face
%  total depth.
%
% Received and Returned;
% HTDepX = Hux [NUMINCX,1].  Total water depth at x-faces.
% HTDepY = Hvy [NUMINCY,1].  Total Water depth at y-faces.
% HTD1 = h1 [NUMNODES,1].   Total water depth at i+1/2,j faces.
% HTD2 = h2 [NUMNODES,1].   Total water depth at i,j-1/2 faces.
% HTD3 = h3 [NUMNODES,1].   Total water depth at i-1/2,j faces.
% HTD4 = h4 [NUMNODES,1].   Total water depth at i,j+1/2 faces.
% BHTx = BHux [NUMINCX,1].  Multiplier for wet x-face. 1.0 = wet. 0.0 =
%           dry.
% BHTy = BHvy [NUMINYX,1].  Multiplier for wet y-face. 1.0 = wet. 0.0 =
%           dry.

global EtaNew EtaXM1 EtaXP1 EtaYM1 EtaYP1 PRECH HCUTOFF HNODE HUX HVY
global NUMINCX NUMINCY NUMNODES PREC TDEPDIRCBC UXP1 UXM1 VYM1 VYP1
global HuxOld HvyOld

% local variables.
BTemp = zeros(NUMNODES,1);             % Boolean calculation variable.
Calc1 = double(zeros(NUMNODES,1));     % Calc. variable.
BTempX = zeros(NUMINCX,1);             % Boolean calculation variable.
CalcX1 = double(zeros(NUMINCX,1));     % Calc variable.
CalcX2 = double(zeros(NUMINCX,1));     % Calc variable.
BTempY = zeros(NUMINCY,1);             % Boolean calculation variable.
CalcY1 = double(zeros(NUMINCY,1));     % Calc variable.
CalcY2 = double(zeros(NUMINCY,1));     % Calc variable.
CentDepth = double(zeros(NUMNODES,1)); % Depth at volume center used to determine
FsX = double(zeros(NUMINCX,1));        % Free surface elevation for depth calc.
FsY = double(zeros(NUMINCY,1));        % Free surface elevation for depth calc.
HDepth1 = double(zeros(NUMINCX,1));    % Depth calculated using free surf from (i,j+1);
HDepth3 = double(zeros(NUMINCX,1));    % Depth calc. using fs from (i,j-1);
HDepth2 = double(zeros(NUMINCY,1));    % Depth calculated using free surf from (i+1,j);
HDepth4 = double(zeros(NUMINCY,1));    % Depth calc. using fs from (i,j);
MnFsX = double(zeros(NUMINCX,1));      % Smaller free surface across x-face.
MxFsX = double(zeros(NUMINCX,1));      % Larger free surface across x-face.
MnFsY = double(zeros(NUMINCY,1));      % Smaller free surface across y-face.
MxFsY = double(zeros(NUMINCY,1));      % Larger free surface across y-face.

% Save to old vectors.
HuxOld = HTDepX;
HvyOld = HTDepY;
% First calculate the volume center depth.
CentDepth = HNODE + EtaNew;
BTemp = ((CentDepth - double(0.0)) > PRECH);
Calc1 = double(BTemp);
CentDepth = Calc1.*CentDepth;
% Now do x-faces.
% First allocate volume center depth to x-face vectors.
HDepth3 = CentDepth(UXM1);
HDepth1 = CentDepth(UXP1);
% Next calculate the maximum free surface elevation across the 
%  x-faces.
BTempX = ((EtaXP1 - EtaXM1) > -PREC);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
DepthX = (CalcX1.*HDepth1) + (CalcX2.*HDepth3);
MxFsX = (CalcX1.*EtaXP1) + (CalcX2.*EtaXM1);
MnFsX = (CalcX1.*EtaXM1) + (CalcX2.*EtaXP1);
BTempX = ((DepthX - double(0.0)) > PRECH);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
FsX = (CalcX1.*MxFsX) + (CalcX2.*MnFsX);
% Finally use this maximum free surface elevation and the 
%  face topographic elevation to determine the face total
%  water depth.
HTDepX = CalcX1.*(HUX + FsX);
% Ensure that depth is greater or equal to zero.
BTempX = ((HTDepX - double(0.0)) > HCUTOFF);
BHTx = double(BTempX);
HTDepX = BHTx.*HTDepX;
% Do Y-faces.
% First allocate volume center depth to x-face vectors.
HDepth2 = CentDepth(VYM1);
HDepth4 = CentDepth(VYP1);
% Next calculate the maximum free surface elevation across the 
%  x-faces.
BTempY = ((EtaYP1 - EtaYM1) > -PREC);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
DepthY = (CalcY1.*HDepth4) + (CalcY2.*HDepth2);
MxFsY = (CalcY1.*EtaYP1) + (CalcY2.*EtaYM1);
MnFsY = (CalcY1.*EtaYM1) + (CalcY2.*EtaYP1);
BTempY = ((DepthY - double(0.0)) > PRECH);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
FsY = (CalcY1.*MxFsY) + (CalcY2.*MnFsY);
% Finally use this maximum free surface elevation and the 
%  face topographic elevation to determine the face total
%  water depth.
HTDepY = CalcY1.*(HVY + FsY);
% Ensure that depth is greater or equal to zero.
BTempY = ((HTDepY - double(0.0)) > HCUTOFF);
BHTy = double(BTempY);
HTDepY = BHTy.*HTDepY;

% Now set total depth for prescribed depth boundaries.
if (TDEPDIRCBC == 1)
   [HTDepX,HTDepY] = tdEPdIRCbc(HTDepX,HTDepY);
end
% Now create the vectors representing total water depth at each volume.
[HTD1,HTD2,HTD3,HTD4] = aLLOCfACE(HTDepX,HTDepY,HTD1,HTD2,HTD3,HTD4);

clear BTemp BTempX BTempY Calc1 CalcX1 CalcX2 CalcY1 CalcY2 CentDepth;
clear FsX FsY HDepth1 HDepth3 HDepth2 HDepth4 MnFsX MxFsX MnFsY MxFsY;
return;
%EOF