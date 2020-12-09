function [Chezy,ChezyX,ChezyY] = cHEZYcALC(MnTotal,nXFace,nYFace,Chezy,ChezyX,ChezyY)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cHEZYcALC calculates the Chezy friction factor given a value of 
% Manning's n for each volume.  These values of Manning's N are stored
% in MN. Upwinded values of n are employed to calculate the Chezy
% friction factor for volume faces.
%
% Recieve:
% MnTotal = MN [NUMNODES,1] or Manning's n.
% nXFace [NUMINCX,1] = MnX or Manning's n at x-faces by upwinding.
% nYFace [NUMINCY,1] = MnY or Manning's n at y-faces by upwinding.
% Chezy = Cz [NUMNODES,1] = volume center Chezy friction factor.
% ChezyX = CzX [NUMINCX,1] = x-face Chezy friction factor.
% ChezyY = CzY [NUMINCY,1] = y-face Chezy friction factor.
%
% Return:
% Chezy = Cz [NUMNODES,1] = volume center Chezy friction factor.
% ChezyX = CzX [NUMINCX,1] = x-face Chezy friction factor.
% ChezyY = CzY [NUMINCY,1] = y-face Chezy friction factor.

global H Hux Hvy NUMNODES NUMINCX NUMINCY PREC

% returned variables.
BTemp = zeros(NUMNODES,1);             % Boolean to avoid divide by zero.
BTempX = zeros(NUMINCX,1);             % Boolean to avoid divide by zero.
BTempY = zeros(NUMINCY,1);             % Boolean to avoid divide by zero.
Calc1 = double(zeros(NUMNODES,1));     % Calc. variable.
Calc2 = double(zeros(NUMNODES,1));     % Calc. variable.
CalcX1 = double(zeros(NUMINCX,1));     % Calc. variable.
CalcX2 = double(zeros(NUMINCX,1));     % Calc. variable.
CalcY1 = double(zeros(NUMINCX,1));     % Calc. variable.
CalcY2 = double(zeros(NUMINCX,1));     % Calc. variable.
DTotal = double(zeros(NUMNODES,1));    % Denominator form of MN.
DXFace = double(zeros(NUMINCX,1));     % Denominator form of MnX.
DYFace = double(zeros(NUMINCY,1));     % Denominator form of MnY.
TCalc = double(zeros(NUMNODES,1));     % Temporary calculation variable.
TCalcX = double(zeros(NUMINCX,1));     % Temporary calculation variable.
TCalcY = double(zeros(NUMINCY,1));     % Temporary calculation variable.

% get denominator forms.
BTemp = ((MnTotal - double(0.0)) > PREC);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
BTempX = ((nXFace - double(0.0)) > PREC);
CalcX1 = double(BTempX);
CalcX2 = double(~BTempX);
BTempY = ((nYFace - double(0.0)) > PREC);
CalcY1 = double(BTempY);
CalcY2 = double(~BTempY);
TCalc = (Calc1.*MnTotal) + (Calc2.*PREC);
TCalcX = (CalcX1.*nXFace) + (CalcX2.*PREC);
TCalcY = (CalcY1.*nYFace) + (CalcY2.*PREC);
DTotal = 1./TCalc;
DXFace = 1./TCalcX;
DYFace = 1./TCalcY;
% Calculate Volume center value.
Chezy = (H.^(1/6)).*DTotal;
% Now calculate face Chezy coefficients.
ChezyX = (Hux.^(1/6)).*DXFace;
ChezyY = (Hvy.^(1/6)).*DYFace;

clear MnTotal nXFace nYFace;
clear BTemp BTempX BTempY DTotal DXFace DYFace TCalc TCalcX TCalcY;
clear Calc1 Calc2 CalcX1 CalcX2 CalcY1 CalcY2;
return;
%EOF