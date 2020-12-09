function [C1] = cOLbaDJ(C1,Sizer)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cOLbaDJ adjusts column indexes to be within domain boundaries.
%
% Received
%
% C1 [Sizer,1] = vector of column indexes.
% Sizer = size of index.
%
% Returned:
% 
% C1 after adjustement to make sure indexes are within domain.

global XINC

% local
BTempD = zeros(Sizer,1);    % boolean.
CalcI1 = zeros(Sizer,1);    % Calc variable.
CalcI2 = zeros(Sizer,1);    % Calc variable.

% Calculations.
BTempD = (C1 < 1);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
C1 = CalcI1.*1 + CalcI2.*C1;
BTempD = (C1 > XINC);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
C1 = CalcI1.*XINC + CalcI2.*C1;

clear BTempD CalcI1 CalcI2;
return;
%EOF