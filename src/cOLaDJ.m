function [C1] = cOLaDJ(C1,Sizer)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% cOLaDJ adjusts column indexes to be within domain boundaries.
%
% Received
%
% C1 [Sizer,1] = vector of column indexes.
% Sizer = size of index.
%
% Returned:
% 
% C1 after adjustement to make sure indexes are within domain.

global NUMCOLS

% local
BTempD = zeros(Sizer,1);    % boolean.
CalcI1 = zeros(Sizer,1);    % Int Calc variable.
CalcI2 = zeros(Sizer,1);    % Int Calc variable.

% Calculations.
BTempD = (C1 < 1);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
C1 = CalcI1.*1 + CalcI2.*C1;
BTempD = (C1 > NUMCOLS);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
C1 = CalcI1.*NUMCOLS + CalcI2.*C1;

clear BTempD CalcI1 CalcI2;
return;
%EOF