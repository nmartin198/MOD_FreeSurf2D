function [R1] = rOWaDJ(R1,Sizer)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% rOWaDJ adjusts row indexes to be within domain boundaries.
%
% Received
%
% R1 [Sizer,1] = vector of row indexes.
% Sizer = size of index.
%
% Returned:
% 
% R1 after adjustement to make sure indexes are within domain.

global NUMROWS

% local
BTempD = zeros(Sizer,1);    % boolean.
CalcI1 = zeros(Sizer,1);    % Int Calc variable.
CalcI2 = zeros(Sizer,1);    % Int Calc variable.

% Calculations.
BTempD = (R1 < 1);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
R1 = CalcI1.*1 + CalcI2.*R1;
BTempD = (R1 > NUMROWS);
CalcI1 = BTempD;
CalcI2 = (1 - BTempD);
R1 = CalcI1.*NUMROWS + CalcI2.*R1;

clear BTempD CalcI1 CalcI2;
return;
%EOF