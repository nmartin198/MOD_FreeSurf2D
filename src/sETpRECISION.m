function [prec,sqprec,mindepth] = sETpRECISION(prec,sqprec,mindepth)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% sETprECISION sets numeric precision values for the program.  prec is the
% smallest number to employ.  sqprec is the square root of the smallest
% number.  mindepth is the minimum water depth; the value of HCUTOFF will
% be used unless this value is smaller than sqprec.
%
% Received and returned:
%
% prec [1] = absolute value of smallest real number. (PREC)
% sqprec [1] = square root of prec. (PRECH)
% mindepth [1] = the minimum water depth. (HCUTOFF)

global HCUTOFF

% local variables.
BTemp = 0;                  % Boolean calc. variable.
Calc1 = double(0.0);        % Calc variable.
Calc2 = double(0.0);        % Calc variable.

% set prec to machine precision.
prec = eps;
% set sqprec to the square root of machine precision.
sqprec = eps^0.5;
% set the minimum water depth using the HCUTOFF parameter.
BTemp = ((HCUTOFF - sqprec) > prec);
Calc1 = double(BTemp);
Calc2 = double(~BTemp);
mindepth = (Calc1.*HCUTOFF) + (Calc2.*sqprec);

clear BTemp Calc1 Calc2;
return;
%EOF