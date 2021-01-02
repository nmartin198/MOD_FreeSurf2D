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
    % along with  MOD_FreeSurf2D.  If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % no longer needed
    %clear BTemp Calc1 Calc2;
end
%EOF