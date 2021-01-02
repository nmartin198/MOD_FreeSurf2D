function [xr,yr] = iNTbOUND(xr,yr,Num)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % iNTbOUND adjusts particle locations for boundaries.
    %
    % Received:
    %
    % xr [Num,1] = x-coordinate particle location.
    % yr [Num,1] = y-coordinate particle lcoation.
    % Num = number of indexes.
    %
    % Returned:
    %
    % xr [Num,1] = adjusted x-coordinate particle location.
    % yr [Num,1] = adjusted y-coordinate particle lcoation.
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

    global NUMCOLS NUMROWS PREC XINDEX YINDEX

    % local variables.
    BoundX1 = zeros(Num,1);    % Boolean to determine if on a boundary.
    BoundX3 = zeros(Num,1);    % Boolean to determine if on a boundary.
    BoundY2 = zeros(Num,1);    % Boolean to determine if on a boundary.
    BoundY4 = zeros(Num,1);    % Boolean to determine if on a boundary.
    BTemp = zeros(Num,1);      % Boolean calc variable.
    Calc1 = double(zeros(Num,1)); % Calc variable.
    Calc2 = double(zeros(Num,1)); % Calc variable.

    % calculations.
    BoundX3 = ((xr - XINDEX(1)) > PREC);
    Calc1 = double(BoundX3);
    Calc2 = double(~BoundX3);
    xr = (Calc2.*XINDEX(1)) + (Calc1.*xr);
    BoundX1 = ((XINDEX(NUMCOLS+1) - xr) > PREC);
    Calc1 = double(BoundX1);
    Calc2 = double(~BoundX1);
    xr = (Calc2.*XINDEX(NUMCOLS+1)) + (Calc1.*xr);
    BoundY2 = ((yr - YINDEX(1)) > PREC);
    Calc1 = double(BoundY2);
    Calc2 = double(~BoundY2);
    yr = (Calc2.*YINDEX(1)) + (Calc1.*yr);
    BoundY4 = ((YINDEX(NUMROWS+1) - yr) > PREC);
    Calc1 = double(BoundY4);
    Calc2 = double(~BoundY4);
    yr = (Calc2.*YINDEX(NUMROWS+1)) + (Calc1.*yr);
    % not needed
    %clear BoundX1 BoundX3 BoundY2 BoundY4 BTemp Calc1 Calc2;
    %clear Num;
end
%EOF