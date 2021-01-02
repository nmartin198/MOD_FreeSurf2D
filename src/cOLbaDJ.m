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
end
%EOF