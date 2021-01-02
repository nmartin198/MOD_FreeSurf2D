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

    global NUMROWS

    % local
    BTempD = zeros(Sizer,1);    % bo olean.
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
    % not needed
    %clear BTempD CalcI1 CalcI2;
end
%EOF