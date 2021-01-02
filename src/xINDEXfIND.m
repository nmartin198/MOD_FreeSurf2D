function [Index1,Index2,Index3,Index4,newq] = xINDEXfIND(Col,Row,pq,...
   NumInc)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % xINDEXfIND finds the indices to be employed with u-velocity terms for
    % bilinear interpolation.  For this function to work, the coordinate location
    % for each particle needs to be known.  This function then determines the
    % stencil to be employed in the bilinear calculation.
    %
    %   +       +   
    %  -3-------2---             + = volume boundaries.
    %   +       +                | = bilinear interpolation boundaries 
    %  ++++++++++++++                       
    %   +       +   
    %  -4-------1---         ->  positive directions.
    %   +       +           \|/
    %  ++++++++++++++
    %   +       +
    %
    % Col is the column location of each particle.
    % Row location of each particle.
    % pq is the weight determined by the distance between m and
    %     Y-particle location coordinate.
    % NumInc is the number of indexes for the variable of interest.
    %
    % Index1,Index2,Index3, and Index4 have the location show above.
    % newq holds the horizontal weights to be employed in the u - velocity
    %     bilinear interpolation.
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

    global PRECH XINC

    %local variables.
    BTemp = zeros(NumInc,1);        % Boolean calculatino variable.
    Calc1 = double(zeros(NumInc,1));% Calc variable.
    Calc2 = double(zeros(NumInc,1));% Calc variable.
    CalcI1 = zeros(NumInc,1);       % Int Calc variable.
    CalcI2 = zeros(NumInc,1);       % Int Calc variable.
    Index1 = zeros(NumInc,1);
    Index2 = zeros(NumInc,1);
    Index3 = zeros(NumInc,1);
    Index4 = zeros(NumInc,1);
    TempRow = zeros(NumInc,1);    % Row placeholder.
    TempRow1 = zeros(NumInc,1);   % Row placeholder.
    newq = double(zeros(NumInc,1));

    % Calculations
    BTemp = ((pq - 0.5) > PRECH);
    CalcI1 = BTemp;
    CalcI2 = (1 - BTemp);
    TempRow = (CalcI1.*Row) + (CalcI2.*(Row + 1));
    TempRow = rOWaDJ(TempRow,NumInc);
    TempRow1 = (CalcI1.*(Row - 1)) + (CalcI2.*Row);
    TempRow1 = rOWaDJ(TempRow1,NumInc);
    Index1 = (((TempRow - 1).*XINC) + Col) + 1;
    Index2 = (((TempRow1 - 1).*XINC) + Col) + 1;
    Index4 = (((TempRow - 1).*XINC) + Col);
    Index3 = (((TempRow1 - 1).*XINC) + Col);
    % Adjust q to be used for bilinear interpolation of u velocities.
    Calc1 = double(CalcI1);
    Calc2 = double(CalcI2);
    newq = (Calc1.*(pq - 0.5)) + (Calc2.*(pq + 0.5));
    % not needed
    %clear Col BTemp Calc1 Calc2 CalcI1 CalcI2 NumInc pq Row TempRow TempRow1;
end
%EOF