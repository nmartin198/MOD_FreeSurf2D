function [UVel,VVel] = qiNfLUXbc(UVel,VVel,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % qiNfLUXbc sets up and enforces Dirichlet inflow flux boundaries.
    % Each boundary flux value is set for the inflow face of the boundary volume.
    % The same volume cannot be employed for both an x-face and y-face 
    % Dirichlet inflow flux boundary.
    %
    % UVel [NUMINCX,1] = u or x-face depth-averaged velocity.
    % VVel [NUMINCY,1] = v or y-face depth-averaged velocity.
    % cETime (float) = current elapsed time for time series lookup
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

    global Hux Hvy NUMCOLS PRECH QINXVOL QINXFLUX QINYVOL QINYFLUX

    % First do x-face sources.
    if (QINXVOL(1) ~= 0)
        % Local variables
        NumX = size(QINXVOL,1);             % Number of x-face sources.
        % Calculations
        %  First get volume, row, column, and x-face vector indexes.
        NewNode = QINXVOL;
        Row = ceil(NewNode./NUMCOLS);
        TempCalc1 = mod(NewNode,NUMCOLS);
        BTempC = (TempCalc1 == 0);
        CalcC1 = double(BTempC);
        CalcC2 = double(~BTempC);
        Col = (CalcC1.*NUMCOLS) + (CalcC2.*TempCalc1);
        BTempC = (Col > 1);
        CalcC1 = double(BTempC);
        CalcC2 = double(~BTempC);
        IndexUH = (CalcC1.*((Row-1) + NewNode + 1 ) ) + ...
                  (CalcC2.*((Row-1) + NewNode ) );
        %  Now determine source volume face depth and put into denominator form.
        TempH = Hux(IndexUH);
        BTemp = ((TempH - double(0.0)) > PRECH);
        Calc1 = double(BTemp);
        Calc2 = double(~BTemp);
        TempCalc1 = (Calc1.*TempH) + (Calc2.*PRECH);
        TempH = Calc1.*(1./TempCalc1);
        %  Finally, determine inflow velocities given the specified flux.
        % Because we now have the container, use a loop.
        for iI = 1:NumX
            cIUH = IndexUH(iI);
            cTH = TempH(iI);
            cTQVol = QINXVOL(iI);
            cTSArray = QINXFLUX(cTQVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            tcalcQ1 = cTSinterp * cTH;
            btcQ1 = ((tcalcQ1 - double(0.0)) > PRECH);
            tcQ1 = double(btcQ1);
            tcalcQ1 = tcQ1 * tcalcQ1;
            UVel(cIUH) = (CalcC1(iI) * (-1.0 * tcalcQ1)) + (CalcC2(iI)*tcalcQ1);
        end
        % no longer needed
        %clear NumX BTemp BTempC Calc1 Calc2 CalcC1 CalcC2 Col IndexUH NewNode;
        %clear Row TempCalc1 TempH;
    end
    % Then do y-face sources.
    if (QINYVOL(1) ~= 0)
        % Local variables
        NumY = size(QINYVOL,1);             % Number of y-face sources.
        % Calculations
        %  First get volume, row, and y-face vector indexes.
        NewNode = QINYVOL;
        Row = ceil(NewNode./NUMCOLS);
        BTempR = (Row > 1);
        CalcR1 = double(BTempR);
        CalcR2 = double(~BTempR);
        IndexVH = (CalcR1.*(NewNode + NUMCOLS)) + (CalcR2.*NewNode);
        %  Now determine source volume face depth and put into denominator form.
        TempH = Hvy(IndexVH);
        BTemp = ((TempH - double(0.0)) > PRECH);
        Calc1 = double(BTemp);
        Calc2 = double(~BTemp);
        TempCalc1 = (Calc1.*TempH) + (Calc2.*PRECH);
        TempH = Calc1.*(1./TempCalc1);
        % Finally, determine inflow velocities given the specified flux.
        % Because we now have the container, use a loop.
        for iI = 1:NumY
            cIVH = IndexVH(iI);
            cTH = TempH(iI);
            cTQVol = QINYVOL(iI);
            cTSArray = QINYFLUX(cTQVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            tcalcQ1 = cTSinterp * cTH;
            btcQ1 = ((tcalcQ1 - double(0.0)) > PRECH);
            tcQ1 = double(btcQ1);
            tcalcQ1 = tcQ1 * tcalcQ1;
            VVel(cIVH) = (CalcC1(iI) * (-1.0 * tcalcQ1)) + (CalcC2(iI)*tcalcQ1);
        end
        % not needed
        %clear NumY BTemp BTempR Calc1 Calc2 CalcR1 CalcR2 IndexVH NewNode Row;
        %clear TempCalc1 TempH;
    end

end
%EOF