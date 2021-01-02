function [FS1,FS2,FS3,FS4] = tdEPsYSbc(FS1,FS2,FS3,FS4,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % tdEPsYSbc sets Dirichlet boundaries for total water depth.
    % The value of total water depth is set by adding the known
    % free surface elevation at the boundary to the RHS of the
    % free surface system of equations.
    %
    % FS1 [NUMROWS,1] = free surface value to be added to rhs.
    % FS2 [NUMCOLS,1] = free surface value to be added to rhs.
    % FS3 [NUMROWS,1] = free surface value to be added to rhs.
    % FS4 [NUMCOLS,1] = free surface value to be added to rhs.
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
    % along with  MOD_FreeSurf2D.  If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global DATUM NUMCOLS TDEPDXDEP TDEPDXVOL TDEPDYDEP TDEPDYVOL
    global XINC ZtX ZtY

    % First do x-face sources.
    if (TDEPDXVOL(1) ~= 0)
        NumX = size(TDEPDXVOL,1);           % Number of x-face total depth sources.
        FSX = double(zeros(NumX,1));        % Specified free surface elevation.
        % Get Indexes.
        NewNode = TDEPDXVOL;
        Row = ceil(NewNode./NUMCOLS);
        TempCalc1 = mod(NewNode,NUMCOLS);
        BTempC = (TempCalc1 == 0);
        CalcI1 = BTempC;
        CalcI2 = (1 - BTempC);
        Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
        BTempC = (Col > 1);
        CalcI1 = BTempC;
        CalcI2 = (1 - BTempC);
        IndexX = (CalcI1.*(((Row - 1).*XINC)+Col+1)) + ...
                 (CalcI2.*(((Row - 1).*XINC)+Col));
        IndexF1 = find(BTempC == 1);
        IndexF3 = find(BTempC == 0);
        RowF1 = Row(IndexF1);
        RowF3 = Row(IndexF3);
        % Calculate the specified free surface elevation in normalized form.
        % Because we now have the container, use a loop.
        for iI = 1:NumX
            cIndX = IndexX(iI);
            cTDepVol = TDEPDXVOL(iI);
            cTSArray = TDEPDXDEP(cTDepVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            cFSX = (cTSinterp + ZtX(cIndX,1) - DATUM);
            FSX(iI) = cFSX;
        end
        FS1(RowF1) = FSX(IndexF1);
        FS3(RowF3) = FSX(IndexF3);
        % not needed
        %clear NumX BTempC CalcI1 CalcI2 Col FSX IndexF1 IndexF3 IndexX NewNode
        %clear Row RowF1 RowF3 TempCalc1;
    end
    % Then y-face sources.
    if (TDEPDYVOL(1) ~= 0)
        NumY = size(TDEPDYVOL,1);           % Number of y-face source volumes.
        FSY = double(zeros(NumY,1));        % Specified free surface elevation.
        % Get Indexes.
        NewNode = TDEPDYVOL;
        TempCalc1 = mod(NewNode,NUMCOLS);
        BTempR = (TempCalc1 == 0);
        CalcI1 = BTempR;
        CalcI2 = (1 - BTempR);
        Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
        Row = ceil(NewNode./NUMCOLS);
        BTempR = (Row > 1);
        CalcI1 = BTempR;
        CalcI2 = (1 - BTempR);
        IndexY = (CalcI1.*(NewNode + NUMCOLS)) + (CalcI2.*NewNode);
        IndexF4 = find(BTempR == 1);
        IndexF2 = find(BTempR == 0);
        ColF4 = Col(IndexF4);
        ColF2 = Col(IndexF2);
        % Calculate the specified free surface elevation in normalized form.
        % Because we now have the container, use a loop.
        for iI = 1:NumY
            cIndY = IndexY(iI);
            cTDepVol = TDEPDYVOL(iI);
            cTSArray = TDEPDYDEP(cTDepVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            cFSY = (cTSinterp + ZtY(cIndY,1) - DATUM);
            FSY(iI) = cFSY;
        end
        FS4(ColF4) = FSY(IndexF4);
        FS2(ColF2) = FSY(IndexF2);
        % not needed
        %clear NumY BTempR CalcI1 CalcI2 Col ColF2 ColF4 FSY IndexY IndexF4 IndexF2;
        %clear NewNode Row TempCalc1;
    end
end
%EOF