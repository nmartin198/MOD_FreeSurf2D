function [UVel,VVel] = vELdIRCbc(UVel,VVel,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vELdIRCbc sets specified inflow velocities at the specified
    % simulation domain boundaries.
    %
    % UVel [NUMINCX,1] =  u or the x-face velocities.
    % VVel [NUMINCY,1] =  v or the y-face velocities.
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

    global NUMCOLS VELDXVOL VELDXVEL VELDYVOL VELDYVEL XINC

    % First do x-face sources.
    if (VELDXVOL(1) ~= 0)
        % local variables.
        NumX = size(VELDXVOL,1);            % Number of x-face total depth sources.
        % Begin calculations.
        NewNode = VELDXVOL;   
        Row = ceil(NewNode./NUMCOLS);
        TempCalc1 = mod(NewNode,NUMCOLS);
        BTemp = (TempCalc1 == 0);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        Col = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
        BTemp = (Col > 1);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        IndexUH = (CalcI1.*(((Row - 1).*XINC) + Col + 1)) + ...
                  (CalcI2.*(((Row - 1).*XINC) + Col));
        % Because we now have the container, use a loop.
        for iI = 1:NumX
            cIUH = IndexUH(iI);
            cTDepVol = VELDXVOL(iI);
            cTSArray = VELDXVEL(cTDepVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            UVel(cIUH,1) = cTSinterp;
        end
        % not needed anymore
        %clear NumX BTemp CalcI1 CalcI2 Col IndexUH NewNode Row TempCalc1;
    end
    % Next do y-face sources.
    if (VELDYVOL(1) ~= 0)
        NumY = size(VELDYVOL,1);  % Number of y-face total depth sources.
        % Begin calculations.
        NewNode = VELDYVOL;   
        Row = ceil(NewNode./NUMCOLS);
        BTemp = (Row > 1);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        IndexVH = (CalcI1.*(NewNode + NUMCOLS)) + (CalcI2.*NewNode);
        % Because we now have the container, use a loop.
        for iI = 1:NumY
            cIVH = IndexVH(iI);
            cTDepVol = VELDYVOL(iI);
            cTSArray = VELDYVEL(cTDepVol)';
            cTSinterp = interp1(cTSArray(:,1),cTSArray(:,2),cETime, ...
                                'linear', cTSArray(size(cTSArray,1),2) );
            VVel(cIVH,1) = cTSinterp;
        end
        % not needed
        %clear NumX BTemp CalcI1 CalcI2 IndexVH NewNode Row;
    end

end
%EOF