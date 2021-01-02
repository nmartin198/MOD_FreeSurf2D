function rITh5iLOCS(H5FileP,OutInt)
    %rITh5iLOCS Write ILOC outputs for each time step
    %  H5FileP : HDF5 file name with path
    %  OutInt : current time step index
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
    
    global H UAve VAve NUMILOCS ILOCS
    
    for iI = 1:NUMILOCS
        cNode = ILOCS(iI,1);
        keyStr = sprintf('/Outputs/ILOCs/%d', cNode);
        curRow = single( [ H(cNode) UAve(cNode) VAve(cNode) ] );
        h5write(H5FileP, keyStr, curRow', [ 1 OutInt ], [3 1]);
    end
end
%EOF