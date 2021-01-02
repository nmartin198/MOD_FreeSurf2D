function rITh5oUTPUT(H5FileP,OutInt)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rITh5oUTPUT writes U, V, Hux, Hvy, and H matrices to the HDF5 file
    %   at the specified output interval.
    %
    %   H5FileP : HDF5 file name with path
    %   OutInt: current index for output time
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

    global H Hux Hvy u v NUMNODES NUMINCX NUMINCY
    
    keyStr = sprintf('/Outputs/Out_%d/H', OutInt);
    h5create(H5FileP, keyStr, [1 NUMNODES], 'DataType', 'single' );
    h5write(H5FileP, keyStr, single(H') );
    keyStr = sprintf('/Outputs/Out_%d/Hux', OutInt);
    h5create(H5FileP, keyStr, [1 NUMINCX], 'DataType', 'single' );
    h5write(H5FileP, keyStr, single(Hux') );
    keyStr = sprintf('/Outputs/Out_%d/Hvy', OutInt);
    h5create(H5FileP, keyStr, [1 NUMINCY], 'DataType', 'single' );
    h5write(H5FileP, keyStr, single(Hvy') );
    keyStr = sprintf('/Outputs/Out_%d/u', OutInt);
    h5create(H5FileP, keyStr, [1 NUMINCX], 'DataType', 'single' );
    h5write(H5FileP, keyStr, single(u') );
    keyStr = sprintf('/Outputs/Out_%d/v', OutInt);
    h5create(H5FileP, keyStr, [1 NUMINCY], 'DataType', 'single' );
    h5write(H5FileP, keyStr, single(v') );

end
%EOF