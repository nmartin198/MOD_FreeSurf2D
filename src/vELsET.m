function [UVel,VVel] = vELsET(UVel,VVel)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % vELsET sets the initial velocity values.  If the files BU.txt exists, then
    % UVel will be == BU.txt, else UVel = zeros.  For VVel, the file BV.txt needs
    % to exist.  For both BU.txt and BV.txt, the file format should be ASCII column
    % vector with one velocity value on each line.  The respective lengths should be
    % NUMINCX for BU.txt and NUMINCY for BV.txt.
    %
    % Received and Returned:
    %
    % UVel [NUMINCX,1] = x-face velocity.(u)
    % VVel [NUMINCY,1] = y-face velocity.(v)
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

    % global variables.
    global NUMINCX NUMINCY
    % local variables.
    Ux = 0;
    Vy = 0;

    % Try to open text files.
    Ux = fopen('BU.txt');
    Vy = fopen('BV.txt');

    % If file exists, read.  Else close.
    if (Ux ~= -1)
       UVel = fscanf(Ux,'%f',NUMINCX);
       fclose(Ux);
    end
    % If file exists, read.  Else close.
    if (Vy ~= -1)
       VVel = fscanf(Vy,'%f',NUMINCY);
       fclose(Vy);
    end
    % no longer needed
    %clear Ux Vy;
end
%EOF