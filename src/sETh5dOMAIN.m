function sETh5dOMAIN(H5FileP)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % sETdOMAIN creates the simulation domain from HDF5 file.
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

    global DATUM DX DY HDEPTH LASTROW MN NUMCOLS NUMINCX NUMINCY
    global NUMNODES NUMROWS RINC ROWBEGIN ROWEND ULASTROW UVP1 UVM1 
    global UXM1 UXP1 X1 X3 XP1 XM1 XINC XINDEX Y2 Y4 YINDEX YM1 YP1
    global ZTopo ZtX ZtY VYP1 VYM1 VVP1 VVM1
    global XXeI XYeI YXeI YYeI

    % intialize variables.
    Depth = double(zeros(NUMROWS,NUMCOLS)); % Water depth from input.
    Fid1 = 0;                               % File indentifier.
    Fid2 = 0;                               % File indentifier.
    Fid3 = 0;                               % File indentifier.
    LASTROW = 0;                            % Last volume in second to last row.
    Mann = double(zeros(NUMROWS,NUMCOLS));  % Manning's n from input.
    MaxTopo = double(0.0);                  % Maximum topographic height [m].
    NUMNODES = 0;                           % Number of nodes in the simulation
    NUMINCX = 0;                            % Number of x-index locations.
    NUMINCY = 0;                            % Number of y-index locations.
    RINC = 0;                               % y-face and volume vector increment.
    TempDepth = double(zeros(NUMCOLS,NUMROWS)); % Temporary read in matrix.
    TempMann = double(zeros(NUMCOLS,NUMROWS));  % Temporary read in matrix.
    TempTopo = double(zeros(NUMCOLS,NUMROWS));  % Temporary read in matrix.
    Topo = double(zeros(NUMROWS,NUMCOLS));  % Domain topography from input.
    ULASTROW = 0;                           % Last x-face loc in 2nd to last row.
    XINC = 0;                               % x-face col. vector increment.
    NodesOut = int32(zeros(NUMROWS, NUMCOLS)); % node indexes for output
    XFIndOut = int32(zeros(NUMROWS, NUMCOLS+1)); % x-face indexes for output
    YFIndOut = int32(zeros(NUMROWS+1, NUMCOLS)); % y-face indexes for output
    % allocate.
    NUMNODES = NUMROWS*NUMCOLS;
    LASTROW = (NUMROWS-1)*NUMCOLS;
    NUMINCX = ((NUMROWS)*(NUMCOLS+1));
    NUMINCY = ((NUMROWS+1)*(NUMCOLS));
    RINC = NUMCOLS;
    ULASTROW = ((NUMROWS-1)*(NUMCOLS+1));
    XINC = NUMCOLS+1;
    % Read in from HDF5 file
    TempTopo = h5read( H5FileP, '/Inputs/Topo/' );
    TempDepth = h5read( H5FileP, '/Inputs/H/' );
    TempMann = h5read( H5FileP, '/Inputs/Mann/' );
    Topo = double(TempTopo');
    Depth = double(TempDepth');
    Mann = double(TempMann');
    % Calculations.
    MaxTopo = max(max(Topo));
    if (DATUM == -9999)
       DATUM = MaxTopo + 1.0;       % DATUM to take free surface measurements from.
    end
    % Now modify HDEPTH, ZTopo, and MN to be in vector form.
    % local variables. Also do first node in row and last node in row.
    % initialize.
    HDEPTH = double(zeros(NUMNODES,1)); % Vector format for undisturbed water depth.
    ZTopo = double(zeros(NUMNODES,1));  % Vector format for topographic elevation at
                                        %     volume centers.
    MN = double(zeros(NUMNODES,1));     % Vector format for channel Manning's roughness
                                        %     coefficients.
    ROWBEGIN = zeros(NUMROWS,1); % Vector to hold index of beginning volume in each
                                        %     row of the domain.
    ROWEND = zeros(NUMROWS,1);   % Vector to hold the index of the ending volume in
                                        %     each row of the domain.
    UXM1 = zeros(NUMINCX,1);     % index representing i,j for i+1/2,j
    UXP1 = zeros(NUMINCX,1);     % index representing i+1,j for i+1/2,j
    UVM1 = zeros(NUMINCX,1);     % index representing i-1/2,j for i+1/2,j
    UVP1 = zeros(NUMINCX,1);     % index representing i+3/2,j for i+1/2,j
    VYP1 = zeros(NUMINCY,1);     % index representing i,j+1 for i,j vel.
    VYM1 = zeros(NUMINCY,1);     % index representing i,j-1 for i,j vel.
    VVP1 = zeros(NUMINCY,1);     % index representing i,j+3/2 for i,j+1/2 vel.
    VVM1 = zeros(NUMINCY,1);     % index representing i,j-1/2 for i,j+1/2 vel.
    X1 = zeros(NUMNODES,1);      % index representing i+1/2,j for i,j
    Y2 = zeros(NUMNODES,1);      % index representing i,j-1/2 for i,j
    X3 = zeros(NUMNODES,1);      % index representing i-1/2,j for i,j
    Y4 = zeros(NUMNODES,1);      % index representing i,j+1/2 for i,j
    XM1 = zeros(NUMNODES,1);     % index representing i-1,j for i,j
    XP1 = zeros(NUMNODES,1);     % index representing i+1,j for i,j
    YM1 = zeros(NUMNODES,1);     % index representing i-1,j for i,j
    YP1 = zeros(NUMNODES,1);     % index representing i+1,j for i,j
    XXeI = double(zeros(NUMINCX,1));    % x-coordinate x-face centers.
    XYeI = double(zeros(NUMINCX,1));    % y-coordinate x-face centers.
    XINDEX = double(zeros(1,NUMCOLS+1));% Spatial location of each volume face in the
                                        %   x-direction.
    XTemp = double(zeros(1,NUMCOLS));   % x-coordinates for y-faces.
    YTemp = double(zeros(1,NUMROWS));   % y-coordinates for x-faces.
    YXeI = double(zeros(NUMINCY,1));    % x-coordinate y-face centers.
    YYeI = double(zeros(NUMINCY,1));    % y-coordinate y-face centers.
    YINDEX = double(zeros(1,NUMROWS+1));% Spatial location of each volume face in the 
                                        %   y-direction.
    ZtX = double(zeros(NUMINCX,1));     % X-face topographic elevation.
    ZtY = double(zeros(NUMINCY,1));     % Y-face topographic elevation.

    % assignments.
    XINDEX = 0:DX:NUMCOLS*DX;
    YINDEX = 0:DY:NUMROWS*DY;
    ROWBEGIN = (1:NUMCOLS:(LASTROW+1))';
    ROWEND = (NUMCOLS:NUMCOLS:NUMNODES)';

    [X1,Y2,X3,Y4,XP1,YM1,XM1,YP1,UXP1,UXM1,UVP1,UVM1,VYP1,VYM1,VVP1,VVM1] = gENfACEiNDEX(X1,...
       Y2,X3,Y4,XP1,YM1,XM1,YP1,UXP1,UXM1,UVP1,UVM1,VYP1,VYM1,VVP1,VVM1);

    XTemp = (XINDEX(1:NUMCOLS)+0.5*DX);
    YTemp = (YINDEX(1:NUMROWS)+0.5*DY);
    YXeI(1:NUMCOLS) = XTemp';
    YYeI(1:NUMCOLS) = YINDEX(1).*ones(NUMCOLS,1);
    for i=1:NUMROWS
       ZTopo((ROWBEGIN(i):1:ROWEND(i)),1) = Topo(i,:)';
       HDEPTH((ROWBEGIN(i):1:ROWEND(i)),1) = Depth(i,:)';
       MN((ROWBEGIN(i):1:ROWEND(i)),1) = Mann(i,:)';
       XXeI(((i-1)*XINC)+1:i*XINC) = XINDEX';
       XYeI(((i-1)*XINC)+1:i*XINC) = YTemp(i).*ones(XINC,1);
       YXeI(ROWBEGIN(i)+NUMCOLS:1:ROWEND(i)+NUMCOLS) = XTemp';
       YYeI(ROWBEGIN(i)+NUMCOLS:1:ROWEND(i)+NUMCOLS) = ...
          YINDEX(i+1)*ones(NUMCOLS,1);
    end

    % Now interpolate topography values.
    [ZtX,ZtY] = sETtOPO(ZTopo,ZtX,ZtY);
    
    % output grid definitions for tracking purposes
    % create our grid definitions
    for iI = 1:NUMROWS
        for jJ = 1:NUMCOLS
            curNode = ( ( iI - 1 ) * NUMCOLS ) + jJ;
            NodesOut(iI, jJ) = curNode;
        end
    end
    for iI = 1:NUMROWS
        for jJ = 1:NUMCOLS+1
            curFace = ( ( iI - 1 ) * (NUMCOLS+1) ) + jJ;
            XFIndOut(iI, jJ) = curFace;
        end
    end
    for iI = 1:NUMROWS+1
        for jJ = 1:NUMCOLS
            curFace = ( ( iI - 1 ) * NUMCOLS ) + jJ;
            YFIndOut(iI, jJ) = curFace;
        end
    end
    % write to HDF5 file for checking
    h5create(H5FileP, '/Outputs/GridSpec/Nodes', [NUMCOLS NUMROWS], 'DataType', 'int32' );
    h5write(H5FileP, '/Outputs/GridSpec/Nodes', int32(NodesOut') );
    h5create(H5FileP, '/Outputs/GridSpec/XFaces', [NUMCOLS+1 NUMROWS], 'DataType', 'int32' );
    h5write(H5FileP, '/Outputs/GridSpec/XFaces', int32(XFIndOut') );
    h5create(H5FileP, '/Outputs/GridSpec/YFaces', [NUMCOLS NUMROWS+1], 'DataType', 'int32' );
    h5write(H5FileP, '/Outputs/GridSpec/YFaces', int32(YFIndOut') );
    % no longer needed
    %clear Tx Ty Topo Depth Mann MaxTopo;
    %clear Fid1 Fid2 Fid3 TempDepth TempMann TempTopo;
    %clear YTemp;
end
%EOF