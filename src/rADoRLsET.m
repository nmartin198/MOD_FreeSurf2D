function rADoRLsET
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rADbcsETuP initializes the carryover values for the radiation
    % free surface boundary conditions.  This script is employed for
    % Orlanski (1976) style absorbing free surface elevtion conditions
    % and for the radiation flux boundary conditions.
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

    global EtaNew EtaXP1 EtaXM1 EtaYP1 EtaYM1 NUMCOLS RORLFSXVOL RORLFSYVOL
    global OrlSizeX HoldFSXB1N2 HoldFSXB1N1 HoldFSXB2N1 HoldFSXBN1 HoldFSXBN
    global FSXCols FSXIndexBX FSXIndexB1 FSXIndexB2
    global OrlSizeY HoldFSYB1N2 HoldFSYB1N1 HoldFSYB2N1 HoldFSYBN1 HoldFSYBN
    global FSYRows FSYIndexBY FSYIndexB1 FSYIndexB2

    if (RORLFSXVOL(1) ~= 0)
        % Global variables.
        OrlSizeX = size(RORLFSXVOL,1);
        HoldFSXB1N2 = double(zeros(OrlSizeX,1));    % Hold free surface value 1 volume in from boundary
        %                                      %     From 2 time steps ago.
        HoldFSXB1N1 = double(zeros(OrlSizeX,1));    % Hold free surface value 1 volume in from boundary
        %                                      %     From 1 time step ago.
        HoldFSXB2N1 = double(zeros(OrlSizeX,1));    % Hold free surface value 2 volumes in from boundary
        %                                      %     From 1 time step ago.
        HoldFSXBN1 = double(zeros(OrlSizeX,1));     % Hold free surface value at boundary volume from
        %                                      %     time step n-1.
        HoldFSXBN = double(zeros(OrlSizeX,1));      % Hold free surface value at boundary volume from
        %                                      %     time step n.
        FSXCols = zeros(OrlSizeX,1);        % Column index for location.
        FSXIndexBX = zeros(OrlSizeX,1);      % Volume index for boundary volume.
        FSXIndexB1 = zeros(OrlSizeX,1);     % Volume index for one volume in from boundary.
        FSXIndexB2 = zeros(OrlSizeX,1);     % Volume index for two volumes in from boundary.
        % Local Variables.
        BTemp = zeros(OrlSizeX,1);          % Boolean calc. variable.
        Calc1 = double(zeros(OrlSizeX,1));  % Calc. variable.
        Calc2 = double(zeros(OrlSizeX,1));  % Calc. variable.
        CalcI1 = zeros(OrlSizeX,1);  % Calc. variable.
        CalcI2 = zeros(OrlSizeX,1);  % Calc. variable.
        NewNode = zeros(OrlSizeX,1);        % Volume indexes for source locations.
        Row = zeros(OrlSizeX,1);            % Row index.
        TempCalc1 = double(zeros(OrlSizeX,1));      % Temporary calculation variable.
        % Calculations
        NewNode = RORLFSXVOL;
        Row = ceil(NewNode./NUMCOLS);
        TempCalc1 = mod(NewNode,NUMCOLS);
        BTemp = (TempCalc1 == 0);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        FSXCols = (CalcI1.*NUMCOLS) + (CalcI2.*TempCalc1);
        BTemp = (FSXCols > 1);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        FSXIndexBX = (CalcI1.*((Row - 1)+NewNode+1)) + (CalcI2.*((Row - 1)+NewNode));
        FSXIndexB1 = NewNode;
        FSXIndexB2 = (CalcI1.*(NewNode - 1)) + (CalcI2.*(NewNode + 1));
        Calc1 = double(CalcI1);
        Calc2 = double(CalcI2);
        HoldFSXBN1 = (Calc1.*EtaXP1(FSXIndexBX)) + (Calc2.*EtaXM1(FSXIndexBX));
        HoldFSXBN = HoldFSXBN1;
        HoldFSXB2N1 = EtaNew(FSXIndexB2);
        HoldFSXB1N1 = EtaNew(FSXIndexB1);
        HoldFSXB1N2 = EtaNew(FSXIndexB1);
        % not needed
        %clear BTemp Calc1 Calc2 CalcI1 CalcI2 NewNode Row TempCalc1;
    else
        OrlSizeX = 0;
        HoldFSXB1N2 = double(0.0);
        HoldFSXB1N1 = double(0.0);
        HoldFSXB2N1 = double(0.0);
        HoldFSXBN1 = double(0.0);
        HoldFSXBN = double(0.0);
        FSXCols = 0;
        FSXIndexBX = 0;
        FSXIndexB1 = 0;
        FSXIndexB2 = 0;
    end
    if (RORLFSYVOL(1) ~= 0)
        % Global variables.
        OrlSizeY = size(RORLFSYVOL,1);
        HoldFSYB1N2 = double(zeros(OrlSizeY,1));    % Hold free surface value 1 volume in from boundary
        %                                      %     From 2 time steps ago.
        HoldFSYB1N1 = double(zeros(OrlSizeY,1));    % Hold free surface value 1 volume in from boundary
        %                                      %     From 1 time step ago.
        HoldFSYB2N1 = double(zeros(OrlSizeY,1));    % Hold free surface value 2 volumes in from boundary
        %                                      %     From 1 time step ago.
        HoldFSYBN1 = double(zeros(OrlSizeY,1));     % Hold free surface value at boundary volume from
        %                                      %     time step n-1.
        HoldFSYBN = double(zeros(OrlSizeY,1));      % Hold free surface value at boundary volume from
        %                                      %     time step n.
        FSYRows = zeros(OrlSizeY,1);        % Row index for location.
        FSYIndexBY = zeros(OrlSizeY,1);     % Volume index for boundary volume.
        FSYIndexB1 = zeros(OrlSizeY,1);     % Volume index for one volume in from boundary.
        FSYIndexB2 = zeros(OrlSizeY,1);     % Volume index for two volumes in from boundary.
        % Local variables.
        BTemp = zeros(OrlSizeY,1);
        Calc1 = double(zeros(OrlSizeY,1));  % Calc. variable.
        Calc2 = double(zeros(OrlSizeY,1));  % Calc. variable.
        CalcI1 = zeros(OrlSizeY,1);  % Calc. variable.
        CalcI2 = zeros(OrlSizeY,1);  % Calc. variable.
        NewNode = zeros(OrlSizeY,1);        % Volume indexes for source locations.
        % Calculations.
        NewNode = RORLFSYVOL;
        FSYRows = ceil(NewNode./NUMCOLS);
        BTemp = (FSYRows > 1);
        CalcI1 = BTemp;
        CalcI2 = (1 - BTemp);
        FSYIndexBY = (CalcI1.*(NewNode+NUMCOLS)) + (CalcI2.*NewNode);
        FSYIndexB1 = NewNode;
        FSYIndexB2 = (CalcI1.*(NewNode - NUMCOLS)) + ...
                     (CalcI2.*(NewNode + NUMCOLS));
        Calc1 = double(CalcI1);
        Calc2 = double(CalcI2);
        HoldFSYBN1 = (Calc1.*EtaYP1(FSYIndexBY)) + (Calc2.*EtaYM1(FSYIndexBY));
        HoldFSYBN = HoldFSYBN1;
        HoldFSYB2N1 = EtaNew(FSYIndexB2);
        HoldFSYB1N1 = EtaNew(FSYIndexB1);
        HoldFSYB1N2 = EtaNew(FSYIndexB1);
        % not needed
        %clear BTemp Calc1 Calc2 CalcI1 CalcI2 NewNode;
    else
        OrlSizeY = 0;
        HoldFSYB1N2 = double(0.0);
        HoldFSYB1N1 = double(0.0);
        HoldFSYB2N1 = double(0.0);
        HoldFSYBN1 = double(0.0);
        HoldFSYBN = double(0.0);
        FSYIndexBY = 0;
        FSYIndexB1 = 0;
        FSYIndexB2 = 0;
        FSYRows = 0;
    end

end
%EOF