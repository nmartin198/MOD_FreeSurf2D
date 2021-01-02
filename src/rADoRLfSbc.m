function [FSXP,FSXM,FSYP,FSYM] = rADoRLfSbc(FSXP,FSXM,FSYP,FSYM)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rADoRLfSbc imposes an absorbing radiation boundary condition on prescribed
    % boundaries.  This radiation boundary condition comes from Orlanski (1976).
    % The free surface elevation is imposed immediately outside the simulation
    % domain.
    %
    % FSXP [NUMINCX,1] = EtaXP1 or the free surface value at (i+1,j) for i+1/2,j.
    % FSXM [NUMINCX,1] = EtaXM1 or the free surface value at (i,j) for i+1/2,j.
    % FSYP [NUMINCY,1] = EtaYP1 or the free surface value at (i,j+1) for i,j+1/2.
    % FSYM [NUMINCY,1] = EtaYM1 or the free surface value at (i,j) for i,j+1/2.
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

    global DX DY Eta EtaNew fdt Hux Hvy PRECH PREC RORLFSXVOL
    global RORLFSYVOL OrlSizeX HoldFSXB1N2 HoldFSXB1N1 HoldFSXB2N1
    global HoldFSXBN1 HoldFSXBN FSXCols FSXIndexBX FSXIndexB1 FSXIndexB2
    global OrlSizeY HoldFSYB1N2 HoldFSYB1N1 HoldFSYB2N1 HoldFSYBN1 HoldFSYBN
    global FSYIndexBY FSYIndexB1 FSYIndexB2 FSYRows

    if (RORLFSXVOL(1) ~= 0)
       % Local variables.
       BTempF = zeros(OrlSizeX,1);            %Boolean temporary variable.
       CalcF1 = double(zeros(OrlSizeX,1));    % Calc Variable.
       CalcF2 = double(zeros(OrlSizeX,1));    % Calc Variable.
       CDenom = double(zeros(OrlSizeX,1));    %Denominator of C term.
       Ce = double(zeros(OrlSizeX,1));        %C term.
       Factor = double(0.0);                  %Multiplication factor.
       FreeS = double(zeros(OrlSizeX,1));     %Calculated free surface.
       FSBN1 = double(zeros(OrlSizeX,1));     %Free surface @ boundary, n-1.
       FSB1N = double(zeros(OrlSizeX,1));     %Free surface @ boundary - 1, n.
       FSB1N2 = double(zeros(OrlSizeX,1));    %Free surface @ boundary - 1, n-2.
       FSB2N1 = double(zeros(OrlSizeX,1));    %Free surface @ boundary - 2, n-1.
       Part1 = double(zeros(OrlSizeX,1));     % Half of free surface calculation.
       Part2 = double(zeros(OrlSizeX,1));     % Half of free surface calculation.
       TempCalc1 = double(zeros(OrlSizeX,1)); % Temporary calculation variable.
       TempCalc2 = double(zeros(OrlSizeX,1)); % Temporary calculation variable.
       % Transfer variables to get the correct values for this calculation.
       FSBN1 = HoldFSXBN1;
       HoldFSXBN1 = HoldFSXBN;
       FSB1N = Eta(FSXIndexB1);
       FSB1N2 = HoldFSXB1N2;
       HoldFSXB1N2 = HoldFSXB1N1;
       HoldFSXB1N1 = Eta(FSXIndexB1);
       FSB2N1 = HoldFSXB2N1;
       HoldFSXB2N1 = Eta(FSXIndexB2);
       % Calculate the C term.
       Factor = DX/fdt;
       TempCalc1 = FSB1N + FSB1N2 - FSB2N1;
       BTempF = (abs(TempCalc1 - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       TempCalc2 = (CalcF1.*TempCalc1) + (CalcF2.*PREC);
       CDenom = 1./TempCalc2;
       TempCalc1 = FSB1N - FSB1N2;
       Ce = TempCalc1.*CDenom;
       % Ensure that the C term is within the proper bounds.
       BTempF = ((Ce - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       Ce = CalcF1.*Ce;
       BTempF = ((Factor - Ce) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       Ce = (CalcF1.*Ce) + (CalcF2.*Factor);
       % Calculate the new Free surface value.
       Factor = fdt/DX;
       TempCalc1 = 1 + (Factor.*Ce);
       BTempF = (abs(TempCalc1 - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       TempCalc2 = (CalcF1.*TempCalc1) + (CalcF2.*PREC);
       TempCalc1 = 1./TempCalc2;
       TempCalc2 = 1 - (Factor.*Ce);
       Part1 = TempCalc2.*TempCalc1;
       TempCalc2 = (2*Factor).*Ce;
       Part2 = TempCalc2.*TempCalc1;
       % Ensure that is wet.
       BTempF = ((Hux(FSXIndexBX) - double(0.0)) > PRECH);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       % Get the new free surface elevation.
       FreeS = (CalcF1.*((Part1.*FSBN1) + Part2.*FSB1N)) + ...
          (CalcF2.*EtaNew(FSXIndexB1));
       HoldFSXBN = FreeS;
       BTempF = (FSXCols > 1);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       FSXP(FSXIndexBX) = (CalcF1.*FreeS) + (CalcF2.*FSXP(FSXIndexBX));
       FSXM(FSXIndexBX) = (CalcF1.*FSXM(FSXIndexBX)) + (CalcF2.*FreeS);
       clear BTempF CalcF1 CalcF2 CDenom Ce Factor FreeS FSBN1 FSB1N FSB1N2;
       clear FSB2N1 Part1 Part2 TempCalc1 TempCalc2;
    end

    if (RORLFSYVOL(1) ~= 0)
       % Local variables.
       BTempF = zeros(OrlSizeY,1);            % see above for definitions.
       CalcF1 = double(zeros(OrlSizeX,1));    % Calc Variable.
       CalcF2 = double(zeros(OrlSizeX,1));    % Calc Variable.
       CDenom = double(zeros(OrlSizeY,1));
       Ce = double(zeros(OrlSizeY,1));
       Factor = double(0.0);
       FreeS = double(zeros(OrlSizeY,1));
       FSBN1 = double(zeros(OrlSizeY,1));
       FSB1N = double(zeros(OrlSizeY,1));
       FSB1N2 = double(zeros(OrlSizeY,1));
       FSB2N1 = double(zeros(OrlSizeY,1));
       Part1 = double(zeros(OrlSizeY,1));
       Part2 = double(zeros(OrlSizeY,1));
       TempCalc1 = double(zeros(OrlSizeY,1));
       TempCalc2 = double(zeros(OrlSizeY,1));
       % Transfer variables to get the correct values for this calculation.   
       FSBN1 = HoldFSYBN1;
       HoldFSYBN1 = HoldFSYBN;
       FSB1N = Eta(FSYIndexB1);
       FSB1N2 = HoldFSYB1N2;
       HoldFSYB1N2 = HoldFSYB1N1;
       HoldFSYB1N1 = Eta(FSYIndexB1);
       FSB2N1 = HoldFSYB2N1;
       HoldFSYB2N1 = Eta(FSYIndexB2);
       % Calculate the C term.
       Factor = DY/fdt;
       TempCalc1 = FSB1N + FSB1N2 - FSB2N1;
       BTempF = (abs(TempCalc1 - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       TempCalc2 = (CalcF1.*TempCalc1) + (CalcF2.*PREC);
       CDenom = 1./TempCalc2;
       TempCalc1 = FSB1N - FSB1N2;
       Ce = TempCalc1.*CDenom;
       % Ensure that the C term is within the proper bounds.
       BTempF = ((Ce - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       Ce = CalcF1.*Ce;
       BTempF = ((Factor - Ce) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       Ce = (CalcF1.*Ce) + (CalcF2.*Factor);
       % Calculate the new Free surface value.
       Factor = fdt/DY;
       TempCalc1 = 1 + (Factor.*Ce);
       BTempF = (abs(TempCalc1 - double(0.0)) > PREC);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       TempCalc2 = (CalcF1.*TempCalc1) + (CalcF2.*PREC);
       TempCalc1 = 1./TempCalc2;
       TempCalc2 = 1 - (Factor.*Ce);
       Part1 = TempCalc2.*TempCalc1;
       TempCalc2 = (2*Factor).*Ce;
       Part2 = TempCalc2.*TempCalc1;
       % Ensure that is wet.
       BTempF = ((Hvy(FSYIndexBY) - double(0.0)) > PRECH);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       % Get the new free surface elevation.
       FreeS = (CalcF1.*((Part1.*FSBN1) + Part2.*FSB1N)) + ...
          (CalcF2.*EtaNew(FSYIndexB1));
       HoldFSYBN = FreeS;
       BTempF = (FSYRows > 1);
       CalcF1 = double(BTempF);
       CalcF2 = double(~BTempF);
       FSYP(FSYIndexBY) = (CalcF1.*FreeS) + (CalcF2.*FSYP(FSYIndexBY));
       FSYM(FSYIndexBY) = (CalcF1.*FSYM(FSYIndexBY)) + (CalcF2.*FreeS);
       % not needed
       %clear BTempF CDenom CalcF1 CalcF2 Ce Factor FreeS FSBN1 FSB1N FSB1N2;
       %clear FSB2N1 Part1 Part2 TempCalc1 TempCalc2;
    end
end
%EOF