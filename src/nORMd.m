function [nd1,nd2,nd3,nd4] = nORMd(OrgVal,nd1,nd2,nd3,nd4)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % nORMd creates the normalized denominator values which will be employed
    % to normalize the LHS of the free surface system of equations.  The
    % procedure for normalization comes directly from Casulli and Cheng (1992).
    %
    % Received:
    %
    % OrgVal is the value of d from lhscALC.m
    % nd1 [NUMNODES,1] = normalized denominator for i+1,j.
    % nd2 [NUMNODES,1] = normalized denominator for i,j-1.
    % nd3 [NUMNODES,1] = normalized denominator for i-1,j.
    % nd4 [NUMNODES,1] = normalized denominator for i,j+1.
    %
    % Returned:
    %
    % nd1 [NUMNODES,1] = normalized denominator for i+1,j.
    % nd2 [NUMNODES,1] = normalized denominator for i,j-1.
    % nd3 [NUMNODES,1] = normalized denominator for i-1,j.
    % nd4 [NUMNODES,1] = normalized denominator for i,j+1.
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

    global NUMNODES PREC XP1 XM1 YM1 YP1

    %local variables.
    BTemp = zeros(NUMNODES,1);             % Boolean temp calculation variable.
    Calc1 = double(zeros(NUMNODES,1));     % Calc variable.
    Calc2 = double(zeros(NUMNODES,1));     % Calc variable.
    Temp = double(zeros(NUMNODES,1));      % Temp calculation variable.
    Temp1 = double(zeros(NUMNODES,1));     % Temp calculation variable.
    Temp2 = double(zeros(NUMNODES,1));     % Temp calculation variable.
    Temp3 = double(zeros(NUMNODES,1));     % Temp calculation variable.
    Temp4 = double(zeros(NUMNODES,1));     % Temp calculation variable.
    Tempa = double(zeros(NUMNODES,1));     % Temp calculation variable.
    Tempd = double(zeros(NUMNODES,1));     % Temp calculation variable.

    % Calculations.
    Temp1 = OrgVal(XP1);
    Temp2 = OrgVal(YM1);
    Temp3 = OrgVal(XM1);
    Temp4 = OrgVal(YP1);
    % First do (i+1,j);
    Tempa = OrgVal.*Temp1;
    Temp = sqrt(Tempa);
    BTemp = ((Temp - double(0.0)) > PREC);
    Calc1 = double(BTemp);
    Calc2 = double(~BTemp);
    Tempa = (Calc1.*Temp) + (Calc2.*PREC);
    nd1 = Calc1.*(1./Tempa);
    % Then (i,j-1);
    Tempa = OrgVal.*Temp2;
    Temp = sqrt(Tempa);
    BTemp = ((Temp - double(0.0)) > PREC);
    Calc1 = double(BTemp);
    Calc2 = double(~BTemp);
    Tempa = (Calc1.*Temp) + (Calc2.*PREC);
    nd2 = Calc1.*(1./Tempa);
    % Then (i-1,j);
    Tempa = OrgVal.*Temp3;
    Temp = sqrt(Tempa);
    BTemp = ((Temp - double(0.0)) > PREC);
    Calc1 = double(BTemp);
    Calc2 = double(~BTemp);
    Tempa = (Calc1.*Temp) + (Calc2.*PREC);
    nd3 = Calc1.*(1./Tempa);
    % Finally (i,j+1);
    Tempa = OrgVal.*Temp4;
    Temp = sqrt(Tempa);
    BTemp = ((Temp - double(0.0)) > PREC);
    Calc1 = double(BTemp);
    Calc2 = double(~BTemp);
    Tempa = (Calc1.*Temp) + (Calc2.*PREC);
    nd4 = Calc1.*(1./Tempa);
    % not needed
    %clear BTemp Calc1 Calc2 Temp Temp1 Temp2 Temp3 Temp4 Tempd Tempa;
    %clear OrgVal;
end
%EOF