function [FTrans,rhs] = lhscALC(FTrans,rhs,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % lhscALC creates the matrix for the LHS of the system of equations for the free
    % surface solution.  Also normalizes the values on both the LHS and RHS to 
    % ensure a symmetric positive definite system.
    %
    % rhs = qRHS [NUMNODES,1] not yet normalized.
    % FTrans = EToEta [NUMNODES,1] which is the normalizing factor.
    %
    % Returned:
    % FTrans = EToEta [NUMNODES,1] which is the normalizing factor.
    % also returns sLHS as a global to save storage.
    % rhs = qRHS [NUMNODES,1] which is now in normalized form.
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

    global aXDen aYDen DX DY fdt G Hux Hvy NUMCOLS NUMINCX
    global NUMINCY NUMNODES PREC sLHS
    global THETA

    %local variables.  The d* are for the main diagonal (d) calculations.
    %  The s* will be the auxiliary diagonals.
    BTemp = zeros(NUMNODES,1);             % Boolean temporary for centers.
    Calc1 = double(zeros(NUMNODES,1));     % Calc variable.
    Calc2 = double(zeros(NUMNODES,1));     % Calc variable.
    d = double(zeros(NUMNODES,1));         % not-normalized main diagonal.
    d1 = double(zeros(NUMNODES,1));        % s_{i+1,j}
    d2 = double(zeros(NUMNODES,1));        % s_{i,j-1}
    d3 = double(zeros(NUMNODES,1));        % s_{i-1,j}
    d4 = double(zeros(NUMNODES,1));        % s_{i,j+1}
    dfor1 = double(zeros(NUMNODES,1));     % 1/sqrt(d_{i,j}*d_{i+1,j})
    dfor2 = double(zeros(NUMNODES,1));     % 1/sqrt(d_{i,j}*d_{i,j-1})
    dfor3 = double(zeros(NUMNODES,1));     % 1/sqrt(d_{i,j}*d_{i-1,j})
    dfor4 = double(zeros(NUMNODES,1));     % 1/sqrt(d_{i,j}*d_{i,j+1})
    DSideX = double(zeros(NUMINCX,1));     % d-value for x-faces.
    DSideY = double(zeros(NUMINCY,1));     % d-value for y-faces.
    NMainD = double(ones(NUMNODES,1));     % Normalized main diagonal (is simply 1s).
    s1 = double(zeros(NUMNODES,1));        % a_{i+1,j}
    s3 = double(zeros(NUMNODES,1));        % a_{i-1,j}
    s2 = double(zeros(NUMNODES,1));        % a_{i,j-1}
    s4 = double(zeros(NUMNODES,1));        % a_{i,j+1}
    SqD = double(zeros(NUMNODES,1));       % Sqrt of d.
    Temp = double(ones(NUMNODES,1));       % Temporary calculation variable.
    Tempa = double(zeros(NUMNODES,1));     % Temporary calculation variable.
    TempX = double(zeros(NUMINCX,1));      % Temporary calculation variable.
    TempXa = double(zeros(NUMINCX,1));     % Temporary calculation variable.
    TempXh = double(zeros(NUMINCX,1));     % Temporary calculation variable.
    TempY = double(zeros(NUMINCY,1));      % Temporary calculation variable.
    TempYh = double(zeros(NUMINCY,1));     % Temporary calculation variable.
    XMult = double(0.0);                   % X-face multiplier.
    YMult = double(0.0);                   % Y-face multiplier.

    % Calculations.
    % X-face calculations.
    XMult = G*(THETA^2)*((fdt^2)/(DX^2));
    TempXh = Hux.^2;
    TempX = TempXh.*aXDen;
    DSideX = XMult.*TempX;
    % Y-face calculations.
    YMult = G*(THETA^2)*((fdt^2)/(DY^2));
    TempYh = Hvy.^2;
    TempY = TempYh.*aYDen;
    DSideY = YMult.*TempY;
    % split DSideX and DSideY into d1, d2, d3 ,d4.  The d1,d2,d3,d4 variables are
    %     equivalent to the s variables in Casulli and Cheng (1992).
    [d1,d2,d3,d4] = aLLOCfACE(DSideX,DSideY,d1,d2,d3,d4);
    % Main diagonal calculation.
    % Get the d value which is the un-normalized main diagonal from
    %     Casulli and Cheng (1992).
    d = Temp + d1 + d2 + d3 + d4;
    % adjust for domain boundaries.
    [d1,d2,d3,d4,rhs] = dOMsYSbOUND(d1,d2,d3,d4,rhs,cETime);
    % Now start the normalization process.
    SqD = sqrt(d);
    BTemp = ((SqD - double(0.0)) > PREC);
    Calc1 = double(BTemp);
    Calc2 = double(~BTemp);
    Temp = (Calc1.*SqD) + (Calc2.*PREC);
    Tempa = Calc1.*(1./Temp);
    % Calculate multiplier to get back "real" values from normalized system.
    % The multiplier is just 1/sqrt(d).
    FTrans = 1.*Tempa;
    % Calculate normalized RHS.
    rhs = rhs.*FTrans;
    % next need to map the d vector to values representing the adjacent 
    % nodes.  adjacent nodes are represented by the dfor* vectors.  The
    % dfor1 vector would be 1/sqrt(d_{i,j}*d_{i+1,j}).
    [dfor1,dfor2,dfor3,dfor4] = nORMd(d,dfor1,dfor2,dfor3,dfor4);
    % Calculate normalized vectors which will be normalized diagonals of the LHS.
    % These are the a* values from Casulli and Cheng (1992).
    d1 = d1.*dfor1;
    d2 = d2.*dfor2;
    d3 = d3.*dfor3;
    d4 = d4.*dfor4;
    % d1 is now actually a1.  s1/sqrt(d_{i,j}*d_{i+1,j});
    % shift for the matlab matrix creation routine.
    s1(2:NUMNODES,1) = d1(1:NUMNODES-1,1);
    s1(1,1) = 0.0;
    s2(1:(NUMNODES-NUMCOLS),1) = d2((NUMCOLS+1):NUMNODES,1);
    s2(((NUMNODES-NUMCOLS)+1):NUMNODES,1) = 0.0;
    s3(1:NUMNODES-1,1) = d3(2:NUMNODES,1);
    s3(NUMNODES,1) = 0.0;
    s4((NUMCOLS+1):NUMNODES,1) = d4(1:(NUMNODES-NUMCOLS),1);
    s4(1:NUMCOLS,1) = 0.0;
    % make negative because will be - a*.
    s1 = -1.*s1;
    s2 = -1.*s2;
    s3 = -1.*s3;
    s4 = -1.*s4;
    % now set up the LHS matrix.  will be passed back to the main program for solving.
    dgs = [-NUMCOLS -1 0 1 NUMCOLS];
    sLHS = spdiags([s2 s3 NMainD s1 s4],dgs,NUMNODES,NUMNODES);
    % not needed
    %clear BTemp Calc1 Calc2 d d1 d2 d3 d4 dgs dfor1 dfor2 dfor3 dfor4 DSideX;
    %clear DSideY NMainD s1 s3 s2 s4 SqD Temp Tempa TempX TempXa TempXh;
    %clear TempY TempYh XMult YMult;
end
%EOF