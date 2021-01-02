function [D1,D2,D3,D4,SRhs] = dOMsYSbOUND(D1,D2,D3,D4,SRhs,cETime)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % dOMsYSbOUND sets the domain free surface boundary terms.
    % Open boundaries will have linearized constant slope boundaries and
    % closed boundaries will have linearized constant slope == 0 boundaries.
    % For open boundaries, the actual water depth/free surface value
    % employed may be updated later, if certain boundary conditions are
    % specified.
    %
    % D1 [NUMNODES,1] = d1 or diagonal i+1,j from the system.
    % D2 [NUMNODES,1] = d2 or diagonal i,j-1 from the system.
    % D3 [NUMNODES,1] = d3 or diagonal i-1,j from the system.
    % D4 [NUMNODES,1] = d4 or diagonal i,j+1 from the system.
    % SRhs [NUMNODES,1] = rhs or the system righhand side.
    % cETime (double) = current elapsed time for time series lookup
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

    global Eta Hux Hvy PRECH LASTROW NUMCOLS NUMINCX NUMINCY NUMNODES
    global NUMROWS ROWBEGIN ROWEND TDEPDIRCBC ULASTROW XINC

    % local variables.
    Calc1A = double(zeros(NUMROWS,1)); % Calc form of Side 1.
    Calc1B = double(zeros(NUMROWS,1)); % Calc form of Side 1.
    Calc2A = double(zeros(NUMCOLS,1)); % Calc form of Side 2.
    Calc2B = double(zeros(NUMCOLS,1)); % Calc form of Side 2.
    Calc3A = double(zeros(NUMROWS,1)); % Calc form of Side 3.
    Calc3B = double(zeros(NUMROWS,1)); % Calc form of Side 3.
    Calc4A = double(zeros(NUMCOLS,1)); % Calc form of Side 4.
    Calc4B = double(zeros(NUMCOLS,1)); % Calc form of Side 4.
    CFS1 = double(zeros(NUMROWS,1));  % Constant slope eelvation east side.
    CFS2 = double(zeros(NUMCOLS,1));  % Constant slope eelvation south side.
    CFS3 = double(zeros(NUMROWS,1));  % Constant slope eelvation west side.
    CFS4 = double(zeros(NUMCOLS,1));  % Constant slope eelvation north side.
    Side1 = zeros(NUMROWS,1);         % Boolean for wet east side.
    Side2 = zeros(NUMCOLS,1);         % Boolean for wet south side.
    Side3 = zeros(NUMROWS,1);         % Boolean for wet west side.
    Side4 = zeros(NUMCOLS,1);         % Boolean for wet north side.

    % First determine which volume faces are open on each boundary.
    Side1 = ((Hux(XINC:XINC:NUMINCX) - double(0.0)) > PRECH);
    Calc1A = double(Side1);
    Calc1B = double(~Side1);
    Side2 = ((Hvy(ROWBEGIN(1):1:ROWEND(1)) - double(0.0)) > PRECH);
    Calc2A = double(Side2);
    Calc2B = double(~Side2);
    Side3 = ((Hux(1:XINC:ULASTROW+1) - double(0.0)) > PRECH);
    Calc3A = double(Side3);
    Calc3B = double(~Side3);
    Side4 = ((Hvy(ROWEND(NUMROWS)+1:1:NUMINCY) - double(0.0)) > PRECH);
    Calc4A = double(Side4);
    Calc4B = double(~Side4);
    % Calculate constant slope values.
    CFS1 = (2.*Eta(NUMCOLS:NUMCOLS:NUMNODES)) - ...
       Eta((NUMCOLS-1):NUMCOLS:(NUMNODES-1));
    CFS2 = (2.*Eta(ROWBEGIN(1):1:ROWEND(1))) - Eta(ROWBEGIN(2):1:ROWEND(2));
    CFS3 = (2.*Eta(1:NUMCOLS:LASTROW+1)) - Eta(2:NUMCOLS:LASTROW+2);
    CFS4 = (2.*Eta(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS))) - ...
       Eta(ROWBEGIN(NUMROWS-1):1:ROWEND(NUMROWS-1));
    % Adjust free surface values for specified depth boundaries.
    if (TDEPDIRCBC == 1)
       [CFS1,CFS2,CFS3,CFS4] = tdEPsYSbc(CFS1,CFS2,CFS3,CFS4,cETime);
    end
    % Allocate constant slope values to adjacent volume vectors if boundary face
    %  is wet.  Otherwise set slope to zero.
    %  East side
    SRhs(NUMCOLS:NUMCOLS:NUMNODES) = SRhs(NUMCOLS:NUMCOLS:NUMNODES) + ...
       (D1(NUMCOLS:NUMCOLS:NUMNODES).*((Calc1A.*CFS1) + (Calc1B.*...
       Eta(NUMCOLS:NUMCOLS:NUMNODES))));
    %  South side.
    SRhs(ROWBEGIN(1):1:ROWEND(1)) = SRhs(ROWBEGIN(1):1:ROWEND(1)) + ...
       (D2(ROWBEGIN(1):1:ROWEND(1)).*((Calc2A.*CFS2) + (Calc2B.*...
       Eta(ROWBEGIN(1):1:ROWEND(1)))));
    % West side.
    SRhs(1:NUMCOLS:LASTROW+1) = SRhs(1:NUMCOLS:LASTROW+1) + ...
       (D3(1:NUMCOLS:LASTROW+1).*((Calc3A.*CFS3) + (Calc3B.*...
       Eta(1:NUMCOLS:LASTROW+1))));
    % North Side.
    SRhs(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS)) = SRhs(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS)) + ...
       (D4(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS)).*((Calc4A.*CFS4) + (Calc4B.*...
       Eta(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS)))));
    % zero out off-diagonals.
    D1(NUMCOLS:NUMCOLS:NUMNODES) = 0.0;
    D2(ROWBEGIN(1):1:ROWEND(1)) = 0.0;
    D3(1:NUMCOLS:LASTROW+1) = 0.0;
    D4(ROWBEGIN(NUMROWS):1:ROWEND(NUMROWS)) = 0.0;
    % not needed
    %clear Calc1A Calc1B Calc2A Calc2B Calc3A Calc3B Calc4A Calc4B;
    %clear CFS1 CFS2 CFS3 CFS4 Side1 Side2 Side3 Side4;
end
%EOF