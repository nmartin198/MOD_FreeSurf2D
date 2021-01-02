function FMass = mASSfLUXcALC(FMass)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % mASSfLUXcALC calculates the mass leaving the domain during a time step.
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

    global DX DY fdt Hux Hvy NUMCOLS NUMINCX NUMINCY NUMROWS RHOW
    global ROWBEGIN ROWEND u ULASTROW v XINC

    % local variables.
    Temp1 = double(zeros(NUMROWS,1));
    Temp2 = double(zeros(NUMCOLS,1));
    Temp3 = double(zeros(NUMROWS,1));
    Temp4 = double(zeros(NUMCOLS,1));
    Tempu = double(zeros(NUMROWS,1));
    Tempuh = double(zeros(1,NUMROWS));
    Tempv = double(zeros(NUMCOLS,1));
    Tempvh = double(zeros(1,NUMCOLS));

    % calcs.
    % u3 side.
    Temp3 = u(1:XINC:ULASTROW+1);
    Tempu = fdt.*Temp3;
    Temp3 = Hux(1:XINC:ULASTROW+1)';
    Tempuh = DY.*Temp3;
    FMass = FMass + ((Tempuh*Tempu)*RHOW);
    %u1 side.
    Temp1 = u(NUMCOLS+1:XINC:NUMINCX);
    Tempu = fdt.*Temp1;
    Temp1 = Hux(NUMCOLS+1:XINC:NUMINCX)';
    Tempuh = DY.*Temp1;
    FMass = FMass - ((Tempuh*Tempu)*RHOW);
    % v2 side.
    Temp2 = v(ROWBEGIN(1):1:ROWEND(1));
    Tempv = fdt.*Temp2;
    Temp2 = Hvy(ROWBEGIN(1):1:ROWEND(1))';
    Tempvh = DX.*Temp2;
    FMass = FMass + ((Tempvh*Tempv)*RHOW);
    % v4 side.
    Temp4 = v(ROWEND(NUMROWS)+1:1:NUMINCY);
    Tempv = fdt.*Temp4;
    Temp4 = Hvy(ROWEND(NUMROWS)+1:1:NUMINCY)';
    Tempvh = DX.*Temp4;
    FMass = FMass - ((Tempvh*Tempv)*RHOW);
    % not needed
    %clear Tempuh Tempu Tempvh Tempv;
    %clear Temp1 Temp2 Temp3 Temp4;
end
%EOF