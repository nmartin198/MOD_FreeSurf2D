function mann = tOTmANNINGcALC(mann)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % tOTmANNINGcALC calculates the total Manning roughness coefficient (n)
    % given the bottom roughness and sidewall roughness.  Potentially, extra
    % roughness for floodplains could also be added here.
    %
    % mann = Manning's n calculated to take into account grain resistance and
    %     or bed form resistance.
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

    global MNWAL NUMNODES Sides

    % local variables.
    ChannelEdge = zeros(NUMNODES,1);       % Boolean for if volume is a channel edge.
    Calc1 = double(zeros(NUMNODES,1));     % Calc. variable.
    MnBottom = double(zeros(NUMNODES,1));  % Manning friction factor due to bottom.
    MnWall = double(zeros(NUMNODES,1));    % Manning factor due to walls.

    % calculate Mn for each volume center.  Use Sides to determine how many
    % walls to employ in resistance calculation.  Assuming that the energy
    % slope is constant and divide the hydraulic radius according to resistance
    % components.
    MnBottom = mann;
    % Now calculate walls.
    MnWall = MNWAL.*ones(NUMNODES,1);
    ChannelEdge = ((Sides >= 1) & (Sides < 4));
    Calc1 = double(ChannelEdge);
    MnWall = Calc1.*MnWall;
    mann = MnBottom + MnWall;
    % not needed
    %clear ChannelEdge Calc1 MnBottom MnWall;
end
%EOF