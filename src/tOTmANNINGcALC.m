function mann = tOTmANNINGcALC(mann)
% This script is a component of "MOD_FreeSurf2D: a Matlab surface
%   fluid flow model for rivers and streams."
%   by N. Martin and S. Gorelick (2004)
%
% tOTmANNINGcALC calculates the total Manning roughness coefficient (n)
% given the bottom roughness and sidewall roughness.  Potentially, extra
% roughness for floodplains could also be added here.

% mann = Manning's n calculated to take into account grain resistance and
%     or bed form resistance.

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

clear ChannelEdge Calc1 MnBottom MnWall;
return;
%EOF