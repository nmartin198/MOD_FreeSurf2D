function sETh5Outputs(H5FileP,NumSteps)
    %sETh5Outputs Sets up outputs to HDF5 file
    %   Outputs include output times for arrays by OUTINT as wells as 
    %   outputs by ILOC.
    %
    %   H5FileP : HDF5 file name with path
    %   NumSteps : Number of time steps in simulation
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
    
    global OUTINT FLUID_DT STARTTIME ENDTIME NUMILOCS ILOCS
    % calculate the output interval in days
    OutDur = ( ( OUTINT * FLUID_DT ) * ( 1.0 / ( 24.0 * 60.0 * 60.0 ) ) );
    NumOuts = ceil( NumSteps / OUTINT );
    OutTimes = STARTTIME : OutDur : (STARTTIME + (NumOuts*OutDur));
    TotNumOuts = size(OutTimes,2);
    OutTimes(TotNumOuts) = ENDTIME;
    % Now are ready to write out the time intervals for the outputs
    h5create(H5FileP, '/Outputs/OutTimes', [1 TotNumOuts], 'DataType', 'double' );
    h5write(H5FileP, '/Outputs/OutTimes', double(OutTimes) );
    % Next check for ILOCs 
    if (NUMILOCS > 0)
        DaysInc = FLUID_DT / (24.0 * 60.0 * 60.0);
        StartDays = STARTTIME + DaysInc;
        AllOutTimes = StartDays : DaysInc : (STARTTIME + (NumSteps*DaysInc));
        AllOutTimes(NumSteps) = ENDTIME;
        h5create(H5FileP, '/Outputs/ILOCs/OutTimes', [1 NumSteps], 'DataType', 'double' );
        h5write(H5FileP, '/Outputs/ILOCs/OutTimes', double(AllOutTimes) );
        for iI = 1:NUMILOCS
            cNode = ILOCS(iI,1);
            keyStr = sprintf('/Outputs/ILOCs/%d', cNode);
            h5create(H5FileP, keyStr, [3 NumSteps], 'DataType', 'single' );
        end
        h5writeatt(H5FileP, '/Outputs/ILOCs/', "Column_1", "U_m/s", 'TextEncoding', 'UTF-8');
        h5writeatt(H5FileP, '/Outputs/ILOCs/', "Column_2", "V_m/s", 'TextEncoding', 'UTF-8');
        h5writeatt(H5FileP, '/Outputs/ILOCs/', "Column_3", "H_m", 'TextEncoding', 'UTF-8');
    end 
end
%EOF