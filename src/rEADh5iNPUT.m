function rEADh5iNPUT(H5FileP)
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rEADiNPUT reads the parameters from the HDF5 file that is used for
    %   input and output.
    %  Modified - 12/2020
    % H5FileP : HDF5 file name with path
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
    % along with MOD_FreeSurf2D. If not, see <https://www.gnu.org/licenses/>.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    global CENTRALLATITUDE COROMEGA DATUM FLUID_DT DX DY ENDTIME
    global EPSILON EVIS G GAMMATX GAMMATY HCUTOFF KAPPA MAXCR MAXITER
    global MAXSTEPS MINSTEPS MNWAL NUK NUMCOLS NUMROWS OUTINT 
    global PATHTRAC PREC PRECH PRECOND QINBC QINXFLUX QINXVOL QINYFLUX
    global QINYVOL RADFLUXBC RADORLFSBC RADVELBC RFLUXXFLUX RFLUXYFLUX
    global RORLFSXVOL RORLFSYVOL RVELXVOL RVELYVOL RHOW
    global STARTTIME TDEPDIRCBC TDEPDXDEP TDEPDXVOL TDEPDYDEP TDEPDYVOL
    global THETA UA VA VELDIRCBC VELDXVOL VELDYVOL VELDXVEL
    global VELDYVEL WINDACT NUMILOCS ILOCS
    % Initialize all of the non-string globals and all of the non-source
    %  values.
    PREC = double(0.0);
    PRECH = double(0.0);
    GAMMATX = double(0.0);
    UA = double(0.0);
    GAMMATY = double(0.0);
    VA = double(0.0);
    WINDACT = 0;  % always false
    % local variable.
    % with HDF5 file can just read the attributes from the known/required
    %  location
    STARTTIME = double( h5readatt( H5FileP, '/Inputs/', 'STARTTIME' ) );
    ENDTIME = double( h5readatt( H5FileP, '/Inputs/', 'ENDTIME' ) );
    DATUM = double( h5readatt( H5FileP, '/Inputs/', 'DATUM' ) );
    DX = double( h5readatt( H5FileP, '/Inputs/', 'DX' ) );
    DY = double( h5readatt( H5FileP, '/Inputs/', 'DY' ) );
    NUMROWS = double( h5readatt( H5FileP, '/Inputs/', 'NUMROWS' ) );
    NUMCOLS = double( h5readatt( H5FileP, '/Inputs/', 'NUMCOLS' ) );
    OUTINT = double( h5readatt( H5FileP, '/Inputs/', 'OUTINT' ) );
    FLUID_DT = double( h5readatt( H5FileP, '/Inputs/', 'FLUID_DT' ) );
    THETA = double( h5readatt( H5FileP, '/Inputs/', 'THETA' ) );
    HCUTOFF = double( h5readatt( H5FileP, '/Inputs/', 'HCUTOFF' ) );
    EPSILON = double( h5readatt( H5FileP, '/Inputs/', 'EPSILON' ) );
    MAXITER = double( h5readatt( H5FileP, '/Inputs/', 'MAXITER' ) );
    PRECOND = double( h5readatt( H5FileP, '/Inputs/', 'PRECOND' ) );
    MAXCR = double( h5readatt( H5FileP, '/Inputs/', 'MAXCR' ) );
    PATHTRAC = double( h5readatt( H5FileP, '/Inputs/', 'PATHTRAC' ) );
    MAXSTEPS = double( h5readatt( H5FileP, '/Inputs/', 'MAXSTEPS' ) );
    MINSTEPS = double( h5readatt( H5FileP, '/Inputs/', 'MINSTEPS' ) );
    G = double( h5readatt( H5FileP, '/Inputs/', 'G' ) );
    KAPPA = double( h5readatt( H5FileP, '/Inputs/', 'KAPPA' ) );
    MNWAL = double( h5readatt( H5FileP, '/Inputs/', 'MNWAL' ) );
    NUK = double( h5readatt( H5FileP, '/Inputs/', 'NUK' ) );
    RHOW = double( h5readatt( H5FileP, '/Inputs/', 'RHOW' ) );
    EVIS = double( h5readatt( H5FileP, '/Inputs/', 'EVIS' ) );
    CENTRALLATITUDE = double( h5readatt( H5FileP, '/Inputs/', 'CENTRALLATITUDE' ) );
    COROMEGA = double( h5readatt( H5FileP, '/Inputs/', 'COROMEGA' ) );
    TDEPDIRCBC = double( h5readatt( H5FileP, '/Inputs/', 'TDEPDIRCBC' ) );
    VELDIRCBC = double( h5readatt( H5FileP, '/Inputs/', 'VELDIRCBC' ) );
    QINBC = double( h5readatt( H5FileP, '/Inputs/', 'QINBC' ) );
    RADVELBC = double( h5readatt( H5FileP, '/Inputs/', 'RADVELBC' ) );
    RADORLFSBC = double( h5readatt( H5FileP, '/Inputs/', 'RADORLFSBC' ) );
    RADFLUXBC = double( h5readatt( H5FileP, '/Inputs/', 'RADFLUXBC' ) );
    RFLUXXFLUX = double( h5readatt( H5FileP, '/Inputs/', 'RFLUXXFLUX' ) );
    RFLUXYFLUX = double( h5readatt( H5FileP, '/Inputs/', 'RFLUXYFLUX' ) );
    NUMILOCS = double( h5readatt( H5FileP, '/Inputs/ILocs/', 'Num_ILoc' ) );
    % Now need to read the BC specifications
    % TDEPDIRCBC
    if true(TDEPDIRCBC)
        numTX = h5readatt( H5FileP, '/Inputs/BC/TDEPDIRCBC/Locations/', 'Num_XFace' );
        numTY = h5readatt( H5FileP, '/Inputs/BC/TDEPDIRCBC/Locations/', 'Num_YFace' );
        if ( numTX > 0 )
            TDEPDXVOL = double( h5read( H5FileP, '/Inputs/BC/TDEPDIRCBC/Locations/XLocs/' ) );
            % create our map
            TDEPDXDEP = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cXF = TDEPDXVOL(iI);
                keyStr = sprintf('/Inputs/BC/TDEPDIRCBC/Forcing/X_%d/', cXF);
                cXFTs = h5read( H5FileP, keyStr );
                TDEPDXDEP(cXF) = double( cXFTs );
            end
        else
            TDEPDXVOL = double( zeros(1) );
        end
        if ( numTY > 0 )
            TDEPDYVOL = double( h5read( H5FileP, '/Inputs/BC/TDEPDIRCBC/Locations/YLocs/' ) );
            % create our map
            TDEPDYDEP = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cYF = TDEPDYVOL(iI);
                keyStr = sprintf('/Inputs/BC/TDEPDIRCBC/Forcing/Y_%d/', cYF);
                cYFTs = h5read( H5FileP, keyStr );
                TDEPDYDEP(cXF) = double( cYFTs );
            end
        else
            TDEPDYVOL = double( zeros(1) );
        end
    end
    % VELDIRCBC
    if true(VELDIRCBC)
        numTX = h5readatt( H5FileP, '/Inputs/BC/VELDIRCBC/Locations/', 'Num_XFace' );
        numTY = h5readatt( H5FileP, '/Inputs/BC/VELDIRCBC/Locations/', 'Num_YFace' );
        if ( numTX > 0 )
            VELDXVOL = double( h5read( H5FileP, '/Inputs/BC/VELDIRCBC/Locations/XLocs/' ) );
            % create our map
            VELDXVEL = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cXF = VELDXVOL(iI);
                keyStr = sprintf('/Inputs/BC/VELDIRCBC/Forcing/X_%d/', cXF);
                cXFTs = h5read( H5FileP, keyStr );
                VELDXVEL(cXF) = double( cXFTs );
            end
        else
            VELDXVOL = double( zeros(1) );
        end
        if ( numTY > 0 )
            VELDYVOL = double( h5read( H5FileP, '/Inputs/BC/VELDIRCBC/Locations/YLocs/' ) );
            % create our map
            VELDYVEL = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cYF = VELDYVOL(iI);
                keyStr = sprintf('/Inputs/BC/VELDIRCBC/Forcing/Y_%d/', cYF);
                cYFTs = h5read( H5FileP, keyStr );
                VELDYVEL(cXF) = double( cYFTs );
            end
        else
            VELDYVOL = double( zeros(1) );
        end
    end
    % QINBC
    if true(QINBC)
        numTX = h5readatt( H5FileP, '/Inputs/BC/QINBC/Locations/', 'Num_XFace' );
        numTY = h5readatt( H5FileP, '/Inputs/BC/QINBC/Locations/', 'Num_YFace' );
        if ( numTX > 0 )
            QINXVOL = double( h5read( H5FileP, '/Inputs/BC/QINBC/Locations/XLocs/' ) );
            % create our map
            QINXFLUX = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cXF = QINXVOL(iI);
                keyStr = sprintf('/Inputs/BC/QINBC/Forcing/X_%d/', cXF);
                cXFTs = h5read( H5FileP, keyStr );
                QINXFLUX(cXF) = double( cXFTs );
            end
        else
            QINXVOL = double( zeros(1) );
        end
        if ( numTY > 0 )
            QINYVOL = double( h5read( H5FileP, '/Inputs/BC/QINBC/Locations/YLocs/' ) );
            % create our map
            QINYFLUX = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            for iI = 1:numTX
                cYF = QINYVOL(iI);
                keyStr = sprintf('/Inputs/BC/QINBC/Forcing/Y_%d/', cYF);
                cYFTs = h5read( H5FileP, keyStr );
                QINYFLUX(cXF) = double( cYFTs );
            end
        else
            QINYVOL = double( zeros(1) );
        end
    end
    % RADVELBC
    if true(RADVELBC)
        numTX = h5readatt( H5FileP, '/Inputs/BC/RADVELBC/Locations/', 'Num_XFace' );
        numTY = h5readatt( H5FileP, '/Inputs/BC/RADVELBC/Locations/', 'Num_YFace' );
        if ( numTX > 0 )
            RVELXVOL = double( h5read( H5FileP, '/Inputs/BC/RADVELBC/Locations/XLocs/' ) );
        else
            RVELXVOL = double( zeros(1) );
        end
        if ( numTY > 0 )
            RVELYVOL = double( h5read( H5FileP, '/Inputs/BC/RADVELBC/Locations/YLocs/' ) );
        else
            RVELYVOL = double( zeros(1) );
        end
    end
    % RADORLFSBC
    if true(RADORLFSBC)
        numTX = h5readatt( H5FileP, '/Inputs/BC/RADORLFSBC/Locations/', 'Num_XFace' );
        numTY = h5readatt( H5FileP, '/Inputs/BC/RADORLFSBC/Locations/', 'Num_YFace' );
        if ( numTX > 0 )
            RORLFSXVOL = double( h5read( H5FileP, '/Inputs/BC/RADORLFSBC/Locations/XLocs/' ) );
        else
            RORLFSXVOL = double( zeros(1) );
        end
        if ( numTY > 0 )
            RORLFSYVOL = double( h5read( H5FileP, '/Inputs/BC/RADORLFSBC/Locations/YLocs/' ) );
        else
            RORLFSYVOL = double( zeros(1) );
        end
    end
    % read the ILOCS if needed
    if (NUMILOCS > 0)
        ILOCS = double( h5read( H5FileP, '/Inputs/ILocs/ILocList/' ) );
    else
        ILOCS = double( zeros(1) );
    end
end
%EOF