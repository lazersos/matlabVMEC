function data = read_openadas22(filename)
%READ_OPENADAS22(filename) Reads an OpenADAS22 file
%   This subroutine reads an OpenADAS22 file as obtained from the website
%   http://open.adas.ac.uk/.  The routine returns a structure with the
%   following elements:
%       itv:         Target Ion Charge 
%       svref:       Stopping Coef. at ref energy, temp, density [cm^-3]
%       tsym:        Target Ion Symbol
%       nbe:         Number of beam energies
%       ntdense:     Number of target temperatures
%       ttref:       Reference target temperature [eV]
%       be:          Beam energies [eV/amu]
%       tdense:      Target Densities [cm^-3]
%       sved:        Beam stopping coef. [cm^3/s] (nbe,ntdense)
%       nttemp:      Number of target temperatures [eV]
%       beref:       Reference beam energy [eV/amu]
%       tdref:       Referecne target density [cm^-3]
%       ttemp:       Target temperatures [eV]
%       svt:         Stopping Coef. at ref energy and density [cm^3/s]
%
% Example usage
%      data=read_openadas22('bme10#h_h1.dat');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

fid = fopen(filename,'r');
if (fid <0)
    disp(' ******************************** xxdata_22 ERROR ***********************************');
    disp(['  AN ERROR OCCURRED READING FROM FILE.  FILE: ' filename]);
    disp(' ******************************* PROGRAM TERMINATED *********************************');
end
line=fgetl(fid);
temp=sscanf(line,'%5i %*7c %e');
data.itz = temp(1);
data.svref = temp(2);
temp=sscanf(line,'%*29c %2c');
data.tsym = strtrim(temp);
fgetl(fid); % dummy
line = fgetl(fid);
temp=sscanf(line,'%i %i %*6c %f');
data.nbe = temp(1);
data.ntdens = temp(2);
data.ttref  = temp(3);
fgetl(fid); % dummy
data.be=fscanf(fid,'%f',data.nbe);
data.tdens=fscanf(fid,'%f',data.ntdens);
fgetl(fid); % dummy
fgetl(fid); % dummy
data.sved = fscanf(fid,'%f',[data.ntdens data.nbe]);
fgetl(fid); % dummy
fgetl(fid); % dummy
line=fgetl(fid);
temp=sscanf(line,'%i %*6c %f %*6c %f');
data.nttemp=temp(1);
data.beref=temp(2);
data.tdref=temp(3);
fgetl(fid); % dummy
data.ttemp=fscanf(fid,'%f',data.nttemp);
fgetl(fid); % dummy
fgetl(fid); % dummy
data.svt=fscanf(fid,'%f',data.nttemp);
fclose(fid);
end


%C-----------------------------------------------------------------------
%C
%C  ***************** FORTRAN77 SUBROUTINE: xxdata_21 *******************
%C
%C  PURPOSE: TO READ DATA FROM AN EFFECTIVE BEAM STOPPING DATA SET.
%C           (ADAS FORMAT ADF21).
%C
%C  CALLING PROGRAM: SBMS / ADAS304
%C
%C  SUBROUTINE:
%C
%C  INPUT : (I*4)  IUNIT     = UNIT TO WHICH DATA SET IS CONNECTED.
%C  INPUT : (I*4)  MXBE      = MAXIMUM NUMBER OF BEAM ENERGIES WHICH CAN
%C                             BE READ.
%C  INPUT : (I*4)  MXTD      = MAXIMUM NUMBER OF TARGET DENSITIES WHICH
%C                             CAN BE READ.
%C  INPUT : (I*4)  MXTT      = MAXIMUM NUMBER OF TARGET TEMPERATURES
%C                             WHICH CAN BE READ.
%C  INPUT : (C*80) DSNIN     = NAME OF FILE TO BE READ.
%C  OUTPUT: (I*4)  ITZ       = TARGET ION CHARGE.
%C  OUTPUT: (C*2)  TSYM      = TARGET ION ELEMENT SYMBOL.
%C  OUTPUT: (R*8)  BEREF     = REFERENCE BEAM ENERGY.
%C                             UNITS: EV/AMU
%C  OUTPUT: (R*8)  TDREF     = REFERENCE TARGET DENSITY.
%C                             UNITS: CM-3
%C  OUTPUT: (R*8)  TTREF     = REFERENCE TARGET TEMPERATURE.
%C                             UNITS: EV
%C  OUTPUT: (R*8)  SVREF     = STOPPING COEFFT. AT REFERENCE BEAM ENERGY,
%C                             TARGET DENSITY AND TEMPERATURE.
%C                             UNITS: CM3 S-1
%C  OUTPUT: (I*4)  NBE       = NUMBER OF BEAM ENERGIES.
%C  OUTPUT: (R*8)  BE()      = BEAM ENERGIES.
%C                             UNITS: EV/AMU
%C                             DIMENSION: MXBE
%C  OUTPUT: (I*4)  NTDENS    = NUMBER OF TARGET DENSITIES.
%C  OUTPUT: (R*8)  TDENS()   = TARGET DENSITIES.
%C                             UNITS: CM-3
%C                             DIMENSION: MXTD
%C  OUTPUT: (I*4)  NTTEMP    = NUMBER OF TARGET TEMPERATURES.
%C  OUTPUT: (R*8)  TTEMP()   = TARGET TEMPERATURES.
%C                             UNITS: EV
%C                             DIMENSION: MXTT
%C  OUTPUT: (R*8)  SVT()     = STOPPING COEFFT. AT REFERENCE BEAM ENERGY
%C                             AND TARGET DENSITY.
%C                             UNITS: CM3 S-1
%C                             DIMENSION: MXTT
%C  OUTPUT: (R*8)  SVED(,)   = STOPPING COEFFT. AT REFERENCE TARGET
%C                             TEMPERATURE.
%C                             UNITS: CM3 S-1
%C                             1ST DIMENSION: MXBE
%C                             2ND DIMENSION: MXTD
%C
%C          (I*4)  I         = ARRAY / LOOP INDEX.
%C          (I*4)  J         = ARRAY INDEX.
%C
%C          (C*80) LINE      = TEXT LINE IN DATA SET.
%C
%C ROUTINES:
%C          ROUTINE    SOURCE    BRIEF DESCRIPTION
%C          -------------------------------------------------------------
%C          I4UNIT     ADAS      RETURNS UNIT NO. FOR OUTPUT OF MESSAGES.
%C
%C AUTHOR:  JONATHAN NASH  (TESSELLA SUPPORT SERVICES PLC)
%C          K1/0/87
%C          JET EXT. 5183
%C
%C DATE:    07/12/93
%C
%C UNIX-IDL PORT:
%C
%C VERSION: 1.1                          DATE: 16-11-95
%C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
%C               - FIRST VERSION
%C
%C-----------------------------------------------------------------------
%C
%C NOTES: Copied from c4data.for. This is v1.1 of xxdata_21.
%C 
%C
%C VERSION  : 1.1                          
%C DATE     : 06-02-2004
%C MODIFIED : Martin O'Mullane
%C              - First version
%C              - File unit is closed within the subroutine.
%C
