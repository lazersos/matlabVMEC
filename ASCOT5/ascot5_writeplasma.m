function ascot5_writeplasma(a5file,znum,anum,charge,mass,rho,ne,ni,te,ti,varargin)
%ASCOT5_WRITEPLASMA Writes the ASCOT5 plasma_1D field to an HDF file
%   The ASCOT5_WRITEPLASMA function writes the ASCOT5 plasma_1D group to an
%   HDF5 file.  The code takes the integer charge state, integer Atomic
%   mass, integer charge number, species mass, the rho coordiante, electron
%   density, ion density, electron temperature, and ion temperature as
%   inputs.  
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

eps0= 8.8541878128E-12;
ec  = 1.60217662E-19;
is1ds=0;
lfix=0;

% Handle varargin
if nargin > 10
    i = 1;
    while (i <= length(varargin))
        switch varargin{i}
            case{'1DS'}
                is1ds=1;
            case{'fix_profs'}
                lfix=1;
        end
        i=i+1;
    end
end

% Do a check
%if lcheck_profs
%    tempe = ne.*ec.*ec./(te.*ec);
%    tempi = ni.*charge.*charge.*ec.*ec./(ti.*ec);
%    debye = sqrt(eps0./(tempe+tempi));
%end

% Fix profs
if lfix
    disp('  Correcting minimum value of profs');
    te(te<10) = 10;
    ti(ti<10) = 10;
    ne(ne<1E18) = 1E18;
    ni(ni<1D18) = 1E18;
end

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Setup path
mod = '1D';
if is1ds
    mod = '1DS';
end
path = ['/plasma/plasma_' mod '_' num2str(id,'%10.10i')];

% Handle existing file
if isfile(a5file)
    disp(['  ' a5file ' exists, adding ' path ' to file.']);
end

% Create datasets
h5create(a5file,[path '/nrho'],1,'Datatype','int32');
h5create(a5file,[path '/nion'],1,'Datatype','int32');
h5create(a5file,[path '/znum'],length(znum),'Datatype','int32');
h5create(a5file,[path '/anum'],length(anum),'Datatype','int32');
h5create(a5file,[path '/charge'],length(charge),'Datatype','int32');
h5create(a5file,[path '/mass'],length(mass));
h5create(a5file,[path '/rho'],length(rho));
h5create(a5file,[path '/edensity'],length(ne));
h5create(a5file,[path '/idensity'],length(ni));
h5create(a5file,[path '/etemperature'],length(te));
h5create(a5file,[path '/itemperature'],length(ti));
if is1ds
    h5create(a5file,[path '/rhomin'],1);
    h5create(a5file,[path '/rhomax'],1);
end

% Write Attributes
h5writeatt(a5file,'/plasma','active',id,'TextEncoding','system');
h5writeatt(a5file,path,'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
h5writeatt(a5file,path,'description','Created by ascot5_writeplasma in Matlab.','TextEncoding','system');

% Write varaibles
h5write(a5file,[path '/nrho'],int32(length(rho)));
h5write(a5file,[path '/nion'],int32(1));
h5write(a5file,[path '/znum'],int32(znum));
h5write(a5file,[path '/anum'],int32(anum));
h5write(a5file,[path '/charge'],int32(charge));
h5write(a5file,[path '/mass'],mass);
h5write(a5file,[path '/rho'],rho);
h5write(a5file,[path '/edensity'],ne);
h5write(a5file,[path '/idensity'],ni);
h5write(a5file,[path '/etemperature'],te);
h5write(a5file,[path '/itemperature'],ti);
if is1ds
    h5write(a5file,[path '/rhomin'],min(rho));
    h5write(a5file,[path '/rhomax'],max(rho));
end

return;
end

