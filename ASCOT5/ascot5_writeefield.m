function ascot5_writeefield(a5file,rho,dvdrho,reff,varargin)
%ASCOT5_WRITEEFIELD Writes the ASCOT5 efield_1D field to an HDF file
%   The ASCOT5_WRITEEFIELD function writes the ASCOT5 efield_1D group to an
%   HDF5 file.  The code takes rho, dv/drho, and reff as inputs.  
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

eps0= 8.8541878128E-12;
ec  = 1.60217662E-19;
is1ds=1;
lfix=0;

% Handle varargin
if nargin > 10
    i = 1;
    while (i <= length(varargin))
        switch varargin{i}
            case{'1DS'}
                is1ds=1;
        end
        i=i+1;
    end
end

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Setup path
mod = '1D';
if is1ds
    mod = '1DS';
end
path = ['/efield/E_' mod '_' num2str(id,'%10.10i')];

% Handle existing file
if isfile(a5file)
    disp(['  ' a5file ' exists, adding ' path ' to file.']);
end

% Create datasets
h5create(a5file,[path '/nrho'],1,'Datatype','int32');
h5create(a5file,[path '/dvdrho'],length(dvdrho));
h5create(a5file,[path '/reff'],1);
if is1ds
    h5create(a5file,[path '/rhomin'],1);
    h5create(a5file,[path '/rhomax'],1);
end

% Write Attributes
h5writeatt(a5file,'/efield','active',id,'TextEncoding','system');
h5writeatt(a5file,path,'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
h5writeatt(a5file,path,'description','Created by ascot5_writeefield in Matlab.','TextEncoding','system');

% Write varaibles
h5write(a5file,[path '/nrho'],int32(length(rho)));
h5write(a5file,[path '/dvdrho'],dvdrho);
h5write(a5file,[path '/reff'],reff);
if is1ds
    h5write(a5file,[path '/rhomin'],min(rho));
    h5write(a5file,[path '/rhomax'],max(rho));
end

return;
end

