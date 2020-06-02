function ascot5_1dsto1d(a5file,plasmaid)
%ASCOT5_1DSTO1D Copies 1DS data to 1D
%   The ASCOT5_1DSTO1D function copies a 1DS plasma ID to a 1D plasma id
%   in an ASCOT5 file setting 1D as active.  It takes a ASCOT5 filename and
%   a plasma id number as input.  If an empty array is passed to plasmaid
%   then the currently active ID is used.
%   
%   Example:
%       ascot5_1dsto1d('ascot5.h5',[]);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00


% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Use active
if isempty(plasmaid)
    plasmaid=h5readatt(a5file,'/plasma','active');
    disp(['  Using plasma: ' plasmaid]);
end

%Pull Profile information
path = ['/plasma/plasma_1DS_' num2str(plasmaid,'%10.10i')];
try
    rho = h5read(a5file,[path '/rho']);
    nion = h5read(a5file,[path '/nion']);
    nrho = h5read(a5file,[path '/nrho']);
    znum = h5read(a5file,[path '/znum']);
    anum = h5read(a5file,[path '/anum']);
    charge = h5read(a5file,[path '/charge']);
    mass = h5read(a5file,[path '/mass']);
    ne = h5read(a5file,[path '/edensity']);
    te = h5read(a5file,[path '/etemperature']);
    ni = h5read(a5file,[path '/idensity']);
    ti = h5read(a5file,[path '/itemperature']);
catch
    disp(['ERROR: Could not find plasma: ' num2str(plasmaid,'%10.10i')]);
    return;
end

% Write to new file
ascot5_writeplasma(a5file,znum,anum,charge,mass,rho,ne,ni,te,ti)

% %WRITE to new path
% newid = num2str(round(rand*1E10),'%10.10i');
% path = ['/plasma/plasma_1D_' num2str(newid,'%10.10i')];
% 
% % Create datasets
% h5create(a5file,[path '/nrho'],1,'Datatype','int32');
% h5create(a5file,[path '/nion'],1,'Datatype','int32');
% h5create(a5file,[path '/znum'],length(znum),'Datatype','int32');
% h5create(a5file,[path '/anum'],length(anum),'Datatype','int32');
% h5create(a5file,[path '/charge'],length(charge),'Datatype','int32');
% h5create(a5file,[path '/mass'],length(mass));
% h5create(a5file,[path '/rho'],length(rho));
% h5create(a5file,[path '/edensity'],length(ne));
% h5create(a5file,[path '/idensity'],length(ni));
% h5create(a5file,[path '/etemperature'],length(te));
% h5create(a5file,[path '/itemperature'],length(ti));
% 
% % Write Attributes
% h5writeatt(a5file,'/plasma','active',newid,'TextEncoding','system');
% h5writeatt(a5file,path,'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
% h5writeatt(a5file,path,'description',...
%     ['Copied from plasma_1DS_' num2str(plasmaid,'%10.10i') ' in MATLAB'],'TextEncoding','system');
% 
% % Write varaibles
% h5write(a5file,[path '/nrho'],int32(nrho));
% h5write(a5file,[path '/nion'],int32(nion));
% h5write(a5file,[path '/znum'],int32(znum));
% h5write(a5file,[path '/anum'],int32(anum));
% h5write(a5file,[path '/charge'],int32(charge));
% h5write(a5file,[path '/mass'],mass);
% h5write(a5file,[path '/rho'],rho);
% h5write(a5file,[path '/edensity'],ne);
% h5write(a5file,[path '/idensity'],ni);
% h5write(a5file,[path '/etemperature'],te);
% h5write(a5file,[path '/itemperature'],ti);

return;


end

