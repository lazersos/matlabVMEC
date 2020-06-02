function ascot5_writebfield_sts(filename,raxis,phiaxis,zaxis,br,bphi,bz,psi,phiedge,nfp,r0,z0)
%ASCOT5_WRITEBFIELD_STS Writes B_STS data to an hdf5 file
%   The ASCOT5_WRITEBFIELD_STS routine writes B_STS_ data to an HDF5 file
%   in the style read by the ASCOT5 code.  If the file exists the
%   /bfield/B_STS_ dataset is added to the file and set as active.
%   Please note that toroidal grids supplied to this routine should not
%   inlucde the symmetric point at 2*pi/nfp.  So if the grid goes from
%   [0,2*pi/nfp] then supply [0,nphi-1] toroidal slices to the routine for
%   br, bphi, bz, psi, phiaxis, r0 and z0.  The r0 and z0 are the axis
%   locations in each toroidal slice.  Also note that psi should be
%   supplied in terms of normalized toroidal flux.
%
%   Example:
%       ascot5_writemarker_fl('test.h5',raxis,phiaxis,zaxis,br,bphi,bz,psi,phiedge,nfp,r0,z0);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

%Constants

% Check data
phi_end=phiaxis(end);
if mod(pi*2,phi_end) == 0
    disp('WARNING: phiaxis(end) = 2*pi/nfp.')
    disp('         ASCOT GRID requires phi to extend to nfpphi-1');
    return;
end
phimin=rad2deg(phiaxis(1));
phimax=rad2deg(phiaxis(end));

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Handle existing file
if isfile(filename)
    disp(['  ' filename ' exists, adding /bfield/B_STS_' id ' to file.']);
end

% Create datasets
h5create(filename,['/bfield/B_STS_' id '/b_nr'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/b_nphi'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/b_nz'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/psi_nr'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/psi_nphi'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/psi_nz'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/axis_nphi'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/toroidalPeriods'],1,'Datatype','int32');
h5create(filename,['/bfield/B_STS_' id '/b_rmin'],1);
h5create(filename,['/bfield/B_STS_' id '/b_rmax'],1);
h5create(filename,['/bfield/B_STS_' id '/b_phimin'],1);
h5create(filename,['/bfield/B_STS_' id '/b_phimax'],1);
h5create(filename,['/bfield/B_STS_' id '/b_zmin'],1);
h5create(filename,['/bfield/B_STS_' id '/b_zmax'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_rmin'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_rmax'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_phimin'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_phimax'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_zmin'],1);
h5create(filename,['/bfield/B_STS_' id '/psi_zmax'],1);
h5create(filename,['/bfield/B_STS_' id '/axis_phimin'],1);
h5create(filename,['/bfield/B_STS_' id '/axis_phimax'],1);
h5create(filename,['/bfield/B_STS_' id '/psi0'],1);
h5create(filename,['/bfield/B_STS_' id '/psi1'],1);
h5create(filename,['/bfield/B_STS_' id '/axisr'],size(r0));
h5create(filename,['/bfield/B_STS_' id '/axisz'],size(z0));
h5create(filename,['/bfield/B_STS_' id '/br'],size(br));
h5create(filename,['/bfield/B_STS_' id '/bphi'],size(bphi));
h5create(filename,['/bfield/B_STS_' id '/bz'],size(bz));
h5create(filename,['/bfield/B_STS_' id '/psi'],size(psi));

% Write Attributes
h5writeatt(filename,'/bfield','active',id,'TextEncoding','system');
h5writeatt(filename,['/bfield/B_STS_' id ],'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
h5writeatt(filename,['/bfield/B_STS_' id ],'description','Written by MATLAB','TextEncoding','system');

% Write Variables
h5write(filename,['/bfield/B_STS_' id '/b_nr'],int32(length(raxis)));
h5write(filename,['/bfield/B_STS_' id '/b_nphi'],int32(length(phiaxis)));
h5write(filename,['/bfield/B_STS_' id '/b_nz'],int32(length(zaxis)));
h5write(filename,['/bfield/B_STS_' id '/psi_nr'],int32(length(raxis)));
h5write(filename,['/bfield/B_STS_' id '/psi_nphi'],int32(length(phiaxis)));
h5write(filename,['/bfield/B_STS_' id '/psi_nz'],int32(length(zaxis)));
h5write(filename,['/bfield/B_STS_' id '/axis_nphi'],int32(length(phiaxis)));
h5write(filename,['/bfield/B_STS_' id '/toroidalPeriods'],int32(nfp));
h5write(filename,['/bfield/B_STS_' id '/b_rmin'],raxis(1));
h5write(filename,['/bfield/B_STS_' id '/b_rmax'],raxis(end));
h5write(filename,['/bfield/B_STS_' id '/b_phimin'],phimin);
h5write(filename,['/bfield/B_STS_' id '/b_phimax'],phimax);
h5write(filename,['/bfield/B_STS_' id '/b_zmin'],zaxis(1));
h5write(filename,['/bfield/B_STS_' id '/b_zmax'],zaxis(end));
h5write(filename,['/bfield/B_STS_' id '/psi_rmin'],raxis(1));
h5write(filename,['/bfield/B_STS_' id '/psi_rmax'],raxis(end));
h5write(filename,['/bfield/B_STS_' id '/psi_phimin'],phimin);
h5write(filename,['/bfield/B_STS_' id '/psi_phimax'],phimax);
h5write(filename,['/bfield/B_STS_' id '/psi_zmin'],zaxis(1));
h5write(filename,['/bfield/B_STS_' id '/psi_zmax'],zaxis(end));
h5write(filename,['/bfield/B_STS_' id '/psi0'],-1.0E-3);
h5write(filename,['/bfield/B_STS_' id '/psi1'],phiedge);
h5write(filename,['/bfield/B_STS_' id '/axis_phimin'],phimin);
h5write(filename,['/bfield/B_STS_' id '/axis_phimax'],phimax);
h5write(filename,['/bfield/B_STS_' id '/axisr'],r0);
h5write(filename,['/bfield/B_STS_' id '/axisz'],z0);
h5write(filename,['/bfield/B_STS_' id '/br'],br);
h5write(filename,['/bfield/B_STS_' id '/bphi'],bphi);
h5write(filename,['/bfield/B_STS_' id '/bz'],bz);
h5write(filename,['/bfield/B_STS_' id '/psi'],psi.*phiedge);


end


