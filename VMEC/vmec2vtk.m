function vmec2vtk(vmec_data,nu,nv)
%VMEC2VTK Routine outputs VMEC |B| on VTK format.
%   The VMEC2VTK routine outputs the VMEC |B| data in the VTK format.  The
%   routine takes a VMEC data structure as returned by READ_VMEC, the
%   number of poloidal points (nu), and number of toroidal points (nv).
%   The code transforms the equilibrium over a field period.
%
% Example usage
%      vmec_data=read_vmec('wout_test.nc');
%      vmec2vtk(vmec_data,64,32);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Grid helpers
nu1=nu-1;
nv1=nv-1;
ns = vmec_data.ns;
ns1 = ns-1;
du = 1./nu1;
dv = 1./nv1;
ds = 1./ns1;

% Poloidal/Toroidal Grid
theta = (0:du:1).*2.*pi;
zeta  = (0:dv:1).*2.*pi./vmec_data.nfp;

% Transform quantities
r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
b = cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
if vmec_data.iasym==1
    r = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    b = b + sfunct(theta,zeta,vmec_data.bmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
end

% Create x and y
zeta3d = shiftdim(repmat(zeta',[1 ns nu]),1);
x = r.*cos(zeta3d);
y = r.*sin(zeta3d);

% Create vtk file
filename='vmec.vtk';
fid = fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,['VMEC by vmec2vtk\n']);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i float\n',ns*nu*nv);
out = [x(:) y(:) z(:)]';
fprintf(fid,'%7.4f %7.4f %7.4f\n',out);
ncells = ns1*nu1*nv1;
fprintf(fid,'CELLS %i %i\n',ncells,9*ncells);
modb = zeros(1,ncells);
cells = zeros(9,ncells);
cells(1,:) = 8; % 8 points
nsu = ns*nu;
m=1;
for k = 1:nv1
    for j = 1:nu1
        for i = 1:ns1
            n = (k-1)*nsu+(j-1)*ns+i-1;
            cells(2:9,m) = [n; n+1; n+1+nsu; n+nsu; n+ns; n+ns+1; n+ns+1+nsu; n+ns+nsu];
            modb(m) = mean(b([i i+1],[j j+1],[k k+1]),'all');
            m = m+1;
        end
    end
end
fprintf(fid,'%i %i %i %i %i %i %i %i %i\n',cells);
fprintf(fid,'CELL_TYPES %i\n',ncells);
cell_types = ones(1,ncells).*12; % VTK Hexahedron
fprintf(fid,'%i\n',cell_types);
fprintf(fid,'CELL_DATA %i\n',ncells);
fprintf(fid,'SCALARS mod(B) float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%7.4f\n',modb);
fclose(fid);



return;
end

