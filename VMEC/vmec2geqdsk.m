function vmec2geqdsk( varargin )
%VMEC2GEQDSK Writes a WOUT data structure to GEQDSK File
% This funciton writes a GEQDSK file from a VMEC wout file.  The user must
% pass a VMEC data structure.  Options include passing of a vessel data
% structure and file name to output the data.  Note that axisymmetry is
% assumed non-axisymmetric equilibria will have their cross section output
% at the phi=0 plane.  
%
% Example usage
%      vmec_data=read_vmec('wout.test');     % Reads VMEC wout file
%      ves_data = read_vessel('vessel.dat'); % Read a vessel file
%      vmec2geqdsk(vmec_data,ves_data,'gtemp.0000');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0


vmec_data=[];
ves_data=[];
for i=1:nargin-1
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch upper(varargin{i}.datatype)
                case 'WOUT'
                    vmec_data=varargin{i};
                case 'VESSEL'
                    ves_data=varargin{i};
            end
        end
    end
    if ischar(varargin{i})
        filename = varargin{i};
    end
end

if isempty(vmec_data)
    return;
end

filename='junk.0000';
nx     = vmec_data.ns;
nz     = vmec_data.ns;
ntheta = 360;
theta = 0:2*pi/(ntheta-1):2*pi;
zeta  = 0;

% Transform the data
r=cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z=sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r=r+sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z=z+cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
phin = vmec_data.phi./vmec_data.phi(end);
for i=1:vmec_data.ns
    chirz(i,1:ntheta,1) = vmec_data.chif(i);
end
r2 = reshape(r,[prod(size(r)) 1]);
z2 = reshape(z,[prod(size(z)) 1]);
chirz = reshape(chirz,[prod(size(chirz)) 1]);
Fchi = TriScatteredInterp(r2,z2,chirz);
r1   = 0.9*vmec_data.rminsurf;
r2   = 1.1*vmec_data.rmaxsurf;
z1   = 1.1*min(z(end,:,1));
z2   = 1.1*max(z(end,:,1));
[xgrid,zgrid]=meshgrid(r1:(r2-r1)./(nx-1):r2,z1:(z2-z1)./(nz-1):z2);
psixz = Fchi(xgrid,zgrid);
psixz(isnan(psixz)) = 0.0;

% Write the file
fid = fopen(filename,'w');
fprintf(fid,'%-10s','  VMEC');
fprintf(fid,'%10s',datestr(now,'mm/dd/yyyy'));
fprintf(fid,'%10s','  #000001');
fprintf(fid,'%10s','  0000ms');
fprintf(fid,'%4d %4d %4d\n',1,nx,nz);
fprintf(fid,'% 15.9E',vmec_data.rmaxsurf.*1.1);
fprintf(fid,'% 15.9E',vmec_data.zmaxsurf.*1.1);
fprintf(fid,'% 15.9E',vmec_data.rminsurf.*0.9);
fprintf(fid,'% 15.9E',vmec_data.rmaxsurf);
fprintf(fid,'% 15.9E\n',0.0);

fprintf(fid,'% 15.9E',r(1,1,1));
fprintf(fid,'% 15.9E',z(1,1,1));
fprintf(fid,'% 15.9E',vmec_data.chif(1));
fprintf(fid,'% 15.9E',vmec_data.chif(end));
fprintf(fid,'% 15.9E\n',vmec_data.b0);

fprintf(fid,'% 15.9E',vmec_data.Itor);
fprintf(fid,'% 15.9E',vmec_data.phi(end));
fprintf(fid,'% 15.9E',0.0);
fprintf(fid,'% 15.9E',r(1,1,1));
fprintf(fid,'% 15.9E\n',0.0);

fprintf(fid,'% 15.9E',z(1,1,1));
fprintf(fid,'% 15.9E',0.0);
fprintf(fid,'% 15.9E',vmec_data.phi(end));
fprintf(fid,'% 15.9E',r(end,1,1));
fprintf(fid,'% 15.9E\n',z(end,1,1));

fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',pchip(phin,vmec_data.phi,0:1/(nx-1):1)); %sf
fprintf(fid,'\n');
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',pchip(phin,vmec_data.presf,0:1/(nx-1):1)); %sp
fprintf(fid,'\n');
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',pchip(phin,vmec_data.presf,0:1/(nx-1):1)); %sffp
fprintf(fid,'\n');
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',pchip(phin,gradient(vmec_data.presf),0:1/(nx-1):1)); %spp
fprintf(fid,'\n');


fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',psixz); %psirz
fprintf(fid,'\n');
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',pchip(phin,1./vmec_data.iotaf,0:1/(nx-1):1)); %sf
fprintf(fid,'\n');

fprintf(fid,'% 15.9E',ntheta);
fprintf(fid,'% 15.9E\n',ntheta);
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',r(end,:,1)); %rbdry
fprintf(fid,'\n');
fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',z(end,:,1)); %zbdry
fprintf(fid,'\n');
if isempty(ves_data)
    fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',r(end,:,1)); %rlim
    fprintf(fid,'\n');
    fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',r(end,:,1)); %zlim
    fprintf(fid,'\n');
else
    dex = ves_data.coords(4,:)==1;
    r_v = ves_data.coords(1,dex);
    z_v = ves_data.coords(2,dex);
    fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',r_v); %rlim
    fprintf(fid,'\n');
    fprintf(fid,'% 15.9E% 15.9E% 15.9E% 15.9E% 15.9E\n',z_v); %zlim
    fprintf(fid,'\n');
end

fclose(fid);


end

