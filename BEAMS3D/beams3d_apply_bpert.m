function [ br,bphi,bz ] = beams3d_apply_bpert(filename_in,filename_out,fluxi0, ni,mi,varargin)
%beams3d_apply_bpert(filename_in,filename_out,fluxi0, ni,mi,varargin)
%applies a helical perturbation to the magnetic field in filename_in and
%writes the result to filename_out. Since only the magnetic field
%components are written in filename_out, this can also be a FIELDLINES
%restart file. The relative amplitude (to the maximum field found in the
%arrays) fluxi0, as well as the toroidal and poloidal mode numbers ni and
%mi can be specified as scalars or arrays, in which case all components are
%added up. Plotting of the perturbation amplitudes is available.
%
%The resulting magnetic field is calculated in the following way:
%The magnetic field in filename_in is scaled helically within the plasma
%according to fluxi0, ni and mi. The curl is calculated from the result,
%and this is taken as the perturbation field for B = B0 + B1.
%   Options:
%       'plot':        Plots the helical perturbation(s)
%
%   Usage:
%   beams3d_apply_bpert('beams3d_38581_3350_b3d.h5','fieldlines_38581_3350_b3d.h5',0.3,1,2);
%
%   Created by: D. Kulla (david.kulla@ipp.mpg.de)
%   Version:    1.0
%   Date:       22/11/30

lplot = 0;
lsave = 0;

% Handle varargin
if nargin > 5
    i = 1;
    while i < nargin -4
        switch varargin{i}
            case 'plot'
                lplot=1;
            case 'save'
                lsave = 1;
        end
        i=i+1;
    end
end

br = h5read(filename_in,'/B_R');
bphi = h5read(filename_in,'/B_PHI');
bz = h5read(filename_in,'/B_Z');
r = h5read(filename_in,'/raxis');
phi = h5read(filename_in,'/phiaxis');
z = h5read(filename_in,'/zaxis');
sarr = h5read(filename_in,'/S_ARR');
uarr = h5read(filename_in,'/U_ARR');
rhoarr = sqrt(sarr);
%[r0, phi0, z0] = meshgrid(r,phi,z);
[rg, phig, zg] = ndgrid(r,phi,z);
xg = rg .* cos(phig);
yg = rg .* sin(phig);
%bx = br .* cos(bphi);
%by = br.* sin(bphi);

fluxphi = zeros(size(rhoarr));
for i = numel(fluxi0)
fluxir = fluxi0(i) * (rhoarr.^2) .*  (1- rhoarr).^2;
fluxphi = fluxphi + fluxir .* cos(mi.*uarr + ni .* phig);
end
fluxphi(sarr>1)=0;

brc = br.* fluxphi;
bphic = bphi.* fluxphi;
bzc = bz.* fluxphi;

% %%
% figure
% colors = parula(20);
% for i = 1:20
% p = patch(isosurface(r0,phi0,z0,permute(fluxphi,[2,1,3]),i/10));
% isonormals(r0,phi0,z0,permute(fluxphi,[2,1,3]),p)
% p.FaceColor = colors(i,:);
% p.EdgeColor = 'none';
% hold on
% end
% view(3)
% axis equal
% camlight
% lighting phong
% xlabel('R')
% ylabel('PHI')
% zlabel('Z')
%
% %%
if lplot
    for j = numel(fluxi0)
        figure
        colors = parula(20);
        for i = 1:20
            p = patch(isosurface(xg,yg,zg,fluxphi,i/20*max(fluxphi,[],'all')));
            %isonormals(xg,yg,zg,fluxphi,p)
            p.FaceColor = colors(i,:);
            p.EdgeColor = 'none';
            hold on
        end
        view(3)
        axis equal
        camlight
        lighting phong
        xlabel('R')
        ylabel('PHI')
        zlabel('Z')
    end
end
%%

%[curlr,curlphi,curlz,cav] = curl(xg,yg,zg,br,bphi,bz);
%[curlr,curlphi,curlz,cav] = curl(rg,phig,zg,br,bphi,bz);
[curlr,curlphi,curlz,~] = curl(brc,bphic,bzc);
curlr(isnan(curlr)|isinf(curlr)|sarr>1) = 0;
curlphi(isnan(curlphi)|isinf(curlphi)|sarr>1) = 0;
curlz(isnan(curlz)|isinf(curlz)|sarr>1) = 0;


modb=sqrt(br.^2+bphi.^2+bz.^2);

% modcurl=sqrt(curlr.^2+curlphi.^2+curlz.^2);
% [~,I] = max(modcurl,[],'all');
% 
% fac = modb(I) ./ modcurl(I);


% figure
% %quiver3(xg,yg,zg,curlr,curlphi,curlz);
% %quiver3(rg,phig,zg,curlr,curlphi,curlz);
% quiver(squeeze(rg(:,1,:)),squeeze(zg(:,1,:)),squeeze(curlr(:,9,:)),squeeze(curlz(:,9,:)));
% %pixplot(squeeze(curlphi(:,:,50)))
%
% pixplot(squeeze(curlr(:,:,50))./squeeze(br(:,:,50)))
% caxis([-1 1])

if lplot
figure
subplot(1,2,1)
pixplot(squeeze(modb(:,1,:)))
title('Before Pert')
axis equal
end

br = br + curlr;
bphi = bphi + curlphi;
bz = bz + curlz;

if lplot
subplot(1,2,2)
modb=sqrt(br.^2+bphi.^2+bz.^2);
pixplot(squeeze(modb(:,1,:)))
title('After Pert')
axis equal
end

if lsave
rbphi = h5read(filename_out,'/B_PHI');
if sum(size(rbphi)-size(bphi))~=0
delete_hdf5_group(filename_out,'/B_R');
h5create(filename_out,'/B_R',size(br));
delete_hdf5_group(filename_out,'/B_PHI');
h5create(filename_out,'/B_PHI',size(br));
delete_hdf5_group(filename_out,'/B_Z');
h5create(filename_out,'/B_Z',size(br));
delete_hdf5_group(filename_out,'/raxis');
h5create(filename_out,'/raxis',size(r));
delete_hdf5_group(filename_out,'/phiaxis');
h5create(filename_out,'/phiaxis',size(phi));
delete_hdf5_group(filename_out,'/zaxis');
h5create(filename_out,'/zaxis',size(z));
delete_hdf5_group(filename_out,'/nr');
h5create(filename_out,'/nr',1);
delete_hdf5_group(filename_out,'/nphi');
h5create(filename_out,'/nphi',1);
delete_hdf5_group(filename_out,'/nz');
h5create(filename_out,'/nz',1);
end
%h5create(filename_out,'/B_R',size(br));
h5write(filename_out,'/B_R',br)
h5write(filename_out,'/B_PHI',bphi)
h5write(filename_out,'/B_Z',bz)
h5write(filename_out,'/raxis',r)
h5write(filename_out,'/phiaxis',phi)
h5write(filename_out,'/zaxis',z)
h5write(filename_out,'/nr',numel(r))
h5write(filename_out,'/nphi',numel(phi))
h5write(filename_out,'/nz',numel(z))


end
end

function delete_hdf5_group(filename, group_name)
fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,group_name,'H5P_DEFAULT');
H5F.close(fid);
end
