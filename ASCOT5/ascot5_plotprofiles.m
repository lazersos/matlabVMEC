function fig = ascot5_plotprofiles(a5file,bfieldid,plasmaid,varargin)
%ASCOT5_PLOTPROFILES Plots profile summary information
%   The ASCOT5_PLOTPROFILES function plots the profile data for a given
%   ASCOT5 file.  The code takes a filename, magnetic field id, and plasma
%   id as inputs.  If empty arrays are passed to either of the ids, the
%   active group id is used for each.  
%
%   Example:
%       a5file='ascot5_test.h5';
%       bid=0838288192;
%       pid=0396210459;
%       fig=ascot5_plotwall(a5file,bid,pid);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Defaults
fig=[];
nphi=1;

% Handle varargin
if nargin > 3
    i=1;
    while i <= nargin-3
        switch varargin{i}
            case{'nphi'}
                i=i+1;
                nphi=varargin{i};
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Use active
if isempty(bfieldid)
    bfieldid=h5readatt(a5file,'/bfield','active');
    disp(['  Using bfield: ' bfieldid]);
end
if isempty(plasmaid)
    plasmaid=h5readatt(a5file,'/plasma','active');
    disp(['  Using plasma: ' plasmaid]);
end

%Pull Grid information
path = ['/bfield/B_STS_' num2str(bfieldid,'%10.10i')];
try
    psi = h5read(a5file,[path '/psi']);
    psi_r0 = h5read(a5file,[path '/psi_rmin']);
    psi_r1 = h5read(a5file,[path '/psi_rmax']);
    psi_z0 = h5read(a5file,[path '/psi_zmin']);
    psi_z1 = h5read(a5file,[path '/psi_zmax']);
    psi_p0 = h5read(a5file,[path '/psi_phimin']);
    psi_p1 = h5read(a5file,[path '/psi_phimax']);
    psi0   = h5read(a5file,[path '/psi0']);
    psi1   = h5read(a5file,[path '/psi1']);
    rax   = h5read(a5file,[path '/axisr']);
    zax   = h5read(a5file,[path '/axisz']);
catch
    disp(['ERROR: Could not find bfield: ' num2str(bfieldid,'%10.10i')]);
    return;
end

% Check nphi consistency
if nphi > size(psi,2)
    fprintf('ERROR: nphi = %i > nphi_max = %i\n',nphi,size(psi,2));
    return;
end

%Pull Profile information (could be 1D or 1DS)
path = ['/plasma/plasma_1D_' num2str(plasmaid,'%10.10i')];
try
    rho = h5read(a5file,[path '/rho']);
catch
    path = ['/plasma/plasma_1DS_' num2str(plasmaid,'%10.10i')];
    try
        rho = h5read(a5file,[path '/rho']);
    catch
        disp(['ERROR: Could not find plasma: ' num2str(plasmaid,'%10.10i')]);
        return;
    end
end
ne = h5read(a5file,[path '/edensity']);
te = h5read(a5file,[path '/etemperature']);
ni = h5read(a5file,[path '/idensity']);
ti = h5read(a5file,[path '/itemperature']);

% Plot rho
fig = figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
subplot(2,4,[1 2 5 6]);
np = size(psi,2);
f = squeeze(psi(:,nphi,:));
f(f<psi0) = psi0;
f = sqrt(f./psi1);
nr = size(f,1);
nz = size(f,2);
xaxis = psi_r0:(psi_r1-psi_r0)./(nr-1):psi_r1;
xtick = psi_r0:(psi_r1-psi_r0)./(4):psi_r1;
yaxis = psi_z0:(psi_z1-psi_z0)./(nz-1):psi_z1;
ytick = psi_z0:(psi_z1-psi_z0)./(4):psi_z1;
paxis = psi_p0:(psi_p1-psi_p0)./(np-1):psi_p1;
pixplot(xaxis,yaxis,f);
hold on;
[~,hac]=contour(xaxis,yaxis,f',[1 1]);
set(hac,'LineWidth',4,'Color','white');
[~,hac]=contour(xaxis,yaxis,f',[0.5 0.5]);
set(hac,'LineWidth',2,'Color','white');
plot(rax(nphi),zax(nphi),'+w','MarkerSize',18);
set(gca,'FontSize',24,'XTick',xtick,'YTick',ytick,...
    'XTickLabelMode','auto','YTickLabelMode','auto');
xlabel('R [m]');
ylabel('Z [m]');
title('ASCOT5 Grid');
colormap parula;
ha = colorbar('Location','north');
set(ha,'Color','black');
ylabel(ha,'Rho');
text(min(xlim)+0.1.*range(xlim),min(ylim)+0.1.*range(ylim),...
    ['\phi= ' num2str(round(paxis(nphi)),'%i') '^o'],'FontSize',24);

%Plot ne
subplot(2,4,[3 4]);
plot(rho,ne./1E19,'k','LineWidth',2);
hold on;
plot(rho,ni./1E19,'--k','LineWidth',2);
set(gca,'FontSize',24);
xlabel('Rho');
ylabel('n x10^{19} [m^{-3}]');
xlim([0 max(rho)]);
legend('Elec.','Ion');

%Plot te
subplot(2,4,[7 8]);
plot(rho,te./1E3,'k','LineWidth',2);
hold on;
plot(rho,ti./1E3,'--k','LineWidth',2);
set(gca,'FontSize',24);
xlabel('Rho');
ylabel('T [keV]');
xlim([0 max(rho)]);

end

