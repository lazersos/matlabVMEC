function data=beams3d_bes(beam_data,geofile,varargin)
%BEAMS3D_BES(beam_data,geofile) Calculates simulated BES data from BEAMS3D
%   The BEAMS3D_BES routine calculates simulated BES data using data as
%   read by the READ_BEAMS3D routine, and a text geometry file.
% 
%   Optional Arguments
%       'plots'     : Create plots
%       'beams'     : Downselect beamlines considered
%                       data=beams3d_bes(beam_data,file,'beams',4:6);
%       'geo_dex'   : Select which viewing cords to use 6:2:38 AEA default
%                       data=beams3d_bes(beam_data,file,'geo_dex',6:2:38);
%
% Example usage
%      beam_data = read_beams3d('beams3d_test.h5');
%      data=beams3d_bes(beam_data,'QSK-aea21-losCoords-op12b-design-cad.txt');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% Helpers
ec = 1.60217662E-19;
lplot = 0;
geo_dex = 6:2:38; % Use for downselecting to actual channels of ILS Red AEA
beam_dex = 1:6; % Use to downselect beams

% Handle varargin
if nargin > 2
    i=1;
    while i < nargin-2
        switch varargin{i}
            case 'plots'
                lplot=1;
            case 'beams'
                i=i+1;
                beam_dex=varargin{i};
            case 'geo_dex'
                i=i+1;
                geo_dex=varargin{i};
        end
        i=i+1;
    end
end

% Read geo file
fid = fopen(geofile);
if fid < 0
    disp(['  Error opening: ' geofile]);
end
geo =[];
header = fgetl(fid);
while ~feof(fid)
    geo = [geo; fscanf(fid,'%*s %f, %f, %f, %f, %f',[1 6])]; 
end
fclose(fid);
if ~isempty(geo_dex)
    geo = geo(geo_dex,:);
end

% Plot geo
if lplot
    figure('Position',[1 1 1024 768],'Color','white');
    plot3(geo(:,1),geo(:,2),geo(:,3),'ok');
    hold on;
    quiver3(geo(:,1),geo(:,2),geo(:,3),geo(:,4),geo(:,5),geo(:,6),500,'k');
    axis equal;
end

% Extract information from BEAMS3D
NEUT   = beam_data.neut_lines(2,:);
dex    = NEUT == 0;
X_BEAM = beam_data.X_lines(2,dex);
Y_BEAM = beam_data.Y_lines(2,dex);
Z_BEAM = beam_data.Z_lines(2,dex);
E_BEAM = (0.5/ec).*beam_data.mass'.*beam_data.vll_lines(1,:).^2;
E_BEAM = E_BEAM(dex);
B_BEAM = beam_data.Beam(dex)';

% Downselect beams
ntotal = [];
dex = zeros(1,length(B_BEAM));
for i = beam_dex(1:end)
    dex = dex + (B_BEAM == i);
end
dex = dex > 0;
X_BEAM = X_BEAM(dex);
Y_BEAM = Y_BEAM(dex);
Z_BEAM = Z_BEAM(dex);
E_BEAM = E_BEAM(dex);
B_BEAM = B_BEAM(dex);

if lplot
    plot3(X_BEAM,Y_BEAM,Z_BEAM,'.r','MarkerSize',0.1);
end

% Construct path and count particles
nsteps = 64;
d_lim = 0.01;
bes = [];
bes_sig=[];
rplot=[];
vol=[];
eplot=zeros(size(geo,1),max(beam_dex));
rplot=zeros(size(geo,1),max(beam_dex));
bes=zeros(size(geo,1),max(beam_dex));
for i = 1:size(geo,1)
    x = geo(i,1)-X_BEAM;
    y = geo(i,2)-Y_BEAM;
    z = geo(i,3)-Z_BEAM;
    d1 = min(min(min(sqrt(x.^2+y.^2+z.^2))))*0.95;
    d2 = max(max(max(sqrt(x.^2+y.^2+z.^2))))*1.05;
    dl = (d2-d1)./nsteps;
    dex = zeros(1,length(X_BEAM));
    for l = d1:dl:d2
        x = geo(i,1)+l.*geo(i,4);
        y = geo(i,2)+l.*geo(i,5);
        z = geo(i,3)+l.*geo(i,6);
        dx = (X_BEAM-x);
        dy = (Y_BEAM-y);
        dz = (Z_BEAM-z);
        d  = sqrt(dx.^2+dy.^2+dz.^2);
        dex = or(dex,d < d_lim);
    end
        
    dx = []; dy = []; dz = []; d = [];
    % Now calc BES data
    for j = beam_dex
        dex2 = and(dex,B_BEAM==j);
        x = X_BEAM(dex2);
        y = Y_BEAM(dex2);
        z = Z_BEAM(dex2);
        if lplot
            plot3(x,y,z,'.b','MarkerSize',1.0);
        end
        if all(~dex2)
            continue
        end 
        energy = E_BEAM(dex2);
        eplot(i,j) = mean(energy);
        r = sqrt(x.^2+y.^2);
        rplot(i,j) = mean(r);
        bes(i,j) = sum(dex2);
    end
end
% Now switch to beam density
bes = cumsum(bes,1);
ntotal = 100.*max(bes,[],1)./(100-beam_data.Shinethrough');
bes = ntotal-bes;

if lplot
    for i = beam_dex
        figure;
        shine = beam_data.Shinethrough(i).*1E-2;
        r = rplot(:,i);
        plot(rplot(:,i),sum(bes(:,i))/(1-shine)-cumsum(bes(:,i)),'ok');
        set(gcf,'Position',[1 1 1024 768],'Color','white');
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('BES Simulated');
        title(['BEAMS3D Beam Emission (E = ' num2str(mean(eplot(:,i)).*1E-3,'%5.2f') ' keV)']);
        xlim([min(r(r>0)) max(r(r>0))]);
    end
end

data.bes = bes;
data.r   = rplot;
data.Energy = eplot;
data.bes_sig = bes_sig;
data.bes_total = ntotal;

return;

end

