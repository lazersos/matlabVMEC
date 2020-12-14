function fig = ascot5_walloverview(a5file,wallid,runid,bid,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Helpers and Defaults
fig=[];
amu = 1.66053906660E-27;
l3d = 0;

% Handle varargin
if nargin > 4
    i=1;
    while i <= nargin-4
        switch varargin{i}
            case{'3d'}
                l3d=1;
        end
        i = i + 1;
    end
end


% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Pull Wall
path = ['/wall/wall_3D_' num2str(wallid,'%10.10i')];
try
    x1x2x3 = h5read(a5file,[path '/x1x2x3']);
catch
    disp(['ERROR: Could not find wall: ' num2str(wallid,'%10.10i')]);
    return;
end
y1y2y3 = h5read(a5file,[path '/y1y2y3']);
z1z2z3 = h5read(a5file,[path '/z1z2z3']);
p=patch(x1x2x3,y1y2y3,z1z2z3);
[faces, vertex]=reducepatch(p,1.0);


% Pull magnetic axis
path = ['/bfield/B_STS_' num2str(bid,'%10.10i')];
try
    axisr = h5read(a5file,[path '/axisr']);
catch
    disp(['ERROR: Could not find bfield: ' num2str(bid,'%10.10i')]);
    return;
end
axisz = h5read(a5file,[path '/axisz']);
axisp0 = h5read(a5file,[path '/axis_phimin']);
axisp1 = h5read(a5file,[path '/axis_phimax']);
axisp = axisp0:(axisp1-axisp0)/(length(axisr)-1):axisp1;
axisp = deg2rad(axisp);
pmax  = deg2rad(axisp1);

% Pull run information
path = ['/results/run_' num2str(runid,'%10.10i') '/endstate'];
try
    walltile = h5read(a5file,[path '/walltile']);
catch
    disp(['ERROR: Could not result: ' num2str(runid,'%10.10i')]);
    return;
end
weight = h5read(a5file,[path '/weight']);
vr = h5read(a5file,[path '/vr']);
vphi = h5read(a5file,[path '/vphi']);
vz = h5read(a5file,[path '/vz']);
mass = h5read(a5file,[path '/mass']).*amu; %in amu
v2 = vr.*vr+vphi.*vphi+vz.*vz;
q  = 0.5.*mass.*v2.*weight;

% Calc Area
V0=[x1x2x3(2,:)-x1x2x3(1,:);y1y2y3(2,:)-y1y2y3(1,:);z1z2z3(2,:)-z1z2z3(1,:)];
V1=[x1x2x3(3,:)-x1x2x3(1,:);y1y2y3(3,:)-y1y2y3(1,:);z1z2z3(3,:)-z1z2z3(1,:)];
F = cross(V0,V1);
A = 0.5.*sqrt(sum(F.*F));

% Calculate Mean position of each triangle
xm = mean(x1x2x3);
ym = mean(y1y2y3);
zm = mean(z1z2z3);
pm = atan2(ym,xm);
%pm(pm<0) = pm(pm<0) + 2*pi;
rm = sqrt(xm.*xm+ym.*ym);
r0m   = pchip(axisp,(axisr),mod(pm,pmax));
z0m   = pchip(axisp,(axisz),mod(pm,pmax));
tm    = atan2(zm-z0m,rm-r0m);


% Count hits
wall_strikes = zeros(1,length(xm));
%nhits=[];
mask = unique(walltile)+1;
%for i = mask'
%    if (i==0), continue; end
%    dex = walltile==i;
%    nhits = [nhits; sum(dex)];
%end
%wall_strikes(mask) = nhits;

%Construct wall_load using accumarray
walltile(walltile==0) = size(x1x2x3,2) + 1; %set 0's from before to extra value to be able to truncate (accumarray needs positive integers)
% Convert to heatflux
wall_load = accumarray(walltile,q); %Sum all entries in q that have the same value in walltile, put result at the position given by walltile
wall_load = wall_load(1:end-1)'; %truncate "0's"
wall_load = wall_load./A;

%temp = ones(1,length(pm))+10;
%scatter(pm,tm,temp,wall_load);

% Adjsut theta to peel from outboard
tm(tm>=0) = abs(tm(tm>=0)-pi);
tm(tm<0) = -(tm(tm<0)+pi);

% Now divide up into modules
nfp = round(2*pi/pmax);
factor = pi/nfp; mdex={};
p1 = -factor;
p2 = factor;
for i = 1:nfp
    mdex{i} = and(pm>=p1,pm<p2);
    p1 = p2;
    p2 = p2 + 2*factor;
    if p2>pi
        p2 = -pi + (p2 - pi);
    end
    if p1 >= pi
        p1 = -pi + p2 - pi;
    end
end
%mdex{1} = and(pm>=-factor,pm<factor);
%mdex{2} = and(pm>=factor,pm<3*factor);
%mdex{3} = and(pm>=3*factor,pm<pi);
%mdex{4} = and(pm>=-pi,pm<-3*factor);
%mdex{5} = and(pm>=-3*factor,pm<-factor);
fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
for i=1:nfp
    subplot(1,5,i);
    x_temp = pm(mdex{i});
    y_temp = tm(mdex{i});
    c_temp = wall_load(mdex{i})./1E6;
    dex_save = 1:length(c_temp);
    dex = c_temp > 1E-3;
    dex_save = dex_save(dex);
    temp = ones(1,length(x_temp))+10;
    scatter(x_temp(dex),y_temp(dex),temp(dex),c_temp(dex),'fill');
    colormap jet;
    caxis([0 1]);
    ylim([-1 1].*pi);
    title(['Module ' num2str(i,'%i')]);
    xlabel('\phi');
    ylabel('\theta');
    outdex=overplot_largest(x_temp(dex),y_temp(dex),c_temp(dex));
    if (l3d)
        x_3d = x1x2x3(:,mdex{i});
        y_3d = y1y2y3(:,mdex{i});
        z_3d = z1z2z3(:,mdex{i});
        x0_3d = r0m(mdex{i}).*cos(pm(mdex{i}));
        y0_3d = r0m(mdex{i}).*sin(pm(mdex{i}));
        z0_3d = z0m(mdex{i});
        for j = 1: length(outdex)
            tdex = dex_save(outdex(j));
            fig3d=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            patch(x_3d,y_3d,z_3d,c_temp,'EdgeColor','none');
            campos([x0_3d(tdex),y0_3d(tdex),z0_3d(tdex)]);
            camtarget([mean(x_3d(:,tdex)),mean(y_3d(:,tdex)),mean(z_3d(:,tdex))]);
            patch(x_3d(:,tdex),y_3d(:,tdex),z_3d(:,tdex),c_temp(:,tdex),'EdgeColor','white');
            camva(60);
            camproj perspective;
            title(['Module ' num2str(i,'%i') ' Spot #' num2str(j,'%i')]);
            colormap jet;
            caxis([0 c_temp(tdex)]);
            drawnow;
        end
    end
end
ha = colorbar('Location','EastOutside');
set(ha,'FontSize',24);
ylabel(ha,'Heat Flux [MW/m^2]');


end

