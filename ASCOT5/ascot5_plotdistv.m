function fig = ascot5_plotdistv(a5file,runid)
%ASCOT5_POLOTDISTV plots the Vpara/Vperp distribution function in 2D
%   The ASCOT5_PLOTDISTV function plots the parallel and perpedicular
%   velocity distribuiton function from as ASCOT5 run.  It utilizes the
%   DIST5D distribution to make the plot. A handle to the figure graphic is
%   returned.
%
%   Example:
%       a5file='ascot5_W7X_20180821_012_5100_h8.h5';
%       id=0614695199;
%       fig=ascot5_plotdistv(a5file,id);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

fig=[];
amu=1.66054E-27;
% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

if isempty(runid)
    runid=h5readatt(a5file,'/results','active');
    disp(['  Using runid: ' runid]);
end


path = ['/results/run_' num2str(runid,'%10.10i') '/distrho5d/'];
try
    dist = h5read(a5file,[path '/ordinate']);
catch
    disp(['ERROR: Could not find run number, dist5d, or ordinate: ' num2str(runid,'%10.10i')]);
    return;
end
x1 = h5read(a5file,[path '/abscissa_vec_07']); %charge
x2 = h5read(a5file,[path '/abscissa_vec_06']); %time
x3 = h5read(a5file,[path '/abscissa_vec_05']); %vperp
x4 = h5read(a5file,[path '/abscissa_vec_04']); %vpara
x5 = h5read(a5file,[path '/abscissa_vec_03']); %Z
x6 = h5read(a5file,[path '/abscissa_vec_02']); %phi [deg]
x7 = h5read(a5file,[path '/abscissa_vec_01']); %R


path = ['/results/run_' num2str(runid,'%10.10i') '/inistate/'];
mass = h5read(a5file,[path '/mass']); %mass
mass = mean(mass).*amu;


% Make factors
dx1 = (x1(end)-x1(1))./length(x1);
dx2 = (x2(end)-x2(1))./length(x2);
dx3 = (x3(end)-x3(1))./length(x3);
dx4 = (x4(end)-x4(1))./length(x4);
dx5 = (x5(end)-x5(1))./length(x5);
dx6 = (x6(end)-x6(1))./length(x6);
dx7 = (x7(end)-x7(1))./length(x7);


% Make dist
dist=squeeze(sum(dist,[5 6 7]))';
%factor = dx1.*dx2.*dx5.*dx6.*dx7;
factor = dx1.*dx5.*dx6.*dx7;
x = x4./(1E6.*mass);
y = x3./(1E6.*mass);
xtick = unique(round(x));
ytick = unique(round(y));
xticklabel={};
yticklabel={};
for i=1:length(xtick)
    xticklabel = [xticklabel num2str(xtick(i),'%i')];
end
for i=1:length(ytick)
    yticklabel = [yticklabel num2str(ytick(i),'%i')];
end

% Make plot
fig = figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
pixplot(x,y,dist./factor);
colormap jet;
set(gca,'FontSize',24);
xlabel('Parallel Velocity (V_{||}) x1000 [km/s]');
ylabel('Perpendicular Velocity (V_{\perp}) x1000 [km/s]');
set(gca,'XTick',xtick,'YTick',ytick,'XTickLabel',xticklabel,'YTickLabel',yticklabel);
title('ASCOT5 Distribution Function');
colorbar;


end

