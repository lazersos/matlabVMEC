function plot_stellopt(data,varargin)
%PLOT_STELLOPT(data) plots stellopt data
%   This function plots the data read from a stellopt output file.  It
%   plots the data according to the datatype field set by read_stellopt.
%
%   Example:
%       out_data=read_stellopt('output.test');
%       press_data=read_stellopt('p_prof.test.min');
%
%   See also read_stellopt and make_stellopt_movie.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/18/11


vmec_data=[];
if ~isfield(data,'datatype')
    disp('ERROR: data not a known datatype!');
    return
end

plottype = 'stellopt_output';
floop_data=[];
segrog_data=[];
ves_data=[];

for i=1:nargin-1
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch upper(varargin{i}.datatype)
                case 'WOUT'
                    vmec_data=varargin{i};
                case 'FLUXLOOP'
                    floop_data=varargin{i};
                case 'VESSEL'
                    ves_data=varargin{i};
            end
        end
    end
    if isstr(varargin{i})
        plottype = varargin{i};
    end
end
        
switch data.datatype
    case {'stellopt_pressure','stellopt_neteti','stellopt_magdiag',...
            'stellopt_mse','stellopt_linne','stellopt_jacobian'}
        plottype=data.datatype;
    case {'stellopt_new','stellopt_output'}
    otherwise
        disp(['ERROR: Unsuported datatype ' data.datatype '!']);
end

switch plottype
    case 'stellopt_output'
        plot_stellopt_output(data);
    case 'stellopt_pressure'
        plot_stellopt_pressure(data);
    case 'stellopt_neteti'
        plot_stellopt_neteti(data);
    case 'stellopt_magdiag'
        plot_stellopt_magdiag(data);
    case 'stellopt_magnetics'
        plot_stellopt_magnetics(data,vmec_data,floop_data,segrog_data,ves_data);
    case 'stellopt_mse'
        plot_stellopt_mse(data,vmec_data);
    case 'stellopt_linne'
        plot_stellopt_linene(data,vmec_data);
    case 'stellopt_te'
        plot_stellopt_te(data,vmec_data);
    case 'stellopt_ne'
        plot_stellopt_ne(data,vmec_data);
    case 'stellopt_kinetic'
        plot_stellopt_kinetic(data,vmec_data);
    case 'stellopt_mse_new'
        plot_stellopt_mse_new(data,vmec_data);
    case 'stellopt_kinetic_grad'
        plot_stellopt_kinetic_grad(data,vmec_data);
    case 'stellopt_jacobian'
        plot_stellopt_jacobian(data);
    otherwise
        disp(['ERROR: Unsuported plottype ' plottype '!']);
end
return
end

function plot_stellopt_linene(data,vmec_data)
ne_norm = 1E19;
theta = 0:2*pi./359:2*pi;
zeta  = data.NELINE_PHI0(1,1);
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
subplot(2,2,[1 3]);
ns = vmec_data.ns;
plot(r(ns,:,1),z(ns,:,1),'r','LineWidth',2);
hold on;
dex = round(2:ns/10:ns);
plot(r(1,1,1),z(1,1,1),'+k','LineWidth',2);
plot(r(dex,:,1)',z(dex,:,1)','k','LineWidth',2);
j = size(data.NELINE_R0,1);
for i = 1: size(data.NELINE_R0,2)
    r2 = [data.NELINE_R0(j,i) data.NELINE_R1(j,i)];
    z2 = [data.NELINE_Z0(j,i) data.NELINE_Z1(j,i)];
    plot(r2,z2,'k','LineWidth',2);
end
axis equal;
xlim([3.2 4.2]);
ylim([-1.5 1.5]);
xlabel('R [m]');
ylabel('Z [m]');
hold off;
subplot(2,2,2);
errorbar(data.NELINE_R0(j,:),data.NELINE_target(j,:)./ne_norm,data.NELINE_sigma(j,:)./ne_norm,'ok','LineWidth',2,'MarkerSize',18);
hold on;
plot(data.NELINE_R0(j,:),data.NELINE_equil(j,:)./ne_norm,'xr','LineWidth',2,'MarkerSize',18);
plot([1 1].*r(1,1,1),ylim,'k','LineWidth',2);
r2 = [1 1].*min(r(ns,:,1));
r2 = [r2; [1 1].*max(r(ns,:,1))];
z2 = [ylim; ylim];
plot(r2',z2','r','LineWidth',2);
xlim([3.2 4.2]);
ylabel('Line Ne [10^{19} m^{-2}]');
legend('Data','Recon','Axis','Edge');
hold off;
subplot(2,2,4);
bar(data.NELINE_R0(j,:),sqrt(data.NELINE_chisq(j,:)),'r');
xlim([3.2 4.2]);
xlabel('R [m]');
ylabel('Error Norm');


end


function plot_stellopt_magnetics(data,vmec_data,floop_data,segrog_data,ves_data)
num_pannel=0;
if isfield(data,'B_PROBES_X')
    num_pannel = num_pannel+1;
end
if isfield(data,'FLUXLOOPS_target')
    num_pannel = num_pannel+1;
end
if isfield(data,'SEGROG_target')
    num_pannel = num_pannel+1;
end
num_vert = 2;
if isempty(vmec_data)
    num_vert=1;
end
if num_vert ==2
    theta = 0:2*pi./179:2*pi;
    zeta  = 0:2*pi./179:2*pi;
    rv     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    bv     = cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm,vmec_data.xn);
    if vmec_data.iasym
        rv     = rv + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv     = zv + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
        bv     = bv + sfunct(theta,zeta,vmec_data.bmns,vmec_data.xm,vmec_data.xn);
    end
end
i=1;
iter = length(data.iter);
if isfield(data,'B_PROBES_X')
    subplot(num_vert,num_pannel,i);
    ha=isotoro(rv,zv,zeta,vmec_data.ns,bv);
    set(ha,'FaceAlpha',1./3.);
    hold on;
    if ~isempty(ves_data)
        hv = plot_vessel(ves_data,'phi3d',0.0);
        set(hv,'LineWidth',3.0,'Color','black');
    end
    title('');
    axis off;
    view([25 37]);
    zoom(1.8);
    colorbar off;
    vals = data.B_PROBES_equil(iter,:);
    target = data.B_PROBES_target(iter,:);
    sigma = data.B_PROBES_sigma(iter,:);
    target(sigma >= 1.0E10) = vals(sigma >= 1.0E10);
    %valg = 0.025;
    %valy = 0.050;
    %delta = abs(vals-target)./abs(target);
    valg = 0.5;
    valy = 1.0;
    delta = abs(vals-target)./abs(sigma);
    x = data.B_PROBES_X(iter,delta <= valg);
    y = data.B_PROBES_Y(iter,delta <= valg);
    z = data.B_PROBES_Z(iter,delta <= valg);
    plot3(x,y,z,'+g','LineWidth',4.0,'MarkerSize',18);
    x = data.B_PROBES_X(iter,delta > valg);
    y = data.B_PROBES_Y(iter,delta > valg);
    z = data.B_PROBES_Z(iter,delta > valg);
    plot3(x,y,z,'+m','LineWidth',4.0,'MarkerSize',18);
    x = data.B_PROBES_X(iter,delta > valy);
    y = data.B_PROBES_Y(iter,delta > valy);
    z = data.B_PROBES_Z(iter,delta > valy);
    plot3(x,y,z,'+r','LineWidth',4.0,'MarkerSize',18);
    x = data.B_PROBES_X(iter,sigma >= 1.0E10);
    y = data.B_PROBES_Y(iter,sigma >= 1.0E10);
    z = data.B_PROBES_Z(iter,sigma >= 1.0E10);
    plot3(x,y,z,'+k','LineWidth',4.0,'MarkerSize',18);
    subplot(num_vert,num_pannel,i+num_pannel);
    sigma2=sigma;
    sigma(sigma>1.0E10) = max(vals);
    errorbar(data.B_PROBES_target(iter,:),sigma,'ok');
    ylim(1.1.*[min(vals) max(vals)]);
    hold on;
    dex = find(delta <= valg);
    plot(dex,vals(delta <= valg),'xg','LineWidth',4.0,'MarkerSize',18);
    dex = find(delta > valg);
    plot(dex,vals(delta > valg),'xm','LineWidth',4.0,'MarkerSize',18);
    dex = find(delta > valy);
    plot(dex,vals(delta > valy),'xr','LineWidth',4.0,'MarkerSize',18);
    dex = find(sigma2 >= 1.0E10);
    plot(dex,vals(sigma2 >= 1.0E10),'xk','LineWidth',4.0,'MarkerSize',18);
    i=i+1;
    xlabel('Probe Index');
    ylabel('B [T]')
    xlim([0 size(data.B_PROBES_target,2)+1]);
end
if isfield(data,'FLUXLOOPS_target')
    subplot(num_vert,num_pannel,i);
    ha = isotoro(rv,zv,zeta,vmec_data.ns,bv);
    set(ha,'FaceAlpha',1./3.);
    hold on;
    if ~isempty(ves_data)
        hv = plot_vessel(ves_data,'phi3d',0.0);
        set(hv,'LineWidth',3.0,'Color','black');
    end
    title('');
    axis off;
    view([25 37]);
    zoom(1.8);
    colorbar off;
    vals = data.FLUXLOOPS_equil(iter,:);
    target = data.FLUXLOOPS_target(iter,:);
    sigma = data.FLUXLOOPS_sigma(iter,:);
    target(sigma >= 1.0E10) = vals(sigma >= 1.0E10);
    %valg = 0.025;
    %valy = 0.050;
    %delta = abs(vals-target)./abs(target);
    valg = 0.5;
    valy = 1.0;
    delta = abs(vals-target)./abs(sigma);
    if ~isempty(floop_data)
        %delta = abs(vals-target)./abs(target);
        loops = floop_data.loops(delta <= valg);
        for k=1:length(loops)
            plot3(loops{k}(1,:),loops{k}(2,:),loops{k}(3,:),'g','LineWidth',4.0,'MarkerSize',18);
        end
        loops = floop_data.loops(delta > valg);
        for k=1:length(loops)
            plot3(loops{k}(1,:),loops{k}(2,:),loops{k}(3,:),'m','LineWidth',4.0,'MarkerSize',18);
        end
        loops = floop_data.loops(delta > valy);
        for k=1:length(loops)
            plot3(loops{k}(1,:),loops{k}(2,:),loops{k}(3,:),'r','LineWidth',4.0,'MarkerSize',18);
        end
        loops = floop_data.loops(sigma >= 1.0E10);
        for k=1:length(loops)
            plot3(loops{k}(1,:),loops{k}(2,:),loops{k}(3,:),'k','LineWidth',4.0,'MarkerSize',18);
        end
    end
    subplot(num_vert,num_pannel,i+num_pannel);
    sigma2=sigma;
    sigma(sigma>1.0E10) = max(vals);
    errorbar(data.FLUXLOOPS_target(iter,:),sigma,'ok');
    ylim(1.1.*[min(vals) max(vals)]);
    hold on;
    dex = find(delta <= valg);
    plot(dex,vals(delta <= valg),'xg','LineWidth',4.0,'MarkerSize',18);
    dex = find(delta > valg);
    plot(dex,vals(delta > valg),'xm','LineWidth',4.0,'MarkerSize',18);
    dex = find(delta > valy);
    plot(dex,vals(delta > valy),'xr','LineWidth',4.0,'MarkerSize',18);
    dex = find(sigma2 >= 1.0E10);
    plot(dex,vals(sigma2 >= 1.0E10),'xk','LineWidth',4.0,'MarkerSize',18);
    i=i+1;ylabel('Flux [Wb]');xlabel('Loop Index');
    xlim([0 size(data.FLUXLOOPS_target,2)+1]);
end

end

function plot_stellopt_mse_new(data,vmec_data)
theta = 0:2*pi./179:2*pi;
zeta  = 0:2*pi./179:2*pi;
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
j = size(data.MSE_R,1);
subplot(3,2,[1 3 5]);
ns = vmec_data.ns;
h=isotoro(r,z,zeta,ns);
set(h,'FaceAlpha',0.25);
camlight left;
hold on;
plot3(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),squeeze(z(1,1,:)),'k','LineWidth',2);
x = data.MSE_R(j,:).*cos(data.MSE_PHI(j,:));
y = data.MSE_R(j,:).*sin(data.MSE_PHI(j,:));
z = data.MSE_Z(j,:);
plot3(x,y,z,'+k','LineWidth',2);
hold on;
subplot(3,2,2);
errorbar(data.MSE_R(j,:),data.MSE_target(j,:),data.MSE_sigma(j,:),'xk','LineWidth',2,'MarkerSize',18);
hold on;
plot(data.MSE_R(j,:),data.MSE_equil(j,:),'or','LineWidth',2,'MarkerSize',18);
ylabel('Polarization');
xlabel('R [m]');
hold off;
subplot(3,2,4);
plot(vmec_data.phi./vmec_data.phi(ns),vmec_data.jcurv./1000,'k','LineWidth',2,'MarkerSize',18);
ylabel('J_TOR [kA]');
xlabel('Normalized Flux');
xlim([0 1.1]);
subplot(3,2,6);
plot(vmec_data.phi./vmec_data.phi(ns),1./vmec_data.iotaf,'k','LineWidth',2,'MarkerSize',18);
ylabel('q');
xlabel('Normalized Flux');
xlim([0 1.1]);
end

function plot_stellopt_kinetic(data,vmec_data)
theta = 0:2*pi./359:2*pi;
zeta  = mean(data.TE_PHI(1,:));
%zeta  = 0.0;
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
subplot(3,2,[1 3 5]);
ns = vmec_data.ns;
plot(r(ns,:,1),z(ns,:,1),'r','LineWidth',2);
hold on;
dex = round(2:ns/10:ns);
plot(r(1,1,1),z(1,1,1),'+k','LineWidth',2);
plot(r(dex,:,1)',z(dex,:,1)','k','LineWidth',2);
j = size(data.TE_R,1);
plot(data.TE_R(j,:),data.TE_Z(j,:),'xk','LineWidth',2,'MarkerSize',18);
plot(data.NE_R(j,:),data.NE_Z(j,:),'ok','LineWidth',2,'MarkerSize',18);
if (isfield(data,'TI_R'))
    j = size(data.TI_R,1);
    plot(data.TI_R(j,:),data.TI_Z(j,:),'sk','LineWidth',2,'MarkerSize',18);
end
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
hold off;
subplot(3,2,2);
errorbar(data.TE_S(j,:),data.TE_target(j,:)./1000,data.TE_sigma(j,:)./1000,'xk','LineWidth',2,'MarkerSize',18);
hold on;
plot(data.TE_S(j,:),data.TE_equil(j,:)./1000,'r','LineWidth',2,'MarkerSize',18);
if (isfield(data,'TI_R'))
    hold on;
    errorbar(data.TI_S(j,:),data.TI_target(j,:)./1000,data.TI_sigma(j,:)./1000,'sk','LineWidth',2,'MarkerSize',18);
end
xlim([0 1.1]);
ylim([0 1.2*max(data.TE_target(j,:)./1000)]);
ylabel('T [keV]');
hold off;
subplot(3,2,4);
errorbar(data.NE_S(j,:),data.NE_target(j,:),data.NE_sigma(j,:),'ok','LineWidth',2,'MarkerSize',18);
hold on;
plot(data.NE_S(j,:),data.NE_equil(j,:),'r','LineWidth',2,'MarkerSize',18);
xlim([0 1.1]);
xlabel('');
ylabel('Ne (norm)');
subplot(3,2,6);
plot(vmec_data.phi./vmec_data.phi(vmec_data.ns),vmec_data.presf./1000,'k','LineWidth',2,'MarkerSize',18);
hold on;
pgrad = [0.0 diff(vmec_data.presf)]./(1./vmec_data.ns)./1000;
pgrad(pgrad > 0.0) = 0.0;
pgrad = max(vmec_data.presf/1000)*pgrad / max(abs(pgrad));
plot(vmec_data.phi./vmec_data.phi(vmec_data.ns),abs(pgrad),'g','LineWidth',2,'MarkerSize',18);
ylabel('P [kPa]');
xlabel('Normalized Flux');
xlim([0 1.1]);
end

function plot_stellopt_kinetic_grad(data,vmec_data)
theta = 0:2*pi./359:2*pi;
zeta  = 0.0;
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
subplot(3,2,[1 3 5]);
ns = vmec_data.ns;
plot(r(ns,:,1),z(ns,:,1),'r','LineWidth',2);
hold on;
dex = round(2:ns/10:ns);
plot(r(1,1,1),z(1,1,1),'+k','LineWidth',2);
plot(r(dex,:,1)',z(dex,:,1)','k','LineWidth',2);
j = size(data.TE_R,1);
plot(data.TE_R(j,:),data.TE_Z(j,:),'xk','LineWidth',2,'MarkerSize',18);
plot(data.NE_R(j,:),data.NE_Z(j,:),'ok','LineWidth',2,'MarkerSize',18);
if (isfield(data,'TI_R'))
    j = size(data.TI_R,1);
    plot(data.TI_R(j,:),data.TI_Z(j,:),'sk','LineWidth',2,'MarkerSize',18);
end
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
hold off;
subplot(3,2,2);
phi = vmec_data.phi./vmec_data.phi(vmec_data.ns);
q   = 1./vmec_data.iotaf;
s = data.TE_S(j,:);
f = data.TE_target(j,:)./1000;
e = data.TE_sigma(j,:)./1000;
[s dex] = sort(s);
f = f(dex); f = f(s<=1.0);
e = e(dex); e = e(s<=1.0);
s = s(s<=1.0);
plot(s,f,'ob');
p = pchip([-0.01 s 1.01],[f(1) f 0.0]);
pgrad = [0.0 diff(ppval(p,phi))./diff(phi)];
pgrad(pgrad > 0) = 0.0;
pgrad = abs(pgrad);
pgrad = max(f)*pgrad./max(pgrad);
hold on;
plot(q, ppval(p,phi),'b');
plot(q, abs(pgrad) ,'b--');
if (isfield(data,'TI_R'))
    s = data.TI_S(j,:);
    f = data.TI_target(j,:)./1000;
    e = data.TI_sigma(j,:)./1000;
    [s dex] = sort(s);
    f = f(dex); f = f(s<=1.0);
    e = e(dex); e = e(s<=1.0);
    s = s(s<=1.0);
    plot(s,f,'or');
    p = pchip([-0.01 s 1.01],[f(1) f 0.0]);
    pgrad = [0.0 diff(ppval(p,phi))./diff(phi)];
    pgrad(pgrad > 0) = 0.0;
    pgrad = abs(pgrad);
    pgrad = max(f)*pgrad./max(pgrad);
    hold on;
    plot(q, ppval(p,phi),'r');
    plot(q, abs(pgrad) ,'r--');
end
xlim([min(q) max(q)*1.1]);
%ylim([0 1.2*max(data.TE_target(j,:)./1000)]);
ylabel('T [keV]');
hold off;
subplot(3,2,4);
s = data.NE_S(j,:);
f = data.NE_target(j,:);
e = data.NE_sigma(j,:);
[s dex] = sort(s);
f = f(dex); f = f(s<=1.0);
e = e(dex); e = e(s<=1.0);
s = s(s<=1.0);
plot(s,f,'ob');
p = pchip([-0.01 s 1.01],[f(1) f 0.0]);
pgrad = [0.0 diff(ppval(p,phi))./diff(phi)];
pgrad(pgrad > 0) = 0.0;
pgrad = abs(pgrad);
pgrad = max(f)*pgrad./max(pgrad);
hold on;
plot(q, ppval(p,phi),'b');
plot(q, abs(pgrad) ,'b--');
xlim([min(q) max(q)*1.1]);
xlabel('');
ylabel('Ne (norm)');
subplot(3,2,6);
plot(1./vmec_data.iotaf,vmec_data.presf./1000,'k','LineWidth',2,'MarkerSize',18);
hold on;
pgrad = [0.0 diff(vmec_data.presf)]./(1./vmec_data.ns)./1000;
pgrad(pgrad > 0.0) = 0.0;
pgrad = max(vmec_data.presf/1000)*pgrad / max(abs(pgrad));
plot(1./vmec_data.iotaf,abs(pgrad),'g','LineWidth',2,'MarkerSize',18);
ylabel('P [kPa]');
xlabel('Normalized Flux');
xlim([min(q) max(q)*1.1]);
%xlim([0 1.1]);
end

function plot_stellopt_te(data,vmec_data)
theta = 0:2*pi./359:2*pi;
zeta  = data.TE_PHI(1,1);
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
subplot(2,2,[1 3]);
ns = vmec_data.ns;
plot(r(ns,:,1),z(ns,:,1),'r','LineWidth',2);
hold on;
dex = round(2:ns/10:ns);
plot(r(1,1,1),z(1,1,1),'+k','LineWidth',2);
plot(r(dex,:,1)',z(dex,:,1)','k','LineWidth',2);
j = size(data.TE_R,1);
plot(data.TE_R(j,:),data.TE_Z(j,:),'+k','LineWidth',2,'MarkerSize',18);
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
hold off;
xtemp=xlim;
subplot(2,2,2);
errorbar(data.TE_R(j,:),data.TE_target(j,:)./1000,data.TE_sigma(j,:)./1000,'ok','LineWidth',2,'MarkerSize',18);
hold on;
ylim([0 1.2*max(data.TE_target(j,:)./1000)]);
plot(data.TE_R(j,:),data.TE_equil(j,:)./1000,'xr','LineWidth',2,'MarkerSize',18);
plot([1 1].*r(1,1,1),ylim,'k','LineWidth',2);
r2 = [1 1].*min(r(ns,:,1));
r2 = [r2; [1 1].*max(r(ns,:,1))];
z2 = [ylim; ylim];
plot(r2',z2','r','LineWidth',2);
ylabel('Te [keV]');
legend('Data','Recon','Axis','Edge');
hold off;
xlim(xtemp);
subplot(2,2,4);
bar(data.TE_R(j,:),sqrt(data.TE_chisq(j,:)),'r');
xlim(xtemp);
xlabel('R [m]');
ylabel('Error Norm');
end

function plot_stellopt_ne(data,vmec_data)
theta = 0:2*pi./359:2*pi;
zeta  = data.TE_PHI(1,1);
r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z     = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r     = r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z     = z + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end
subplot(2,2,[1 3]);
ns = vmec_data.ns;
plot(r(ns,:,1),z(ns,:,1),'r','LineWidth',2);
hold on;
dex = round(2:ns/10:ns);
plot(r(1,1,1),z(1,1,1),'+k','LineWidth',2);
plot(r(dex,:,1)',z(dex,:,1)','k','LineWidth',2);
j = size(data.NE_R,1);
plot(data.NE_R(j,:),data.NE_Z(j,:),'+k','LineWidth',2,'MarkerSize',18);
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
hold off;
xtemp=xlim;
subplot(2,2,2);
plot(data.NE_R(j,:),data.NE_target(j,:),'ok','LineWidth',2,'MarkerSize',18); % Plot because normalized
hold on;
ylim([0 1.2*max(data.NE_target(j,:))]);
plot(data.NE_R(j,:),data.NE_equil(j,:),'xr','LineWidth',2,'MarkerSize',18);
plot([1 1].*r(1,1,1),ylim,'k','LineWidth',2);
r2 = [1 1].*min(r(ns,:,1));
r2 = [r2; [1 1].*max(r(ns,:,1))];
z2 = [ylim; ylim];
plot(r2',z2','r','LineWidth',2);
ylabel('Ne (normalized)');
legend('Data','Recon','Axis','Edge');
hold off;
xlim(xtemp);
subplot(2,2,4);
bar(data.NE_R(j,:),sqrt(data.NE_chisq(j,:)),'r');
xlim(xtemp);
xlabel('R [m]');
ylabel('Error Norm');
end

function plot_stellopt_output(data)
% Plot chi^2
axh=subplot(2,1,1);
set(axh,'NextPlot','add');
colormap('lines');
colors=get(gca,'ColorOrder');
names=['Total';data.names(:)];
plot(axh,data.func_eval,data.chi_total,'o','Color','black',...
    'MarkerSize',10);
for i=1:data.numels
    name=strcat('chisq_',data.struc_names{i});
    plot(axh,data.func_eval,data.(name),'.','Color',colors(i,:));
end
set(gca,'Yscale','log');
ylabel('\chi^2');
legend(names,'Location','NorthEastOutside');
ylim([0.01 max(data.chi_total)]);
% Plot RMS ERROR
subplot(2,1,2);
hold on
for i=1:data.numels
    name=strcat('rms_',data.struc_names{i});
    plot(data.func_eval,data.(name),'Color',colors(i,:));
end
set(gca,'Yscale','log');
hold off
ylabel('ERROR');
return
end

function plot_stellopt_pressure(data)
% Pressure plots
% LHD p=p0*(1-phi)^2 (phi=sqrt(s))
subplot(3,1,[1 2]);
%Find VMEC DOMAIN
dex1=find(data.data(:,4)<1,1,'first');
dex2=find(data.data(:,4)<1,1,'last');
x1=data.data(dex1,1);
x2=data.data(dex2,1);
y1=data.norm*max(data.data(:,7))*.01;
sigmas=data.data(:,9);
np_prof_total=length(sigmas);
%sigmas=sigmas./sqrt(np_prof_total).*data.norm;
sigmas=sigmas.*data.norm;
%Plot
hp1=errorbar(data.data(:,1),data.data(:,7).*data.norm,sigmas,'ok');
%hp1=plot(data.data(:,1),data.data(:,5),'ok');
hold on
fill([x1 x2 x2 x1],[0 0 y1 y1],'red');
%plot(data.data(dex1:dex2,1),vmec_line,'r');
hp2=plot(data.data(:,1),data.data(:,8),'k','LineWidth',2);
flux=sqrt(data.data(:,4));
p_lhd=max(data.data((flux<1),8)).*(1-flux.^2);
p_lhd(flux>1)=0.0;
%hp3=plot(data.data(:,1),p_lhd,'blue');
text((x2-x1)/3+x1,y1*5,'VMEC Domain');
hold off
title('Pressure Profile');
%xlabel('R [m]');
ylabel('Pressure [Pa]');
ylim([0 max([max(data.data(:,7).*data.norm) max(data.data(:,8))])]);
%legend([hp1,hp2,hp3],'Experiment','Simulation','LHD Profile (1-\Phi^2)');
legend([hp1,hp2],'Experiment','Simulation');
subplot(3,1,3);
wegt=(data.data(:,8)-data.data(:,7).*data.norm)./sigmas;
reddex=abs(wegt)>4.0;
yeldex=(abs(wegt)>1.0).*(abs(wegt)<=4.0);
gredex=abs(wegt)<=1.0;
if (max(reddex) > 0)
    bar(data.data(:,1),wegt.*reddex,'r','EdgeColor','red');
    hold on;
end
if (max(yeldex) > 0)
    bar(data.data(:,1),wegt.*yeldex,'y','EdgeColor','yellow');
    hold on;
end
if (max(gredex) > 0)
    bar(data.data(:,1),wegt.*gredex,'g','EdgeColor','green');
    hold on;
end
plot(data.data(:,1),data.data(:,1).*0,'k'); % Black Axis Line
hold off;
ylim([-5 5]);
title('Weighted Deviation')
xlabel('R[m]');
ylabel('(p_{sim}-p_{exp})/\sigma');
return
end

function plot_stellopt_neteti(data)
% Thomson plots
subplot(3,1,1)
if isfield(data,'ne_data')
    title('Profile Diagnostics');
    mask=data.ne_data(:,8)~=0.0;
    hp1=errorbar(data.ne_data(mask,4),data.ne_data(mask,5),data.ne_data(mask,7),'xk');
    ylabel('n_e [m^-3]');
    hold on
    hp2=plot(data.ne_data(mask,4),data.ne_data(mask,6),'ro');
    ylim([0 1.2*max([data.ne_data(mask,5); data.ne_data(mask,6)])]);
    xlim([0 1.2]);
end
if isfield(data,'te_data')
    subplot(3,1,2)
    mask=data.te_data(:,8)~=0.0;
    hp1=errorbar(data.te_data(mask,4),data.te_data(mask,5),data.te_data(mask,7),'xk');
    ylabel('T_e [eV]');
    hold on
    hp2=plot(data.te_data(mask,4),data.te_data(mask,6),'ro');
    ylim([0 1.2*max([data.te_data(mask,5); data.te_data(mask,6)])]);
    xlim([0 1.2]);
end
if isfield(data,'ti_data')
    subplot(3,1,3)
    mask=data.ti_data(:,8)~=0.0;
    hp1=errorbar(data.ti_data(mask,4),data.ti_data(mask,5),data.ti_data(mask,7),'xk');
    ylabel('T_i [eV]');
    xlabel('Normalized Toroidal Flux');
    hold on
    hp2=plot(data.ti_data(mask,4),data.ti_data(mask,6),'ro');
    ylim([0 1.2*max([data.ti_data(mask,5); data.ti_data(mask,6)])]);
    xlim([0 1.2]);
end
end

function plot_stellopt_magdiag(data)
% Magnetic Diagnostics Plot

% Plot simulation values against actual values
subplot(3,1,[1 2]);
hp1=errorbar(1:size(data.data,1),data.data(:,2),data.data(:,3),'kx');  % Measured
set(gca,'XTick',1:size(data.data,1));
%set(gca,'XTickLabel',data.names);
hold on
plot(1:size(data.data,1),data.data(:,1),'or');
plot(0:(size(data.data,1)+1),zeros(1,size(data.data,1)+2),'r--');
title('Magnetic Diagnostics');
ylabel('Signals');
xlim([0 size(data.data,1)+1]);
hold off
subplot(3,1,3);
% Plot residuals
dex1 = abs(data.data(:,4)) <= 1;
dex2 = abs(data.data(:,4)) > 1;
dex3 = abs(data.data(:,4)) <= 3;
dex4 = abs(data.data(:,4)) > 3;
bar(1:size(data.data,1),abs(data.data(:,4)).*dex1,'g');
hold on
bar(1:size(data.data,1),abs(data.data(:,4)).*dex2.*dex3,'y');
bar(1:size(data.data,1),abs(data.data(:,4)).*dex4,'r');
hold off
ylabel('Deviation');
set(gca,'XTick',1:size(data.data,1));
%set(gca,'XTickLabel',data.names);
xlim([0 size(data.data,1)+1]);
return
end

function plot_stellopt_mse(data,vmec_data)
% MSE Plot
% LHD
geo_x_min=-2.0; geo_x_max=0.0; geo_y_min=-5.0; geo_y_max=-2.5;
phimin = pi;
phimax = 3.*pi/4.;
nphi=90;
% DIII-D
geo_x_min=-3.0; geo_x_max=3.0; geo_y_min=-3.0; geo_y_max=3.0;
phimin = 0.0;
phimax = 2.*pi;
nphi = 360.;
% Make Chord Plot
mse_r=data.data(:,1);
mse_z=data.data(:,2);
mse_phi=data.data(:,3);
mse_s=data.data(:,4);
mse_pol=data.data(:,5);
mse_vac=data.data(:,6);
mse_vmec=data.data(:,7);
mse_sigma=data.data(:,8);
mse_wegt=(mse_pol-mse_vmec)./mse_sigma;
if ~isempty(vmec_data)
    subplot(2,2,[1]);
    %phimin=(round(min(mse_phi)/(0.5*pi))-1)*pi/2;
    %phimax=phimin+pi/2;
    phi=phimin:(phimax-phimin)./(nphi-1):phimax;
    theta=[0 pi+pi/10.];
    rho=round(2:(vmec_data.ns-2)/9:vmec_data.ns);
    r=cfunct(theta,phi,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    cosphi=cos(phi);
    sinphi=sin(phi);
    % Now we need to find zeros of function
    for u=1:nphi
        for s=1:length(rho)
            f = @(x) sfunct(x,phi(u),vmec_data.zmns(:,rho(s)),vmec_data.xm,vmec_data.xn);
            theta1=fzero(f,[-pi/2 pi/2]);
            theta2=fzero(f,[pi/2 3*pi/2]);
            r_temp=cfunct([theta1 theta2],phi(u),vmec_data.rmnc(:,rho(s)),vmec_data.xm,vmec_data.xn);
            r(rho(s),1,u)=r_temp(1);
            r(rho(s),2,u)=r_temp(2);
        end
    end
    % Now make plots
    plot(squeeze(r(1,1,:)).*cosphi(:),squeeze(r(1,1,:)).*sinphi(:),'r--');
    hold on
    for i=1:length(rho)-1
        plot(squeeze(r(rho(i),1,:)).*cosphi(:),squeeze(r(rho(i),1,:)).*sinphi(:),'r');
        plot(squeeze(r(rho(i),2,:)).*cosphi(:),squeeze(r(rho(i),2,:)).*sinphi(:),'r');
    end
    i=length(rho);
    plot(squeeze(r(rho(i),1,:)).*cosphi(:),squeeze(r(rho(i),1,:)).*sinphi(:),'r','LineWidth',2);
    plot(squeeze(r(rho(i),2,:)).*cosphi(:),squeeze(r(rho(i),2,:)).*sinphi(:),'r','LineWidth',2);
    plot3(mse_r.*cos(mse_phi),...
        mse_r.*sin(mse_phi),...
        mse_z,'ok');
    %axis equal
    xlim([geo_x_min geo_x_max]);
    ylim([geo_y_min geo_y_max]);
    %axis tight
    xlabel('X [m]');
    ylabel('Y [m]');
    title('MSE Geometry');
    subplot(2,2,3);
    [ax,handle_j,handle_iota]=plotyy(vmec_data.phi./vmec_data.phi(vmec_data.ns),vmec_data.jcurv./1000,vmec_data.phi./vmec_data.phi(vmec_data.ns),vmec_data.iotaf);
    if strcmp(strtrim(vmec_data.pcurrtype),'akima_spline_ip') || strcmp(strtrim(vmec_data.pcurrtype),'akima_spline')
        hold on
        %plotyy(vmec_data.phi,vmec_data.jcurv./1000,vmec_data.phi,vmec_data.iotaf);
        plot(vmec_data.acauxs(vmec_data.acauxs>0),vmec_data.acauxs(vmec_data.acauxs>0).*0.0,'+k')
        plot(vmec_data.acauxs(1),[0.0],'+k')
        dex = (mse_s > 1./vmec_data.ns);
        plot(mse_s(dex),mse_s(dex).*0.0,'ob');
        dex = mse_s < 1./vmec_data.ns;
        plot(mse_s(dex),mse_s(dex).*0.0,'or');
        hold off
    end
    title('Current/Iota Profiles');
    xlabel('Normalized Toroidal Flux');
    ylabel(ax(1),'Current Profile [kA]');
    ylabel(ax(2),'Rotational Transform');
    axis tight;
    axes(ax(2));
    axis tight;
    %set(handle_j,'YColor','black');
    %set(handle_iota,'YColor','black');
    subplot(2,2,2);
else
    subplot(2,1,1);
end
dex = find(mse_sigma < 1.e10);
mse_precer=mse_sigma(dex)./mse_pol(dex);
mse_pol_deg=180.*mse_pol(dex)./pi;
mse_sigma_deg=mse_precer.*mse_pol_deg;
hp1=errorbar(data.data(dex,4),mse_pol_deg,mse_sigma_deg,'kx');  % Measured
%hp1=errorbar(data.data(dex,4),180*mse_pol(dex)/pi,180*mse_sigma(dex)/pi/10.,'kx');  % Measured
%hp1=errorbar(1:size(data.data,1),180*mse_pol/pi,180*mse_sigma/pi,'kx');  % Measured
%set(gca,'XTick',1:size(data.data,1));
hold on
title('MSE Polarization');
ylabel('Degrees');
%xlim([0 size(data.data,1)+1]);
xlim([0 1]);
if ~isempty(vmec_data)
    mse_dex=data.data(:,4);
    dex = mse_s > 1./vmec_data.ns;
    if ~isempty(dex), plot(mse_dex(dex),180*mse_vmec(dex)/pi,'ob'); end
    dex = mse_s < 1./vmec_data.ns;
    if ~isempty(dex), plot(mse_dex(dex),180*mse_vmec(dex)/pi,'or'); end
    dex = mse_s > 1;
    if ~isempty(dex), plot(mse_dex(dex),180*mse_vmec(dex)/pi,'or'); end
    subplot(2,2,4);
else
    plot(1:size(data.data,1),180*mse_vmec/pi,'or');
    subplot(2,1,2);
end
% Plot residuals
%mse_error=abs(data.data(:,5)-data.data(:,6))./data.data(:,5);
%bar(1:size(data.data,1),mse_error*100.);
dev = (data.data(:,5)-data.data(:,7))./data.data(:,9);
dex_5 = abs(dev) >= 5.0;
dex_2 = and(abs(dev) < 5.0, abs(dev) > 1.0);
dex_1 = abs(dev) <= 1.0;
x_temp=1:size(data.data,1);
bar(x_temp,dev.*dex_5,'r','EdgeColor','none');
hold on;
bar(x_temp,dev.*dex_2,'y','EdgeColor','none');
bar(x_temp,dev.*dex_1,'r','EdgeColor','none');
%ylabel('% Difference');
ylabel('Deviation');
xlabel('MSE Datapoint');
set(gca,'XTick',1:size(data.data,1));
xlim([0 size(data.data,1)+1]);
ylim([-5 5]);
return
end

function plot_stellopt_jacobian(data)
%Jacobian Plot
pixplot(data.jac);
% Sort x ticks
dex=1;
temp_name=data.var_name{1};
var_label=zeros(1,data.nvars);
for i=1:data.nvars
    if strcmp(temp_name,data.var_name{i})
        var_label(i)=dex;
    else
        dex=dex+1;
        var_label(i)=dex;
        temp_name=data.var_name{i};
    end
end
xtick=cell(1,data.nvars);
hold on;
for i=1:max(var_label)
    offset=sum(var_label<i);
    plot([1 1].*offset+1,[0 data.mtargets+1],'k');
    dex=round(sum(var_label==i)./2);
    %xtick{offset+dex}=data.var_name{offset+dex};
    text(offset+dex+0.5,-0.5,strjust(deblank(data.var_name{offset+dex}),'right'),'Rotation',90,...
        'HorizontalAlignment','right','Interpreter','none');
end
set(gca,'XTick',1:data.nvars+0.5);
set(gca,'XTickLabel',xtick);
% Sort Y ticks
dex=1;
temp_name=data.target_name{1};
target_label=zeros(1,data.mtargets);
for i=1:data.mtargets
    if strcmp(temp_name,data.target_name{i})
        target_label(i)=dex;
    else
        dex=dex+1;
        target_label(i)=dex;
        temp_name=data.target_name{i};
    end
end
ytick=cell(1,data.mtargets);
hold on;
for i=1:max(target_label)
    offset=sum(target_label<i);
    plot([0 data.nvars+1],[1 1].*offset+1,'k');
    dex=round(sum(target_label==i)./2);
    ytick{offset+dex}=data.target_name{offset+dex};
end
set(gca,'YTick',1:data.mtargets);
set(gca,'YTickLabel',ytick);
caxis([-1 1].*abs(min(caxis)));
colorbar;
title('STELLOPT Jacobian');
return
end

