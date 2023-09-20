function [I,s,t]=vmec_thrift(vmec_data,varargin)
%VMEC_THRIFT Computes the evolution the current profile
%   The VMEC_THRIFT routine calculates the current profile evolution
%   assuming no feedback on the equilbrium.

% Constants
ec=1.60217663E-19;
Temax = 5000; % only used if no Te provided

%Initializations
t = [];
te_prof=[];
ne_prof=[];
zeff_prof=[];
jdotb_prof=[];
lplot = 1;

% Input arguments
if nargin > 1
    i = 1;
    while i < length(varargin)
        switch varargin{i}
            case 't'
                i = i + 1;
                t = varargin{i};
            case 'te'
                i = i + 1;
                te_x = varargin{i};
                i = i + 1;
                te_t = varargin{i};
                i = i + 1;
                te_prof = varargin{i};
            case 'ne'
                i = i + 1;
                ne_x = varargin{i};
                i = i + 1;
                ne_t = varargin{i};
                i = i + 1;
                ne_prof = varargin{i};
            case 'zeff'
                i = i + 1;
                zeff_x = varargin{i};
                i = i + 1;
                zeff_t = varargin{i};
                i = i + 1;
                zeff_prof = varargin{i};
            case 'jdotb'
                i = i + 1;
                jdotb_x = varargin{i};
                i = i + 1;
                jdotb_t = varargin{i};
                i = i + 1;
                jdotb_prof = varargin{i};
        end
        i=i+1;
    end
end

% time grid
if isempty(t)
    t=0:0.05:10.0;
end

% Calculate some helpers
[s11,s12,s21,s22]=calc_susceptance(vmec_data);
nu = 64;
nv = 36.*vmec_data.nfp;
theta = 0:2*pi/(nu-1):2*pi;
zeta = 0:2*pi/(nv-1):2*pi;
g  = cfunct(theta,zeta,vmec_data.gmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
b  = cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
if vmec_data.iasym
     g =  g + sfunct(theta,zeta,vmec_data.gmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
     b =  b + sfunct(theta,zeta,vmec_data.bmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
end
Bsqav = -trapz(theta,trapz(zeta,(b.*b).*g,3),2)./(vmec_data.vp'.*4.*pi.*pi);
Bav = -trapz(theta,trapz(zeta,(b).*g,3),2)./(vmec_data.vp'.*4.*pi.*pi);
phia = vmec_data.phi(end);

% Electron Temperatures
if isempty(te_prof)
    te_x = 0:0.01:1;
    te_t = t;
    [X1,X2]=ndgrid(te_x,te_t);
    s_temp = linspace(0,1,vmec_data.ns);
    f_temp = @(x) pchip(s_temp,Temax.*vmec_data.presf./max(vmec_data.presf),x);
    te_prof = f_temp(X1);
else
    [X1,X2]=ndgrid(te_x,te_t);
end
F_te = griddedInterpolant(X1,X2,te_prof,'makima','nearest');

% Electron Density
if isempty(ne_prof)
    ne_x = 0:0.01:1;
    ne_t = t;
    [X1,X2]=ndgrid(ne_x,ne_t);
    s_temp = linspace(0,1,vmec_data.ns);
    f_temp = @(x) pchip(s_temp,ec.*Temax.*vmec_data.presf./max(vmec_data.presf),x);
    f2_temp = @(x) pchip(s_temp,vmec_data.presf./f_temp(s_temp),x);
    ne_prof = f2_temp(X1);
else
    [X1,X2]=ndgrid(ne_x,ne_t);
end
F_ne = griddedInterpolant(X1,X2,ne_prof,'makima','nearest');

% Z-Effective
if isempty(zeff_prof)
    zeff_x = 0:0.01:1;
    zeff_t = t;
    [X1,X2]=ndgrid(zeff_x,zeff_t);
    zeff_prof = X1.*0.0+1.0;
else
    [X1,X2]=ndgrid(zeff_x,zeff_t);
end
F_zeff = griddedInterpolant(X1,X2,zeff_prof,'makima','nearest');

% J.B
if isempty(jdotb_prof)
    jdotb_x = 0:0.01:1;
    jdotb_t = t;
    [X1,X2]=ndgrid(jdotb_x,jdotb_t);
    s_temp = linspace(0,1,vmec_data.ns);
    %f_eps = @(x) sqrt(x).*vmec_data.Aminor./vmec_data.Rmajor;
    %f_phip = @(x) pchip(s_temp,vmec_data.phipf,x).*sqrt(x).*2;
    %f_p = @(x) pchip(linspace(-1,1,2.*vmec_data.ns-1),[vmec_data.presf(end:-1:2) vmec_data.presf],x);
    %f_pprime = @(x) (f_p(x+0.0005)-f_p(x-0.0005)).*1000;
    %f_j = @(x) sqrt(f_eps(x)).*f_pprime(x).*vmec_data.Rmajor./f_phip(x);
    f_j = @(x) pchip(s_temp,vmec_data.jcurv,x);
    f_jdotb = @(x) f_j(x).*pchip(s_temp,Bav,x);
    jdotb_prof = f_jdotb(X1);
else
    [X1,X2]=ndgrid(jdotb_x,jdotb_t);
end
F_jdotb = griddedInterpolant(X1,X2,jdotb_prof,'makima','nearest');

% Flux grid
s = 0:1./(vmec_data.ns-1):1;
s2 = -1:1./(vmec_data.ns-1):1;
rho = sqrt(s);

% Calculate pprime
p2 = [vmec_data.presf(end:-1:2) vmec_data.presf];
h = 1./(vmec_data.ns-1);
pprime = (pchip(s2,p2,s+h)-pchip(s2,p2,s-h))./2.*h;


clog = @(te,ne) real(24-log( sqrt(ne.*1E-6)./(te) ));
etaperp = @(te,ne,zeff) 1.0313621201E-04.*zeff.*clog(ne,te).*max(te,14).^-1.5;
F = @(zeff) (1+1.198.*zeff+0.222.*zeff.*zeff)./(1+2.966.*zeff+0.752.*zeff.*zeff);
etapara = @(te,ne,zeff) F(zeff).*etaperp(te,ne,zeff);
etafunc = @(x,t) etapara(F_te(x,t),F_ne(x,t),F_zeff(x,t)).*F(F_zeff(x,t));

% Neoclassical resistivity
% https://crppwww.epfl.ch/~sauter/neoclassical/
[~,ft] = vmec_fraction(vmec_data);
F_trapped = @(x) pchip(s,ft,x);
[X1,X2]=ndgrid(s,t);
frac_trapped = F_trapped(X1);
ne_temp = F_ne(X1,X2);
te_temp = F_te(X1,X2);
Ze_temp = F_zeff(X1,X2);
[signeo, sigspitzer] = sigmaneo(frac_trapped,ne_temp,max(te_temp,14));
F_etaneo = griddedInterpolant(X1,X2,1./abs(signeo),'makima','nearest');

% if (lplot)
%     fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardcopy','off');
% end




% Spline helpers
s11_spl = pchip(s,abs(s11));
vp_spl = pchip(s,vmec_data.vp);
pprime_spl = pchip(s,pprime);
Bsqav_spl = pchip(s,Bsqav);

% Extra argument helpers
thriftpde2 = @(x,t,u,dudx) thriftpde(x,t,u,dudx,s11_spl,vp_spl,pprime_spl,Bsqav_spl,F_jdotb,F_etaneo,phia);
thriftic2 = @(x) thriftic(x);
thriftbc2 = @(xl,ul,xr,ur,t) thriftbc(xl,ul,xr,ur,t);
utoI = @(u) u.*phia./(pi.*4E-7);

% Calculate solution
m = 0;
sol = pdepe(m,thriftpde2,thriftic2,thriftbc2,s,t);
I = utoI(sol);

% Calculate derivative quantities
x2 = [-s(2) s 2.*s(end)-s(end-1)];
x2_h = 0.5.*(x2(2:end)+x2(1:end-1));
dIds = zeros(size(I));
iota = zeros(size(I));
j = zeros(size(I));
A = pi.*s*vmec_data.Aminor.*vmec_data.Aminor;
dAds = pi.*vmec_data.Aminor.*vmec_data.Aminor;
for i = 1:length(t)
    I2 = [I(i,2) I(i,:) 2.*I(i,end)-I(i,end-1)];
    I_spl = pchip(x2,I2);
    dIds(i,:) = diff(ppval(I_spl,x2_h))./(s(2)-s(1));
    j(i,:) = dIds(i,:)./dAds;
    iota(i,:) = (sol(i,:)-s12')./abs(s11)';
end

% Finish
if lplot
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardcopy','off');
    subplot(2,2,1);
    surf(s,t,I,'LineStyle','none');
    set(gca,'FontSize',24);
    xlabel('Norm. Tor. Flux (s)');
    ylabel('Time [s]');
    zlabel('I [A]');
    view(3);
    subplot(2,2,2);
    plot(t,I(:,end),'k','LineWidth',4);
    set(gca,'FontSize',24);
    xlabel('Time [s]');
    ylabel('Total Current [A]');
    x = t'; y = I(:,end);
    g = fittype('a-b.*exp(-x./c)');
    f0 = fit(x,y,g,'StartPoint',[I(end,end),I(end,end),20])
    hold on;
    plot(x,f0(x),'--r');
    legend('THRIFT','A-B*(-EXP/C)','Location','southeast')
    subplot(2,2,3);
    surf(s,t,j,'LineStyle','none');
    set(gca,'FontSize',24);
    xlabel('Norm. Tor. Flux (s)');
    ylabel('Time [s]');
    zlabel('j [A/m^2]');
    view(3);
    subplot(2,2,4);
    px = round(linspace(2,length(s),4));
    px(px>length(s)) = length(s);
    plot(t,iota(:,px),'LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Time [s]');
    ylabel('Rot. Trans. \iota');
end

end

function [c,f,s] = thriftpde(x,t,u,dudx,s11_spl,vp_spl,pprime_spl,Bsqav_spl,jdotbfunc,etafunc,phia)
mu0 = pi*4E-7;
s11_local = ppval(s11_spl,x);
vp_local = ppval(vp_spl,x);
pprime_local = ppval(pprime_spl,x);
Bsqav_local = ppval(Bsqav_spl,x);
jdotb_local = jdotbfunc(x,t);
eta_local=etafunc(x,t);
c = phia.*phia./s11_local;
f = (pprime_local.*u + Bsqav_local.*dudx./mu0 - jdotb_local).*vp_local.*eta_local;
s = 0;
end

function u0 = thriftic(x)
u0 = 0.0;
end

function [pl,ql,pr,qr] = thriftbc(xl,ul,xr,ur,t)
pl = ul;
ql = 0.0;
pr = 0;
qr = 1.0;
end


function [signeo, sigspitzer] = sigmaneo(ft,ne,te,varargin)
%
%   [signeo, sigspitzer, nuestar] = sigmaneo(ft,ne,te,varargin)
%
% Compute neoclassical conductivity using formula in Ref.1: O. Sauter et al, Phys. Plasmas 6 (1999) 2834.
%
% dimensions from size(ft) in input
%
% Outputs:
%    signeo : neoclassical conductivity at each value of the input parameters
%    sigsptz: Spitzer conductivity
%
% Inputs:
%  Required:
%    ft()   : trapped fraction (Note, can use ftav.m with formula in Lin-Liu et al, Phys. Plasmas 2 (1995) 1666.)
%    ne()   : local electron density
%    Te()   : Local electron temperature
%  Optionals:
%    varargin{1}: zeff: effective charge (default: 3)
%    varargin{2}: nuestar: electron collisionality (default set to 0.)
%               can be calculated from nustar.m
%
% examples of calls:
%
% [signeo,sigsptz]=sigmaneo(ft,ne,te) (defaults : zeff=3, nuestar=0)
% [signeo,sigsptz]=sigmaneo(ft,ne,te,2.5.*ft./ft) (defaults : nuestar=0)
% [signeo,sigsptz]=sigmaneo(ft,ne,te,zeff,nuestar)
%
%

nuestar=zeros(size(ft));
zeff=3. .* ones(size(ft));

nargeff=nargin-3;
if nargeff > 0
  if ~isempty(varargin{1}); zeff=varargin{1}; end
end
if nargeff > 1
  if ~isempty(varargin{2}); nuestar=varargin{2}; end
end

NZ = 0.58 + 0.74 ./ (0.76 + zeff);
lnLam = 17.*ones(size(ft));
ii=find(ne>0 & te>0);
if length(ii)>0
  lnLam = 31.3 - log(sqrt(ne)./te);
end

sigspitzer = 1.9012E+04 .* te.^1.5 ./ zeff ./ NZ ./ lnLam;

ft33eff = ft ./ (1. + (0.55-0.1.*ft).*sqrt(nuestar) + 0.45.*(1.-ft).*nuestar./zeff.^1.5);
signeo = sigspitzer .* (1. - ft33eff.*(1.+0.36./zeff - ft33eff.*(0.59./zeff - 0.23./zeff.*ft33eff)));
end