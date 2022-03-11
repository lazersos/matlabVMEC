function data = eqdisk2vmec( varargin )
%EQDSK2VMEC Converts EQDSK files to VMEC INDATA namelists
%   The EQDSK2VMEC routine takes an EQDSK file and computes a VMEC INDATA
%   namelist.  Plots showing the fitting routines used are also made.  
%   This is a work in progress

if (nargin == 1)
    if isstr(varargin{1})
        filename = varargin{1};
        data = eqdisk2vmec_file(varargin{1});
    end
elseif (nargin >= 2)
    i = 1;
    while i <= nargin
        switch varargin{i}
            case{'NSTX'}
                i=i+1;
                shot_num = varargin{i};
                i=i+1;
                time = varargin{i};
                data = eqdisk2vmec_NSTX(shot_num,time);
                filename = ['NSTX' num2str(shot_num,'%6.6d') '_' num2str(time*1000,'%4.4d')];
            case{'DIIID'}
        end
        i=i+1;
    end
end
% Now output vmec data
vmec_input = vmec_namelist_init('indata');
vmec_input.delt = 1.0;
vmec_input.niter = 20000;
vmec_input.tcon0 = 1.0;
vmec_input.ns_array = [16 32 64 128];
vmec_input.ftol_array = [1e-30  1e-30  1e-30  1e-12];
vmec_input.niter_array = [1000  2000  4000  20000];
vmec_input.lasym = 1;
vmec_input.nfp =1;
vmec_input.mpol = max(size(data.refou));
vmec_input.ntor = 0;
vmec_input.nzeta = 1;
vmec_input.nstep = 200;
vmec_input.ntheta = 2*vmec_input.mpol+6;
vmec_input.phiedge = data.phiedge;
vmec_input.lfreeb = 0;
vmec_input.mgrid_file = '';
vmec_input.extcur =[0 0 0];
vmec_input.nvacskip = 6;
vmec_input.gamma = 0.0;
vmec_input.bloat = 1.0;
vmec_input=rmfield(vmec_input,'bcrit');
vmec_input.spres_ped = 1.0;
vmec_input.pres_scale = 1.0;
vmec_input.pmass_type = 'akima_spline';
vmec_input.am = data.am;
vmec_input.am_aux_s = data.am_aux_s;
vmec_input.am_aux_f = data.am_aux_f;
vmec_input.pcurr_type = 'akima_spline_ip';
vmec_input.curtor = data.curtor;
vmec_input.ac = data.ac;
vmec_input.ac_aux_s = data.ac_aux_s;
vmec_input.ac_aux_f = data.ac_aux_f;
vmec_input.piota_type = 'akima_spline';
vmec_input.ai = data.ai;
vmec_input.ai_aux_s = data.ai_aux_s;
vmec_input.ai_aux_f = data.ai_aux_f;
vmec_input.raxis_cc = data.raxis;
vmec_input.raxis_cs = 0.0;
vmec_input.zaxis_cc = data.zaxis;
vmec_input.zaxis_cs = 0.0;
vmec_input.rbc = data.refou';
vmec_input.zbs = data.zefou';
vmec_input.rbs = data.refou2';
vmec_input.zbc = data.zefou2';
vmec_input=rmfield(vmec_input,'raxis');
vmec_input=rmfield(vmec_input,'zaxis');
dex = strfind(filename,'.');
if ~isempty(dex)
    filename = filename(1:(dex(1)-1));
end
new_file = ['input.' filename];
write_vmec_input(new_file,vmec_input);
new_file = [filename '.fig'];
saveas(gcf,new_file);




end

function  data2 = eqdisk2vmec_file(filename)
data=read_efit(filename);
subplot(3,2,[1 3 5]);
contourf(data.xgrid,data.zgrid,data.psixz',100); axis equal; axis tight; 
hold on; 
plot(data.xbndry,data.zbndry,'r','LineWidth',2); 
plot(data.xaxis,data.zaxis,'+r','LineWidth',2); 
plot(data.xlim,data.zlim,'k','LineWidth',2); hold off;
xlabel('R'); ylabel('Z'); title('GEQDSK');
mpol = 12;
ntor = 0;
R = data.xbndry;
Z = data.zbndry;
np = length(R);
rmaj = mean(R);
Rminor = R - rmaj;
theta = atan2(Z,Rminor);
theta(theta<0) = theta(theta<0) + 2*pi;
[theta, dex] =sort(theta);
R = R(dex);
Z = Z(dex);
Rminor = Rminor(dex);
rho = sqrt(Rminor.^2+Z.^2);
nu = np;
nv = 1;
refou=zeros(mpol+1,2*ntor+1);
zefou=zeros(mpol+1,2*ntor+1);
refou2=zeros(mpol+1,2*ntor+1);
zefou2=zeros(mpol+1,2*ntor+1);
cosu=zeros(nu,mpol+1);
cosv=zeros(nv,2*ntor+1);
sinu=zeros(nu,mpol+1);
sinv=zeros(nv,2*ntor+1);
alu=2*pi/nu;
for i=1:nu
    for j=1:mpol+1
        m=j-1;
        cosu(i,j)=cos(m*(i-1)*alu);
        sinu(i,j)=sin(m*(i-1)*alu);
    end
end
alv=2*pi/nv;
for i=1:nv
    for j=1:2*ntor+1
        n=j-ntor-1;
        cosv(i,j)=cos(n*(i-1)*alv);
        sinv(i,j)=sin(n*(i-1)*alv);
    end
end
% Our Amplitude is 1 everywhere since we do a full period
% Our Kernels is
%    read_VMEC arrays
%    cos(mu+nv) = cos(mu)cos(nv)-sin(mu)sin(nv)
%    sin(mu+nv) = sin(mu)cos(nv)+cos(mu)sin(nv)
fnuv=zeros(1,mpol+1);
fnuv(1)=1./(nu*nv);
for i=2:mpol+1
    fnuv(i)=2*fnuv(1);
end
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                refou(m1,n1)=refou(m1,n1)+R(i)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou(m1,n1)=zefou(m1,n1)+Z(i)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                refou2(m1,n1)=refou2(m1,n1)+R(i)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou2(m1,n1)=zefou2(m1,n1)+Z(i)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end
data2.refou = refou;
data2.zefou = zefou;
data2.refou2 = refou2;
data2.zefou2 = zefou2;
% Now handle profiles
pflux=data.psiaxis:(data.psilim-data.psiaxis)/(length(data.sf)-1):data.psilim;
pflux_norm=pflux-pflux(1);
pflux_norm=pflux_norm./max(pflux_norm);
jdotb=data.spp+data.sffp;
press=data.sp;
qprof=data.qpsi;
iotaf=1./qprof;
% dphi = q * dpsi
q_spl=pchip(pflux',qprof);
fun = @(x) ppval(q_spl,x);
for i=1:length(pflux)
    tflux(i) = integral(fun,pflux(1),pflux(i));
end
phiedge = tflux(end).*2.*pi.*sign(data.btor);
tflux = tflux./tflux(length(tflux));
s = 0:1/99:1;
am_aux_s = s;
am_aux_f = pchip(tflux,press,s);
p = polyfit(am_aux_s,am_aux_f,9);
p(length(p)) = am_aux_f(1);
p = [-sum(p(1:10))+am_aux_f(100) p];
am = fliplr(p);
subplot(3,2,2); plot(tflux,press,'k'); title('Pressure [Pa]'); hold on; plot(am_aux_s,am_aux_f,'or'); plot(am_aux_s,polyval(p,am_aux_s),'b'); legend('EQDSK','SPLINE','POLY'); hold off
ac_aux_s = s;
ac_aux_f = pchip(tflux,jdotb,s);
p = polyfit(ac_aux_s,ac_aux_f,9);
p(length(p)) = ac_aux_f(1);
p = [-sum(p(1:10))+ac_aux_f(100) p];
ac = fliplr(p);
subplot(3,2,4); plot(tflux,jdotb,'k'); title('<J.B>/<B/R> [A]'); hold on; plot(ac_aux_s,ac_aux_f,'or'); plot(ac_aux_s,polyval(p,ac_aux_s),'b'); hold off
ai_aux_s = s;
ai_aux_f = pchip(tflux,iotaf,s);
p = polyfit(ai_aux_s,ai_aux_f,9);
p(length(p)) = ai_aux_f(1);
p = [-sum(p(1:10))+ai_aux_f(100) p];
ai = fliplr(p);
subplot(3,2,6); plot(tflux,qprof,'k'); title('q'); hold on; plot(ai_aux_s,1./ai_aux_f,'or'); plot(ai_aux_s,1./polyval(p,ai_aux_s),'b'); hold off
data2.phiedge = phiedge;
data2.curtor = data.totcur;
data2.raxis = data.xaxis;
data2.zaxis = data.zaxis;
data2.am = am;
data2.am_aux_s = am_aux_s;
data2.am_aux_f = am_aux_f;
data2.ac = ac;
data2.ac_aux_s = ac_aux_s;
data2.ac_aux_f = ac_aux_f;
data2.ai = ai;
data2.ai_aux_s = ai_aux_s;
data2.ai_aux_f = ai_aux_f;
% Check the data
mn = 1;
for m1 = 1: mpol + 1
    for n1=1: 2*ntor +1
        xm(mn) = m1-1;
        xn(mn) = n1-ntor-1;
        rmnc(mn,1) = refou(m1,n1);
        rmns(mn,1) = refou2(m1,n1);
        zmnc(mn,1) = zefou2(m1,n1);
        zmns(mn,1) = zefou(m1,n1);
        mn = mn + 1;
    end
end
ntheta=360;
theta=0:2*pi/(ntheta-1):2*pi;
zeta = 0;
r=cfunct(theta,zeta,rmnc,xm,xn)+sfunct(theta,zeta,rmns,xm,xn);
z=cfunct(theta,zeta,zmnc,xm,xn)+sfunct(theta,zeta,zmns,xm,xn);
subplot(3,2,[1 3 5]);
hold on
plot(r,z,'b');
hold off

return
end

function data = eqdisk2vmec_NSTX(shotnum,time)
% Open MDSPlus connection
import MdsPlus.*;   % So we can call MdsPlus java calls
server = 'skylark.pppl.gov';  % Server you wish to connect to
server_port = 8501;        % Server port
mds_server=MdsPlus(server,server_port);  % Open the connection
mds_server.OpenTree('NSTX',shotnum);  % Open the Tree
% Get Data
type = '\EFIT03';
try
    efit_nbdry = mds_server.Value([type '::NBDRY']).Double;
catch
    type = '\EFIT02';
    try
        efit_nbdry = mds_server.Value([type '::NBDRY']).Double;
    catch
        type = '\EFIT01';
    end
end
efit_nbdry = mds_server.Value([type '::NBDRY']).Double;
efit_nmass = mds_server.Value([type '::NMASS']).Double;
efit_gtime = mds_server.Value([type '::GTIME']).Double;
efit_zbdry = mds_server.Value([type '::ZBDRY']).Double;
efit_rbdry = mds_server.Value([type '::RBDRY']).Double;
efit_psirz = mds_server.Value([type '::PSIRZ']).Double;
efit_rlim  = mds_server.Value([type '::RLIM']).Double;
efit_zlim  = mds_server.Value([type '::ZLIM']).Double;
efit_ipmhd  = mds_server.Value([type '::IPMHD']).Double;
efit_rgrid = mds_server.Value([type '::R']).Double;
efit_zgrid = mds_server.Value([type '::Z']).Double;
efit_ffp = mds_server.Value([type '::FFPRIM']).Double;
efit_pp = mds_server.Value([type '::PPRIME']).Double;
efit_pres = mds_server.Value([type '::PRES']).Double;
efit_fpol = mds_server.Value([type '::FPOL']).Double;
efit_qpsi = mds_server.Value([type '::QPSI']).Double;
efit_psin = mds_server.Value([type '::PSIN']).Double;
efit_psiaxis = mds_server.Value([type '::PSI0']).Double;
efit_psilim = mds_server.Value([type '::PSIBDY']).Double;
efit_psi = mds_server.Value([type '::PSI']).Double;
efit_rax = mds_server.Value([type '::R0']).Double;
efit_zax = mds_server.Value([type '::Z0']).Double;
efit_btorvac = mds_server.Value([type '::BT0VAC']).Double;
efit_betap = mds_server.Value([type '::BETAP']).Double;
efit_betat = mds_server.Value([type '::BETAT']).Double;
% Adjust the data shapes
nbdry = round(max(efit_nbdry));
ntime = length(efit_gtime) ;
nrho  = length(efit_ffp)/ntime;
efit_ffp = reshape(efit_ffp, [nrho ntime]);
efit_pres = reshape(efit_pres, [nrho ntime]);
efit_fpol = reshape(efit_fpol, [nrho ntime]);
efit_qpsi = reshape(efit_qpsi, [nrho ntime]);
efit_psin = reshape(efit_psin, [nrho ntime]);
efit_psi = reshape(efit_psi, [nrho ntime]);
%efit_psiaxis = reshape(efit_psiaxis, [nrho ntime]);
%efit_psilim = reshape(efit_pslim, [nrho ntime]);
efit_pp = reshape(efit_pp, [nrho ntime]);
efit_rgrid = reshape(efit_rgrid, [nrho ntime]);
efit_zgrid = reshape(efit_zgrid, [nrho ntime]);
efit_psirz = reshape(efit_psirz, [nrho nrho ntime]);
i1 = 1;
for i=1:length(efit_nbdry)
    i2=i1+nbdry-1;
    efit_r(1:nbdry,i) = efit_rbdry(i1:i2);
    efit_z(1:nbdry,i) = efit_zbdry(i1:i2);
    efit_rlim(1:nbdry,i) = efit_rlim(i1:i2);
    efit_zlim(1:nbdry,i) = efit_zlim(i1:i2);
    i1 = i2 + 1;
end
% Filter Pressure
dex = length(efit_pres(:,1));
for i=2:length(efit_pres(1,:))-1
    if (efit_pres(dex,i) < 0.0)
        efit_pres(:,i) = efit_pres(:,i-1);
    end
end
% Now get the values
%for j=1:size(efit_r,1)
%    rbdry(j)=pchip(efit_gtime,efit_r(j,:),time);
%    zbdry(j)=pchip(efit_gtime,efit_z(j,:),time);
%end
data.xbndry   =pchip(efit_gtime,efit_r,time);
data.zbndry   =pchip(efit_gtime,efit_z,time);
data.xlim    =pchip(efit_gtime,efit_rlim,time);
data.zlim    =pchip(efit_gtime,efit_zlim,time);
data.xaxis   =pchip(efit_gtime,efit_rax,time);
data.zaxis   =pchip(efit_gtime,efit_zax,time);
data.psiaxis =pchip(efit_gtime,efit_psiaxis,time);
data.psilim =pchip(efit_gtime,efit_psilim,time);
data.sp    =pchip(efit_gtime,efit_pres,time);
data.spp      =pchip(efit_gtime,efit_pp,time);
data.sffp      =pchip(efit_gtime,efit_ffp,time);
data.qpsi       =pchip(efit_gtime,efit_qpsi,time);
data.totcur     =pchip(efit_gtime,efit_ipmhd,time);
data.beta     =pchip(efit_gtime,efit_betat,time);
data.betap     =pchip(efit_gtime,efit_betap,time);
psi     =pchip(efit_gtime,efit_psi,time);
data.xgrid   =pchip(efit_gtime,efit_rgrid,time);
data.zgrid   =pchip(efit_gtime,efit_zgrid,time);
data.psixz   =pchip(efit_gtime,efit_psirz,time);
% Clean up the data
dex = data.xbndry >= 0.1;
data.xbndry = data.xbndry(dex);
data.zbndry = data.zbndry(dex);
dex = data.xlim ~= 0;
data.xlim = data.xlim(dex);
data.zlim = data.zlim(dex);
mds_server.CloseTree('NSTX',shotnum);
mds_server.disconnect(server);
% Plot the data
subplot(3,2,[1 3 5]);
contourf(data.xgrid,data.zgrid,data.psixz',100); axis equal; axis tight; 
hold on; 
plot(data.xbndry,data.zbndry,'r','LineWidth',2); 
plot(data.xaxis,data.zaxis,'+r','LineWidth',2); 
plot(data.xlim,data.zlim,'k','LineWidth',2); hold off;
xlabel('R'); ylabel('Z'); title(['\' type '  GEQDSK']);
mpol = 12;
ntor = 0;
R = data.xbndry;
Z = data.zbndry;
np = length(R);
rmaj = mean(R);
Rminor = R - rmaj;
theta = atan2(Z,Rminor);
theta(theta<0) = theta(theta<0) + 2*pi;
[theta, dex] =sort(theta);
R = R(dex);
Z = Z(dex);
Rminor = Rminor(dex);
rho = sqrt(Rminor.^2+Z.^2);
nu = np;
nv = 1;
refou=zeros(mpol+1,2*ntor+1);
zefou=zeros(mpol+1,2*ntor+1);
refou2=zeros(mpol+1,2*ntor+1);
zefou2=zeros(mpol+1,2*ntor+1);
cosu=zeros(nu,mpol+1);
cosv=zeros(nv,2*ntor+1);
sinu=zeros(nu,mpol+1);
sinv=zeros(nv,2*ntor+1);
alu=2*pi/nu;
for i=1:nu
    for j=1:mpol+1
        m=j-1;
        cosu(i,j)=cos(m*(i-1)*alu);
        sinu(i,j)=sin(m*(i-1)*alu);
    end
end
alv=2*pi/nv;
for i=1:nv
    for j=1:2*ntor+1
        n=j-ntor-1;
        cosv(i,j)=cos(n*(i-1)*alv);
        sinv(i,j)=sin(n*(i-1)*alv);
    end
end
% Our Amplitude is 1 everywhere since we do a full period
% Our Kernels is
%    read_VMEC arrays
%    cos(mu+nv) = cos(mu)cos(nv)-sin(mu)sin(nv)
%    sin(mu+nv) = sin(mu)cos(nv)+cos(mu)sin(nv)
fnuv=zeros(1,mpol+1);
fnuv(1)=1./(nu*nv);
for i=2:mpol+1;
    fnuv(i)=2*fnuv(1);
end
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                refou(m1,n1)=refou(m1,n1)+R(i)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou(m1,n1)=zefou(m1,n1)+Z(i)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                refou2(m1,n1)=refou2(m1,n1)+R(i)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou2(m1,n1)=zefou2(m1,n1)+Z(i)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end
data.refou = refou;
data.zefou = zefou;
data.refou2 = refou2;
data.zefou2 = zefou2;
% Create the profiles
pflux=data.psiaxis:(data.psilim-data.psiaxis)/(length(data.sp)-1):data.psilim;
pflux_norm=pflux-pflux(1);
pflux_norm=pflux_norm./max(pflux_norm);
jdotb=data.spp+data.sffp;
press=data.sp;
qprof=data.qpsi;
iotaf=1./qprof;
q_spl=pchip(pflux',qprof);
fun = @(x) ppval(q_spl,x);
for i=1:length(pflux)
    tflux(i) = integral(fun,pflux(1),pflux(i));
end
phiedge = tflux(length(pflux))*2*pi;
tflux = tflux./tflux(length(tflux));
s = 0:1/99:1;
am_aux_s = s;
am_aux_f = pchip(tflux,press,s);
p = polyfit(am_aux_s,am_aux_f,9);
p(length(p)) = am_aux_f(1);
p = [-sum(p(1:10))+am_aux_f(100) p];
am = fliplr(p);
subplot(3,2,2); plot(tflux,press,'k'); title('Pressure [Pa]'); hold on; plot(am_aux_s,am_aux_f,'or'); plot(am_aux_s,polyval(p,am_aux_s),'b'); legend('EQDSK','SPLINE','POLY'); hold off
ac_aux_s = s;
ac_aux_f = pchip(tflux,jdotb,s);
p = polyfit(ac_aux_s,ac_aux_f,9);
p(length(p)) = ac_aux_f(1);
p = [-sum(p(1:10))+ac_aux_f(100) p];
ac = fliplr(p);
subplot(3,2,4); plot(tflux,jdotb,'k'); title('<J.B>/<B/R> [A]'); hold on; plot(ac_aux_s,ac_aux_f,'or'); plot(ac_aux_s,polyval(p,ac_aux_s),'b'); hold off
ai_aux_s = s;
ai_aux_f = pchip(tflux,iotaf,s);
p = polyfit(ai_aux_s,ai_aux_f,9);
p(length(p)) = ai_aux_f(1);
p = [-sum(p(1:10))+ai_aux_f(100) p];
ai = fliplr(p);
subplot(3,2,6); plot(tflux,qprof,'k'); title('q'); hold on; plot(ai_aux_s,1./ai_aux_f,'or'); plot(ai_aux_s,1./polyval(p,ai_aux_s),'b'); hold off
data.phiedge = phiedge;
data.curtor = data.totcur;
data.raxis = data.xaxis;
data.zaxis = data.zaxis;
data.am = am;
data.am_aux_s = am_aux_s;
data.am_aux_f = am_aux_f;
data.ac = ac;
data.ac_aux_s = ac_aux_s;
data.ac_aux_f = ac_aux_f;
data.ai = ai;
data.ai_aux_s = ai_aux_s;
data.ai_aux_f = ai_aux_f;
% Check the data
mn = 1;
for m1 = 1: mpol + 1
    for n1=1: 2*ntor +1
        xm(mn) = m1-1;
        xn(mn) = n1-ntor-1;
        rmnc(mn,1) = refou(m1,n1);
        rmns(mn,1) = refou2(m1,n1);
        zmnc(mn,1) = zefou2(m1,n1);
        zmns(mn,1) = zefou(m1,n1);
        mn = mn + 1;
    end
end
ntheta=360;
theta=0:2*pi/(ntheta-1):2*pi;
zeta = 0;
r=cfunct(theta,zeta,rmnc,xm,xn)+sfunct(theta,zeta,rmns,xm,xn);
z=cfunct(theta,zeta,zmnc,xm,xn)+sfunct(theta,zeta,zmns,xm,xn);
subplot(3,2,[1 3 5]);
hold on
plot(r,z,'b');
hold off
return
end

