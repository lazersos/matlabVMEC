function [rho, Qi, Qe] = ascot5_calc_heat(data)
%ASCOT5_CALC_HEAT(data) Calculates heating profile from ASCOT5 data.
%   The ASCOT5_CALC_HEAT function calculates the ion and electron heating
%   using the BEAMS3D slowing down operator.  Units returned are in
%   [W/PHInorm] where PHInorm is normalized flux.  To get into units of
%   [W/m] the values should be multiplied by the VMEC VP array.  The value
%   rho is the normalized minor radius rho=sqrt(s).
%
% Example usage
%      data=read_hdf5('ascot5.h5');     % Reads VMEC wout file
%      [rho,Qi,Qe]=ascot5_calc_heat(data); % Reads VMEC mercier file
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00


rundex=1;

coul_log_type=2; %1: BEAMS; 2: ASCOT5

% Constants
me = 9.10938356D-31;
ec = 1.60217662E-19;
amu = 1.66053906660E-27;
e0 = 8.8542e-12;
hbar = 1.0546e-34;
Zplasma = 1;

% Extract run info
temp = fieldnames(data.plasma);
plasmaidstr = temp{rundex};
temp = fieldnames(data.results);
resultsidstr = temp{rundex};

% Extract Profile Stuffstuff
pte=data.plasma.(plasmaidstr).etemperature;
pne=data.plasma.(plasmaidstr).edensity;
pti=data.plasma.(plasmaidstr).itemperature;
pni=data.plasma.(plasmaidstr).idensity;
prho = data.plasma.(plasmaidstr).rho;
myZ = double(data.plasma.(plasmaidstr).znum);
mass = double(data.plasma.(plasmaidstr).mass).*amu;

% Extract distribution stuff
rho = data.results.(resultsidstr).distrho5d.abscissa_vec_01;
theta = data.results.(resultsidstr).distrho5d.abscissa_vec_02;
vpara = data.results.(resultsidstr).distrho5d.abscissa_vec_04;
vperp = data.results.(resultsidstr).distrho5d.abscissa_vec_05;
dist_rho = data.results.(resultsidstr).distrho5d.ordinate;
phi = data.results.(resultsidstr).distrho5d.abscissa_vec_03;
time = data.results.(resultsidstr).distrho5d.abscissa_vec_06;
ntheta = double(data.results.(resultsidstr).distrho5d.abscissa_nbin_02);
nphi  = double(data.results.(resultsidstr).distrho5d.abscissa_nbin_03);
ntime = double(data.results.(resultsidstr).distrho5d.abscissa_nbin_06);
dtheta = range(theta)./ntheta;
dphi  = range(phi)./nphi;
dtime = range(time)./ntime;
rho = (rho(1:end-1)+rho(2:end)).*0.5;
vpara = (vpara(1:end-1)+vpara(2:end)).*0.5;
vperp = (vperp(1:end-1)+vperp(2:end)).*0.5;
vtotal = sqrt(vpara'.^2.+vperp.^2);

% Create distribution in vll/vperp
fv=squeeze(sum(dist_rho,[1 2 5 6])).*dtheta.*dphi.*dtime;

% Spline to rho
te = pchip(prho,pte,rho);
ne = pchip(prho,pne,rho);
ti = pchip(prho,pti,rho);
ni = pchip(prho,pni,rho);
te(rho>1) = 0;
ne(rho>1) = 0;
te3 = te.^3;

% Calc Coulomb log
coulomb_log=[];
if coul_log_type == 1
    dex= te < 10.*myZ.*myZ;
    coulomb_log(dex) = 23 - log(myZ.*sqrt(ne(dex).*1E-6./te3(dex)));
    dex = ~dex;
    coulomb_log(dex) = 24 - log(sqrt(ne(dex).*1E-6)./te(dex));
    v_crit = ((0.75*sqrt(pi.*mass./me).*mass./mass).^(1./3.)).*sqrt(2.*te.*ec./mass);
    vcrit_cube = v_crit.^3;
    tau_spit = 3.777183E41.*mass.*sqrt(te3)./(ne.*myZ.*myZ.*coulomb_log');
    C1 = 1./tau_spit;
    C2 = vcrit_cube./tau_spit;
    C1(rho>1) = 0;
    C2(rho>1) = 0;
    f=squeeze(sum(dist_rho,[1 2 5 6])).*dtheta.*dphi.*dtime;
    Ei = zeros(size(f));
    Ee = zeros(size(f));
    for i=1:size(f,3)
        Ee(:,:,i) = C1(i).*f(:,:,i).*vtotal.*vtotal.*mass;
        Ei(:,:,i) = C2(i).*f(:,:,i).*vtotal./(vtotal.*vtotal).*mass;
    end
    Qe = squeeze(sum(Ee,[1 2]));
    Qi = squeeze(sum(Ei,[1 2]));
elseif coul_log_type ==2
    ma = amu;
    qa = myZ.*ec;
    mb_i = mass;
    mb_e = me;
    qb_i = Zplasma.*ec;
    qb_e = -ec;
    tej = te.*ec;
    tij = ti.*ec;
    % Calc Debye Length
    debye_i=ni.*qb_i.*qb_i./tij;
    debye_e=ne.*qb_e.*qb_e./tej;
    debye = sqrt(e0./(debye_i+debye_e));
    mr_i = ma.*mb_i/(ma+mb_i);
    mr_e = ma.*mb_e/(ma+mb_e);
    % We need to loop over surface
    for k=1:length(rho)
        vbar_i = vtotal.*vtotal + 2*tij(k)./mb_i;
        vbar_e = vtotal.*vtotal + 2*tej(k)./mb_e;
        bcl_i = abs(qa.*qb_i./(4.*pi.*mr_i.*vbar_i));
        bcl_e = abs(qa.*qb_e./(4.*pi.*mr_e.*vbar_e));
        bqm_i = abs(hbar./(2.*mr_i.*sqrt(vbar_i)));
        bqm_e = abs(hbar./(2.*mr_e.*sqrt(vbar_e)));
        itemp = bcl_i.*0.0;
        etemp = bcl_e.*0.0;
        dex = bcl_i>bqm_i;
        itemp(dex) = log(debye(k) ./ bcl_i(dex));
        itemp(~dex) = log(debye(k) ./ bqm_i(~dex));
        dex = bcl_e>bqm_e;
        etemp(dex) = log(debye(k) ./ bcl_e(dex));
        etemp(~dex) = log(debye(k) ./ bqm_e(~dex));
        %clog_i(k,:,:) = itemp;
        %clog_e(k,:,:) = etemp;
        tau_i  = 3.777183E41.*mass.*sqrt(te3(k))./(ne(k).*myZ.*myZ.*itemp);
        tau_e  = 3.777183E41.*mass.*sqrt(te3(k))./(ne(k).*myZ.*myZ.*etemp);
        v_crit = ((0.75*sqrt(pi.*mass./me).*mass./mass).^(1./3.)).*sqrt(2.*te(k).*ec./mass);
        C1_e   = 1./tau_e;
        C1_i   = 1./tau_i;
        C2_e   = v_crit.*v_crit.*v_crit./tau_e;
        C2_i   = v_crit.*v_crit.*v_crit./tau_i;
        ftemp  = squeeze(fv(:,:,k));
        Qe(k)         = squeeze(sum(C1_e.*ftemp.*vtotal.*vtotal.*mass,[1 2]));
        Qe(k)         = Qe(k) + squeeze(sum(C1_i.*ftemp.*vtotal.*vtotal.*mass,[1 2]));
        Qi(k)         = squeeze(sum(C2_e.*ftemp.*vtotal.*mass./(vtotal.*vtotal),[1 2]));
        Qi(k)         = Qi(k) + squeeze(sum(C2_i.*ftemp.*vtotal.*mass./(vtotal.*vtotal),[1 2]));
    end
    
end

% Calc Coefficients
end

