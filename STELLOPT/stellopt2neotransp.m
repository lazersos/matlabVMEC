function stellopt2neotransp(ext,varargin)
%STELLOPT2NEOTRASP Calculates radial electric field using NEOTRANSP
%   The STELLOPT2NEOTRANSP subroutine uses the NEOTRANSP library to
%   calculate the radial electric field from a STELLOPT reconstruction run.
%   The input extension of a run is passed to the routine and the 'tprof'
%   and 'wout' files are read.  The options are as follows:
%       list:       List all available equilID's
%       equilid:    Specify specific equilID (w7x-sc1 is default)
%
%   Example
%       stellopt2neotransp('test');
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00


equilID='w7x-sc1'; % Standard

% Handle varargin
if nargin > 1
    i=1;
    while i < nargin
        switch varargin{i}
            case{'list','LIST'}
                database=w7xdatabase();
                list(database)
                return;
            case{'equilid','equilID','EQUILID'}
                i=i+1;
                equilID=varargin{i};
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

% Get data
tprof=importdata(['tprof.' ext]);
tprof=tprof.data;
wout_nc_file = ['wout_' ext '.nc'];
wout_txt_file = ['wout.' ext];
wout_txt2_file = ['wout.' ext '.txt'];
if isfile(wout_nc_file)
    vmec_data=read_vmec(wout_nc_file);
elseif isfile(wout_txt_file)
    vmec_data=read_vmec(wout_txt_file);
elseif isfile(wout_txt2_file)
    vmec_data=read_vmec(wout_txt2_file);
else
    disp('ERROR: Cannot file wout file');
    return;
end
    
B00axis= vmec_data.b0;
phiedge = vmec_data.phi(end);

rho = sqrt(tprof(:,1));
ne  = tprof(:,2).*1E-20;
te  = tprof(:,3).*1E-3;
ti  = tprof(:,4).*1E-3;
Z   = tprof(:,5);
ni  = ne./Z;
h2  = 0.025;
rho2 = 0.1:h2:0.95;

% Electrons
Prof(1).Z=-1;
p_spl = spline(rho,te);
pp_spl = spline_deriv(p_spl,1);
Prof(1).TkeV=ppval(p_spl,rho2);
Prof(1).dTkeVdrho=ppval(pp_spl,rho2);
p_spl = spline(rho,ne);
pp_spl = spline_deriv(p_spl,1);
Prof(1).n20=ppval(p_spl,rho2);
Prof(1).dn20drho=ppval(pp_spl,rho2);

% Ions
Prof(2).Z=1;
Prof(2).m_au=1;
p_spl = spline(rho,ti);
pp_spl = spline_deriv(p_spl,1);
Prof(2).TkeV=ppval(p_spl,rho2);
Prof(2).dTkeVdrho=ppval(pp_spl,rho2);
p_spl = spline(rho,ni);
pp_spl = spline_deriv(p_spl,1);
Prof(2).n20=ppval(p_spl,rho2);
Prof(2).dn20drho=ppval(pp_spl,rho2);

Opt.makeErsearchplot=0; %set to 0 to turn off, or to the figure number
Opt.roots='i&e';          %can be 'i', 'e', 'i&e' 

Transp = neotransp(equilID,rho2,Prof,B00axis,Opt);

% Now we need to fit the data
factor = 1E-3;
dphi1 = Transp.dPhidskV(1,:).*1E3;
dphi2 = Transp.dPhidskV(2,:).*1E3;

% Handle single root
dphi1(isnan(dphi1)) = dphi2(isnan(dphi1));
dphi2(isnan(dphi2)) = dphi1(isnan(dphi2));

rho  = Transp.rho;
s    = rho.*rho;
fig = figure('Position',[ 1 1 1024 768],'Color','white');
subplot(1,2,1);
plot(s,factor.*dphi1./phiedge,'ob');
hold on;
plot(s,factor.*dphi2./phiedge,'or');
xlim([0 1]);
set(gca,'FontSize',24);
xlabel('Norm. Toridal Flux (s)');
ylabel('dV/ds [kV/Wb]');
title('NEOTRANSP E');
dphim=0.5*(dphi1+dphi2);
s = s(~isnan(dphim));
dphi1 = dphi1(~isnan(dphim));
dphi2 = dphi2(~isnan(dphim));
dphim = dphim(~isnan(dphim));
plot(s,factor.*dphim./phiedge,'ok');
pp=polyfit([0 s 1],[0 dphim 0], 6);
h = 0.01;
s2 = 0:h:1;
plot(s2,factor.*polyval(pp,s2)./phiedge,'g');
legend('Root1','Root2','Avg','Polyfit');
subplot(1,2,2);
plot(s,factor.*cumsum(dphi1),'b');
hold on;
plot(s,factor.*cumsum(dphi2),'r');
plot(s2,factor.*cumsum(polyval(pp,s2)),'ok');
xlim([0 1]);
set(gca,'FontSize',24);
xlabel('Norm. Toridal Flux (s)');
ylabel('e-static Potential [kV]');
title('NEOTRANSP Potential');
legend('Root1','Root2','Polyfit');


s3 = 0:0.05:1;
disp(['  POT_AUX_S = ' num2str(s3,' %20.10E')]);
disp(['  POT_AUX_F = ' num2str(cumsum(polyval(pp,s3)),' %20.10E')]);

end

