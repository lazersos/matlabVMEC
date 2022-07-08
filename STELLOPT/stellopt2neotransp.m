function stellopt2neotransp(ext,varargin)
%STELLOPT2NEOTRASP Calculates radial electric field using NEOTRANSP
%   The STELLOPT2NEOTRANSP subroutine uses the NEOTRANSP library to
%   calculate the radial electric field from a STELLOPT reconstruction run.
%   The input extension of a run is passed to the routine and the 'tprof'
%   and 'wout' files are read.  The options are as follows:
%       list:       List all available equilID's
%       equilid:    Specify specific equilID (w7x-sc1 is default)
%       b0:         Specify field on axis (default is to use WOUT value)
%       D:          Deutrium main ion mass (H default)
%       T:          Tritium main ion mass (H default)
%       He:         Helium main ion mass (H default)
%
%   Example
%       stellopt2neotransp('test');
%       stellopt2neotransp([],'list'); % List equil database
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.00


equilID='w7x-sc1'; % Standard
B00axis=[];
mion=1;

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
            case{'b0'}
                i=i+1;
                B00axis=varargin{i};
            case{'D'}
                mion=2;
            case{'T'}
                mion=3;
            case{'He'}
                mion=4;
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

if isempty(B00axis)
    B00axis= vmec_data.b0;
end
phiedge = vmec_data.phi(end); %Only used for plotting

rho = sqrt(tprof(:,1));
ne  = tprof(:,2).*1E-20;
te  = tprof(:,3).*1E-3;
ti  = tprof(:,4).*1E-3;
Z   = tprof(:,5);
ni  = ne./Z;
h2  = 0.01;
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
Prof(2).m_au=mion;
p_spl = spline(rho,ti);
pp_spl = spline_deriv(p_spl,1);
Prof(2).TkeV=ppval(p_spl,rho2);
Prof(2).dTkeVdrho=ppval(pp_spl,rho2);
p_spl = spline(rho,ni);
pp_spl = spline_deriv(p_spl,1);
Prof(2).n20=ppval(p_spl,rho2);
Prof(2).dn20drho=ppval(pp_spl,rho2);

Opt.makeErsearchplot=0; %set to 0 to turn off, or to the figure number
Opt.roots='i/e';          %can be 'i', 'e', 'i&e' 

Transp = neotransp(equilID,rho2,Prof,B00axis,Opt);

rho  = Transp.rho;
s    = Transp.s;
% Now we need to fit the data
dphi = Transp.dPhidskV(1,:).*1E3; %kV to V
dphi(isnan(dphi)) = 0;
sh = 0.5.*(s(2:end)+s(1:end-1)); %half grid
f =(dphi(2:end)+dphi(1:end-1)).*0.5; %half grid
phi = cumsum(f.*diff(s))+dphi(1);
pp = pchip(sqrt(sh),phi);
h3 = 0.05;
s3 = 0:h3:1;
rho3=sqrt(s3);
% Make a plot
plottransp(Transp); set(gcf,'Position',[1 1 1024 768],'Color','white');
drawnow;
subplot(2,3,3); hold on; f=ppval(pp,rho3); x=0.5.*(s3(2:end)+s3(1:end-1)); dphids=-diff(f)./diff(s3); Er=dphids.*sqrt(x).*2./Transp.minorradiusW7AS./1E3; plot(sqrt(x),Er,'or');
saveas(gcf,['transp_' ext '.fig']);
saveas(gcf,['transp_' ext '.png']);
close(gcf);

h3 = 0.05;
s3 = 0:h3:1;
disp(['  POT_AUX_S = ' num2str(s3,' %20.10E')]);
%disp(['  POT_AUX_F = ' num2str(cumsum(polyval(pp,s3)),' %20.10E')]);
disp(['  POT_AUX_F = ' num2str(cumsum(ppval(pp,s3)).*h3,' %20.10E')]);

end

