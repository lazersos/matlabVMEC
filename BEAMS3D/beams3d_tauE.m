function [tauE, nttau] = beams3d_tauE(beam_data,varargin)
%BEAMS3D_TAUE Calculates energy confinement and tripple product
%   The BEAMS3D_TAUE function calculates the thermal confinement time and
%   tripple product for a BEAMS3D run.  It takes a BEAMS3D data structure
%   as read by READ_BEAMS3D as input.  It returns the thermal confinement
%   time in seconds and the tripple product in units of m^-3*eV*s.
%
% Example usage
%      beam_data=read_beams3d('beams3d_test.h5');  % Reads BEAMS3D HDF5 file
%      [tauE, nttau]=beams3d_tauE(beam_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Defaults
ec = 1.60217662E-19;
tau_e=[];
nttau=[];
beamdex=[];

% Handle varargin
if ~isempty(varargin)
    i=1;
    while i<=numel(varargin)
        switch varargin{i}
            case 'beamdex'
                i=i+1;
                beamdex=varargin{i};
        end
        i=i+1;
    end
end

% Handle beamdex
if isempty(beamdex)
    beamdex = 1:double(beam_data.nbeams);
end

% Extract P
p = ec.*beam_data.NE.*(beam_data.TE+beam_data.TI./beam_data.ZEFF_ARR);

% Grid helpers
dr = range(beam_data.raxis)./double(beam_data.nr);
dz = range(beam_data.zaxis)./double(beam_data.nz);
dp = range(beam_data.phiaxis)./double(beam_data.nphi);
nfp = round(2*pi/max(beam_data.phiaxis));

% Make volume
dV=repmat(beam_data.raxis.*dr.*dp.*dz,[1 beam_data.nphi beam_data.nz]);

% Calc Ekin
E=sum(p.*dV,'all').*nfp;

% Call table for total power
table = beams3d_powertable(beam_data,'quiet');
totals = sum(table(:,beamdex),2);
pnbi = sum(totals([4 5 8]));

% Calc tauE
tauE=E./pnbi;

temp = [];
for i=1:beam_data.nphi
    s=beam_data.S_ARR(:,i,:);
    n=beam_data.NE(:,i,:)./beam_data.ZEFF_ARR(:,i,:);
    t=beam_data.TI(:,i,:);
    %smin = min(min(s));
    %j=find(s==smin);
    j=s<0.01; % rho 10%
    temp(i) = mean(n(j).*t(j));
end
nttau = mean(temp).*tauE;
return
end

