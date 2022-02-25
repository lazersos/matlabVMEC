function dist_out = beams3d_normdist(beam_data,vmec_data)
%BEAMS3D_NORMDIST Fixed normalization of distribution to s^3/m^6
%   The BEAMS3D_NORMDIST routine correct the missing 1/dV factor present in
%   the distribution function of VMEC based BEAMS3D runs.  It requires the
%   user to provide both a VMEC data structure as returned by READ_VMEC and
%   a BEAMS3D data strucutre as returned by READ_BEAMS3D.  The function
%   returns the correctly normalized distribution function in units of
%   s^3/m^6.
%
%   Example Usage
%       beam_data=read_beams3d('beams3d_test.h5');
%       vmec_data=read_vmec('wout_test.nc');
%       dist = beams3d_normdist(beam_data,vmec_data);
%
%   Maintained by:  Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0


dist_out = [];

% Set dist out the whatever is in the BEAMS3D file
dist_out=beam_data.dist_prof;

% Probably should check version in future
if beam_data.VERSION > 999
    disp(['Distribution already in correct units.'])
    dist_out=beam_data.dist_prof;
    return;
end

%Figure out grid
edges = 0:1./beam_data.ns_prof1:1;
rhob = 0.5.*(edges(1:end-1)+edges(2:end));
sb = rhob.*rhob;
edges = 0:1./beam_data.ns_prof2:1;
ub = pi.*(edges(1:end-1)+edges(2:end));
edges = 0:1./beam_data.ns_prof3:1;
phib = pi.*(edges(1:end-1)+edges(2:end));
dsv = 1./(vmec_data.ns-1);
sv = 0:dsv:1;
rhov = sqrt(sv);

% Create Interpolating functions (odd modes in rho)
i1 = floor(sb./dsv)+1;
i1 = max(min(i1,vmec_data.ns-1),2);
i2 = i1+1;
gmnc = zeros(vmec_data.mnmax_nyq,beam_data.ns_prof1);
x = (sb-(i1-1)*dsv)./dsv;
for mn = 1:vmec_data.mnmax_nyq
    f1 = vmec_data.gmnc(mn,i1);
    f2 = vmec_data.gmnc(mn,i2);
    if ~(mod(vmec_data.xm_nyq(mn),2)==0) % odd
        f1 = f1.*rhob./rhov(i1);
        f2 = f2.*rhob./rhov(i2);
    end
    %plot([f1;f2]'); pause(1);
    gmnc(mn,:) = f1.*(1.0-x)+f2.*x;
end
if vmec_data.iasym
    gmns = zeros(vmec_data.mnmax_nyq,beam_data.ns_prof1);
    for mn = 1:vmec_data.mnmax_nyq
        x = (sb-(i1-1)*dsv)./dsv;
        f1 = vmec_data.gmns(mn,i1);
        f2 = vmec_data.gmns(mn,i2);
        if ~(mod(vmec_data.xm_nyq(mn),2)==0) % odd
            f1 = f1.*rhob./rhov(i1);
            f2 = f2.*rhob./rhov(i1);
        end
        gmns(mn,:) = f1.*(1-x)+f2.*x;
    end
end

% Calculate Jacobian
sqrtg = cfunct(ub,phib,gmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
if vmec_data.iasym
    sqrtg = sqrtg + sfunct(ub,phib,gmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
end

% Adjust from dV/ds to dV/drho dV/drho = dV/ds * ds/drho = 2*rho*dVds
sqrtg =2.*sqrtg.*repmat(rhob',[1 beam_data.ns_prof2 beam_data.ns_prof3]);

% Adjust shape
sqrtg = repmat(sqrtg,[1 1 1 beam_data.ns_prof4 beam_data.ns_prof5 beam_data.nbeams]);
sqrtg = permute(sqrtg,[6 1 2 3 4 5]);

% Normalize
dist_out = dist_out./sqrtg;

end