function [ p_spline j_spline beta betai lp lj k iote] = pies_profiles(varargin)
%pies_profiles Converts VMEC presf and jcurv arrays into splines over phi
%   Detailed explanation goes here

offset=0;
i=1;
vmec_data = [];
use_poly = 0;
mu0=4*pi*1.e-7;

% Handle varargin
while i <= nargin-offset
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'wout'
                    vmec_data=varargin{i};
            end
        end
    elseif ischar(varargin{i})
        switch varargin{i}
            case 'poly'
                use_poly=1;
        end
    end
    i=i+1;
end


% Get the profiles
phiedge = vmec_data.phi(vmec_data.ns);
scoord = vmec_data.phi./phiedge;
pressf = vmec_data.presf;
iprime = vmec_data.jcurv./max(abs(vmec_data.jcurv));
beta = mu0;
betai = vmec_data.Itor*mu0/(trapz(scoord,iprime)*phiedge);
iote  = vmec_data.Itor*mu0/(2*pi);
p_spline = pchip(scoord,pressf);
j_spline = pchip(scoord,iprime);

% Flip the coefs for PIES
p_spline.coefs = flipdim(p_spline.coefs,2);
j_spline.coefs = flipdim(j_spline.coefs,2);
lp=int32(p_spline.pieces)+1;          % Number of breakpoints for p
lj=int32(j_spline.pieces)+1;          % Number of breakpoints for j
k=int32(p_spline.order);          % Order of spline (must be 4)


return



end

