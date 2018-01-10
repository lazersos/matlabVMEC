function [deltamn,xm_d,xn_d] = convert_boundary_PG(rmnc,zmns,xm,xn)
%CONVERT_BOUNDARY_PG(rmnc,zmns,xm,xn) Convert to Garabedian coordiantes.
% This funciton reads takes arrays rmnc zmns with size MNMAX,NS and
% array xm and xn size MNMAX and converts them into the Garabedian
% represenation.  It returns a DELTAMN, XM_D, and XN_D arrays.
%
% Example usage
%      vmec_data=read_vmec('wout.test');     % Reads VMEC wout file
%      [deltamn, xm_d, xn_d] =convert_boundary_PG(vmec_data.rmnc,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
%
% See also unique_boundary_PG and read_vmec.
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% Helpers
nfp  = min(xn(xn>0));
nmax = max(xn)/nfp;
mmax = max(xm);
mmax2 = mmax+1;
ns = size(rmnc,2);
xn2=xn./nfp;
nm00=and(xn2==0,xm==0);
nm01=and(xn2==0,xm==1);
r00 = rmnc(nm00,:);
rnorm = 2./(rmnc(nm01,:)+zmns(nm01,:));

for nm = 1:length(xm)
    rmnc2(nm,:) = rmnc(nm,:).*rnorm;
    zmns2(nm,:) = zmns(nm,:).*rnorm;
end

n1 = 2*nmax+1;
n2 = 2*mmax2+1;
deltamn=zeros(n1*n2,ns);
xm_d=zeros(n1*n2,1);
xn_d=zeros(n1*n2,1);
nm =1;
for m=-mmax2:mmax2
    for n= -nmax:nmax
        xm_d(nm) = m;
        xn_d(nm) = n;
        nm = nm + 1;
    end
end
for m = 0:mmax
    for n=-nmax:nmax
        nm1 = and(xn2==n,xm==m);
        if isempty(find(nm1)), continue; end;
        nm2 = and(xn_d==n,xm_d==m+1);
        nm3 = and(xn_d==-n,xm_d==-m+1);
        try
        deltamn(nm2,:) = 0.5*(rmnc2(nm1,:)-zmns2(nm1,:)) + deltamn(nm2,:);
        deltamn(nm3,:) = 0.5*(rmnc2(nm1,:)+zmns2(nm1,:)) + deltamn(nm3,:);
        catch
            return;
        end
    end
end

nm00=and(xn_d==0,xm_d==0);
deltamn(nm00,:) = r00; % Stores R00 in DELTA(0,0)
return

end

