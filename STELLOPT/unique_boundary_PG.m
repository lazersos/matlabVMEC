function [rmnc, zmns, xm, xn] = unique_boundary_PG(deltamn,xm_d,xn_d)
%UNIQUE_BOUNDARY_PG Converts the Garabedian coordinates to VMEC.
% This funciton reads takes arrays DELTAMN with size MNMAX,NS and
% array XM_D and XN_D size MNMAX and converts them into the VMEC
% represenation.  It returns a RMNC, ZMNS, XM, XN.
%
% Example usage
%      vmec_data=read_vmec('wout.test');     % Reads VMEC wout file
%      [deltamn, xm_d, xn_d] =convert_boundary_PG(vmec_data.rmnc,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
%      [rmnc,zmns,xm,xn] = unique_boundary_PG(deltamn,xm_d,xn_d);
%
% See also convert_boundary_PG and read_vmec.
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% Helpers
s=size(deltamn);
mnmax = s(1);
ns   = s(2);
nmin = min(xn_d);
mmin = min(xm_d);
nmax = max(xn_d);
mmax = max(xm_d);
nm00 = and(xn_d==0,xm_d==0);
nm01 = and(xn_d==0,xm_d==1);

% Extract normalization
rnorm = deltamn(nm00,:)./deltamn(nm01,:);
r00   = deltamn(nm00,:);
deltamn(nm00,:) = 1;

% Create rmnc zmns
n1 = nmax-nmin+1;
n2 = max(abs([mmin mmax])+2);
mnmax = n1*n2;
rmnc = zeros(mnmax,ns);
zmns = zeros(mnmax,ns);
mn = 1;
for m = 0:mmax+1
    for n=nmin:nmax
        xm(mn) = m;
        xn(mn) = n;
        mn = mn + 1;
    end
end


% Negative m
for m=mmin:0
    for n=nmin:nmax
        mn1 = and(xn==-n,xm==(abs(m)+1));
        mn2 = and(xn_d==n,xm_d==m);
        rmnc(mn1,:) = rmnc(mn1,:) + deltamn(mn2,:);
        zmns(mn1,:) = zmns(mn1,:) + deltamn(mn2,:);
    end
end
% Positive m
for m=1:mmax
    for n=nmin:nmax
        mn1 = and(xn==n,xm==(m-1));
        mn2 = and(xn_d==n,xm_d==m);
        rmnc(mn1,:) = rmnc(mn1,:) + deltamn(mn2,:);
        zmns(mn1,:) = zmns(mn1,:) - deltamn(mn2,:);
    end
end
% m=0
for n=1:nmax
    mn1 = and(xn==n,xm==0);
    mn2 = and(xn==-n,xm==0);
    rmnc(mn1,:) = rmnc(mn1,:) + rmnc(mn2,:);
    rmnc(mn2,:) = 0;
    zmns(mn1,:) = zmns(mn1,:) - zmns(mn2,:);
    zmns(mn2,:) = 0;
end
mn1 = and(xn==0,xm==0);
zmns(mn1,:) = 0;

% Normalize
for i=1:mnmax
    rmnc(i,:) = rmnc(i,:).*rnorm(1,:);
    zmns(i,:) = zmns(i,:).*rnorm(1,:);
end
mn1 = and(xn==0,xm==0);
rmnc(mn1,:) = r00;


end

