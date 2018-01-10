function stellopt_space( ext )
%STELLOPT_SPACE(EXT) Plots the variation of equilibira
%   Detailed explanation goes here

xvec=read_stellopt('xvec.dat');
vmec_data=read_vmec(['wout_' ext '.00000.nc']);
[deltamn0,xm_d0,xn_d0]=convert_boundary_PG(vmec_data.rmnc,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
deltamn0=deltamn0(:,end);
theta = 0:2*pi./127:2*pi;

rmnc0 = vmec_data.rmnc;
zmns0 = vmec_data.zmns;
xm0   = vmec_data.xm;
xn0   = vmec_data.xn;
ns0   = vmec_data.ns;

% Sort the var names
fid = fopen('var_labels','r');
n=fscanf(fid,'%d',1);
fgetl(fid);
mn = 1;
dex = zeros(n,1);
for i=1:n
    line=strtrim(fgetl(fid));
    switch(line(1:5))
        case 'DELTA'
            dex(i) = 1;
            xn_d(mn)=sscanf(line(7:10),'%d',1);
            xm_d(mn)=sscanf(line(12:15),'%d',1);
            mn = mn + 1;
    end
end

% Extract just the Fourier harmonics
x = xvec.x(dex==1,:);

delta_save = deltamn0;
for i = 1:size(x,2)
    %if xvec.chisq(i) >= 1E12; continue; end
    deltamn0 = delta_save;
    deltamn = x(:,i);
    for mn = 1:numel(xn_d)
        m = xm_d(mn);
        n = xn_d(mn);
        mn1 = and(xn_d0==n,xm_d0==m);
        deltamn0(mn1) = deltamn(mn);
    end
    [rmnc,zmns,xm,xn] = unique_boundary_PG(deltamn0,xm_d0,xn_d0);
    r = cfunct(theta,0,rmnc,xm,xn);
    z = sfunct(theta,0,zmns,xm,xn);
    hold on;
    linestyle = 'k';
    linewidth = 0.5;
    if (i==1), linestyle='r'; linewidth=2.0; end
    plot(r,z,linestyle,'LineWidth',linewidth);
    hold off;
end

end

