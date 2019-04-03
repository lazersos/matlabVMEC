function fieldlines_fit_surf( data, isurf )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lplots=1;
mpol=10;
ntheta=2.^nextpow2(mpol*8);
theta = 0:2*pi./(ntheta-1):2*pi;
rpts = data.R_lines(isurf,:);
zpts = data.Z_lines(isurf,:);
npoinc = data.npoinc;
nsteps = data.nsteps;

%Filter
zpts = zpts(rpts>0);
rpts = rpts(rpts>0);

% Setup Arrays
rf = zeros(ntheta,npoinc);
zf = zeros(ntheta,npoinc);

for n=1:npoinc
    rtemp = rpts(n:npoinc:end);
    ztemp = zpts(n:npoinc:end);
    dr    = rtemp-rtemp(1);
    dz    = ztemp-ztemp(1);
    dl    = dr.*dr+dz.*dz;
    dl(1) = max(dl);
    [l(1),dex] = min(dl);
    r2    = rtemp(2);
    z2    = ztemp(2);
    rtemp(2) = rtemp(dex);
    ztemp(2) = ztemp(dex);
    rtemp(dex) = r2;
    ztemp(dex) = z2;
    for j=3:length(rtemp)
        dr    = rtemp-rtemp(j-1);
        dz    = ztemp-ztemp(j-1);
        dl    = dr.*dr+dz.*dz;
        dl(1:j-1) = max(dl);
        [l(j-1),dex] = min(dl);
        r2    = rtemp(j);
        z2    = ztemp(j);
        rtemp(j) = rtemp(dex);
        ztemp(j) = ztemp(dex);
        rtemp(dex) = r2;
        ztemp(dex) = z2;
    end
    % Define theta-like vector
    l = cumsum(dl);
    l = 2.*pi.*l./max(l);
    % Put on equidistant grid
    r2 = pchip(l,rtemp,theta);
    z2 = pchip(l,ztemp,theta);
    % Rescale
    r2 = r2 - mean(r2);
    z2 = z2 - mean(z2);
    %[r_scale,dex] = max(r2);
    %r2 = r2-r_scale;
    %r2 = [r2(dex:end) r2(1:dex-1)];
    %
    rf(:,n) = r2;
    zf(:,n) = z2;
    % Plot 1
    if lplots
        subplot(1,2,1);
        plot(r2,z2,'ok'); hold on;
    end
    % Fit R/Z
    x0 = [zeros(1,mpol-1) zeros(1,mpol-1) 0];
    x0(1) = mean(r2);
    x = lsqcurvefit(@(x,x0) f_ctemp(x,x0,mpol),x0,[theta theta],[r2 z2]);
    rmn = x(1:mpol);
    zmn(2:mpol) = x(mpol+1:2*mpol-1);
    phase = x(end);
    theta2 = theta - phase;
    % FFT R
%     temp = fftshift(fft(r2));
%     j = length(temp)/2 + 1;
%     dex1=j:j+mpol;
%     dex2=j-1:-1:j-mpol;
%     rmn  = real(temp(dex1));
%     rmn(2:end)  = rmn(2:end)+real(temp(dex2));
%     rmn = rmn ./ length(temp);
%     % FFT Z
%     temp = fftshift(fft(z2));
%     j = length(temp)/2 + 1;
%     dex1=j:j+mpol;
%     dex2=j-1:-1:j-mpol;
%     zmn  = imag(temp(dex1));
%     zmn(2:end)  = zmn(2:end)-imag(temp(dex2));
%     zmn = zmn ./ length(temp);
%     theta2 = theta;
    % Back Transformation
    xm   = 0:mpol-1;
    xn   = xm.*0;
    ro   = cfunct(theta2,0,rmn',xm,xn);
    zo   = sfunct(theta2,0,zmn',xm,xn);
    if lplots
        plot(ro,zo,'r'); hold off; drawnow;
        subplot(1,2,2);
        plot(theta2,r2,'ok'); hold on; plot(theta2,ro,'b');
        plot(theta2,z2,'+k');plot(theta2,zo,'g');
        hold off;
        drawnow; pause(1.);
    end
end



return
end


function [F,J]=f_ctemp(x,xdata,mpol)
phase = x(end);
nout = length(xdata)/2;
%phase = mod(x(end),2*pi);
F =zeros(1, length(xdata));
m = 0;
i = 1;
for k=1:mpol
    F(1:nout) = F(1:nout)+x(i).*cos((xdata(1:nout)-phase).*m);
    m = m+1;
    i = i+1;
end
m = 1;
for k=1:mpol-1
    F(nout+1:2*nout) = F(nout+1:2*nout)+x(i).*sin((xdata(nout+1:2*nout)-phase).*m);
    m = m+1;
    i = i+1;
end
return
end

