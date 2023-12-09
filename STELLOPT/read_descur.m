function data = read_descur(filename,varargin)
%READ_DESCUR Read DESCUR files
%   The READ_DESCUR subroutine reads DESCUR output into a structure. This
%   structure can be used for initalizing a VMEC indata namelist structure.
%   The rouine can also optionally plot the DESCUR data.
%   Currently only supports reading plotout files.
%
%   Example usage
%
%       data = read_descur('plotout'); % just read
%       data = read_descur('plotout','plot'); % make a plot too
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

lplot=0;
data=[];

% Handle inputs
if (nargin > 1)
    j=1;
    while j<=numel(varargin)
        if ischar(varargin{j})
            switch(varargin{j})
                case {'plot'}
                    lplot = 1;
            end
        end
        j=j+1;
    end
end


fid=fopen(filename,'r');
if contains(filename,'plotout')
    head = fgetl(fid);
    temp=sscanf(head,'%i',7);
    data.mpol = temp(1);
    data.ntheta = temp(2);
    data.nphi = temp(3);
    data.mpol1 = temp(4);
    data.nphi2 = temp(5);
    data.nfp = temp(6);
    data.mpnt = temp(7);
    temp = fscanf(fid,'%e %e',[2 data.nphi.*data.ntheta]);
    data.rin = reshape(temp(1,:),[data.ntheta,data.nphi]);
    data.zin = reshape(temp(2,:),[data.ntheta,data.nphi]);
    temp = fscanf(fid,'%e %e %e %e',[4 data.mpol.*(2.*data.nphi2+1)]);
    data.rmnc = temp(1,:)';
    data.zmns = temp(2,:)';
    data.rmns = temp(3,:)';
    data.zmnc = temp(4,:)';
    i=1;
    for m = 0:data.mpol1
        for n = -data.nphi2:data.nphi2
            data.xm(i) = m;
            data.xn(i) = -n.*data.nfp;
            i = i + 1;
        end
    end
elseif contains(filename,'outcurve')
end

fclose(fid);

if lplot
    figure('Position',[1 1 1024*2 768],'Color','white','InvertHardCopy','off');
    zeta = [0 pi/2 pi]./data.nfp;
    dex1  = round([1 data.nphi/4 data.nphi/2])-1;
    theta = deg2rad(0:360);
    for i = 1:3
        subplot(1,3,i);
        dex = (1:data.ntheta) + dex1(i).*data.ntheta;
        plot(data.rin(dex),data.zin(dex),'.k');
        hold on;
        r = cfunct(theta,zeta(i),data.rmnc,data.xm,data.xn);
        r = r + sfunct(theta,zeta(i),data.rmns,data.xm,data.xn);
        z = sfunct(theta,zeta(i),data.zmns,data.xm,data.xn);
        z = z + cfunct(theta,zeta(i),data.zmnc,data.xm,data.xn);
        plot(squeeze(r),squeeze(z),'r');
    end
end
end