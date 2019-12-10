function data= read_bbnbi_injector(varargin)
%READ_BBNBI_INJECTOR Reads the BBNBI 'injector' file.
%The READ_BBNBI_INJECTOR function reads a BBNBI 'injector' file.  Called
%with no options it attempts to read the 'injector' file from the current
%directory.
%
% Example usage
%      bbnbi_data = read_bbnbi_injector; % Read from directory
%      bbnbi_data = read_bbnbi_injector('test/injector'); % Read from path
%      bbnbi_data = read_bbnbi_injector('plot'); % Make a plot
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.00
data=[];

lplot = 0;
filename='./injector';
if nargin > 0
    for i=1:nargin
        switch varargin{i}
            case {'plot'}
                lplot=1;
            otherwise
                if ~isempty(strfind(varargin{i},'/injector'))
                    filename=varargin{i};
                end
        end
    end
end

if ~isfile(filename)
    disp([' No injector file found at:' filename]);
    return;
end

fid = fopen(filename,'r');
line = fgetl(fid);
line = fgetl(fid);
data.ninj_id=sscanf(line,'%i %*s',1);
line = fgetl(fid);
data.npin_id=sscanf(line,'%i %*s',1);
for j = 1:data.npin_id
    line = fgetl(fid);
    line = fgetl(fid);
    data.weight_id(j)=sscanf(line,'%i %*s',1);
    line = fgetl(fid);
    data.energy_id(j)=sscanf(line,'%i %*s',1);
    line = fgetl(fid);
    temp=sscanf(line,'%d %*s',2);
    data.horz_mis(j) = temp(1);
    data.vert_mis(j) = temp(2);
    line = fgetl(fid);
    data.nbeamlet(j)=sscanf(line,'%i %*s',1);
    for i = 1: data.nbeamlet(j)
        line = fgetl(fid);
        line = fgetl(fid);
        temp=sscanf(line,'%i %i %*s',2);
        data.disp_id(j,i) = temp(1);
        data.anum_id(j,i) = temp(2);
        line = strrep(fgetl(fid),'D','e');
        temp=sscanf(line,'%g %g %g',3);
        data.x(j,i) = temp(1);
        data.y(j,i) = temp(2);
        data.z(j,i) = temp(3);
        line = strrep(fgetl(fid),'D','e');
        temp=sscanf(line,'%g %g %g',3);
        data.theta(j,i) = temp(1);
        data.phi(j,i) = temp(2);
        data.length(j,i) = temp(3);
    end
end

data.nx = data.length.*sin(data.theta).*cos(data.phi);
data.ny = data.length.*sin(data.theta).*sin(data.phi);
data.nz = data.length.*cos(data.theta);

fclose(fid);

if (lplot)
    plot3(data.x',data.y',data.z','s');
    hold on;
    axis equal;
    for j = 1:data.npin_id
        mx = mean(data.x(j,:),2);
        my = mean(data.y(j,:),2);
        mz = mean(data.z(j,:),2);
        dx = data.x(j,:)-mx;
        dy = data.y(j,:)-my;
        dz = data.z(j,:)-mz;
        d  = sqrt(dx.*dx+dy.*dy+dz.*dz);
        [~,i]=min(d);
        quiver3(mx,my,mz,data.nx(j,i),data.ny(j,i),data.nz(j,i));
    end
end

return;
end

