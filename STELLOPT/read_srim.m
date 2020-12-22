function data = read_srim(filename)
%READ_SRIM Reads the output of the SRIM package
%   The READ_SRIM funciton reads the output of the Standard Ranging Into
%   Matter package.  Current it supports reading the RANGE_3D.txt file.
%
% Example usage
%      data=read_srim('RANGE_3D.txt');     % Reads VMEC wout file
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.00

data=[];
if contains(filename,'RANGE_3D.txt')
    fid = fopen(filename,'r');
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    line = fgetl(fid);
    dex = strfind(line,'=')+1;
    data.ion = strtrim(line(dex(1):dex(1)+10));
    data.ion_mass = sscanf(line(dex(2):end),'%f');
    line = fgetl(fid);
    dex = strfind(line,'=')+1;
    data.ion_energy = sscanf(line(dex(1):end),'%e').*1E3;
    line = fgetl(fid);
    dex = strfind(line,'=')+1;
    data.ion_angle = sscanf(line(dex(1):end),'%f').*1E3;
    fgetl(fid);
    line = fgetl(fid);
    lcount = 0;
    while ~contains(line,'-------')
        lfirst = 1;
        while contains(line,'Layer')
            if lfirst
                lfirst = 0;
                lcount = lcount + 1;
                data.layer{lcount}.name = line(14:end);
                data.num_layers = lcount;
            elseif contains(line,'Depth')
                dex = strfind(line,'=')+1;
                data.layer{lcount}.depth = sscanf(line(dex:end),'%e').*1E-10; % Ang to m
            elseif contains(line,'Density')
                dex = strfind(line,'=')+1;
                data.layer{lcount}.density = sscanf(line(dex:end),'%e').*1000; % g/cm^3 to kg/m^3
            end
            line = fgetl(fid);
        end
        line = fgetl(fid);
    end
    data.range3d=fscanf(fid,'%d %e %e %e',[4 inf]);
    data.range3d(2:4,:) = data.range3d(2:4,:).*1E-10;
    fclose(fid);
    data.units='kg,m,eV';
else
    disp(['  Unknown filetype: ' filename]);
end
return
end

