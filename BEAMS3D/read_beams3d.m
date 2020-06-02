function data = read_beams3d(filename)
%READ_BEAMS3D Reads the HDF5 file created by BEAMS3D
% This funciton reads the fieldlines file and returns the data from the
% file in a structure.
%
% Example usage
%      data=read_beams3d('fieldline_test.h5');  % Reads BEAMS3D HDF5 file
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

if (strcmp(filename(end-1:end),'h5'))
    data = read_hdf5(filename);
    if ~isstruct(data)
        disp('ERROR: File not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
        return;
    end
    data.datatype='BEAMS3D';
    data.X_lines=data.R_lines.*cos(data.PHI_lines);
    data.Y_lines=data.R_lines.*sin(data.PHI_lines);
    data.phiend = data.PHI_lines(data.npoinc,:);
    
    if ndims(data.X_lines) > 2
        data.X_lines=squeeze(data.X_lines(:,1,:));
        data.Y_lines=squeeze(data.Y_lines(:,1,:));
        data.Z_lines=squeeze(data.Z_lines(:,1,:));
    end
    data.lwall=0;
    if isfield(data,'wall_strikes')
        data.wall_strikes = double(data.wall_strikes);
        data.lwall=1;
    end
    % Catch some old format issues
    if isfield(data,'ns_prof1')
        data.ns_prof = data.ns_prof1;
    end
elseif (strcmp(filename(1:12),'beams3d_diag'))
    
    fid=fopen(filename,'r');
    line=fgetl(fid);
    if (strfind(line,'BEAMLINES'))  % New multi-beam file format
        data.nbeams = sscanf(line,'%*s %i',1);
        for i=1:data.nbeams
            line = fgetl(fid);
            line = fgetl(fid);
            data.energy(i) = sscanf(line,'%*i %g %*g %*g',1);
            data.charge(i) = sscanf(line,'%*i %*g %g %*g',1);
            data.mass(i)   = sscanf(line,'%*i %*g %*g %g',1);
            line = fgetl(fid);
            line = fgetl(fid);
            data.npart(i)=sscanf(line,'%i %*i %*d',1);
            data.nlost(i)=sscanf(line,'%*i %i %*d',1);
            data.lost_frac(i)=sscanf(line,'%*i %*i %e',1);
            data.t_end(i)=sscanf(line,'%*i %*i %*e %e',1);
            line=fgetl(fid);
            line=fgetl(fid);
            line=fgetl(fid);
            data.vllaxis=sscanf(line(4:end),'%e');
            line=fgetl(fid);
            data.vlldist(i,:,:)=fscanf(fid,'%d',[length(data.vllaxis) inf]);
            line=fgetl(fid);
            line=fgetl(fid);
        end
    else  % Old format
        line=fgetl(fid);
        data.npart=sscanf(line,'%i %*i %*d',1);
        data.nlost=sscanf(line,'%*i %i %*d',1);
        data.lost_frac=sscanf(line,'%*i %*i %e',1);
        data.t_end=sscanf(line,'%*i %*i %*e %e',1);
        line=fgetl(fid);
        line=fgetl(fid);
        line=fgetl(fid);
        data.vllaxis=sscanf(line(4:end),'%e');
        line=fgetl(fid);
        data.vlldist=fscanf(fid,'%d',[length(data.vllaxis) inf]);
    end
    fclose(fid);
end

end

