function data = read_beams3d(filename, varargin)
%READ_BEAMS3D Reads the HDF5 file created by BEAMS3D
% This funciton reads the fieldlines file and returns the data from the
% file in a structure.
%
% Example usage
%      data=read_beams3d('beams3d_test.h5');  % Reads BEAMS3D HDF5 file
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Handle varargin
if ~isempty(varargin)
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'beams'
            
        end
        i=i+1;
    end
end

if (strcmp(filename(end-1:end),'h5'))
    data = read_hdf5(filename);
    if ~isstruct(data)
        disp('ERROR: File not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
        return;
    end
    data.datatype='BEAMS3D';
    if ~isfield(data,'R_lines')
        data.npoinc = 0;
        disp('WARNING: Incomplete file no particle data!');
        return;
    end
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
    if data.VERSION < 3
        disp('Old version be careful with distribution function');
        % Fix non-double values
        for k={'ns_prof1','ns_prof2','ns_prof3','ns_prof4','ns_prof5','Beam',...
                'end_state','neut_lines','wall_faces'}
            if isfield(data,k{1})
                data.(k{1})=double(data.(k{1}));
            end
        end
        % Catch some old format issues
        if isfield(data,'ns_prof1')
            data.ns_prof = data.ns_prof1;
            h = 1./data.ns_prof1;
            data.rho = h/2.:h:(1-h/2);
        end
        if isfield(data,'dist_prof') && ~isfield(data,'dist2d_prof')
            % This is a patch for 2018a still used by IPP-HGW because of
            % reasons
            data.dist2d_prof = squeeze(sum(sum(sum(data.dist_prof,4),3),2));
            %data.dist2d_prof = squeeze(sum(data.dist_prof,[2 3 4]));
        end
        % Fix quantities multiplied not multiplied by drho
        if data.VERSION<2.7
            drho = 1./double(data.ns_prof1);
            data.ndot_prof   = data.ndot_prof ./drho;
            data.epower_prof = data.epower_prof./drho;
            data.ipower_prof = data.ipower_prof./drho;
            data.j_prof      = data.j_prof./drho;
            if isfield(data,'dense_prof')
                data.dense_prof  = data.dense_prof./drho;
            end
        end
        % Make the 5D Axis variables
        data.dist_rhoaxis=(double(1:data.ns_prof1)-0.5)./(data.ns_prof1);
        data.dist_uaxis=2.*pi.*(double(1:data.ns_prof2)-0.5)./(data.ns_prof2);
        data.dist_paxis=2.*pi.*(double(1:data.ns_prof3)-0.5)./(data.ns_prof3);
        d4=data.partvmax.*2./(data.ns_prof4);
        d5=data.partvmax./data.ns_prof5;
        data.dist_Vaxis= (-data.partvmax+d4.*0.5):d4:(data.partvmax-d4.*0.5);
        data.dist_Waxis=data.partvmax.*(double(1:data.ns_prof5)-0.5)./(data.ns_prof5);
    else
        % Fix non-double values
        for k={'ndist1','ndist2','ndist3','ndist3','ndist1','Beam',...
                'end_state','neut_lines','wall_faces','npoinc'}
            data.(k{1})=double(data.(k{1}));
            data.ns_prof1 = size(data.ipower_prof,2);
        end
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

