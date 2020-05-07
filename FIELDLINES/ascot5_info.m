function ascot5_info(a5file)
%ASCOT5_INFO Outputs summary of groups in ASCOT5 file.
%   The ASCOT5_INFO subroutine prints the group strcuture of an ASCOT5
%   file.
%
%   Example:
%       a5file='ascot5_W7X_20180821_012_5100_h8.h5';
%       ascot5_info(a5file);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

root = h5info(a5file);
groups = root.Groups;
disp([a5file ' contains:']);
for i = 1 : length(groups)
    disp(['   ' groups(i).Name]);
    for j = 1 : length(groups(i).Groups)
        disp(['      ' groups(i).Groups(j).Name]);
        if strfind(groups(i).Groups(j).Name,'run_')
            try
                name = [groups(i).Groups(j).Name '/endstate/endcond'];
                data = h5read(a5file,name);
                disp(['          NPART_TOTAL:   ' num2str(length(data),'%i') ]);
                disp(['          NPART_MAX_SIM: ' num2str(sum(data==1),'%i') ]);
                disp(['          NPART_MIN_E:   ' num2str(sum(data==2),'%i') ]);
                disp(['          NPART_THERM:   ' num2str(sum(data==4),'%i') ]);
                disp(['          NPART_WALL:    ' num2str(sum(data==8),'%i') ]);
                disp(['          NPART_MAX_RHO: ' num2str(sum(data==10),'%i') ]);
                disp(['          NPART_MIN_RHO: ' num2str(sum(data==20),'%i') ]);
                disp(['          NPART_MAX_WAL: ' num2str(sum(data==100),'%i') ]);
            catch
            end
        end
    end
end


end

