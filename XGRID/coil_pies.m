function coil_pies(coildata,extcur)
%COIL_PIES(coildata,extcur) Creates a PIES coil_data file.
%   This function creates a PIES coil_data file from a coildata structure
%   and an extcur array.
%
%   Usage:
%   coil_data=read_coils('coils.test');
%   extcur=[1.2e4 1.2e4 1.2e4 -3.5e3 1.1e6 -2.5e4];
%   coil_pies(coil_data,extcur);
%
%   See also read_coils, plot_coils.
%
%   Written by: S. Lazerson (lazerson@pppl.gov)
%   Verion:     1.0
%   Date:       12/21/10


if nargin==2
    % First sort the data by current group
    coildata.vert=sortrows(coildata.vert',5)';
    % Now check extcur array size
    nextcur=max(size(extcur));
    ngroups=max(coildata.vert(5,:));
    if ngroups==nextcur
        for i=1:ngroups
            test1=coildata.vert(5,:)==i;
            test2=coildata.vert(4,:) ~= 0.0;
            coildata.vert(4,logical(test1.*test2))=extcur(i);
        end
    else
        disp(['ERROR: ngroups=' num2str(ngroups)...
            '  nextcur=',num2str(nextcur)]);
        return
    end
    fid=fopen('coil_data','w+');
    %for i=1:max(size(coildata.vert))
        fprintf(fid,'%+20.10e  %+20.10e  %+20.10e  %+20.10e\n',coildata.vert(1:4,:));
    %end
    fclose(fid);
else
    disp('ERROR: Usage coil_pies(coildata,extcur)');
end
return
end

