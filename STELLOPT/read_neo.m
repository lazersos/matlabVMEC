function data = read_neo( filename )
%READ_NEO Reads NEO output files.
% This funciton reads the NEO output files.  This includes the neo_out,
% neolog, and bootstrap current file.  The value of these files is returned
% as an array.  See the VMECwiki for details.
%  http://vmecwiki.pppl.wikispaces.net/NEO
%
% Example usage
%      neo_log=('neolog.test');     % Reads NEO log file
%      neo_data=('neo_out.test');     % Reads NEO out file
%      neo_current=('neocurr.test');     % Reads NEO current file.
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

fid=fopen(filename,'r');
if strfind(filename,'neo_out')
    len = length(fgetl(fid));
    frewind(fid);
    if len <64
        data=fscanf(fid,'%e %e %e',[3 inf]);
    elseif len<120
        data=fscanf(fid,'%i %e %e %e %e %e',[6 inf]);
    else
        data=fscanf(fid,'%i %e %e %e %e %e %e %e %e %e %e %e %e',[13 inf]);
    end
elseif strfind(filename,'neolog')
    data=fscanf(fid,' %i %i %i %i %i %i %i %e\n',[8 inf]);
else
    try
        data=fscanf(fid,'%i %e %e %e %e %e',[6 inf]);
        
    catch
        disp('This was not a reconginzed file type!');
        data=-1;
    end
end
fclose(fid);

