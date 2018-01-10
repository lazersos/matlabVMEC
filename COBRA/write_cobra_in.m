function write_cobra_in(vmec_ext,ntheta,nzeta,ns)
%WRITE_CORBRA_IN(vmec_ext,ntheta,nzeta) Creates a in_cobra.ext file.
%   The write_cobra_in function outputs a in_cobra.ext file in the current
%   directory.  This is done so that the cobra code can be run.  
%
%   Example:
%       write_cobra_in('test',10,10,99);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           8/22/13

fid = fopen(['in_cobra.' vmec_ext],'w');
fprintf(fid,'%d\n',nzeta);
fprintf(fid,'%d\n',ntheta);
fprintf(fid,'%s\n',vmec_ext);
zeta = 0:2*pi/(nzeta):2*pi;
theta =  0:2*pi/(ntheta):2*pi;
for i=1:nzeta
    fprintf(fid,' %10.5f ',zeta(i));
end
fprintf(fid,'\n');
for j=1:ntheta
    fprintf(fid,' %10.5f ',theta(j));
end
fprintf(fid,'\n');
fprintf(fid,'%d\n',ns-1);
fprintf(fid,' %d ',2:ns);
fprintf(fid,'\n');
fclose(fid);

end

