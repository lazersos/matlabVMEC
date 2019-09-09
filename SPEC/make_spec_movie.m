function  make_spec_movie( filename )
%MAKE_SPEC_MOVIE(filename) Creates an iteration movie for a SPEC run
%   MAKE_SPEC_MOVIE(filename) Creates an iteration by iteration plot of
%   the SPEC calculation corresponding to filename.  Here filename is the
%   name of the run (part before .spec).  It reads the .filename.iterations
%   file for the surface, force balance, and energy information.
%
% Example usage
%      make_spec_movie('spec_calc');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% NOTES:
%   10/7/13     First working version uses fixed boundary non-stellarator
%               symmetric version of SPEC code.
theta = 0:2*pi/90:2*pi;
zeta  = 0;
spacer_format='float64';
int_format='int';
float_format='float64';
machine_format='n';
['.' filename '.iterations'];
try
fid=fopen(['.' filename '.iterations'],'r',machine_format);
catch
    return
end
temp=fread(fid,1,spacer_format);
mnmax=fread(fid,1,int_format);
nvol=fread(fid,1,int_format);
narr=mnmax*(nvol+1);
temp=fread(fid,1,spacer_format);
temp=fread(fid,1,spacer_format);
xm=real(fread(fid,mnmax,int_format))';
temp=fread(fid,1,spacer_format);
temp=fread(fid,1,spacer_format);
xn=real(fread(fid,mnmax,int_format))';
temp=fread(fid,1,spacer_format);
j=1;
while (~feof(fid))
    try
        temp=fread(fid,1,spacer_format);
        wflag(j)=fread(fid,1,int_format);
        iflag(j)=fread(fid,1,int_format);
        energy(j) = fread(fid,1,float_format);
        rflag(j)=fread(fid,1,float_format);
        temp=fread(fid,1,spacer_format);
        temp=fread(fid,1,spacer_format);
        rbc = fread(fid,narr,float_format);
        temp=fread(fid,1,spacer_format);
        temp=fread(fid,1,spacer_format);
        zbs = fread(fid,narr,float_format);
        temp=fread(fid,1,spacer_format);
        temp=fread(fid,1,spacer_format);
        rbs = fread(fid,narr,float_format);
        temp=fread(fid,1,spacer_format);
        temp=fread(fid,1,spacer_format);
        zbc = fread(fid,narr,float_format);
        temp=fread(fid,1,spacer_format);
        rbc = reshape(rbc,[mnmax nvol+1]);
        zbs = reshape(zbs,[mnmax nvol+1]);
        rbs = reshape(rbs,[mnmax nvol+1]);
        zbc = reshape(zbc,[mnmax nvol+1]);
        r{j}   = cfunct(theta,zeta,rbc,xm,xn) + sfunct(theta,zeta,rbs,xm,xn);
        z{j}   = cfunct(theta,zeta,zbc,xm,xn) + sfunct(theta,zeta,zbs,xm,xn);
        j=j+1;
    catch
        break
    end
end
fclose(fid);

% Do some vetting
dex = iflag == 0;
dex2 =circshift(dex,[1 -1]);
rflag(dex) = rflag(dex2);
energy(dex) = energy(dex2);
dex = energy > 10.*mean(energy);
energy(dex) = mean(energy);

subplot(1,2,2);
for j=1:length(r)
    if (rflag(j) == 0); continue; end;
    subplot(1,2,1);
    plot(r{1}',z{1}','k')
    hold on;
    plot(r{1}(1,1),z{1}(1,1),'ok')
    plot(r{j}',z{j}','r');
    plot(r{j}(1,1),z{j}(1,1),'+r')
    %plot(r{j}(:,[1 15 30 45 60 75]),z{j}(:,[1 15 30 45 60 75]),'r+');
    hold off;
    title('Flux Surfaces');
    xlabel('R [m]');
    ylabel('Z [m]');
    axis equal;
    subplot(1,2,2);
    %[ax,h1,h2]=plotyy(1:length(r),energy,1:length(r),log10(rflag));
    %xlim(ax(1),[1 length(r)]);
    %xlim(ax(2),[1 length(r)]);
    %ylabel(ax(1),'Energy');
    %ylabel(ax(2),'log(Force Err)');
    plot(1:length(r),log10(rflag),'k');
    ylabel('log(Force Err.)');
    xlim([1 length(r)])
    hold on;
    plot([j j],ylim,'r');
    hold off;
    title('Energy / Force Balance');
    xlabel('Iteration');
    pause(0.01);
end


end

