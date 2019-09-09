function get_ipecopt_current( val ,nn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



file=dir('ipecopt.*');
ipec_data=read_stellopt(file.name);
xdata=read_stellopt('xvec.dat');
if isempty(val)
    iter_min=ipec_data.iter(end);
    dex=find(xdata.iter == iter_min);
end

% Just the xvec
disp(['COIL_AMP(1) = ' num2str(xdata.x(1,dex),' %20.10E ')]);
disp(['COIL_PHASE(1) = ' num2str(xdata.x(2,dex),' %20.10E ')]);
disp(['COIL_AMP(2) = ' num2str(xdata.x(3,dex),' %20.10E ')]);
disp(['COIL_PHASE(2) = ' num2str(xdata.x(4,dex),' %20.10E ')]);
disp(['COIL_AMP(3) = ' num2str(xdata.x(5,dex),' %20.10E ')]);
disp(['COIL_PHASE(3) = ' num2str(xdata.x(6,dex),' %20.10E ')]);


% FULL
disp('coil_name(1)="ppu"');
disp('coil_name(2)="ppl"');
disp('coil_name(3)="rwmef"');
disp(['coil_cur(1,1:6) = ' num2str(xdata.x(1,dex).*cos(xdata.x(2,dex) + nn*2*pi*(0:11)./12),' %20.10E ')]);
disp(['coil_cur(2,1:6) = ' num2str(xdata.x(3,dex).*cos(xdata.x(4,dex) + nn*2*pi*(0:11)./12),' %20.10E ')]);
disp(['coil_cur(3,1:6) = ' num2str(xdata.x(5,dex).*cos(xdata.x(6,dex) + nn*2*pi*(0:5)./6),' %20.10E ')]);

% Partial
disp('coil_name(1)="ppu"');
disp('coil_name(2)="ppl"');
disp('coil_name(3)="rwmef"');
disp(['coil_cur(1,1:6) = ' num2str(xdata.x(1,dex).*cos(xdata.x(2,dex) + nn*2*pi*(0:5)./6),' %20.10E ')]);
disp(['coil_cur(2,1:6) = ' num2str(xdata.x(3,dex).*cos(xdata.x(4,dex) + nn*2*pi*(0:5)./6),' %20.10E ')]);
disp(['coil_cur(3,1:6) = ' num2str(xdata.x(5,dex).*cos(xdata.x(6,dex) + nn*2*pi*(0:5)./6),' %20.10E ')]);


end

