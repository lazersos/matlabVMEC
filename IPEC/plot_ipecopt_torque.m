function plot_ipecopt_torque
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%plot(ipec_data.NTV_s(1,:),ipec_data.NTV_target(1,:),'r','LineWidth',2.0)
figure('Position',[1 1 1024 768],'Color','white');
hold on;

file=dir('pent_*');
for i=1:length(file)
    pent=read_ipec(file(i).name);
    plot(pent.psi,pent.T_phi./pent.dvdpsi,'Color',[.31,.31,.31]);
end
plot(pent.psi,pent.T_phi./pent.dvdpsi,'k','LineWidth',2.0);
hold off;
set(gca,'FontSize',18);
xlabel('Rho');
ylabel('Torque Density [N/m^{-2}]')
yr=ylim;
dy=diff(yr)/20;
y1=max(yr);
text(0.75,y1-1*dy,['Total Torque: ' num2str(pent.T_phi_total) 'Nm'],'FontSize',18);

try
    file=dir('ipecopt.*');
    ipec_data=read_stellopt(file.name);
    xdata=read_stellopt('xvec.dat');
    iter_min=ipec_data.iter(end);
    dex=find(xdata.iter == iter_min);
    ppu_amp=xdata.x(1,dex); ppu_phase=xdata.x(2,dex);
    ppl_amp=xdata.x(3,dex); ppl_phase=xdata.x(4,dex);
    rwm_amp=xdata.x(5,dex).*2; rwm_phase=xdata.x(6,dex);
    text(0.05,y1-1*dy,['NCC Upper: ' num2str(ppu_amp/1000) ' kA-t @ ' num2str(180*mod(ppu_phase,2*pi)/pi) '°'],'FontSize',18);
    text(0.05,y1-2*dy,['RWM : ' num2str(rwm_amp/1000) ' kA-t @ ' num2str(180*mod(rwm_phase,2*pi)/pi) '°'],'FontSize',18);
    text(0.05,y1-3*dy,['NCC Lower: ' num2str(ppl_amp/1000) ' kA-t @ ' num2str(180*mod(ppl_phase,2*pi)/pi) '°'],'FontSize',18);
end


end

