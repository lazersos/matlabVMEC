function vmec_metrics( vmec_data )
%VMEC_METRICS(vmec_data)  Plots metric elements of the vmec_data structure.
%   VMEC_METRICS(vmec_data) Creates two plots of the metric elements based
%   on the VMEC rmnc and zmns arrays.  
%
% Example usage
%      data=read_vmec('wout.test');     % Reads VMEC wout file
%      vmec_metrics(vmec_data);         % Makes a plot
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% NOTES:
%   5/27/14      Created!
ns=vmec_data.ns;
for i=1:ns, dRdu(:,i)=-vmec_data.xm'.*vmec_data.rmnc(:,i); end;
for i=1:ns, dRdv(:,i)=-vmec_data.xn'.*vmec_data.rmnc(:,i); end;
for i=1:ns, dZdu(:,i)=vmec_data.xm'.*vmec_data.zmns(:,i); end;
for i=1:ns, dZdv(:,i)=vmec_data.xn'.*vmec_data.zmns(:,i); end;
Rspl=pchip(vmec_data.phi/vmec_data.phi(end),vmec_data.rmnc);
Zspl=pchip(vmec_data.phi/vmec_data.phi(end),vmec_data.zmns);
dZdsspl=spline_deriv(Zspl,1);
dRdsspl=spline_deriv(Rspl,1);
dRds=ppval(dRdsspl,vmec_data.phi/vmec_data.phi(end));
dZds=ppval(dZdsspl,vmec_data.phi/vmec_data.phi(end));
gss=dRds.^2+dZds.^2;
guu=dRdu.^2+dZdu.^2;
gvv=dRdv.^2+dZdv.^2;
gsu=dRds.*dRdu+dZds.*dZdu;
gsv=dRds.*dRdv+dZds.*dZdv;
guv=dRdu.*dRdv+dZdu.*dZdv;
rho=vmec_data.phi./vmec_data.phi(end);
fig=figure('Position',[1 1 1920 1080],'Color','white');
subplot(3,3,1);
plot(rho,gss');
set(gca,'FontSize',12); xlabel('s'); ylabel('gss');
subplot(3,3,2);
plot(rho,gsu');
set(gca,'FontSize',12); xlabel('s'); ylabel('gsu');
subplot(3,3,3);
plot(rho,gsv');
set(gca,'FontSize',12); xlabel('s'); ylabel('gsv');
subplot(3,3,5);
plot(rho,guu');
set(gca,'FontSize',12); xlabel('s'); ylabel('guu');
subplot(3,3,6);
plot(rho,guv');
set(gca,'FontSize',12); xlabel('s'); ylabel('guv');
subplot(3,3,9);
plot(rho,gvv');
set(gca,'FontSize',12); xlabel('s'); ylabel('gvv');
fig2=figure('Position',[1 1 1920 1080],'Color','white');
subplot(2,3,1);
plot(rho,dRds);
set(gca,'FontSize',12); xlabel('s'); ylabel('dR/ds');
subplot(2,3,2);
plot(rho,dRdu);
set(gca,'FontSize',12); xlabel('s'); ylabel('dR/du');
subplot(2,3,3);
plot(rho,dRdv);
set(gca,'FontSize',12); xlabel('s'); ylabel('dR/dv');
subplot(2,3,4);
plot(rho,dZds);
set(gca,'FontSize',12); xlabel('s'); ylabel('dZ/ds');
subplot(2,3,5);
plot(rho,dZdu);
set(gca,'FontSize',12); xlabel('s'); ylabel('dZ/du');
subplot(2,3,6);
plot(rho,dZdv);
set(gca,'FontSize',12); xlabel('s'); ylabel('dZ/dv');
end

