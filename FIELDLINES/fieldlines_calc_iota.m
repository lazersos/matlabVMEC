function [rho, iota, iota_err] = fieldlines_calc_iota(lines)
%FIELDILNES_INPUT Generates the FIELDLINES input namelist.
%   Detailed explanation goes here
x = lines.R_lines - lines.R_lines(1,:);
y = lines.Z_lines - lines.Z_lines(1,:);
%plot(x(1:10:end,:),y(1:10:end,:),'k.','MarkerSize',0.1)
theta = atan2(y,x);
dtheta = diff(theta,[],2); %axis=0?
dtheta(dtheta<-pi) = dtheta(dtheta<-pi)+2*pi;
dtheta(dtheta>pi) = dtheta(dtheta>pi)-2*pi;
theta = abs(cumsum(dtheta,2));
rho = mean(sqrt(x.^2+y.^2),2);
iota = zeros(size(rho));
iota_err = zeros(size(rho));
for i = 1:lines.nlines
    p = polyfit(lines.PHI_lines(i,1:lines.nsteps-1),theta(i,:),1);
%     f = polyval(p,lines.PHI_lines(i,1:lines.nsteps-1));
%     hold off
%     plot(lines.PHI_lines(i,1:lines.nsteps-1),theta(i,:),'o');
%     hold on
%     plot(lines.PHI_lines(i,1:lines.nsteps-1),f,'-');
%p = polyfit(lines.PHI_lines(i,:),theta(i,:),1);
    iota(i) = p(1);
    %iota_err(i) = sqrt(mu(2));
end
iota(1) = 2*iota(2)-iota(3);
figure
plot(iota);
end

