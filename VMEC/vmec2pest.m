function u = vmec2pest(s,theta,zeta,lmns,xm,xn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dth = 1.0;
n1  = 0;
th  = theta;
th1 = theta;
while (abs(dth) >= 1.0E-05 && n1< 500)
    lam = sfunct(th,zeta,lmns(:,s),xm,xn);
    dlam = cfunct(th,zeta,lmns(:,s).*xm',xm,xn);
    dth = -(th + lam - th1) / (1+dlam);
    th  = th + 0.5*dth;
    n1 = n1 + 1;
end
u = th;
end

