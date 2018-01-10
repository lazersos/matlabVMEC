function jll = calc_jll( vmec_data,theta,zeta )
%CALC_JLL(vmec_data,theta,zeta) Calculates the parallel current density.
% This funciton takes a VMEC data structure (as read by read_vmec) and
% theta/zeta arrays as input and outputs the parallel current density.
%
% Example usage
%      theta=0:2*pi/359:2*pi;
%      zeta=0:2*pi/63:2*pi;
%      data=read_vmec('wout.test');        % Reads VMEC wout file
%      jll=calc_jll(vmec_data,theta,zeta); % Calculate the current
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

b=cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm,vmec_data.xn);
g=cfunct(theta,zeta,vmec_data.gmnc,vmec_data.xm,vmec_data.xn);
bu=cfunct(theta,zeta,vmec_data.bsubumnc,vmec_data.xm,vmec_data.xn);
bv=cfunct(theta,zeta,vmec_data.bsubvmnc,vmec_data.xm,vmec_data.xn);
ju=cfunct(theta,zeta,vmec_data.currumnc,vmec_data.xm,vmec_data.xn);
jv=cfunct(theta,zeta,vmec_data.currvmnc,vmec_data.xm,vmec_data.xn);

if (vmec_data.iasym)
    b=b+sfunct(theta,zeta,vmec_data.bmns,vmec_data.xm,vmec_data.xn);
    g=g+sfunct(theta,zeta,vmec_data.gmns,vmec_data.xm,vmec_data.xn);
    bu=bu+sfunct(theta,zeta,vmec_data.bsubumns,vmec_data.xm,vmec_data.xn);
    bv=bv+sfunct(theta,zeta,vmec_data.bsubvmns,vmec_data.xm,vmec_data.xn);
    ju=ju+sfunct(theta,zeta,vmec_data.currumns,vmec_data.xm,vmec_data.xn);
    jv=jv+sfunct(theta,zeta,vmec_data.currvmns,vmec_data.xm,vmec_data.xn);
end

jll = (bu.*ju+bv.*jv)./(g.*b);

return;
end

