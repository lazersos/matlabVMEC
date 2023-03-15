function [R_lines,PHI_lines,Z_lines] = fieldlines_follow(in_data,start_loc,phi_extent,poinc_loc,grid_extent)
%FIELDLINES_FOLLOW Summary of this function goes here
%   Detailed explanation goes here
nstart = size(start_loc,1);
[r,phi,z] = ndgrid(in_data.raxis,in_data.phiaxis,in_data.zaxis);
maxphi=max(in_data.phiaxis);
dRdphi = r.*in_data.B_R ./ in_data.B_PHI;
dZdphi = r.*in_data.B_Z ./ in_data.B_PHI;

%Setup derivatives
dR_F = griddedInterpolant(r,phi,z,dRdphi,'cubic');
dZ_F = griddedInterpolant(r,phi,z,dZdphi,'cubic');


%Function to be integrated (q is a vector of length 2n with n entries for R
%and Z respectively.
    function dqdphi = f(phi,q)
        dqdphi = zeros(nstart*2,1);
        repphi=repmat(mod(phi,maxphi),nstart,1); %Periodic BC
        dqdphi(1:nstart) = dR_F(q(1:nstart),repphi,q(nstart+1:end));
        dqdphi(nstart+1:end)= dZ_F(q(1:nstart),repphi,q(nstart+1:end));
        dex=[q(1:nstart)<grid_extent(1) |...
            q(1:nstart)>grid_extent(2) ; ...
            q(nstart+1:end)<grid_extent(3) |...
            q(nstart+1:end)>grid_extent(4) ];
        dqdphi(dex)=0;
    end


    function [value,isterminal,direction] = events(t,~)
        value = mod(t-poinc_loc,maxphi)-pi;
        isterminal=0;
        direction=-1;
    end

options = odeset('RelTol',1e-5,'Events',@events); %'vectorized'
% options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
%    'Refine',refine);
%parpool(4)
% R_lines = zeros(nstart,numel(phi_extent));
% Z_lines = zeros(nstart,numel(phi_extent));
% PHI_lines = zeros(nstart,numel(phi_extent));
% R_lines(:,1)=start_loc(:,1);
% Z_lines(:,1)=start_loc(:,2);
% PHI_lines(:,1)=phi_extent;
starts = reshape(start_loc,1,[]);
[phi,y,phie,ye,ie] = ode89(@f,[phi_extent(1) phi_extent(2)],starts,options);
R_lines=ye(:,1:nstart)';
Z_lines=ye(:,nstart+1:end)';
PHI_lines=repmat(phie(:),1,nstart)';
%     R_lines(:,i)=y(end,1:nstart)';
%     Z_lines(:,i)=y(end,nstart+1:end)';
%     PHI_lines(:,i)=phi(end);

end



%
% function [R_lines,Z_lines] = fieldlines_follow(in_data,start_loc,phi_extent)
% %FIELDLINES_FOLLOW Summary of this function goes here
% %   Detailed explanation goes here
%
% [r,phi,z] = ndgrid(in_data.raxis,in_data.phiaxis,in_data.zaxis);
%
% dRdphi = r.*in_data.B_R ./ in_data.B_PHI;
% dZdphi = r.*in_data.B_Z ./ in_data.B_PHI;
%
% dR_F = griddedInterpolant(r,phi,z,dRdphi,'cubic');
% dZ_F = griddedInterpolant(r,phi,z,dZdphi,'cubic');
%
%
%     function dqdphi = f(phi,q)
%         dqdphi = [0; 0];
%         dqdphi(1) = dR_F(q(1),phi,q(2));
%         dqdphi(2) = dZ_F(q(1),phi,q(2));
%     end
%
% % options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,...
% %    'Refine',refine);
% [t,y] = ode45(@f,[phi_extent phi_extent+2*pi],start_loc);
% R_lines=y(:,1);
% Z_lines=y(:,2);
% %    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
% %       'MaxStep',t(nt)-t(1));
% end