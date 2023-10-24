function [srhoxi,beam_data] = beams3d_calc_srhoxi(beam_data,ind,npitch,varargin)
%beams3d_calc_srhoxi calculates a birth map in radial (rho) and pitch (xi)
%coordinates to be used as input for slowing down calculation.
lplot=0;
if strcmp(varargin,'plot')
    lplot=1;
end
if ~isfield(beam_data,'vperp')
beam_data.vperp=beams3d_calc_vperp(beam_data);
end
if ~isfield(beam_data,'pitch')
beam_data.pitch=beam_data.vll_lines./sqrt(beam_data.vll_lines.^2+beam_data.vperp.^2);
end
if ~isfield(beam_data,'birth')
beam_data.birth = beams3d_calc_depo(beam_data,'depo_only');
end

[s, ~, dVds] = beams3d_volume(beam_data);
x=0:1.0/beam_data.ns_prof1:1;
%c=(x(1:end-1)+x(2:end)).*0.5;
%Match beams3d_profiles here
nc = 128;
dc = 1./(nc);
rhoedges = 0:dc:1;
c=(rhoedges(1:end-1)+rhoedges(2:end)).*0.5;

dex=beam_data.end_state==0|beam_data.end_state==1; %Only particles that arent promptly lost
f=sqrt(beam_data.S_lines(ind,dex));

%rhoedges=[0 beam_data.rho(1:end-1)+diff(beam_data.rho)/2, 1, 1.5]; %correct rho edges
pitchedges=linspace(-1.,1.,npitch);

[discrho, ~]=discretize(f,rhoedges);
[discpitch, ~]=discretize(beam_data.pitch(ind,dex),pitchedges);
discpitch(isnan(discpitch))=numel(pitchedges)+1;
discrho(isnan(discrho))=numel(rhoedges);

tmp=repmat(beam_data.Weight(dex)./numel(ind),1,numel(ind))'; %repeat and normalize
beam=repmat(beam_data.Beam(dex),1,numel(ind))';

srhoxi=accumarray([discrho(:),discpitch(:),beam(:)],tmp(:),[numel(rhoedges), numel(pitchedges)+1, max(beam_data.Beam)]);
srhoxi=srhoxi(1:end-1,1:end-1,:); %Truncate ends from nan and normalize to density

v=pchip([0 sqrt(s)],[0 dVds.*2.*sqrt(s)],c);
srhoxi=numel(pitchedges)/2*numel(rhoedges)*srhoxi./v'; %Volume normalization
tmp1 = trapz(pitchedges,squeeze(sum(srhoxi,3)),2);
tmp2=sum(beam_data.birth,1);
tmp2=interp1(beam_data.dist_rhoaxis,tmp2,c,'linear','extrap');
tmp3=tmp2'./tmp1;
srhoxi=srhoxi./tmp3;
tmp4 = trapz(pitchedges,squeeze(sum(srhoxi,3)),2);
beam_data.srhoxi=srhoxi;
beam_data.V=v;
beam_data.c=c;
if lplot
pixplot(c,pitchedges,squeeze(sum(srhoxi,3)))
colorbar
figure
hold on
plot(c,tmp1,'DisplayName','Integral from S(rho,xi)')
plot(c,tmp4,'DisplayName','Renorm S(rho,xi)')
plot(c,tmp2,'--','DisplayName','From calc depo')
legend()
end
end