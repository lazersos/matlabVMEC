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
beam_data.birth = beams3d_calc_depo(beam_data);
end
[s, ~, dVds] = beams3d_volume(beam_data);
x=0:1.0/beam_data.ns_prof1:1;
c=(x(1:end-1)+x(2:end)).*0.5;
%%
%ind=3:1000;
dex=beam_data.end_state==0|beam_data.end_state==1; %Only particles that arent promptly lost
f=sqrt(beam_data.S_lines(ind,dex));

rhoedges=[0 beam_data.rho(1:end-1)+diff(beam_data.rho)/2, 1, 1.5]; %correct rho edges
pitchedges=linspace(-1.,1.,npitch);

[discrho, ~]=discretize(f,rhoedges);
[discpitch, ~]=discretize(beam_data.pitch(ind,dex),pitchedges);
discpitch(isnan(discpitch))=numel(pitchedges)+1;
discrho(isnan(discrho))=numel(rhoedges);

tmp=repmat(beam_data.Weight(dex)./numel(ind),1,numel(ind))'; %repeat and normalize
beam=repmat(beam_data.Beam(dex),1,numel(ind))';

srhoxi=accumarray([discrho(:),discpitch(:),beam(:)],tmp(:),[numel(rhoedges), numel(pitchedges)+1, max(beam_data.Beam)]);
srhoxi=srhoxi(1:end-2,1:end-1,:); %Truncate ends from nan and normalize to density

%%
v=pchip([0 sqrt(s)],[0 dVds.*2.*sqrt(s)],c);
srhoxi=numel(pitchedges)/2*numel(rhoedges)*srhoxi./v'; %Volume normalization
beam_data.srhoxi=srhoxi;
beam_data.V=v;
beam_data.c=c;
if lplot
pixplot(beam_data.rho,pitchedges,squeeze(sum(srhoxi,3)))
colorbar
figure
plot(beam_data.rho,trapz(pitchedges,squeeze(sum(srhoxi,3)),2),'DisplayName','Integral from S(rho,xi)')
%plot(beam_data.rho,sum(birth,2))
hold on
plot(beam_data.rho,sum(beam_data.birth,1),'DisplayName','From calc depo')
legend()
end
end