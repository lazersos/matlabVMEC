function varargout = beams3d_profiles(beam_data)
%function [s, ne, te, ti, zeff] = beams3d_profiles(beam_data)
%BEASM3D_PROFILES Extracts profile data from BEAMS3D Data
%   The BEAMS3D_PROFILES function calculates the radial kinetic profiles
%   (in flux) from the background grid data.  The returned profiles are on
%   a 128 point grid.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [s, ne, te, ti, zeff] = beams3d_profiles(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.5


ns = 128;
ds = 1./(ns-1);
s = 0:ds:1;
ne=[]; te=[]; ti=[]; zeff=[]; ni=[];

% Extract
[C, IA, ~] = unique(beam_data.S_ARR);
te_temp = beam_data.TE(IA);
ne_temp = beam_data.NE(IA);
ti_temp = beam_data.TI(IA);
zeff_temp = beam_data.ZEFF_ARR(IA);
if isfield(beam_data,'NI')
    nni=size(beam_data.NI,1);
    for i=1:nni
        ni_temp(i,:)=beam_data.NI(i,IA);
    end
else
    ni_temp=[];
end

% Make mirror
C=[-flipud(C); C];
ne_temp = [flipud(ne_temp); ne_temp];
te_temp = [flipud(te_temp); te_temp];
ti_temp = [flipud(ti_temp); ti_temp];
zeff_temp = [flipud(zeff_temp); zeff_temp];
if ~isempty(ni_temp)
    ni_temp = [flipud(ni_temp'); ni_temp']';
end

% spline
[C, IA, ~] = unique(C);
ne = pchip(C,ne_temp(IA),s);
te = pchip(C,te_temp(IA),s);
ti = pchip(C,ti_temp(IA),s);
zeff = pchip(C,zeff_temp(IA),s);
if ~isempty(ni_temp)
    for i=1:nni
        ni(i,:) = pchip(C,ni_temp(i,IA),s);
    end
end

% Output 
varargout{1}=s;
varargout{2}=ne;
varargout{3}=te;
varargout{4}=ti;
varargout{5}=zeff;
if ~isempty(ni), varargout{6}=ni; end

end

