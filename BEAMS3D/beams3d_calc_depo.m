function birth = beams3d_calc_depo(beam_data, varargin)
%BEAMS3D_CALC_DEPO Calculated deposition
%   The BEAMS3D_CALC_DEPO routine calculates the radial birth profile for a
%   given BEAMS3D run.  Outputs in units of [part/m^3/s].
%
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       birth = beams3d_calc_depo(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0

nrho = [];
ldepo_only=0;

% Handle varargin
if ~isempty(varargin)
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'nrho'
                nrho = varargin{i+1};
            case 'depo_only'
                ldepo_only=1;
                dexd=(beam_data.end_state==0|beam_data.end_state==1);
        end
        i=i+1;
    end
end

% Setup rho array
if isempty(nrho)
    nrho = beam_data.ns_prof1;
end
edges = 0:1.0/nrho:1;
rho_out = 0.5.*(edges(1:end-1)+edges(2:end));

% Get Volume elements
[s,~,dVds]=beams3d_volume(beam_data,'ns',nrho);
dVdrho = pchip(sqrt(s),2.*sqrt(s).*dVds,rho_out);

% Define starting points
dexs = 1;
if beam_data.lbeam
    dexs=2;
end

% Calc births
birth = zeros(beam_data.nbeams,length(rho_out));
rho_lines = ones(1,beam_data.nparticles);
for i=1:beam_data.nbeams
    dexb = beam_data.Beam == i;
    if ldepo_only
        dexb=and(dexb,dexd);
    end
    rho_lines(:) = 1.5;
    rho_lines(dexb) = sqrt(beam_data.S_lines(dexs,dexb));
    for j=1:nrho
        dexi = and(rho_lines(:)>=edges(j),rho_lines(:)<edges(j+1));
        birth(i,j) = sum(beam_data.Weight(dexi));
    end
    birth(i,:) = nrho*birth(i,:)./dVdrho;
end

end

