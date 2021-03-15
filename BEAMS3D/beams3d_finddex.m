function dex = beams3d_finddex(beam_data,varargin)
%BEAMS3D_FINDDEX Returns indices of specific particles
%   The BEAMS3D_FINDDEX function creates an array which provides specific
%   indicies for each particle.  In general it is used to subsample a
%   specific population of particles at a specific timeslice.  It takes a
%   beams3d data structure as returned by read_beams3d and an optional
%   population type.
%   Options:
%       'orbit_birth':      Initial condition for any orbiting particle.
%       'orbit_last':       Last timepoint for orbiting particle.
%       'therm_birth':      Initial condition for any thermalized particle.
%       'therm_end':        Point at which particle thermalizes.
%       'wall_birth':       Initial condition for wall strikes.
%       'wall_hit':         Point of wall collision
%       'shine_birth':      Initial condition of shinethorugh particles.
%       'shine_hit':        Location of wall strikes for shinethrough.
%       'port_birth':       Initial condition for port strikes.
%       'port_hit':         Location of port strikes.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       dex = beam3d_finddex(beam_data,'wall_hit'); % Wall hits
%
%   Created by: S. Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.00

%Defaults
dex=[];
eorbit=0; etherm=1; ewall=2; eshine=3; eport=4;
llast = 0;
etarg = eorbit;

% Handle varargin
if ~isempty(varargin)
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'orbit_birth'
                etarg=eorbit;
                llast=0;
            case 'orbit_last'
                etarg=eorbit;
                llast=1;
            case 'therm_birth'
                etarg=etherm;
                llast=0;
            case 'therm_last'
                etarg=etherm;
                llast=1;
            case 'wall_birth'
                etarg=ewall;
                llast=0;
            case 'wall_hit'
                etarg=ewall;
                llast=1;
            case 'shine_birth'
                etarg=eshine;
                llast=0;
            case 'shine_hit'
                etarg=eshine;
                llast=1;
            case 'port_birth'
                etarg=eport;
                llast=0;
            case 'port_hit'
                etarg=eport;
                llast=1;
        end
        i=i+1;
    end
end

end_state=double(beam_data.end_state');

mask=~(end_state==etarg);

% Handle a hitonly run
if size(beam_data.R_lines,1) == 3
    disp('WARNING:  Possible hit_only run!');
    disp('      dex=1 Before hit');
    disp('      dex=2 Wall hit or last point');
    disp('      dex=3 Point beyond wall');
    dex = ones(1,beam_data.nparticles);
    % The in this case index=2 is always the last point even for orbiting
    % particles.
    if llast
        dex = dex + 1;
    end
    %return;
else    
    if llast
        dex=zeros(1,beam_data.nparticles);
        for i=1:beam_data.nparticles
            dex(i) = min([find(beam_data.R_lines(:,i)==0,1,'first')-1 beam_data.npoinc+1]);
        end
    else
        dex = ones(1,beam_data.nparticles);
    end
end
dex(mask) = 0;

return;

end

