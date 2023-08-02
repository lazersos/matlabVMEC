function out = beams3d_from_fidasim(filename,vmec_data)
    % Append '_equilibrium.h5' to the filename if not already present
    [~, ~, ext] = fileparts(filename);
    if ~strcmp(ext, '.h5')
        filename = [filename, '_equilibrium.h5'];
    end
    
    % Read the contents of the HDF5 file
    data = read_hdf5(filename);

    out.s=data.plasma.profiles.rho.^2';
    out.ne=data.plasma.profiles.dene'*1e6;
    out.te=data.plasma.profiles.te'*1e3;
    out.ti=data.plasma.profiles.ti'*1e3;
    out.zeff=data.plasma.profiles.zeff';
species='D';
        if strcmp(species,'H')
        NI_AUX_M=[  1.6735576929E-27 1.795208784901260e-26];
    elseif strcmp(species,'D')
        NI_AUX_M=[  3.344325680132399e-27 1.795208784901260e-26];
    end
    NI_AUX_Z  = [1 5];%D and B
    NI_AUX  = [NI_AUX_M; NI_AUX_Z];
    ni_temp=out.ne;
    %n_fi_timepoint=spline(transp.X(:,transp_t_ind),transp.BDENS(:,transp_t_ind)*1e6,rho_b3d);
    %ni_temp = ni_temp-n_fi_timepoint; %Fast ion displacement of thermal ions
    c_B=(out.zeff-NI_AUX_Z(1))/(NI_AUX_Z(2).^2-NI_AUX_Z(2)*NI_AUX_Z(1));
    out.ni = [        ni_temp.*(1-NI_AUX_Z(2).*c_B)/NI_AUX_Z(1);
        ni_temp.*c_B];
%        beams3d_beamnamelist(vmec_data,energyfrac,power,r_beam,p_beam,z_beam,div_beam,'pfrac',powerfrac,...
%         'beam_dex',beam_dex,'per_beam','beamlet',species,...
%         'TE', Te_in, 'NE', ne_in,'TI',Ti_in, 'ZEFF', Zeff_in, 'POT',Pot_in,'t_end', 2.0e-1,'NPOINC', 5, 'NPARTICLES_START',4096,...
%         'nr', 101,'nz', 101,'nphi', 16, 'ZMEAN', ZMEAN, 'NI', NI_AUX, NI_IN,'FREE', fidasim_inputs,'file',fname);
end