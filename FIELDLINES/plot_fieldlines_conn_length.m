function out_handle = plot_fieldlines_conn_length(ext_name)
%PLOT_FIELDLINES_CONN_LENGTH Makes a connection length plot
%   The PLOT_FIELDLINES_CONN_LENGTH routine generates a connection length
%   plot from two fieldlines HDF5 files with extensions <ext_name> and
%   <ext_name>_rev.  It also requires in input file with FIELDLINES_INPUT
%   fortran namelist.  
%
%   Usage:
%       plot_fieldlines_conn_length('test');
%
%       Here input.test, fieldlines_test.h5, and fieldlines_test_rev.h5
%       must all be in the same directory.
%
%   Created by: S. Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.0
%   Date:       03/04/20


out_handle = [];
marker_size = 36;

% Read data
infile = ['input.' ext_name];
file1 = ['fieldlines_' ext_name '.h5'];
file2 = ['fieldlines_' ext_name '_rev.h5'];
try
    input = read_namelist(infile,'FIELDLINES_INPUT');
catch
    disp('Cannot find input file.');
    disp(['   File: ' infile]);
    return;
end
try
    data1 = read_fieldlines(file1);
catch
    disp('Cannot find forward file.');
    disp(['   File: ' file1]);
    return;
end
try
    data2 = read_fieldlines(file2);
catch
    disp('Cannot find reverse file.');
    disp(['   File: ' file2]);
    return;
end

% Extract grid
R_plt= input.r_start;
Z_plt = input.z_start;

% Extract distance
L_plt = log10(data1.L_lines+data2.L_lines);
S_plt = 0.*L_plt+marker_size;

fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
scatter(R_plt,Z_plt,S_plt,L_plt,'filled');
set(gca,'FontSize',24);
xlabel('R [m]');
ylabel('Z [m]');
title(['Connection Length (' strrep(ext_name,'_','\_') ')']);
colormap jet;
ha_cbar=colorbar;
set(ha_cbar,'FontSize',24);
ylabel(ha_cbar,'Log_{10}(Conn. Len.)');
axis equal;
out_handle=fig;

end

