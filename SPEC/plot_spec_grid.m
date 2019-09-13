function hplot = plot_spec_grid(data,nzeta)
%PLOT_SPEC_GRID Plots the SPEC grid.
%   The PLOT_SPEC_GRID function plots the finite element grid stored in a
%   SPEC data structure as read by read_spec.  
%
% Example usage
%      data=read_spec('spec_test.h5');  % Reads SPEC HDF5 file
%      plot_spec_grid(data);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0


nvol=size(data.R_grid,2);


% A note about the grid
% data.R_grid has two dimesions the first is ntheta*nphi*nr the second is
% nvol.  The first dimension is vectorized.  Order of these arrays is
% ntheta, nphi, nr.
% Plot axis radial segments
% Create the initial radial offset
offset = (nzeta-1)*data.ntheta_grid +1;
offset_rho=data.ntheta_grid*data.nzeta_grid;
plot(data.R_grid(offset,1),data.Z_grid(offset,1),'+');
nr = data.grid_num(1);  % Get number of subgrids
% From axis to first grid
for j=1:data.ntheta_grid
    hold on
    x1 = data.R_grid(offset,1);
    x2 = data.R_grid(offset+offset_rho+j,1);
    y1 = data.Z_grid(offset,1);
    y2 = data.Z_grid(offset+offset_rho+j,1);
    plot([x1 x2],[y1 y2],'k');
    hold off
end
% From 1st grid to last subgrid
for k=2:nr-1
    for j=1:data.ntheta_grid
        hold on
        x1 = data.R_grid(offset+offset_rho*(k-1)+j,1);
        x2 = data.R_grid(offset+offset_rho*(k-1)+j+1,1);
        y1 = data.Z_grid(offset+offset_rho*(k-1)+j,1);
        y2 = data.Z_grid(offset+offset_rho*(k-1)+j+1,1);
        plot([x1 x2],[y1 y2],'k');
        x1 = data.R_grid(offset+offset_rho*(k-1)+j,1);
        x2 = data.R_grid(offset+offset_rho*k+j,1);
        y1 = data.Z_grid(offset+offset_rho*(k-1)+j,1);
        y2 = data.Z_grid(offset+offset_rho*k+j,1);
        plot([x1 x2],[y1 y2],'k');
        hold off
    end
end
% Finish rest of subgrids
for l=2:data.nvol_grid
    nr = data.grid_num(l);
    for k=1:nr-1
        for j=1:data.ntheta_grid
            hold on
            x1 = data.R_grid(offset+offset_rho*(k-1)+j,l);
            x2 = data.R_grid(offset+offset_rho*(k-1)+j+1,l);
            y1 = data.Z_grid(offset+offset_rho*(k-1)+j,l);
            y2 = data.Z_grid(offset+offset_rho*(k-1)+j+1,l);
            plot([x1 x2],[y1 y2],'k');
            x1 = data.R_grid(offset+offset_rho*(k-1)+j,l);
            x2 = data.R_grid(offset+offset_rho*k+j,l);
            y1 = data.Z_grid(offset+offset_rho*(k-1)+j,l);
            y2 = data.Z_grid(offset+offset_rho*k+j,l);
            plot([x1 x2],[y1 y2],'k');
            hold off
        end
    end
end
% Do last line
k=nr;
for j=1:data.ntheta_grid
    hold on
    x1 = data.R_grid(offset+offset_rho*(k-1)+j,l);
    x2 = data.R_grid(offset+offset_rho*(k-1)+j+1,l);
    y1 = data.Z_grid(offset+offset_rho*(k-1)+j,l);
    y2 = data.Z_grid(offset+offset_rho*(k-1)+j+1,l);
    plot([x1 x2],[y1 y2],'k');
    hold off
end
axis equal
hplot=gca;
return
end