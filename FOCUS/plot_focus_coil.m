function plot_focus_coil( ext_name )
%PLOT_FOCUS_COIL(ext_name) Plots the FOCUS coil with flux surface Bnormal.
%   Detailed explanation goes here

coil_name=[ext_name '.coils'];
data_name=['focus_' ext_name '.h5'];

coil_data=read_coils(coil_name);
data=read_focus(data_name);

x=reshape(data.xsurf,[1 numel(data.xsurf)]);
y=reshape(data.ysurf,[1 numel(data.ysurf)]);
z=reshape(data.zsurf,[1 numel(data.zsurf)]);
b=reshape(data.Bn,[1 numel(data.Bn)]);
s=ones(1,numel(data.Bn));

plot_coils(coil_data,'field_period');
hold on;
scatter3(x,y,z,s,b);
colormap hot;
colorbar;
end

