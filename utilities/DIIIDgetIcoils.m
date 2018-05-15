function DIIIDgetIcoils( shot_num,time )
%DIIIDgetIcoils( shot_num,time ) Summary of this function goes here
%   DIIIDgetIcoils( shot_num,time ) prints I-coil currents for a given
%   timesliece in [ms].

server = 'atlas.gat.com';  % Server you wish to connect to
server_port = 8000;        % Server port
tree = 'd3d';              % Tree you wish to open
import MdsPlus.*;   % So we can call MdsPlus java calls
mds_server=MdsPlus(server,server_port);  % Open the connection
mds_server.OpenTree(tree,shot_num);  % Open the Tree
%% I Coils
icoil_names={'IU330' 'IU270' 'IU210' 'IU150' 'IU90' 'IU30'...
             'IL330' 'IL270' 'IL210' 'IL150' 'IL90' 'IL30'};
for i =1:length(icoil_names);
    coil_str=['PTDATA("' upper(icoil_names{i}) '")'];
    coil_time_str=['DIM_OF(PTDATA("' upper(icoil_names{i}) '"))'];
    coil_name=upper(icoil_names{i});
    coil.(coil_name)=mds_server.Value(coil_str).Double;
    coil.([coil_name '_time'])=mds_server.Value(coil_time_str).Double;
end
    extcur(22)=pchip(coil.IU330_time,coil.IU330,time);
    extcur(23)=pchip(coil.IU270_time,coil.IU270,time);
    extcur(24)=pchip(coil.IU210_time,coil.IU210,time);
    extcur(25)=pchip(coil.IU150_time,coil.IU150,time);
    extcur(26)=pchip(coil.IU90_time,coil.IU90,time);
    extcur(27)=pchip(coil.IU30_time,coil.IU30,time);
    extcur(28)=pchip(coil.IL330_time,coil.IL330,time);
    extcur(29)=pchip(coil.IL270_time,coil.IL270,time);
    extcur(30)=pchip(coil.IL210_time,coil.IL210,time);
    extcur(31)=pchip(coil.IL150_time,coil.IL150,time);
    extcur(32)=pchip(coil.IL90_time,coil.IL90,time);
    extcur(33)=pchip(coil.IL30_time,coil.IL30,time);
    
j=1;
for i=22:33
    fprintf('  EXTCUR(%2.2d) = %20.10E    !%s\n',i,extcur(i),icoil_names{j});
    j=j+1;
end

end

