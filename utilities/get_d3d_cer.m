function [r phi z ti v sig_ti t] = get_d3d_cer(shot_num,mds_obj,cer_type)
%GET_D3D_CER(shotnum,mds_obj) Retreives CER data for a given shot number
%   Detailed explanation goes here

geom_file=[];
if (shot_num <= 105000)
elseif (shot_num <= 108000)
elseif (shot_num <= 111600)
elseif (shot_num <= 116500)
elseif (shot_num <= 123499)
elseif (shot_num <= 134999)
    geom_file='/Volumes/slazerso/Sims/DIIID/probes/cer_123500_134999.txt';
elseif (shot_num <= 141799)
elseif (shot_num <= 142999)
    geom_file='/Volumes/slazerso/Sims/DIIID/probes/cer_141800_142999.txt';
    %geom_file='/Users/slazerso/Desktop/cer_141800_142999.txt';
elseif (shot_num <= 1477699)
    geom_file='/Volumes/slazerso/Sims/DIIID/probes/cer_143000_147799.txt';
else
end

fid=fopen(geom_file,'r');
header=fgetl(fid);
geom=fscanf(fid,'%*s %f %f/%f %f %f %f %f %f %f',[9 inf]);
fseek(fid,0,-1);
header=fgetl(fid);
temp_names=fscanf(fid,'%s %*f %*f/%*f %*f %*f %*f %*f %*f %*f');
fclose(fid);

% Sort names
i=1; i1=1; i2=2;
while (i2 <= length(temp_names))
    switch char(temp_names(i1))
        case {'R','V','T','B','M'}
            if isempty(str2num(char(temp_names(i2))))
                names{i}=char(temp_names(i1:i2-1))';
                i=i+1;
                i1=i2;
                i2=i1+1;
            elseif (i2 == length(temp_names))
                names{i}=char(temp_names(i1:i2))';
                i1=i2;
                i2=i1+1;
            else
                i2=i2+1;
            end
    end
end

r=geom(1,:)./100;
% CER coordinates are engineering coordinates
phi=pi*(360-geom(9,:))./180; 
z=(geom(4,:)+geom(5,:))./200;

switch upper(cer_type)
    case 'USER'
        type_str='CERAUTO';
    case 'QUICK'
        type_str='CERQUICK';
    case 'FIT'
        type_str='CERFIT';
end


ti=[];
r=[];
z=[];
phi=[];
sig_ti=[];
view_type='TANGENTIAL';
for i=1:31
    cer_r_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':R'];
    cer_z_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':Z'];
    cer_phi_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':VIEW_PHI'];
    cer_t_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TIME'];
    cer_ti_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TEMP'];
    cer_rot_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':ROTC'];
    cer_sig_ti_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TEMP_ERR'];
    try
        t{i}=mds_obj.Value(cer_t_str).Double;
        ti{i}=mds_obj.Value(cer_ti_str).Double;
        sig_ti{i}=mds_obj.Value(cer_sig_ti_str).Double;
        temp=mds_obj.Value(cer_r_str).Double;
        r(i)=temp(1);
        temp=mds_obj.Value(cer_z_str).Double;
        z(i)=temp(1);
        temp=mds_obj.Value(cer_phi_str).Double;
        phi(i)=temp(1);
        v{i}=mds_obj.Value(cer_rot_str).Double;
    catch
        v{i} = 0.0;
        t{i} =0.0;
        ti{i}=[];
        sig_ti{i}=1.0E30;
    end
end
view_type='VERTICAL';
for i=1:24
    cer_r_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':R'];
    cer_z_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':Z'];
    cer_phi_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':VIEW_PHI'];
    cer_t_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TIME'];
    cer_ti_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TEMP'];
    cer_rot_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':ROTC'];
    cer_sig_ti_str=['\D3D::TOP.IONS.CER.' type_str '.' view_type '.CHANNEL' num2str(i,'%2.2d') ':TEMP_ERR'];
    try
        t{i+31}=mds_obj.Value(cer_t_str).Double;
        ti{i+31}=mds_obj.Value(cer_ti_str).Double;
        sig_ti{i+31}=mds_obj.Value(cer_sig_ti_str).Double;
        temp=mds_obj.Value(cer_r_str).Double;
        r(i+31)=temp(1);
        temp=mds_obj.Value(cer_z_str).Double;
        z(i+31)=temp(1);
        temp=mds_obj.Value(cer_phi_str).Double;
        phi(i+31)=temp(1);
        v{i+31}=mds_obj.Value(cer_rot_str).Double;
    catch
        v{i+31} = 0.0;
        t{i+31} =0.0;
        ti{i+31}=[];
        sig_ti{i+31}=1.0E30;
    end
end

phi=pi*(360-phi)./180;


end

