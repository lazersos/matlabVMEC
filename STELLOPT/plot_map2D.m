function [ output_args ] = plot_map2D( map,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dex=[];
for i=1:length(map.target_name)
    if strcmp(map.target_name{i},type)
        dex=[dex i];
    end
end

data=squeeze(sum(map.fval(dex,:,:))).^2;
max_val=max(max(data));
exc_dex=find(data == max_val);
min_val=min(min(data));
data2=data;
data2(exc_dex)=nan;
pixplot(data2);



end

