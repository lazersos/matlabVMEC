function spec_ideal_select(vmec_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gam=(1+sqrt(5))*.5;
flux_spl=pchip(vmec_data.iotaf,vmec_data.phi./vmec_data.phi(vmec_data.ns));
subplot(2,1,1);
plot(vmec_data.iotaf,vmec_data.presf,'k');
iota_max = max(vmec_data.iotaf);
axis tight;
subplot(2,1,2);
plot(vmec_data.iotaf,[abs(diff(vmec_data.presf)) 0.0],'k');
axis tight;
colors={'r','g','b','m','k','c','y','g'};
outarr=[];
nobel_save=[];
for i=8:-1:1
    subplot(2,1,1);
    hold on;
    [farey,p,q]=farey_tree(i);
    if (iota_max > 1)
        p = [p p+1];
        q = [q q];
        farey = p./q;
    end
    len=length(farey);
    plot([farey; farey],repmat(ylim',[1 len]),colors{i},'LineWidth',2);
    nobel = [];
    nobel(1:len-1) = (p(1:len-1) + gam.*p(2:len)) ./ (q(1:len-1) + gam.*q(2:len));
    nobel_save=[nobel_save nobel];
    plot([nobel; nobel],repmat(ylim',[1 len-1]),'Color',colors{i},'LineWidth',0.5);
    xlim([min(vmec_data.iotaf) max(vmec_data.iotaf)]);
    subplot(2,1,2);
    hold on;
    plot([farey; farey],repmat(ylim',[1 len]),colors{i},'LineWidth',2);
    plot([nobel; nobel],repmat(ylim',[1 len-1]),'Color',colors{i},'LineWidth',0.5);
    xlim([min(vmec_data.iotaf) max(vmec_data.iotaf)]);
    %subplot(2,2,[2 4]);
    %textarr={};
    %for j=1:len-1
    %    textarr=[textarr;  [num2str(p(j),'%3i') '+ k*' num2str(p(j+1),'%3i')]];
    %    textarr=[textarr;  ['------------ = ' num2str(nobel(j),'%8.2E')]];
    %    textarr=[textarr;  [num2str(q(j),'%3i') '+ k*' num2str(q(j+1),'%3i')]];
    %    textarr=[textarr;  ' '];
    %    outarr =[outarr; nobel(j) ppval(flux_spl,nobel(j))];
    %end
    %text(0.25*(i-1),0.5,textarr,'Color',colors{i});
    %axis tight;
end
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
subplot(2,1,1);
while but == 1
    [xi,yi,but] = ginput(1);
    iota_target=xi;
    [val, dex]=min(abs(nobel_save - iota_target));
    flux_val = ppval(flux_spl,nobel_save(dex));
    disp(['iota = ' num2str(nobel_save(dex),'%20.10E') '   flux = ' num2str(flux_val,'%20.10E') ]);
end

end

