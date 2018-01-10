function [ output_args ] = plot_terpsichore( data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nlarge=5;

subplot(1,2,1);
sxi = (0:data.nvi)./data.ni;
[sortedX,sortingIndices]=sort(data.ximax,'descend');
leg_text={};
for i=1:nlarge
    hold on;
    dex = sortingIndices(i);
    plot(sxi,data.xi(dex,:)','LineWidth',2);
    leg_text{i} = ['m/n = ' num2str(data.ms(dex),'%2i') '/' num2str(data.ns(dex),'%2i')];
end
set(gca,'FontSize',24);
axis tight;
xlabel('Normalized Flux');
ylabel('\xi');
title('Largest \xi Modes');
legend(leg_text,'FontSize',12);
set(gcf,'Position',[1 -100 1024 768],'Color','white')
sxi = (0:data.nvi-1)./data.ni;
subplot(1,2,2);
[sortedX,sortingIndices]=sort(data.etamax,'descend');
leg_text={};
for i=1:nlarge
    hold on;
    dex = sortingIndices(i);
    plot(sxi,data.eta(dex,:)','LineWidth',2);
    leg_text{i} = ['m/n = ' num2str(data.ms(dex),'%2i') '/' num2str(data.ns(dex),'%2i')];
end
set(gca,'FontSize',24);
axis tight;
xlabel('Normalized Flux');
ylabel('\eta');
title('Largest \eta Modes');
legend(leg_text,'FontSize',12);

end

