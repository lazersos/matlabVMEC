function plot_iota_rat
%PLOT_IOTA_RAT Makes a plot of n=1..3 rationals

xlim([0 2]);
ylim([0 3]);
hold on;
%n=3
for m=1:30
    plot([1 1].*(3/m),[0 1].*(1/m),'g','LineWidth',2);
    text((3/m),(1/m),['m = ' num2str(m,'%i')]);
end
%n=2
for m=1:20
    plot([1 1].*(2/m),[0 1].*(1/m),'r','LineWidth',2);
    text((2/m),(1/m),['m = ' num2str(m,'%i')]);
end
%n=1
for m=1:10
    plot([1 1].*(1/m),[0 1].*(1/m),'b','LineWidth',2);
    text((1/m),(1/m),['m = ' num2str(m,'%i')]);
end
end

