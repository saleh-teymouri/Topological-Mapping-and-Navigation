clear all
clc

X = {'0%'; '1%'; '3%'; '5%'; '7%'};
Y1 = [88825, 74446, 69098, 60987, 77678];             %L457
Y2 = [174235, 139423, 111036, 109364, 162366];        %BC

plot(Y1);
hold on
plot(Y2);

ax = gca;
ax.YAxis.Exponent = 0;
set(gca,'xtick',[1:5],'xticklabel',X)

ylim([20e3 20e4])            %L457: ylim([0 9.75e4])         BC: ylim([0 18.75e4])
xlim([0.375 5.625])

grid on
text(1:length(Y1),Y1,num2str(Y1'),'vert','bottom','horiz','center');
text(1:length(Y2),Y2,num2str(Y2'),'vert','bottom','horiz','center');
legend('Fisrt Environment', 'Second Environment')
xlabel('HIW contribution in constructing 2-simplices')
ylabel('Observations')