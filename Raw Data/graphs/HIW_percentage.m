clear all
clc

X = {'0%'; '1%'; '3%'; '5%'; '7%'};
Y1 = [88825, 74446, 69098, 60987, 77678];             %L457
Y2 = [174235, 139423, 111036, 109364, 162366];        %BC
Y = Y1;
bar(Y);
xticklabels(X);
ax = gca;
ax.YAxis.Exponent = 0;

ylim([0 9.75e4])            %L457: ylim([0 9.75e4])         BC: ylim([0 18.75e4])
xlim([0.375 5.625])

title('Results for the first complex environment')      %L457: first        BC: second
xlabel('HIW contribution in constructing 2-simplices')
ylabel('Observations')

text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center');
grid on;