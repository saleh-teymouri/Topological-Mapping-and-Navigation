%data = [168807.69, 137266.2, 126609.41, 116463.67, 135985.83, 179540.28, 201191.72 ,205160.18, 256996.44, 286902.42]; %BC
data = [83839.84, 72713.85, 70304.33, 62986.46, 68613.13, 77708.93, 83833.70, 91367.27, 115450.7, 131065.87]; %L457
bin = 0:0.001:0.01;
histogram('BinEdges',bin,'BinCounts',data)
% Fancy up the graph.
grid on;
xlim([0, 0.01]);
%ylim([0, 350000])
xlabel('Ratio Bins', 'FontSize', 14);
ylabel('Number of Observations', 'FontSize', 14);
title('First Environment', 'FontSize', 14);