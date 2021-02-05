clear all
clc

A = load('LPA_Minimize.txt');
X = A(:,1);
Y = A(:,2);

bar(X,Y,0.5)
text(X, Y, num2str(Y,'%0.0f'),'HorizontalAlignment','center','VerticalAlignment','bottom')

xlabel('\DeltaS Values')
ylabel('Number of Landmarks')

xlim([-1 19])
ylim([600 2600])

grid on