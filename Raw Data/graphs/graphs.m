clear all
clc

X = 0:1000:12000;
Y = [7.97923, 4.19968, 2.07385, 1.30597, 0.510349, 0.361925, 0.215054, 0.195012, 0.159537, 0.159633, 0.182039, 0.160982, 0.0877289];

Yi = smooth(Y);
plot(X,Yi,'b');

xlabel('time step')
ylabel('growth rate')

grid on