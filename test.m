clc;
clear ALL;
addpath('src');

x1 = [1 5 5 1 1];
y1 = [1 1 5 5 1];
%figure;
plot(x1,y1);
hold on;
x2 = [0 7];
y2 = [0 6];
plot(x2,y2);
[xi, yi,ii] = polyxpoly(x1,y1,x2,y2,'unique') ;

plot(xi,yi,'*r');
