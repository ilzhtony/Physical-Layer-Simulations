close all; clear; clc;

N = 100000;
NP = 1000;
x = mygaussrnd(1, N);

[p, s] = ksdensity(x, 'npoints', NP);
plot(s, p);
