clear; close all; clc;

theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
load data1
[dmodel, perf] = dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb);

addpath('../lhs')

% X = lhsdesign_modified(1000, [0 0], [100 100], 10);

[YX, MSE] = predictor(X, dmodel);