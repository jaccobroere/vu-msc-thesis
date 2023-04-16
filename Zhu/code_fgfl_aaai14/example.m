%% This is an example code to use fast_gfl.m
clear all;
close;

%% load 
addpath('./GFL');  % add GFR

% generate a 2d-grid signal as beta0
d1 = 10;
d2 = 10;
d = d1*d2;
beta0_2d = zeros(d1,d2);
randNum=1;
beta0_2d(d1/2-d1/5:d1/2+1+d1/5,d1/2-d1/5:d1/2+1) = 0.8;
beta0_2d(d1/2-d1/5:d1/2+1+d1/5,d1/2:d1/2+1+d1/5) = 0.2;
beta0 = reshape(beta0_2d,d,1);
figure,
subplot(1,2,1);
imshow(beta0_2d);
title('beta0');

% generate the correspoding graph-structure
[nE E_weight E_in E_out] = calc_w_2d(d1,d2);
Graph{1} = nE;                  
Graph{2} = E_weight;            
Graph{3} = E_in;                  
Graph{4} = E_out;

% generate samples
N = 80;
rng(randNum);
X=randn(N,d);       % the data matrix
rng(randNum);
noise=randn(N,1)*0.05;           % 0.05 is a ratio to control noise level
y = X*beta0 + noise;     % the response

%% generalized fused lasso 
rho1 = 0.1;
rho2 = 0.5;
opts.tol = 10^-5;   % tolerance.
opts.maxIter = 1000; % maximum iteration number of optimization.

tic
[beta, funcVal] = fast_gfl(X, y, Graph, rho1, rho2, opts);
t = toc;
e = norm(beta0-beta);
beta_2d = reshape(beta,d1,d2);
subplot(1,2,2);
imshow(beta_2d);
title('recovered beta by GFL');





