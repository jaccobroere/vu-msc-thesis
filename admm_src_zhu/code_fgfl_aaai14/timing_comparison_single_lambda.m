%% This is a comparison %%
clear all;
close;

%% load code 
addpath('./GFL');  % add GFR


%% load data
p = 3300;
gamma = 10;
X = dlmread(strcat('data/data_matrix_p',int2str(p),'.txt'));
y = dlmread(strcat('data/response_p',int2str(p),'.txt'));
in_node = dlmread(strcat('data/in_node_p',int2str(p),'.txt'));
out_node = dlmread(strcat('data/out_node_p',int2str(p),'.txt'));
idx = dlmread(strcat('data/D_idx_p',int2str(p),'.txt'));
jdx = dlmread(strcat('data/D_jdx_p',int2str(p),'.txt'));
val = dlmread(strcat('data/D_val_p',int2str(p),'.txt'));
lambda = dlmread(strcat('data/lambdas_single_p',int2str(p),...
    '_gamma_',int2str(gamma),'.txt'));
D = sparse(idx,jdx,val);
m = size(in_node,1);

Graph{1} = 2*m;                  
Graph{2} = ones(2*m,1);            
Graph{3} = [in_node;out_node];                  
Graph{4} = [out_node;in_node];


% Graph{1} = m;                  
% Graph{2} = ones(m,1);            
% Graph{3} = in_node;                  
% Graph{4} = out_node;


%% run code 
opts.tol = 1e-10;   % tolerance.
opts.maxIter = 1000;

num_grid = length(lambda);
runtime = zeros(1,num_grid);
funVals = zeros(1,num_grid);
for i = 1:num_grid
    tic
    [beta, funcVal] = fast_gfl(X, y, Graph, ...
        lambda(i)*gamma, 1.0/gamma, opts);
    t = toc;
    fprintf('funVal is %6.8f.\n',funcVal(end));
    runtime(i) = t;
    funVals(i) = funcVal(end);
end

%% save results 
dlmwrite(strcat('../software/tmp/fGFL_funVals_p',int2str(p),...
    '_gamma_',int2str(gamma),'.txt'), funVals,'precision',20);
dlmwrite(strcat('../software/tmp/fGFL_timing_p',int2str(p),...
    '_gamma_',int2str(gamma),'.txt'), runtime);

