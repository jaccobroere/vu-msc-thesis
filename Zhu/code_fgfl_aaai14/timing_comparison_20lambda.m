% function input lambda_graph, gamma = 10, and
% optimal function values at lambda_graph
function timing_comparison_20lambda(num_groups,gamma)
%% load code
if ~isdeployed
    addpath('./GFL');
end


%% check input
if isdeployed
    num_groups = str2num(['int64(' num_groups ')']);
    gamma = str2num(['int64(' gamma ')']);
end
%% load data
p = 11 * num_groups;
X = dlmread(strcat('data/data_matrix_p',int2str(p),'.txt'));
y = dlmread(strcat('data/response_p',int2str(p),'.txt'));
in_node = dlmread(strcat('data/in_node_p',int2str(p),'.txt'));
out_node = dlmread(strcat('data/out_node_p',int2str(p),'.txt'));
lambda = dlmread(strcat('data/lambdas_single_p',int2str(p),...
    '_gamma_',int2str(gamma),'.txt'));
opt_funVals = dlmread(strcat('data/opt_funVals_20lambdas_p',...
    int2str(p),'_gamma_',int2str(gamma),'.txt'));
subopt = dlmread('data/subopt.txt');
m = size(in_node,1);

Graph{1} = 2*m;
Graph{2} = ones(2*m,1);
Graph{3} = [in_node;out_node];
Graph{4} = [out_node;in_node];


% Graph{1} = m;
% Graph{2} = ones(m,1);
% Graph{3} = in_node;
% Graph{4} = out_node;


%% run code for a large # of iterations %%
opts.tol = 1e-30;   % tolerance.
opts.maxIter = 5e3;
num_grid = length(lambda);

% store iterations
% iter_fGFL = zeros(0,length(lambda),length(subopt));
% for i = 1:num_grid
%     disp(lambda(i));
%     [beta, funcVal] = fast_gfl_funVal_stopping(X, y, Graph, ...
%         lambda(i)*gamma, 1.0/gamma, opt_funVals(i),EPS, opts);
%     subopt_tmp = (funcVal - opt_funVals(i)) / opt_funVals(i);
%     for j = 1:length(subopt)
%         iter_fGFL(i,j) = find(subopt_tmp <  subopt(j),1,'first');
%     end
% end

runtime = zeros(num_grid,length(subopt));
for i = 1:num_grid
    disp(lambda(i));
    for j = 1:length(subopt)
        tic
        [beta, funcVal] = fast_gfl_funVal_stopping(X, y, Graph, ...
            lambda(i)*double(gamma), 1.0/double(gamma),opt_funVals(i),subopt(j),opts);
        t = toc;
        runtime(i,j) = t;
    end
end

%% save runtime
dlmwrite(strcat('../software/tmp/fGFL_20lambda_timing_p',int2str(p),...
    '_gamma_',int2str(gamma),'.txt'), runtime);

