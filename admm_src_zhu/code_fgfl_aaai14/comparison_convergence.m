% if maxIter is missing, this function reports 
% the function values when run fGFL for a large 
% number of iterations; 
% if maxIter is a vector, this function runs fGFL
% for the iterations specified in maxIter and 
% report their runtimes. 
function comparison_convergence(num_groups,lambda_graph,...
    gamma,maxIter)
%% load code
if ~isdeployed
  addpath('./GFL');
end

%% check input
if isdeployed
     num_groups = str2num(['int64(' num_groups ')']);
     lambda_graph = str2double(lambda_graph);
     gamma = str2double(gamma);
     if nargin == 4 
         maxIter = str2num(['int64(' maxIter ')']);
     end
end
lambda_sparse = lambda_graph * gamma;
p = num_groups * 11;
if nargin == 4
    disp('Running fGFL for a large number of iterations.');
    max_iter = maxIter;
elseif nargin == 3
    disp('Running fGFL for specified number of iterations.');
    max_iter = dlmread(strcat('data/fGFL_iter_seq_p',int2str(p),...
        '_lam1_',num2str(lambda_graph),...
        '_lam2_',num2str(lambda_sparse),'.txt'));
else
    error('please supply num_groups,lambda_graph and gamma');
end

%% load data
% disp(p); disp(lambda_graph); disp(lambda_sparse);
X = dlmread(strcat('data/data_matrix_p',int2str(p),'.txt'));
y = dlmread(strcat('data/response_p',int2str(p),'.txt'));
in_node = dlmread(strcat('data/in_node_p',int2str(p),'.txt'));
out_node = dlmread(strcat('data/out_node_p',int2str(p),'.txt'));
%idx = dlmread(strcat('data/D_idx_p',int2str(p),'.txt'));
%jdx = dlmread(strcat('data/D_jdx_p',int2str(p),'.txt'));
%val = dlmread(strcat('data/D_val_p',int2str(p),'.txt'));
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
ratio = lambda_graph / lambda_sparse;
opts.tol = 1e-30;   % tolerance.
if (length(max_iter) == 1)
    disp('Running fGFL for a large number of iterations.');
    opts.maxIter = max_iter;
    tic
    [beta, funcVal] = fast_gfl(X, y, Graph, ...
        lambda_sparse, ratio, opts);
    t = toc;
    dlmwrite(strcat('../software/tmp/convergence_timing_p',int2str(p),'.txt'),t);
    dlmwrite(strcat('../software/tmp/funVal_p',int2str(p),...
        '_lambda_graph_',num2str(lambda_graph),...
        '_gamma_',num2str(gamma),'.txt'),funcVal,'precision',20);
    dlmwrite(strcat('../software/tmp/beta_glf_p',int2str(p),'.txt'),beta,'precision',20);
else
    disp('Running fGFL for specified number of iterations.');
    t = zeros(1,length(max_iter));
    for i = 1:length(max_iter)
        opts.maxIter = max_iter(i);
        tic
        fast_gfl(X, y, Graph,lambda_sparse, ratio, opts);
        t(i) = toc;
    end
    dlmwrite(strcat('../software/tmp/convergence_timing_p',int2str(p),...
            '_lam1_',num2str(lambda_graph),...
        '_lam2_',num2str(lambda_sparse),'.txt'),t);
end

