function [beta, funcVal] = fast_gfl2(X, y, Graph, rho1, rho2, opts)
%% A fast algirhtm sovleing the generalized fused lasso model
% paper: AAAI-14: " Efficient Generalized Fused Lasso and 
% its Application to the Diagnosis of Alzheimer’s Disease "
% Bo Xin, Yoshinobu Kawahara, Yizhou Wang and Wen Gao
% implemented by Bo Xin 
% 2014-7-30 Peking University

%%
%$ Description
% this code solves: 
% min_{beta}: 0.5*||y-X*beta||_2^2 + rho1 * ||beta||_1 + rho2 * \Sum (||beta_i-beta_j||_1)
    % input: 
        % X: design matrix   \in R^{N*d}
        % y: response vector \in R^N
        % rho1, rho2 are tuning parameters
        % Graph information: 
        % opts: setting for accelrated proximal methods(FISTA) algorihm
    % output:
        % beta: the variables vector 
        % funcVal: function value of each iteration

% note that this code implemented GFL for "0.5*||y-X*beta||_2^2" loss, for general 
% loss terms, the only difference change is to implement the private
% function: gradVal_eval and funcVal_eval correspodingly.
        
%% Code starts here


if nargin < 5
    error('\n Inputs: X, y, rho1, rho2 and graph information should be specified!\n');
end

% Graph eEdge structure
nE = Graph{1};
E_w = Graph{2};
E_in = Graph{3};
E_out = Graph{4};

% initialize a starting point
beta0 = zeros(size(X,2), 1);
beta_z= beta0;
beta_zold = beta0;

% parameters for accelerated proximal methods (e.g. FISTA)
t = 1;
t_old = 0;

iter = 0;
gamma = 1;
gamma_inc = 2; 

bFlag=0; % this flag tests whether the gradient step only changes a little

funcVal = [];

while iter < opts.maxIter
    
    alpha = (t_old - 1) /t;    
    beta_s = (1 + alpha) * beta_z - alpha * beta_zold;

    % compute function value and gradients of the search point
    gbetas  = gradVal_eval(beta_s);
    Fs   = funVal_eval(beta_s);
    
    while true  % backstrap
        vg = beta_s - gbetas/gamma;
        beta_zp = eff_general_flsa(length(vg), vg, rho1/gamma, rho2/gamma, nE, E_in, E_out, E_w);
        l1c_beta_zp = calc_omegagss(beta_zp,rho1,rho2,E_in,E_out,E_w);       
        Fzp = funVal_eval  (beta_zp);        
        delta_beta_zp = beta_zp - beta_s;
        r_sum = norm(delta_beta_zp, 'fro')^2;
        Fzp_gamma = Fs + trace(delta_beta_zp' * gbetas) + gamma/2 * norm(delta_beta_zp)^2;

        
        
        if (r_sum <=1e-20)
            bFlag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        
        if (Fzp <= Fzp_gamma)
            break;
        else
            gamma = gamma * gamma_inc;
        end
    end   

    currentfuncVal = Fzp + l1c_beta_zp;    
    funcVal = cat(1, funcVal, currentfuncVal);
  
    beta_zold = beta_z;
    beta_z = beta_zp;
   
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    % test stop condition.
    if iter>=2
        if (abs( funcVal(end) - funcVal(end-1) ) <=...
               opts.tol* funcVal(end-1))
          break;
        end
    end

    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
    
end

beta = beta_zp;


% private functions
   function [grad_beta] = gradVal_eval(beta)
            grad_beta = [];
            grad_beta = cat(2, grad_beta, X'*X*beta - X'*y);
        
   end

%     function [regFuncVal] = regFuncVal_eval(beta,rho1,rho2,E_in,E_out,E_w)
%         regFuncVal = 0.5 * norm (y - X*beta)^2 + ...
%             rho1*(norm(beta,1) + ...
%         rho2*(sum(E_w.*abs(beta(E_in) - beta(E_out)))/2);
%     end

    function [funcVal] = funVal_eval (beta)
        funcVal = 0.5 * norm (y - X*beta)^2; 
    end
    
    function [omega] = calc_omegagss(beta,rho1,rho2,E_in,E_out,E_w)
           if rho1 == 0
               omega = rho1*norm(beta,1);
           else
               w_in = beta(E_in);
               w_out = beta(E_out);
               omega = rho1*norm(beta,1) + rho2*sum(E_w.*abs(w_out - w_in))/2;
           end 
    end

end