function [w omega] = eff_general_flsa(n, z, lambda1, lambda2, nE, E_in, E_out, E_w)
% this is fast flsa

if lambda1==0 && lambda2==0
    w = z;
    omega = 0;
elseif lambda1 ~=0 && lambda2 == 0
    w = sign(z).* max(abs(z)-lambda1,0);
else 
    [w omega]  = mexEFLSA(n, z, lambda1, lambda2, nE, E_in, E_out, E_w, 0);
    clear mexEFLSA
end
    

    
    
    
    
