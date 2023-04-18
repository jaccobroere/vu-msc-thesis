function [nEdge Edge_weight Edge_in Edge_out] = calc_w_2d(s1,s2)


p = s2*s1;
IND = zeros(p,1);
 for i = 1:p
     IND(i) = i;
 end    
[R C] = ind2sub([s1,s2],IND);

Edge_weight = [];
Edge_in = [];
Edge_out = [];
nEdge = 0;

for i = 1:p
    pI_row = R(i);
    pI_col = C(i);
    
    
    pJ_row = pI_row;
    pJ_col = pI_col - 1;
    if(validpoint(pJ_row,pJ_col))
        j = s1 * (pJ_col-1) + pJ_row;
        tw = 1;
        Edge_in = cat(1,Edge_in,i);
        Edge_out = cat(1,Edge_out,j);
        Edge_weight = cat(1,Edge_weight,tw);
        nEdge = nEdge + 1;
    end

    pJ_row = pI_row - 1;
    pJ_col = pI_col;
    if(validpoint(pJ_row,pJ_col))
        j = s1 * (pJ_col-1) + pJ_row;
        tw = 1;
        Edge_in = cat(1,Edge_in,i);
        Edge_out = cat(1,Edge_out,j);
        Edge_weight = cat(1,Edge_weight,tw);
        nEdge = nEdge + 1;
    end    
    
    pJ_row = pI_row + 1;
    pJ_col = pI_col;
    if(validpoint(pJ_row,pJ_col))
        j = s1 * (pJ_col-1) + pJ_row;
        tw = 1;
        Edge_in = cat(1,Edge_in,i);
        Edge_out = cat(1,Edge_out,j);
        Edge_weight = cat(1,Edge_weight,tw);
        nEdge = nEdge + 1;
    end
      
    pJ_row = pI_row;
    pJ_col = pI_col + 1;
    if(validpoint(pJ_row,pJ_col))
        j = s1 * (pJ_col-1) + pJ_row;
        tw = 1;
        Edge_in = cat(1,Edge_in,i);
        Edge_out = cat(1,Edge_out,j);
        Edge_weight = cat(1,Edge_weight,tw);  
        nEdge = nEdge + 1;
    end 

end

    % private function

     function v = validpoint(pJrow, pJcol)

        v = pJrow>0&&pJrow<=s1&&pJcol>0&&pJcol<=s2;
     end
end




