%This function will assign the derived values for T_l(i,j), T_g(i,j), and
%Q(i,j) into their corresponding location in their respective matrix
%without changing the inlet boundary conditions.
function [T_l,T_g,Q]=assignment(T_l,T_g,Q,X,T_l_in,T_g_in)
k=1;
for i=1:size(T_l,1)   
    for j=1:size(T_l,2)
        if T_l(i,j)==T_l_in
        else
            T_l(i,j)=X(k);
            k=k+1;
        end            
    end
end
for i=1:size(T_g,1) 
    for j=1:size(T_g,2)
        if T_g(i,j)==T_g_in
        else
            T_g(i,j)=X(k);
            k=k+1;
        end            
    end
end
for i=1:size(Q,1) 
    for j=1:size(Q,2)
         Q(i,j)=X(k);
         k=k+1;
         if k>numel(X)
             break
         end
    end
end
     
