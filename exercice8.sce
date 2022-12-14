//trois boucles //pourquoi invalid index ici?
function [C]= matmat3b(A,B)
    [n,m]= size(A);
    [m,r]=size(B);
    C = zeros(n,r)
    for (i=1:n)
        for (j=1:r)
            for(k=1:m)
                C(i,j) = A(i,k)*B(k,j)+C(i,j);
            end
        end    
    end
endfunction

