//1 seule boucle 
function [C]= matmat1b(A,B)
    [n,m]= size(A);
    [m,r]=size(B);
    C = zeros(n,r)
    for (i=1:n)
        for (j=1:n)
                C(i,:) = A(i,:)*B+C(i,:);
                       end
               end    

endfunction
