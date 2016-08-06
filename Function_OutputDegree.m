% output the degree of vertex

function [degree]=Function_OutputDegree(set,N)

degree=zeros(N,1);
for i=1:N
    for j=1:N
        if set(i,j)==1
            degree(i)=degree(i)+1;
        end
    end
end