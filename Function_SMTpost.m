% find SMT of the inputed graph

function [mst_set]=Function_SMTpost(N,prune_set,L,set)

exist=zeros(N,1); % 1 means this vertex will be in the final solution
for i=1:N
    for j=1:N
        if prune_set(i,j)==1
            exist(i)=1;
            break;
        end
    end
end

a=0;
for i=1:N
    if exist(i)==1
        a=a+1;
        ind(a)=i; % vertex a in the MST processing graph is vertex i in the inital graph
    end
end

NN=sum(exist); % the number of vertices in the SMT

newL=zeros(NN);
for i=1:NN
    for j=1:i
        newL(i,j)=L(ind(i),ind(j))*set(ind(i),ind(j));
    end
end
if sum(sum(newL))==0
    mst_set=zeros(N);
else
    newL=sparse(newL);
    [tree, pred]=graphminspantree(newL);
    [a1,a2,s]=find(tree);
    mst_set=zeros(N);
    im=size(a1); im=im(1);
    for i=1:im
        x1=ind(a1(i)); x2=ind(a2(i));
        mst_set(x1,x2)=1; mst_set(x2,x1)=1;
    end
end