function [netPrice]=Function_netPrize(set,L,N,node_weight) % find the net-weight of the subnetwork


netPrice=0; % the net-weight
for i=1:N-1
    for j=i+1:N
        netPrice=netPrice-set(i,j)*L(i,j);
    end
end
[degree]=Function_OutputDegree(set,N);
for i=1:N
    if degree(i)>0
        netPrice=netPrice+node_weight(i);
    end
end