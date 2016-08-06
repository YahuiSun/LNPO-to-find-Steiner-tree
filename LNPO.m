clear all; clc;


load('Example_data'); % load data
savename=['savename']; % the name of the saved data

optimal=1e8; % the optimal solution; use this value if the optimal solution is already known
FElimit=1e6; % the upper bound of the fitness evaluation times
%

% initialization
fe_PO=0; % function evaluation times
fit_PO=zeros(5e7,1); % the best obtained solution
N=node_num; % the number of vertices in the initial graph
success=0; % terminate the program when success=1
Maxw=max(node_weight); % the biggest node weight
terminal_num=sum(Terminal); % the number of terminals in the initial graph
foodpoint=Terminal; % the terminal set; foodpoint(i)=1 if vertex i is a terminal, or foodpoint(i)=0
foodnum=sum(foodpoint); % the number of terminals in the initial graph
edgevalue=ones(N,N); % the parameter alpha(i,j) used in the conductivity update equation
food=zeros(foodnum,1);
r=0;
for i=1:N
    if foodpoint(i)==1
        r=r+1;
        food(r)=i;  %  Terminal r is Vertex food(r)
    end
end
%


% algorithm parameter
I=1e0; % the flux flowing into each source node
kk=50;  % the inner iteration times
cutoff=1e-3;   % cutoff value of edge conductivity
alpha=0.222; % parameter in the evoluation function of conductivity
sigma=1; % parameter in the evoluation function of conductivity
%


Rset=set; % the connectivity matrix of the initial graph
RL=L; % the edge length matrix of the initial graph
Q=zeros(N); % the flux matrix
%
whole=0;
while success==0
    whole=whole+1;   % outer iteration time
    set=Rset; % the connectivity matrix of the initial graph
    L=RL; % the edge length matrix of the initial graph
    D=set; % the initial edge connectivity matrix
    
    exist=ones(N,1); % 1 means a vertex is in the subnetwork, 0 means not
    
    % iteration of D and P
    for k=1:kk  % inner iteration time
        NUM=zeros(N,1);
        % whether a vertex exist
        [degree]=Function_OutputDegree(set,N); % degree of vertex
        for i=1:N
            if degree(i)==0
                exist(i)=0; % 1 means a vertex is in the subnetwork, 0 means not
            end
        end
        vertex=sum(exist);  % number of existing vertex
        
        % unequal possibility 3: choose the sink node and the source nodes
        adjacentlength=zeros(foodnum,1);
        for i=1:foodnum
            for j=1:N
                if set(food(i),j)==1
                    adjacentlength(i)=adjacentlength(i)+L(food(i),j);
                end
            end
        end
        [B,ind]=sort(adjacentlength);   %%%  ascending order,  B is the ordered vector, ind is the intex
        TotalAdjacent=0;
        for i=1:foodnum
            TotalAdjacent=TotalAdjacent+B(i);
        end
        sumad=zeros(foodnum,1);
        for i=1:foodnum
            for j=1:i
                sumad(i)=sumad(i)+B(j);
            end
        end
        random=ceil(rand(1)*TotalAdjacent);
        if random<=sumad(1)
            sink=food(ind(foodnum)); luck=ind(foodnum);
        else
            for i=2:foodnum
                if random>sumad(i-1) && random<=sumad(i)
                    sink=food(ind(foodnum+1-i)); luck=ind(foodnum+1-i);
                end
            end
        end
        
        %  calculate the pressure
        A=zeros(vertex-1);
        g=0;
        for w1=1:N
            if exist(w1)==1 % only calculate the existing vertice
                if w1==sink % neglect the vertex if it's the sink point
                else
                    g=g+1; % row number of A
                    h=0;
                    NUM(w1)=g; %  record the row number's relationship with vertex number, which will be used when using pressure
                    for w3=1:N  % the vertex corresponds to the columes of A
                        if exist(w3)==1  % only calculate the existing point
                            if w3==sink(1)  % neglect the vertex if it's the sink point
                            else
                                h=h+1; % colume number of A
                                if w3==w1 %  the row vertex is the colume vertex
                                    for i=1:N
                                        if set(w1,i)==1
                                            A(g,h)=A(g,h)-D(w1,i)/L(w1,i);
                                        end
                                    end
                                else % the row vertex is not the colume vertex
                                    if set(w1,w3)==1
                                        A(g,h)=D(w1,w3)/L(w1,w3);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        X=zeros(vertex-1,1);
        for i=1:foodnum
            if i==luck
            else
                X(NUM(food(i)))=-I;
            end
        end
        
        P=A\X; % the pressures
        
        % calculate the flux, update the conductivity and cut edges
        for i=1:N
            for j=1:N
                if set(i,j)==1
                    if i==sink
                        Q(i,j)=D(i,j)*(0-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if j==sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-0)/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if i~=sink && j~=sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                end
            end
        end
        
        % move out disconnected part
        [set]=function_prim(set,foodpoint,N);
        %
        
        [TotalPrice]=Function_netPrize(set,L,N,node_weight); % the net-weight of the subnetwork
        
        fe_PO=fe_PO+1; % function evaluation times
        fit_PO(fe_PO)=TotalPrice; % the obtained solution after this step
        
        [bestsofar]=max(fit_PO(1:fe_PO)); % the best obtained solution
        
        if TotalPrice==bestsofar
            MaxSet=set;
            for i=1:N
                for j=1:N
                    if set(i,j)==1
                        edgevalue(i,j)=1.2; % update alpha(i,j) in the conductivity update equation
                    else
                        edgevalue(i,j)=0.8;
                    end
                end
            end
        end
        
        if fe_PO>FElimit % terminate the program when the fitness evaluation times reaches its upper bound.
            success=2;
            break;
        end
        
    end  % inner iteration

    if bestsofar==optimal % terminate the program when the optimal solution is found.
        success=1;
    end
    
end   % the outer iteration

[MaxSet]=Function_SMTpost(N,MaxSet,L,Rset); % find the MST of the best obtained solution

save([savename]);  %  save data

