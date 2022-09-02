%Genetico
function [x,ClasErr]=GeneticAlgorithm(TrainFeat, TrainClas)

% Genetic algorithm that optimized the mse

Npop=1000; % Population
Ngen=100; % Generations
Nmax=10;  % Max number of features
Nfeat=size(TrainFeat,1);

Q=[TrainFeat;ones(1,size(TrainFeat,2))];
A=TrainClas*Q';
B=Q*Q';
if(size(TrainClas,1)==1)
    ClaseTrain=TrainClas;
else
    [~,ClaseTrain]=max(TrainClas);
end

% Initiate population randomly
X=rand(Nfeat,Npop)>0.5;

for g=1:Ngen
    % Population should have a maximum number of features
    aux=X+0.1*rand(size(X));
    aux2=sort(aux,'descend');
    X(aux<(ones(size(aux,1),1)*aux2(Nmax+1,:)))=0;
    
    % Search for identical individuals and modify them
    for n=2:Npop
        while sum((sum(X(:,1:n-1)~=(X(:,n)*ones(1,n-1))))==0)
            ind0=find(X(:,n)==0);
            ind1=find(X(:,n)==1);
            X(ind0(ceil(rand*length(ind0))),n)=1;
            X(ind1(ceil(rand*length(ind1))),n)=0;
        end
    end           
                
    % Evaluate each individual of the population
    for n=1:Npop
        ind=[find(X(:,n));Nfeat+1];
    
        V=A(:,ind)*pinv(B(ind,ind));
        Y=V*Q(ind,:);
            
        if(size(TrainClas,1)==1)
            ClaseEst=Y>0.5;
        else
            [~,ClaseEst]=max(Y);
        end
        Etrain(1,n)=mean(ClaseEst~=ClaseTrain);
    end
                
    PERF=mean(Etrain,1);
    
    % Select 10% to survive
    [aux,ind]=sort(PERF,'ascend');
    ind_parents=ind(1:round(0.1*Npop));
    ind_children=ind(round(0.1*Npop)+1:end);
    
    % Cross-over 
    for n=1:length(ind_children)
        a=ind_parents(ceil(rand*length(ind_parents))); % Select one of the parents
        X(:,ind_children(n))=X(:,a);
        a=ind_parents(ceil(rand*length(ind_parents))); % Select the other parent
        b=rand(Nfeat,1)>0.5;   % Select half of the features fot this parent
        X(b,ind_children(n))=X(b,a);
    end
    
    % Mutate randomly 1% of the genes keeping the best one
    PERFbest=aux(1);
    x=X(:,ind_parents(1));
    ind=rand(size(X))<0.01;
    X(ind)=~X(ind);
    X(:,ind_parents(1))=x;
    
end
aux=mean(Etrain,1);
ClasErr=aux(ind_parents(1));