function [Acc,x] = Classification(data, stress, exp, ndrivers)

% Prepare features for the LSLD
data = table2array(data);
data = normalize(data)';
data = [ones(1,size(data,2));data];

% Define matrices for each class: Classes'=[1 0 0] if stress=-1 => rest
% Classes'=[0 1 0] if stress=0 => hwy y Classes'=[0 0 1] if stress=1 => city
Classes=zeros(3,size(stress,2));
Classes(1,(stress==-1))=1; %Class 1: rest [1 0 0]
Classes(2,(stress==0))=1; %Class 2: Hwy [0 1 0]
Classes(3,(stress==1))=1; %Class 3: city [0 0 1]

% LOO since there are very few drivers
for m=1:ndrivers
    TrainFeat=data(:,exp~=m);
    TrainClas=Classes(:,exp~=m);

    [x(:,m),PerrorTrain(m)]=GeneticAlgorithm(TrainFeat, TrainClas);

    TrainFeat=data(x(:,m),exp~=m);
    TrainClas=Classes(:,exp~=m);

    TestFeat=data(x(:,m),exp==m);
    TestClas=Classes(:,exp==m);

    M=(TrainFeat*TrainFeat');
    V{m}=TrainClas*TrainFeat'*pinv(M);
    Y(:,(exp==m))=V{m}*TestFeat;
end
    
[a,b]=max(Y);
[c,d]=max(Classes);
Acc=mean(b==d);   % Accuracy

x = x(2:end,:);
