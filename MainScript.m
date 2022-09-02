clear
close all

%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

Nsamples=[... %ECG    %EMG    %Foot   %Hand   %HR     %Marker %Resp
              32      128     2       2       1       0       1; ...%drive01
              32      0       2       0       1       1       2; ... %drive02
              16      0       1       1       0       0       1; ... %drive03
              32      0       2       2       1       1       2; ... %drive04
              32      1       2       2       1       1       2; ... %drive05
              32      1       2       2       1       1       2; ... %drive06
              32      1       2       2       1       1       2; ... %drive07
              32      1       2       2       1       1       2; ... %drive08
              32      1       2       2       1       1       2; ... %drive09
              32      1       2       2       1       1       2; ... %drive10
              32      1       2       2       1       1       2; ... %drive11
              23      1       2       2       1       1       2; ... %drive12
              32      1       2       0       1       1       2; ... %drive13
              32      1       2       2       0       1       2; ... %drive14
              31      1       2       2       1       1       2; ... %drive15
              32      1       2       2       1       1       2; ... %drive16
              32      1       2       2       1       1       2; ... %drive17
              32      1       2       2       1       1       2; ... %drive18
              ];

for n=1:18
    P=fopen(sprintf('drive%02.0f.dat',n));  
    DATA{n}=fread(P,[sum(Nsamples(n,:)) 200000],'int16');
    fclose(P);
    ECG{n}=DATA{n}(1:Nsamples(n,1),:);
    EMG{n}=DATA{n}(sum(Nsamples(n,1))+(1:Nsamples(n,2)),:);
    FOOT{n}=DATA{n}(sum(Nsamples(n,1:2))+(1:Nsamples(n,3)),:);
    HAND{n}=DATA{n}(sum(Nsamples(n,1:3))+(1:Nsamples(n,4)),:);
    HR{n}=DATA{n}(sum(Nsamples(n,1:4))+(1:Nsamples(n,5)),:);
    MARKER{n}=DATA{n}(sum(Nsamples(n,1:5))+(1:Nsamples(n,6)),:);
    RESP{n}=DATA{n}(sum(Nsamples(n,1:6))+(1:Nsamples(n,7)),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EXPLORATORY ANALYSIS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep those subjects with all data for all sensors
find(cellfun(@isempty,MARKER)) %drivers 1 and 3. They are discarded since we need marker information
find(cellfun(@isempty,ECG)) 
find(cellfun(@isempty,EMG)) %drivers 2, 3 and 4 have no EMG data 
find(cellfun(@isempty,FOOT))
find(cellfun(@isempty,HAND)) %drivers 2 and 13 have no GSR-HAND data
find(cellfun(@isempty,HR)) %drivers 3 and 14 have no HR data 
find(cellfun(@isempty,RESP)) 

% After visual inspection of markers for remaining subjects, drivers 9, 14,
% 16, 17 and 18 have no appropriate marker signal so they have to be discarded
drivers=[5 6 7 8 10 11 12 15];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EXTRACTION OF MARKERS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
fs=15.5;  % Minimum common sample frequency among signals
N=fs*60; % samples in 1 minute

inip = []; endp = []; %initial and end points
L = []; stress = []; exp = [];
for m = 1:length(drivers)
    marker=MARKER{drivers(m)};
    
    %search for 8 marker values: rest, city, highway, city, highway, city, rest
    ind1=zeros(size(marker));
    ind1(5:end-4)=(marker(5:end-4)-(marker(1:end-8)+marker(9:end))/2);
    % Search for 8 max points, we leave the signal between 0 and 1
    ind2=zeros(size(marker));
    for n=1:8
        [a,b]=max(ind1);
        ind2(b)=1;
        ind1(max(b-4000,1):min(b+4000,length(ind1)))=0;
    end
    
    %Generate labels -1:rest   0:highway   1:city
    label=[-1    1    0    1    0    1    -1];
    ind=find(ind2);
    
    for n=1:7
        ind2(ind(n):ind(n+1)-1)=label(n);
    end

    %divide in slots of 1000 samples
    inip(m)=ind(1); % initial point
    temp=ind(8); % end point
    L(m)=floor((temp-inip(m))/N)*N;
    endp(m)=inip(m)+L(m)-1;

    % Stress classes
    stress=[stress reshape(ind2(1,inip(m):endp(m)),N,L(m)/N)];
    exp=[exp m*ones(1,L(m)/N)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FEATURE EXTRACTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ECG %%%
ecg = []; FeatECG=[];
for m=1:length(drivers)
    fecg=reshape(ECG{drivers(m)}(:,inip(m):endp(m)),N*size(ECG{drivers(m)},1),L(m)/N);
    FeatECG = [FeatECG; FeaturesECG(fecg,size(ECG{drivers(m)},1),fs,N)];
    ecg=[ecg reshape(mean(ECG{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end


%%% HR %%%
hr = [];
for m=1:length(drivers)
    hr=[hr reshape(mean(HR{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end
FeatHR = FeaturesHR(hr,fs,N, FeatECG.meanrr');


%%% Hand %%%
FeatHAND=[]; hand=[];
for m=1:length(drivers)
    fhand=reshape(HAND{drivers(m)}(:,inip(m):endp(m)),N*size(HAND{drivers(m)},1),L(m)/N);
    FeatHAND = [FeatHAND; FeaturesGSR(fhand,size(HAND{drivers(m)},1),fs,N,'hand')];
    hand=[hand reshape(mean(HAND{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end
FeatHAND = [FeatHAND FeaturesGSR(hand,size(HAND{drivers(m)},1), fs,N, 'handstats')];


%%% Foot %%%
FeatFOOT=[]; foot=[];
for m=1:length(drivers)
    ffoot=reshape(FOOT{drivers(m)}(:,inip(m):endp(m)),N*size(FOOT{drivers(m)},1),L(m)/N);
    FeatFOOT = [FeatFOOT; FeaturesGSR(ffoot,size(FOOT{drivers(m)},1),fs,N,'foot')];
    foot=[foot reshape(mean(FOOT{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end
FeatFOOT = [FeatFOOT FeaturesGSR(foot,size(FOOT{drivers(m)},1), fs,N, 'footstats')];


%%% EMG %%%
emg = []; FeatEMG=[];
for m=1:length(drivers)
    femg=reshape(EMG{drivers(m)}(:,inip(m):endp(m)),N*size(EMG{drivers(m)},1),L(m)/N);
    FeatEMG = [FeatEMG; FeaturesEMG(femg, 'emgfeat')];
    emg=[emg reshape(mean(EMG{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end
FeatEMG = [FeatEMG FeaturesEMG(emg, 'emgstats')];


%%% RESPIRATION %%%
resp = []; FeatRESP=[];
for m=1:length(drivers)
    fresp=reshape(RESP{drivers(m)}(:,inip(m):endp(m)),N*size(RESP{drivers(m)},1),L(m)/N);
    FeatRESP = [FeatRESP; FeaturesRESP(fresp, size(RESP{drivers(m)},1), fs, N, 'respfeat')];
    resp=[resp reshape(mean(RESP{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end
FeatRESP = [FeatRESP FeaturesRESP(resp, size(RESP{drivers(m)},1), fs, N, 'respstats')];


%%%%%%%%%%%%%%%%%%%%
%%%%% CLEANING %%%%%
%%%%%%%%%%%%%%%%%%%%

% All features
data=[FeatECG, FeatHR, FeatHAND, FeatFOOT, FeatEMG, FeatRESP];

ind_nan=[]; n_nan=[];
for i=1:size(data,2)
    count=sum(isnan(table2array(data(:,i))));
    if count > 0
        ind_nan=[ind_nan i];
        n_nan=[n_nan count];
    end
end
% There are 20 features with nan values. One features has 605 nan values so
% we discarded. The other ones have 2 or 4 nan values, so we impute their values.
for i=1:length(ind_nan)
    f=table2array(data(:,ind_nan(i)));
    fimpute=median(f(~isnan(f))); %value to impute
    f(isnan(f))=fimpute;
    data(:,ind_nan(i))=array2table(f);
    clear f
end
data(:,ind_nan(1))=[]; %discard feature with 605 nan values

% Classes
stress = mean(stress);
% Round decimals to -1, 0 and 1
stress(stress>0.5)=1;
stress((-0.5<=stress)&(stress<=0.5))=0;
stress(stress<-0.5)=-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CLASSIFICATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

[Accuracy, featselected] = Classification(data, stress, exp, length(drivers));

clc
disp ('***** RESULTS *****')
disp(['Accuracy: ', num2str(Accuracy*100), '%'])
fprintf('\n')

% Features most selected
%bar(sum(featselected'))
count=sum(featselected');
ind = find(count > (length(drivers)/2));
count = count(ind);

display('Features most selected:')
for i = 1:length(ind)
    disp([data.Properties.VariableNames(ind(i)), num2str(count(i)), 'times'])
end
