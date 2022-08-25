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
L = [];
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FEATURE EXTRACTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% HR %%%
hr = [];
for m=1:length(drivers)
    hr=[hr reshape(mean(HR{drivers(m)}(:,inip(m):endp(m)),1),N,L(m)/N)];
end

FeatHR = FeaturesHR(hr,fs,N);





