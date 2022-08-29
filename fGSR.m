function [Sm,Sd,num_resp,tm]=fGSR(s,FS,lim)

% input: gsr signal, number of slots, amplitude limit to restrict max and minimums to detect
% output: impulse magnitudes, durations, number of responses and average
% time of recovery

if(~exist('FS'))
    FS=2;
end

% Amplitude limit. Those below will be discarded
if(~exist('lim'))
    lim=25;
end

% Sample frequency
fs=FS*15.5;

% Filter high frequency noise
lgsr=length(s);
lgsr2=lgsr/2;
t=(1:lgsr)/fs;
% Low pass filter at 4Hz to eliminate high pass noise
[b,a]=ellip(4,0.1,40,1/fs); 
[H,w]=freqz(b,a,lgsr);
sf=filter(b,a,s);

% There is a buzz in the signal so first 350 points are discarded
sf350=sf(350:length(sf));
% Calculate slopes of each point in the signal
sf_pend=diff(sf350);

% Detect sign changes. From + to - is a maximum. Otherwise a minimum
ind=1;
indmax=1;
indmin=1;
while ind<length(sf_pend)
    if sign(sf_pend(ind))>sign(sf_pend(ind+1))
        posmax(indmax)=ind;
        indmax=indmax+1;
    end
    if sign(sf_pend(ind))<sign(sf_pend(ind+1))
        posmin(indmin)=ind;
        indmin=indmin+1;
    end
    ind=ind+1;
end

% If no posmax or posmin is created is cause the signal just goes up or down
% We exist returning 0
if (exist('posmin', 'var') + exist('posmax', 'var')) <2 %either one or both don't exist
    [Sm,Sd,num_resp,tm] = deal(0);
    return
end

% If position of the first detected maximum is lower than the first minimum, it is discarded 
% since we consider the signals to start with a minimum
if posmax(1)<posmin(1)
    for n=1:(length(posmax)-1)
        posmax(n)=posmax(n+1);
    end
    posmax=posmax(1:end-1);
end
% If there is a last remaining posmin with its corresponding posmax it is discarded
if length(posmin)>length(posmax)
    posmin=posmin(1:end-1);
end

ind=1;
while ind<length(posmax)
    if(sf350(posmax(ind))-sf350(posmin(ind)))<=lim
        for n=ind:(length(posmax)-1)
            posmax(n)=posmax(n+1);
            posmin(n)=posmin(n+1);
        end
        % Decrease ind, since I advanced one position, the look would skip a value
        ind=ind-1;
        % Truncate, since when advancing I remove the last one so it won't be twice
        posmax=posmax(1:end-1);
        posmin=posmin(1:end-1);
    end
    ind=ind+1;
end
% return number of responses
num_resp=length(posmax);

% Calculate average time of recovery: time between the response peak and the moment in which
% the response line is the middle of the baseline
tm=0;   % inicialize so it doesn't give an error if there is one or less responses
for i=1:length(posmax)-1
    tm(i)=(posmin(i+1)-posmax(i))/2;
end

% to avoid some errors
if isempty(posmax) || isempty(posmin)
    [Sm,Sd,num_resp,tm] = deal(0);
    return
end

% Calculate vectors with response magnitudes and durations
for ind=1:length(posmax)
    Sm(ind)=sf350(posmax(ind))-sf350(posmin(ind));
    Sd(ind)=posmax(ind)-posmin(ind);
end
