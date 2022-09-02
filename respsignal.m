function [num_resp,Vt,Ti,Te]=respsignal(s)

% Function that calculates the respiratory interval, detecting those maximum moments
% of inhaling and minimum if exhaling.
fs=15.5*2; 

L=ceil(fs); %Number of samples equivalent to the minimum possible period

n=1;
ind=1;   
val=s(ind);
searchm=1; %This value will be 1 when we search for a max and -1 when searching for a minimum

while ind<length(s)-L
    [a,b]=max(searchm*s(ind+0:min(ind+L-1,length(s))));
    if(a~=val)
        val=a;
        ind=b+ind-1;
    else
        P(n)=ind;      
        if length(P)>1 && P(n)==P(n-1)
            break
        end
        searchm=-searchm;
        val=searchm*s(P(n));
        n=n+1;
    end
end
    
INSP=P(1:2:end);
EXP=P(2:2:end);

if (length(INSP)+length(EXP)==2) && INSP==1 && EXP==1
    [num_resp, Vt, Ti, Te]=deal(0);
else
    num_resp=length(EXP);
    
    %Calculate inhaling and exhaling times and volume
    for ind=1:length(EXP)
        Vt(ind)=s(EXP(ind))-s(INSP(ind));
        Ti(ind)=EXP(ind)-INSP(ind);
        if ind<length(EXP)
            Te(ind)=INSP(ind+1)-EXP(ind);
        end
    end
end

% plot(s)
% hold on
% plot(INSP,s(INSP),'ro');
% plot(EXP,s(EXP),'rx');