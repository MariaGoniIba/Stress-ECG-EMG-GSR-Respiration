function [rr,S1,S2,hurst]=intervalRR(ecg,FS)
    % Function that calculates the RR interval for an ecg segment. FS: number of samples by slot,
    % by default 32
    
    if(~exist('FS'))
        FS=32;
    end
    
    %Sample frequency
    fs=FS*15.5;
    
    %Check if peaks go up or down
    MEAN=mean(ecg);
    
    if(mean(ecg>MEAN)>0.5)
        %Flip it
        ecg=-ecg;
    end
    
    L=FS*30;
    
    H=zeros(floor(length(ecg)/L)-2,L/2+1);
    %Calculate a first approximation of the ecg waveform
    n=1;
    for ind=2:floor(length(ecg)/L)-2
        [a,b]=max(ecg(ind*L+(1:L)));
        H(n,:)=ecg(ind*L+b+(floor(-L/4):floor(L/4)));
        n=n+1;
    end
    pattern=mean(H);
    
    y=filter(pattern(end:-1:1),1,ecg);
    
    B=firls(1000,[0 0.004*FS/32 0.005*FS/32 1],[0 0 1 1]);
    y2=filter(B,1,y);
    
    ind=0;
    n=1;
    while ind<length(ecg)-L-500-L/2
        x1=y2(ind+(1:L)+floor(L/4)+500);
        x2=ecg(ind+(1:L));
        [a,b]=max(x1>0.7*max(x1));
        [c,d]=max(x2(b+1:min(b+30,length(x2))));
        r(n)=b+d+ind;
        ind=ind+b+30;
        n=n+1;
    end
    
    r2=r/fs;
    
    rr=diff(r2);
    
    %Delete values out of range
    rr=rr([1 find((rr(2:end)<1.75*rr(1:end-1)))+1]);
    rr=rr([1 find((rr(2:end)>0.75*rr(1:end-1)))+1]);
    
    s1=rr(1:end-1)*sqrt(2)/2+rr(2:end)*sqrt(2)/2;
    s2=rr(1:end-1)*sqrt(2)/2-rr(2:end)*sqrt(2)/2;
    
    S1=std(s1);
    S2=std(s2);
    
    %Hurst coefficient
    for i = 1:20
        m = floor(length(rr)/(2*i));
        for j=1:2*i;
            interval = rr(1+(j-1)*m:(j*m));
            M = mean(interval);
            x = (interval-M);
            V = cumsum(x);
            R(j) = max(V)-min(V);
            S(j) = std(interval);      
        end
        tau(i)=m;
        RS(i)=mean(R./S);
    end;
    %plot(log10(tau),log10(RS),'+');
    % xlabel('log(\tau)','FontSize',12);
    % ylabel('log(R/S)','FontSize',12);
    % hold on;
    q = polyfit(log10(tau),log10(RS),1);
    t = 0.8:0.01:2.7;   %chosen values based on the common values we were taking
    y = q(1)*t+(q(2));
    %hold off
    
    hurst=q(1);