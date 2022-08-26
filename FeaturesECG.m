function FeatECG = FeaturesECG(fecg,sizeecg,fs,N)

    D=fs/N;
    f=(0:N-1)*D;

    for k=1:size(fecg,2)
        [rr,s1(k),s2(k),hurst(k)]=intervalRR(fecg(:,k),sizeecg);
        %More features
        meanrr(k)=mean(rr);
        stdrr(k)=std(rr);
        rsaindex(k)=(max(rr)-min(rr))/mean(rr);
        dbd(k)=max(rr)-min(rr);
        ei(k)=max(rr)/min(rr);
        ie(k)=min(rr)/max(rr);
        sdnn(k)=std(rr);
        [Mo,posMo]=mode(rr);
        loadindex(k)=f(posMo)/(2*Mo*length(rr)*(max(rr)-min(rr)));
        hrvindex(k)=length(rr)/f(posMo);
        more50=sum(abs(([rr 0]-[0 rr]))>=0.05);
        less50=sum(abs(([rr 0]-[0 rr]))<0.05);
        nn50(k)=more50;
        pnn50(k)=(more50-1)/(more50+less50);
        DARR=abs([rr 0]-[0 rr]);
        rmssd(k)=std(DARR);
        mirr(k)=prctile(rr,75)-prctile(rr,25);
        mdarr(k)=median(DARR);
        %Approximate entropy
        P=hist(rr,100);     
        P=P(find(P));
        P=P/sum(P);
        ApEnRr(k)=sum(-P.*log2(P));
    end

    FeatECG=table(s1', s2', hurst', meanrr', stdrr', rsaindex', dbd', ei', ie', sdnn', loadindex', hrvindex', ...
        nn50', pnn50', rmssd', mirr', mdarr', ApEnRr');
    FeatECG.Properties.VariableNames = {'s1', 's2', 'hurst', 'meanrr', 'stdrr', 'rsaindex', 'dbd', 'ei', 'ie', 'sdnn', 'loadindex', 'hrvindex', ...
        'nn50', 'pnn50', 'rmssd', 'mirr', 'mdarr', 'ApEnRr'};
   