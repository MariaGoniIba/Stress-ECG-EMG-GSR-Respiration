function FeatHR = FeaturesHR(hr, fs, N)
    
    %Number of samples in LF,MF and HF
    D=fs/N;
    f=(0:N-1)*D;
    samplesVLF=find(f>0.04,1)-1;  
    samplesLF=find(f>0.08,1)-1;    %number of samples below LF
    samplesMF=find(f>0.15,1)-1;    %no. samples below MF
    samplesHF=find(f>0.5,1)-1;     %no. samples below HF

    samples01=find(f>0.1,1)-1;    
    samples02=find(f>0.2,1)-1;    
    samples03=find(f>0.3,1)-1;    
    samples04=find(f>0.4,1)-1;     

    %Extract features ratioL=lf/hf,ratioM=(lf+mf)/hf
    energ=((abs(fft(hr))).^2)/N;
    vlf=(sum(energ((2:samplesVLF),:)))/N;
    lf=(sum(energ((2:samplesLF),:)))/N;
    mf=(sum(energ((samplesLF+1:samplesMF),:)))/N;
    hf=(sum(energ((samplesMF+1:samplesHF),:)))/N;
    ratioL=lf./hf;
    ratioM=((lf+mf)./hf);
    
    %Other features
    RSA=hf./samplesHF; %integration HF band: 0.15-0.5Hz
    frec01=(sum(energ((2:samples01),:)))/N;  %energy between 0-0.1Hz
    frec12=(sum(energ((samples01+1:samples02),:)))/N; %energy between 0.1-0.2Hz
    frec23=(sum(energ((samples02+1:samples03),:)))/N; %energy between 0.2-0.3Hz
    frec34=(sum(energ((samples03+1:samples04),:)))/N; %energy between 0.3-0.4Hz
    
    FeatHR=[];    
    FeatHR=[FeatHR;mean(hr)];
    FeatHR=[FeatHR;std(hr)];
    FeatHR=[FeatHR;ratioL];
    FeatHR=[FeatHR;ratioM];
    FeatHR=[FeatHR;vlf];
    FeatHR=[FeatHR;lf];
    FeatHR=[FeatHR;mf];
    FeatHR=[FeatHR;hf];
    FeatHR=[FeatHR;log(vlf)];
    FeatHR=[FeatHR;log(lf)];       
    FeatHR=[FeatHR;log(mf)];        
    FeatHR=[FeatHR;log(hf)];        
    FeatHR=[FeatHR;log(lf+mf+hf)];  
    FeatHR=[FeatHR;log(hf./(lf+mf+hf))];    
    FeatHR=[FeatHR;log(lf./hf)];  
    FeatHR=[FeatHR;hf./(hf+lf)];
    FeatHR=[FeatHR;RSA];
    FeatHR=[FeatHR;frec01];
    FeatHR=[FeatHR;frec12];
    FeatHR=[FeatHR;frec23];
    FeatHR=[FeatHR;frec34];