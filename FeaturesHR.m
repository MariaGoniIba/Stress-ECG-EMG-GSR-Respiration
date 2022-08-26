function FeatHR = FeaturesHR(hr, fs, N, meanrr)
    
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
    
    % Create table with features
    FeatHR=table(mean(hr)', std(hr)', ratioL', ratioM', vlf', lf', mf', hf', log(vlf)', log(lf)', log(mf)', log(hf)', ...
        log(lf+mf+hf)', log(hf./(lf+mf+hf))', log(lf./hf)', (hf./(hf+lf))', RSA', frec01', frec12', frec23', frec34', ...
        (sqrt(lf)./meanrr)', (sqrt(hf)./meanrr)', min(hr)',max(hr)',median(hr)',prctile(hr,25)', prctile(hr,75)', ...
        geomean(abs(hr))',harmmean(hr)', mode(hr)', range(hr)', iqr(hr)',diag(cov(hr)),mad(hr)',std(hr,1)',var(hr,1)', skewness(hr)', ...
        kurtosis(hr)');
    FeatHR.Properties.VariableNames = {'meanhr', 'stdhr', 'ratioL', 'ratioM', 'vlf', 'lf', 'mf', 'hf', 'logvlf', 'loglf', ...
        'logmf', 'loghf', 'loglf+mf+hf', 'loghf./lf+mf+hf', 'loglf./hf', 'hf./hf+lf', 'RSA', 'frec01', 'frec12', 'frec23', ...
        'frec34', 'sqrt(lf)/meanrr', 'sqrt(hf)/meanrr', 'minhr','maxhr', 'medianhr','prctilehr25', 'prctilehr75', ...
        'geomeanabshr','harmmeanhr', 'modehr', 'rangehr', 'iqrhr','diagcovhr','madhr','stdhr1','varhr1', 'skewnesshr', ...
        'kurtosishr'};
