function FeatRESP = FeaturesRESP(fresp, sizeresp, fs,N, modality)
    
    if strcmp(modality, 'respfeat')
        [num_respir, frec_resp, Vt_sum, Ti_sum, Te_sum, Ttot_sum, Vt_mean, Ti_mean, Te_mean, Ttot_mean, Vt_mean_std, ...
            Ti_mean_std, Te_mean_std, Ttot_mean_std, Vt_std, Ti_std, Te_std, Ttot_std, Vt_min, Ti_min, Te_min, Ttot_min, ...
            Vt_max, Ti_max, Te_max, Ttot_max, ratio_resp1_sum, ratio_resp2_sum, ratio_resp3, ApEnVt, ApEnTi, ApEnTe, ...
            ApEnTtot, Vt_cov, Ti_cov, Te_cov, Ttot_cov, ratio_resp1_cov, ratio_resp2_cov, ratio_resp3_cov, Vt_ar, Ti_ar, ...
            Te_ar, Ttot_ar, ratio_resp2_ar, ratio_resp1_ar, ratio_resp3_ar, Ti_plot, Te_plot, Ttot_plot]=deal([]);
       
            for k=1:size(fresp,2)
                [num_respir(k),Vt,Ti,Te]=respsignal(fresp(:,k));
                frec_resp(k)=num_respir(k)/((N/fs)*sizeresp);    
                %seconds, time of the slot by number of resp slots
                Vt_sum(k)=sum(Vt);  %sum of total volume
                Ti_sum(k)=sum(Ti);  %sum of inhaling times
                Te_sum(k)=sum(Te);  %sum of exhaling times
                Vt_mean(k)=mean(Vt);  %mean of total volume
                Ti_mean(k)=mean(Ti);  %mean of inhaling times
                Te_mean(k)=mean(Te);  %mean of exhaling times
                Vt_std(k)=std(Vt);  %std of total volume
                Ti_std(k)=std(Ti);  %std of inhaling times
                Te_std(k)=std(Te);  %std of exhaling times
                ratio_resp1=Vt./Ti; %mean inhaling flux
                ratio_resp1_sum(k)=sum(ratio_resp1);
                Ttot=Ti(1:length(Ti)-1)+Te;
                if isempty(Ttot)
                    Ttot=0;
                end
                Ttot_sum(k)=sum(Ttot);
                Ttot_mean(k)=mean(Ttot);
                Ttot_std(k)=std(Ttot);  
                Vt_mean_std(k)=Vt_mean(k)+Vt_std(k);
                Ti_mean_std(k)=Ti_mean(k)+Ti_std(k);
                Te_mean_std(k)=Te_mean(k)+Te_std(k);
                Ttot_mean_std(k)=Ttot_mean(k)+Ttot_std(k);
                Vt_min(k)=min(Vt);  
                Ti_min(k)=min(Ti);
                Te_min(k)=min(Te);
                Ttot_min(k)=min(Ttot);
                Vt_max(k)=max(Vt);  
                Ti_max(k)=max(Ti);
                Te_max(k)=max(Te);
                Ttot_max(k)=max(Ttot);
                ratio_resp2=Ti(1:length(Ti)-1)./Ttot;   
                ratio_resp2_sum(k)=sum(ratio_resp2);
                ratio_resp3(k)=frec_resp(k)/Vt_sum(k);  
                Vt_cov(k)=cov(Vt);  
                Ti_cov(k)=cov(Ti);  
                Te_cov(k)=cov(Te);  
                Ttot_cov(k)=cov(Ttot);  
                ratio_resp1_cov(k)=cov(ratio_resp1);
                ratio_resp2_cov(k)=cov(ratio_resp2);  
                ratio_resp3_cov(k)=cov(ratio_resp3); 
                %Approximate entropy
                P=hist(Vt,1000);       
                P=P(find(P));
                P=P/sum(P);
                ApEnVt(k)=sum(-P.*log2(P));
                P=hist(Ti,1000);   
                P=P(find(P));
                P=P/sum(P);
                ApEnTi(k)=sum(-P.*log2(P));
                P=hist(Te,1000);    
                P=P(find(P));
                P=P/sum(P);
                ApEnTe(k)=sum(-P.*log2(P));
                P=hist(Ttot,1000);  
                P=P(find(P));
                P=P/sum(P);
                ApEnTtot(k)=sum(-P.*log2(P));
                Vt_ar1=lpc(Vt,1);   %AR model 1st order
                if length(Vt_ar1)<2
                    Vt_ar(k)=0;
                else
                    Vt_ar(k)=(Vt_ar1(:,2))';   
                end
                Ti_ar1=lpc(Ti,1);
                if length(Ti_ar1)<2
                    Ti_ar(k)=0;
                else
                    Ti_ar(k)=(Ti_ar1(:,2))';   
                end
                Te_ar1=lpc(Te,1);
                if length(Te_ar1)<2
                    Te_ar(k)=0;
                else
                    Te_ar(k)=(Te_ar1(:,2))';   
                end
                Ttot_ar1=lpc(Ttot,1);
                if length(Ttot_ar1)<2
                    Ttot_ar(k)=0;
                else
                    Ttot_ar(k)=(Ttot_ar1(:,2))';
                end
                if isempty(ratio_resp2)
                    ratio_resp2_ar(k)=0;
                else
                    ratio_resp2_ar1=lpc(ratio_resp2,1);
                    if length(ratio_resp2_ar1)<2
                        ratio_resp2_ar(k)=0;
                    else
                        ratio_resp2_ar(k)=(ratio_resp2_ar1(:,2))';    %for ratio_resp2
                    end
                end
                ratio_resp1_ar1=lpc(ratio_resp1,1);
                if length(ratio_resp1_ar1)<2
                    ratio_resp1_ar(k)=0;
                else
                    ratio_resp1_ar(k)=(ratio_resp1_ar1(:,2))';    %for ratio_resp1
                end
                ratio_resp3_ar1=lpc(ratio_resp3,1);
                if length(ratio_resp3_ar1)<2
                    ratio_resp3_ar(k)=0;
                else
                    ratio_resp3_ar(k)=(ratio_resp3_ar1(:,2))';    %for ratio_resp3
                end
            end
            
            %Create table with features
            FeatRESP=table(num_respir', frec_resp', Vt_sum', Ti_sum', Te_sum', Ttot_sum', Vt_mean', Ti_mean', Te_mean', Ttot_mean', ...
                Vt_mean_std', Ti_mean_std', Te_mean_std', Ttot_mean_std', Vt_std', Ti_std', Te_std', Ttot_std', Vt_min', Ti_min', ...
                Te_min', Ttot_min', Vt_max', Ti_max', Te_max', Ttot_max', ratio_resp1_sum', ratio_resp2_sum', ratio_resp3', ApEnVt', ...
                ApEnTi', ApEnTe', ApEnTtot', Vt_cov', Ti_cov', Te_cov', Ttot_cov', ratio_resp1_cov', ratio_resp2_cov', ratio_resp3_cov', ...
                Vt_ar', Ti_ar', Te_ar', Ttot_ar', ratio_resp1_ar', ratio_resp2_ar', ratio_resp3_ar');
            FeatRESP.Properties.VariableNames = {'num_respir', 'frec_resp', 'Vt_sum', 'Ti_sum', 'Te_sum', 'Ttot_sum', 'Vt_mean', 'Ti_mean', 'Te_mean', 'Ttot_mean', ...
                'Vt_mean_std', 'Ti_mean_std', 'Te_mean_std', 'Ttot_mean_std', 'Vt_std', 'Ti_std', 'Te_std', 'Ttot_std', 'Vt_min', 'Ti_min', ...
                'Te_min', 'Ttot_min', 'Vt_max', 'Ti_max', 'Te_max', 'Ttot_max', 'ratio_resp1_sum', 'ratio_resp2_sum', 'ratio_resp3', 'ApEnVt', ...
                'ApEnTi', 'ApEnTe', 'ApEnTtot', 'Vt_cov', 'Ti_cov', 'Te_cov', 'Ttot_cov', 'ratio_resp1_cov', 'ratio_resp2_cov', 'ratio_resp3_cov', ...
                'Vt_ar', 'Ti_ar', 'Te_ar', 'Ttot_ar', 'ratio_resp1_ar', 'ratio_resp2_ar', 'ratio_resp3_ar'};

    elseif strcmp(modality, 'respstats')
        D=fs/N;
        f=(0:N-1)*D;
        muestrasVLF=find(f>0.04,1)-1;  
        muestrasLF=find(f>0.08,1)-1;    %no. samples below LF
        muestrasMF=find(f>0.15,1)-1;    %no. samples below MF
        muestrasHF=find(f>0.5,1)-1;     %no. samples below HF
        
        %no. samples below 01Hz,12Hz,23Hz y 34Hz
        muestras01=find(f>0.1,1)-1;    
        muestras02=find(f>0.2,1)-1;    
        muestras03=find(f>0.3,1)-1;    
        muestras04=find(f>0.4,1)-1;    

        energ_resp=((abs(fft(fresp))).^2)/N;
        vlf_resp=(sum(energ_resp((2:muestrasVLF),:)))/N;
        lf_resp=(sum(energ_resp((2:muestrasLF),:)))/N;
        mf_resp=(sum(energ_resp((muestrasLF+1:muestrasMF),:)))/N;
        hf_resp=(sum(energ_resp((muestrasMF+1:muestrasHF),:)))/N;
        frec01_resp=(sum(energ_resp((2:muestras01),:)))/N;  %energy between 0-0.1Hz
        frec12_resp=(sum(energ_resp((muestras01+1:muestras02),:)))/N; %energy between 0.1-0.2Hz
        frec23_resp=(sum(energ_resp((muestras02+1:muestras03),:)))/N; %energy between 0.2-0.3Hz
        frec34_resp=(sum(energ_resp((muestras03+1:muestras04),:)))/N; %energy between 0.3-0.4Hz
        ratioL_resp=lf_resp./hf_resp;  
        
        %Recurrence analysis
        Nr=zeros(1,size(fresp,2));
        Ni=zeros(1,size(fresp,2));
        E=zeros(1,size(fresp,2));
        for k=1:size(fresp,2)
            [Nr(k),Ni(k),E(k)]=recurrence_analysis(fresp(:,k));
        end
        
        % Welch periodogram to find maximum peaks and frequencies. To check if it detects correctly the maximum,
        % do a semilogy(y)
        potmax=zeros(1,size(fresp,2));
        pospotmax=zeros(1,size(fresp,2));
        for k=1:size(fresp,2)
            y=pwelch(fresp(:,k),512);
            [potmax(k),pospotmax(k)]=max(y(3:18));
            pospotmax(k)=pospotmax(k)+2;
        end
        
        VTE=std(fresp)./mean(fresp);
        aux=sign(diff(fresp));
        MARIAFEAT=mean(aux(2:end,:)~=aux(1:end-1,:));
        DRESP=diff(fresp);   
        DDRESP=diff(DRESP);
        DDDRESP=diff(DDRESP);

        % Autoregressive model 1st order
        lpc_resp=lpc(fresp,1);
        lpc_resp=(lpc_resp(:,2))';
        % Sample entropy of respiration
        P=hist(fresp,1000);
        P=P./(ones(1000,1)*sum(P));
        ApEnResp=sum(-P.*log2(P+eps));
        % Sample entropy of first derivative of resp
        P=hist(DRESP,1000);
        P=P./(ones(1000,1)*sum(P));
        ApEnDResp=sum(-P.*log2(P+eps));
        % Sample entropy of second derivative of resp
        P=hist(DDRESP,1000);
        P=P./(ones(1000,1)*sum(P));
        ApEnDDResp=sum(-P.*log2(P+eps));
        % Sample entropy of third derivative of resp
        P=hist(DDDRESP,1000);
        P=P./(ones(1000,1)*sum(P));
        ApEnDDDResp=sum(-P.*log2(P+eps));

        FeatRESP=table(mean(fresp)', std(fresp)', VTE', MARIAFEAT', mean(DRESP)', mean(DDRESP)', mean(DDDRESP)', min(fresp)', max(fresp)', ...
            median(fresp)', prctile(fresp,25)', prctile(fresp,75)', geomean(abs(fresp))', harmmean(fresp)', mode(fresp)', range(fresp)', ...
            iqr(fresp)', diag(cov(fresp)), mad(fresp)', std(fresp,1)', var(fresp,1)', skewness(fresp)', kurtosis(fresp)', vlf_resp', ...
            lf_resp', mf_resp', hf_resp', frec01_resp', frec12_resp', frec23_resp', frec34_resp', ratioL_resp',lpc_resp', ApEnResp', ...
            ApEnDResp', ApEnDDResp', ApEnDDDResp', potmax', pospotmax', Nr', Ni', E');
        FeatRESP.Properties.VariableNames = {'mean_resp','std_resp','VTE','MARIAFEAT','meanDRESP','meanDDRESP','meanDDDRESP','min_resp', ...
            'max_resp','median_resp','prctile25_resp','prctile75_resp','geomean_resp','harmmean_resp','mode_resp','range_resp', ...
            'iqr_resp','diag_resp','mad_resp','std1_resp','var1_resp','skewness_resp','kurtosis_resp', 'vlf_resp', 'lf_resp', ...
            'mf_resp','hf_resp','frec01_resp','frec12_resp','frec23_resp','frec34_resp','ratioL_resp', 'lpc_resp', 'ApEnResp', ...
            'ApEnDResp','ApEnDDResp','ApEnDDDResp','potmax','pospotmax','Nr','Ni','E'};
    end