function FeatEMG = FeaturesEMG(femg, modality)
    % Filter not necessary. Some papers apply LPF between 20 and 500 Hz but this is above the current fs
    if strcmp(modality, 'emgfeat')
        [zc, ssc, ApEn, periodog, coeff1, coeff2, coeff3] = deal([]);
        for k=1:size(femg,2)
            s = femg(:,k);
    
            %Hallo los cruces por cero
            zc(k)=sum((s(2:end)>0)~=(s(1:end-1)>0));
            
            %Hallo los cambios de pendiente
            s_pend=diff(s);
            ssc(k)=sum((s_pend(2:end)>0)~=(s_pend(1:end-1)>0));
    
            %Hallo la entropia aproximada
            P=hist(femg(:,k),100);
            P=P(find(P));
            P=P/sum(P);
            ApEn(k)=sum(-P.*log2(P));
            %Hallo el periodograma
            periodog(k)=sum(abs(fft(femg(:,k))).^2);
            %Hallo los coeficientes del modelo AR
            coeff=lpc(femg(:,k),5);
            coeff1(k)=coeff(2); 
            coeff2(k)=coeff(3); 
            coeff3(k)=coeff(4);
    
            clear s s_pend
        end
    
        FeatEMG=table(zc', ssc', ApEn', periodog', coeff1', coeff2', coeff3');
        FeatEMG.Properties.VariableNames = {'zc', 'ssc', 'ApEn_emg', 'periodog_emg', 'coeff1_emg', 'coeff1_em2', 'coeff1_em3'};

    elseif strcmp(modality, 'emgstats')
        %20 cepstral coefficients
        COEFCEPS = [];
        for k=1:size(femg,2)
            [coefceps]=cceps(femg(:,k),20);
            COEFCEPS=[COEFCEPS coefceps];
        end
        COEFCEPS=COEFCEPS';
        COEFCEPS=array2table(COEFCEPS);

        FeatEMG=table(mean(femg)', std(femg)', mean(abs(femg))', moment(femg,2)', moment(femg,3)', moment(femg,4)', moment(femg,5)', ...
            min(femg)', max(femg)', median(femg)', prctile(femg,10)', prctile(femg,25)', prctile(femg,75)', prctile(femg,90)', ...
            geomean(abs(femg))', harmmean(femg)', mode(femg)', range(femg)', iqr(femg)', diag(cov(femg)), mad(femg)', std(femg,1)', ... 
            var(femg,1)', skewness(femg)', kurtosis(femg)');
        FeatEMG.Properties.VariableNames = {'mean_emg', 'std_emg', 'mav_emg', 'moment2_emg', 'moment3_emg', 'moment4_emg', ...
            'moment5_emg', 'min_emg', 'max_emg', 'median_emg', 'prc10_emg', 'prc25_emg', 'prc75_emg', ...
            'prc90_emg', 'geomean_emg', 'harmmean_emg', 'mode_emg', 'range_emg', 'iqr_emg', 'diag_emg', 'mad_emg', 'std1_emg', ...
            'var1_emg', 'skewness_emg', 'kurtosis_emg'};
        FeatEMG=[FeatEMG COEFCEPS];
    end
