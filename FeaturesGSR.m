function FeatGSR = FeaturesGSR(s, sizegsr, fs,N, modality)
    
    if strcmp(modality,'hand')
        for k=1:size(s,2)
            [Sm,Sd,n_resp,tm]=fGSR(s(:,k),2);
            Num_resp(k)=n_resp;
            frec_oc(k)=n_resp/((N/fs)*sizegsr);    
            %seconds. Time of the slot by number of slots
            Sa=0.5*(Sm.*Sd);    %area of the responses
            Sm_sum(k)=sum(Sm);
            Sd_sum(k)=sum(Sd);
            Sa_sum(k)=sum(Sa);
            tm_sum(k)=sum(tm);
            Sclmad(k)=mad(Sm); 
            Sclstd(k)=std(Sm);
        end
    
        FeatGSR=table(Num_resp', frec_oc', Sm_sum', Sd_sum', Sa_sum', tm_sum', Sclmad', Sclstd');
        FeatGSR.Properties.VariableNames = {'Num_resp', 'frec_oc', 'Sm_sum', 'Sd_sum', 'Sa_sum', 'tm_sum', 'Sclmad', 'Sclstd'};
    
    elseif strcmp(modality, 'handstats')
        FeatGSR=table(mean(s)', mad(s)', std(s)', min(s)', max(s)', median(s)', prctile(s,25)', ...
        prctile(s,75)', geomean(abs(s))', harmmean(s)', mode(s)', range(s)', iqr(s)', diag(cov(s)), mad(s)', ...
        std(s,1)', var(s,1)', skewness(s)', kurtosis(s)');
        
        FeatGSR.Properties.VariableNames = {'mean_hand', 'mand_hand', 'std_hand', 'min_hand', 'max_hand', 'median_hand', ...
            'prctile25_hand', 'prctile75_hand', 'geomean_hand', 'harmmean_hand', 'mode_hande', 'range_hand', 'iqr_hand', 'diag_hand', ...
            'mad_hand', 'std1_hand', 'var1_hand', 'skew_hand', 'kurtosis_hand'};

    elseif strcmp(modality,'foot')
        for k=1:size(s,2)
            [Sm,Sd,n_resp,tm]=fGSR(s(:,k),2);
            Num_resp(k)=n_resp;
            frec_oc(k)=n_resp/((N/fs)*sizegsr);    
            %seconds. Time of the slot by number of slots
            Sa=0.5*(Sm.*Sd);    %area of the responses
            Sm_sum(k)=sum(Sm);
            Sd_sum(k)=sum(Sd);
            Sa_sum(k)=sum(Sa);
            tm_sum(k)=sum(tm);
            Sclmad(k)=mad(Sm); 
            Sclstd(k)=std(Sm);
        end
    
        FeatGSR=table(Num_resp', frec_oc', Sm_sum', Sd_sum', Sa_sum', tm_sum', Sclmad', Sclstd');
        FeatGSR.Properties.VariableNames = {'Num_resp', 'frec_oc', 'Sm_sum', 'Sd_sum', 'Sa_sum', 'tm_sum', 'Sclmad', 'Sclstd'};

    elseif strcmp(modality, 'footstats')
        FeatGSR=table(mean(s)', mad(s)', std(s)', min(s)', max(s)', median(s)', prctile(s,25)', ...
        prctile(s,75)', geomean(abs(s))', harmmean(s)', mode(s)', range(s)', iqr(s)', diag(cov(s)), mad(s)', ...
        std(s,1)', var(s,1)', skewness(s)', kurtosis(s)');
        
        FeatGSR.Properties.VariableNames = {'mean_foot', 'mand_foot', 'std_foot', 'min_foot', 'max_foot', 'median_foot', ...
            'prctile25_foot', 'prctile75_foot', 'geomean_foot', 'harmmean_foot', 'mode_foote', 'range_foot', 'iqr_foot', 'diag_foot', ...
            'mad_foot', 'std1_foot', 'var1_foot', 'skew_foot', 'kurtosis_foot'};
    end