function FeatHAND = FeaturesHAND(fhand, sizegsr, fs, N)
    
    for k=1:size(fhand,2)
        [Sm,Sd,n_resp,tm]=fGSR(fhand(:,k),2);
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

    FeatHAND=table(Num_resp', frec_oc', Sm_sum', Sd_sum', Sa_sum', tm_sum', Sclmad', Sclstd');
    FeatHAND.Properties.VariableNames = {'Num_resp', 'frec_oc', 'Sm_sum', 'Sd_sum', 'Sa_sum', 'tm_sum', 'Sclmad', 'Sclstd'};