function [obs loss] = observables(num_odors,Nexc,Ninh,dt,spikecount_E,spikecount_I,g_oed,g_ed,g_ied)
%calculate

% firing rates
FR=zeros(num_odors,Nexc);
FRs=zeros(num_odors,Nexc);
FRI=zeros(num_odors,Ninh);
mfr_o=zeros(num_odors,1);
mfrs_o=zeros(num_odors,1);
mfri_o=zeros(num_odors,1);
for oo=1:num_odors
    times(:,oo)=[1000/dt+3000/dt*(oo-1);2000/dt+3000/dt*(oo-1)]; %2500
    times_s(:,oo)=[1+3000/dt*(oo-1);1000/dt+3000/dt*(oo-1)];
    times_1(:,oo)=[1000+3000*(oo-1);2500+3000*(oo-1)];
end
for i=1:num_odors
    spikeE_temp=sort(spikecount_E((spikecount_E(:,1)>times(1,i)) & (spikecount_E(:,1)<times(2,i)),2));
    spikes_temp=sort(spikecount_E((spikecount_E(:,1)>times_s(1,i)) & (spikecount_E(:,1)<times_s(2,i)),2));
    spikeI_temp=sort(spikecount_I((spikecount_I(:,1)>times(1,i)) & (spikecount_I(:,1)<times(2,i)),2));
    if ~isempty(spikeE_temp)
        [NspikE,EdspikE]=histcounts(spikeE_temp,'BinMethod','integers');
        FR(i,round(EdspikE(1)):floor(EdspikE(length(EdspikE))))=NspikE/(times(2,i)-times(1,i))/dt*1000;
        mfr_o(i)=mean(FR(i,:));
    end
    if ~isempty(spikes_temp)
        [Nspiks,Edspiks]=histcounts(spikes_temp,'BinMethod','integers');
        FRs(i,round(Edspiks(1)):floor(Edspiks(length(Edspiks))))=Nspiks/(times_s(2,i)-times_s(1,i))/dt*1000;
        mfrs_o(i)=mean(FRs(i,:));
    end
    if ~isempty(spikeI_temp)
        [NspikIb,EdspikIb]=histcounts(spikeI_temp,'BinMethod','integers');
        FRI(i,round(EdspikIb(1)):floor(EdspikIb(length(EdspikIb))))=NspikIb/(times(2,i)-times(1,i))/dt*1000;
        mfri_o(i)=mean(FRI(i,:));
%          mfri_o1(i)=mean(FRI(i,1:500));
%           mfri_o2(i)=mean(FRI(i,501:1000));
    end
    clear spikeE_temp spikeI_temp
end
mfr=mean(mfr_o);
std_fr=std(mfr_o);
mfr_s=mean(mfrs_o);
mfr_i=mean(mfri_o);
% mfr_i1=mean(mfri_o1);
% mfr_i2=mean(mfri_o2);

% correlation between odor-evoked activity patterns
for t=1:num_odors
    for tt=1:num_odors
        corr_E(t,tt)=corr(FR(t,:)',FR(tt,:)');
    end
end

CC_1=mean(triu(corr_E,1));
CC_1(CC_1==0)=[];
CC=mean(CC_1);

% gOE, gSyn, percentage of recurrent input
G_aff=g_oed*70;
G_tot=g_ed*70;
G_totI=g_ied*70;
tim=[];
for ave=1:num_odors
    tim=[1+3000*(ave-1):2500+3000*(ave-1)];
    Gaff_av(ave,:,:)=G_aff(tim,:);
    G_tot_av(ave,:,:)=G_tot(tim,:);
    G_totI_av(ave,:,:)=G_totI(tim,:);
    exc_g(ave,:)=mean(G_tot_av(ave,1000:2500,:)-mean(G_tot_av(ave,501:1000,:),2),2);
    inh_g(ave,:)=mean(G_totI_av(ave,1000:2500,:)-mean(G_totI_av(ave,501:1000,:),2));
end

aff_g_norm=mean(mean(mean(Gaff_av(:,1000:2500,:)-mean(G_tot_av(:,501:1000,:),2),2)))/70;
exc_g_norm=mean(mean(exc_g))/70;
inh_g_norm=mean(mean(inh_g))/70;
syn_g=(exc_g_norm+inh_g_norm);
percent_rec=(exc_g_norm-aff_g_norm)/exc_g_norm*100;

obs=[mfr_s mfr mfr_i CC aff_g_norm percent_rec syn_g];
loss=[mfr_s<0.15 0.9<mfr<1.1 -0.15<CC<0.15 0.03<aff_g_norm<0.16 78<percent_rec<90 1<syn_g<9];

end


