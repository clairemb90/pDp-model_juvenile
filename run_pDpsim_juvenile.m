tit=1;

twoorthree=1; %'network architecture: one inhibitory population (1) or two inhibitory populations, one feedback and one feedworward (2) ');
randortuned=input('network type: random (1), E-I assemblies (2), E assemblies (3)');
analys=input('analysis: yes (1) or no (2) ');

pagexls=1; 
if randortuned==1
    ext='_small'; 
    ext1='_small_tboth';
elseif randortuned==2
    ext='_small_tboth';
elseif randortuned==3
    ext='_small_tboth';
    ext1='_small';
else 
    disp('network not available')
end

num_odors=6;
trials=4;
load('odorseq_bo_poisson.mat')

[param,connec,tot]=xlsread('parameters_juvenile.xlsx',pagexls);

tot_mat=length(param);

r_ob=int8(r_olfbs);clear r_olfbs

for mat=1:tot_mat %can be changed if simulation of only subsets of networks    
   
    Ninh=1000;
    
    Istim=ones(length(r_ob),Ninh);
    PA1=0;
    
    struct=1;
   
    if randortuned==1 
        load(strcat(connec{mat},ext)),load(strcat(connec{mat},ext1),'Ee','EnsE')
    elseif randortuned==2         
        load(strcat(connec{mat},ext))
    elseif randortuned==3
        load(strcat(connec{mat},ext))
        load(strcat(connec{mat},ext1),'w_ie')
        w_ie=1.65*w_ie;
    end
    
    Nexc=length(w_ee); 
    Ninh=length(w_ii);
    Nob=1500;

    y=param(mat,:);
    
    dt=0.1;
    
    %neuronal parameters
    tau_exc=30;
    tau_inh=10;
    E_e=0; 
    E_i=-70; 
    g_rest_e=1.35;
    g_rest_i=0.9;
    tau_e=85; 
    tau_i=50;
    tau_ref_e=8/dt;
    tau_ref_i=8/dt;
    tau_w=40;
    ada=1;
   
    E_th=-38; 
    E_th_i=-45;
  
    E_rest=-60; %
    E_rest_ii=-65;
    beka=10;

    %activity neuron
    
        ds=0.1;
    downs=floor(length(r_ob)*ds);
    perc=1;
    NumNe=[1:Nexc/perc];
    NumNi=[1:Ninh/perc];
    V_id=E_rest_ii.*ones(downs,Ninh/perc);
    V_ed=E_rest.*ones(downs,Nexc/perc);
    g_oed=zeros(downs,Nexc/perc);
    g_eed=zeros(downs,Nexc/perc);
    g_ied=zeros(downs,Nexc/perc);
    g_ief=zeros(downs,Nexc/perc);
    g_ieb=zeros(downs,Nexc/perc);
   
    %% solve
    t=1;
    ts=1;
   
    tsc=1;
    tsci=1;
    tsci_b=1;
    VI=E_rest_ii.*ones(1,Ninh);
    VE=E_rest.*ones(1,Nexc);
    GEE=zeros(1,Nexc);
    GIE=zeros(1,Nexc);
    GIEB=zeros(1,Nexc);
    GIEF=zeros(1,Nexc);
    GEI=zeros(1,Ninh);
    GII=zeros(1,Ninh);
    GOE=zeros(1,Nexc);
    GOI=zeros(1,Ninh);
    
    refrac_E=single(zeros(1,Nexc));
    refrac_I=single(zeros(1,Ninh));
    WE=zeros(1,Nexc);
  
    spikecount_E=zeros(1000000,2);
    spikecount_I=zeros(3000000,2);
    
    while t<=length(r_ob)
        
        V_e=VE;
        V_i=VI;
        g_ee=GEE;
        g_ibe=GIEB;
        g_ife=GIEF;
        if twoorthree==1
            g_ie=GIE;
        else
            g_ie=GIEB+GIEF;
        end
        g_ei=GEI;
        g_ii=GII;
        g_oe=GOE;
        g_oi=GOI;
        w_e=WE;
        
        VE=V_e+dt/tau_e.*((E_rest-V_e)+((g_oe+g_ee).*(E_e-V_e)+g_ie.*(E_i-V_e)-w_e)/g_rest_e);
        
        VI=V_i+dt/tau_i.*((E_rest_ii-V_i)+((g_oi+g_ei).*(E_e-V_i)+g_ii.*(E_i-V_i))/g_rest_i);
        
        WE=w_e+dt/tau_w*(ada*(V_e-E_rest)-w_e);
        
        refre=find(refrac_E(1,:)>0);
        refri=find(refrac_I(1,:)>0);
        VE(refre)=E_rest;
        VI(refri)=E_rest_ii;
        neufire=find(VE>E_th);
        neufiri=find(VI>E_th_i);
        
        spikecount_E(tsc:tsc+length(neufire)-1,1)=t;
        spikecount_E(tsc:tsc+length(neufire)-1,2)=neufire;
        tsc=tsc+length(neufire);
        spikecount_I(tsci:tsci+length(neufiri)-1,1)=t;
        spikecount_I(tsci:tsci+length(neufiri)-1,2)=neufiri;
        tsci=tsci+length(neufiri);
        obfire=find(r_ob(t,:)>0);
        refrac_E=refrac_E-1;
        refrac_I=refrac_I-1;
        refrac_E(neufire)=tau_ref_e;
        refrac_I(neufiri)=tau_ref_i;
        T=2;
        WE(1,neufire)=WE(1,neufire)+beka;
       
        GEE=g_ee.*exp(-dt/tau_exc)+y(1).*sum(w_ee(neufire,:),1);
        GEI=g_ei.*exp(-dt/tau_exc)+y(4).*sum(w_ei(neufire,:),1);
        GIE=g_ie.*exp(-dt/tau_inh)+y(2).*sum(w_ie(neufiri,:),1);
        GII=g_ii.*exp(-dt/tau_inh)+y(6).*sum(w_ii(neufiri,:),1);
        
        GOE=g_oe.*exp(-dt/tau_exc)+y(3).*sum(w_oe(obfire,:),1);
        GOI=g_oi.*exp(-dt/tau_exc)+y(5).*sum(w_oi(obfire,:),1);
        if mod(t-1,1/ds)==0 || t==1
            V_id(ts,:)=VI(1,NumNi);
            V_ed(ts,:)=VE(1,NumNe);
        
            g_oed(ts,:)=GOE(1,NumNe);
            g_eed(ts,:)=GEE(1,NumNe);
            g_ied(ts,:)=GIE(1,NumNe);
            if twoorthree==2
                g_ieb(ts,:)=GIEB(1,NumNe);
                g_ief(ts,:)=GIEF(1,NumNe);
            end
            ts=ts+1;
        end
        t=t+1;
        clear refre refri neufire neufiri obfire
    end
    clear refrac_E refrac_I 
    spikecount_E=spikecount_E(1:tsc-1,:);
    spikecount_I=spikecount_I(1:tsci-1,:);
    g_ed=g_oed+g_eed;
    
    %% ANALYSIS
    ACT{tit}=spikecount_E;
    ACT_I{tit}=spikecount_I;
    %gE_save{tit}=g_ed;
   % gI_save{tit}=g_ied;
    %V_save{tit}=V_ed;
    
    if analys==1
        [obs loss] = observables(num_odors,Nexc,Ninh,dt,spikecount_E,spikecount_I,g_oed,g_ed,g_ied);
        obs_mat(mat,:)=obs;
        loss_mat(mat,:)=loss;
        %[res_cot] = cotuning_bo(num_odors,Nexc,dt,g_oed,g_ed,g_ied,EnsE)
        %cot_mat(mat,:)=res_cot;
    end
    tit=tit+1
    clearvars -except mat obs_mat gE_save gI_save V_save twoorthree tit randortuned analys pagexls ext ext1 fact ods num_odors trials r_ob tot_mat PIN ACT ACT_I param connec tot
    %save()
end


%plot timecourse
dt=0.1;
for mat=1:tot_mat
    spikecount_E=ACT{mat};
for i=1:360
    times(:,i)=[200/dt*(i-1)+1;200/dt*i];
    %FR(mat,i,:)=sum(r_olfbs(times(1,i):times(2,i),:))/(times(2,i)-times(1,i))/dt*1000;  
    spikeE_temp=sort(spikecount_E((spikecount_E(:,1)>times(1,i)) & (spikecount_E(:,1)<times(2,i)),2));
     if ~isempty(spikeE_temp)
        [NspikE,EdspikE]=histcounts(spikeE_temp,'BinMethod','integers');
        FR(mat,i,round(EdspikE(1)):floor(EdspikE(length(EdspikE))))=NspikE/(times(2,i)-times(1,i))/dt*1000;      
    end
end
end
% mean networks
FRav=squeeze(squeeze(mean(mean(FR),3)));
% mean 4 trials
figure,
for od=1:6
    for ii=1:18
        FRtrav(od,ii)=mean(FRav([ii+15*(od-1) ii+15*(od-1)+90 ii+15*(od-1)+180]));
    end
    hold on,plot(smooth(squeeze(FRtrav(od,:)),2))
end

%plot correlation
dt=0.1;

for i=1:6*4
    times(:,i)=[1000/dt+3000/dt*(i-1);2500/dt+3000/dt*(i-1)];
    FR(i,:)=sum(r_olfbs(times(1,i):times(2,i),:))/(times(2,i)-times(1,i))/dt*1000;
end

for t=1:6*4
    for tt=1:6*4
        corr_OB(t,tt)=corr(FR(t,:)',FR(tt,:)');
    end
end
figure,imagesc(corr_OB([1 7 13 2 8 14 3 9 15 4 10 16 5 11 17 6 12 18],[1 7 13 2 8 14 3 9 15 4 10 16 5 11 17 6 12 18]))

