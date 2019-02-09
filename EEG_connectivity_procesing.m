% Copyright @ Jimmy Cao at University of Tasmania, Australia;
% Refer to paper: Cao et al (2016). Resting-state EEG power and coherence vary between migraine phases. The journal of headache and pain, 17(1), 102.

clc;
clear;
close all;

HeadFilePath = 'eeglab14_1_2b\plugins\dipfit2.3\standard_BESA\';

file_structure = dir('D:\chronic pain project\data(raw)\CM\');
m=size(file_structure,1);

channel = {  '1'  '2'  '3'  '4'   '5' '6'  '7'  '8'  '9'  '10'  '11' '12'  '13'  '14' '15' '16'  '17' '18' '19' '20'...
    '21' '22' '23' '24' '25' '26' '27' '28' '29'  '30' } ;
color = colormap(jet);

addpath('C:\Users\shouchungchang\Documents\MATLAB\');
addpath('C:\Users\shouchungchang\Documents\MATLAB\eeglab14_1_1b\');

hdmfile = [HeadFilePath 'standard_BESA.mat'];
mrifile = [HeadFilePath 'avg152t1.mat'];
chanfile = [HeadFilePath 'standard-10-5-cap385.elp'];


for ii = 3:m
    
    % clearvars -except file_structure channel color   ii
    close all;

    %Pro-processing
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    % ------------------------------------------------
    location1 = strcat('D:\chronic pain project\data(raw)\CM\',file_structure(ii).name,'\open\');
    location2 = strcat('D:\chronic pain project\data(raw)\CM\',file_structure(ii).name,'\close\');

  

    mkdir('D:\chronic pain project\result\CM\Coh\',file_structure(ii).name);
    mkdir('D:\chronic pain project\result\CM\Pdc\',file_structure(ii).name);
    mkdir('D:\chronic pain project\result\CM\iCoh\',file_structure(ii).name);

    locationsave21=strcat('D:\chronic pain project\result\CM\Coh\',file_structure(ii).name,'\open_Coh');
    locationsave22=strcat('D:\chronic pain project\result\CM\Coh\',file_structure(ii).name,'\close_Coh');
    locationsave23=strcat('D:\chronic pain project\result\CM\Pdc\',file_structure(ii).name,'\open_Pdc');
    locationsave24=strcat('D:\chronic pain project\result\CM\Pdc\',file_structure(ii).name,'\close_Pdc');
    locationsave25=strcat('D:\chronic pain project\result\CM\iCoh\',file_structure(ii).name,'\open_iCoh');
    locationsave26=strcat('D:\chronic pain project\result\CM\iCoh\',file_structure(ii).name,'\close_iCoh');

  


    %location:open
    file_structure2 = dir(location1);
    m2=size(file_structure2,1);
    for i1=3:m2
        location_open{i1-2,:}=strcat(location1,file_structure2(i1).name);
    end

    %location:close
    file_structure3 = dir(location2);
    m3=size(file_structure3,1);
    for i2=3:m3
        location_close{i2-2,:}=strcat(location2,file_structure3(i2).name);
    end

    % open
    EEG = pop_loadcnt( location_open{1,1}, 'dataformat', 'auto','memmapfile', '');
    EEG = pop_editset(EEG, 'setname',  [Set ExtractEpoch{j, 2} rj]);
    eeglab redraw;
    EEG = eeg_checkset( EEG ); eeglab redraw;

    EEG = pop_resample(EEG, 250); eeglab redraw;

    EEG = pop_chanedit(EEG, 'lookup', chanfile);
    EEG = eeg_checkset( EEG ); eeglab redraw;
    EEG = pop_select( EEG,'nochannel',{'EKG(L)' 'EKG(R)'  'A1' 'A2' 'VEOU' 'VEOL'});
    eeglab redraw;

    EEG = pop_eegfilt( EEG, 0, 50,  0, 0);
    EEG = pop_eegfilt( EEG, .5, 0,   0, 0);
    eeglab redraw;
    
    %% reject artifacts by visual inspection (optional)
%     pop_eegplot( EEG, 1, 0, 1);
%     keyboard;
%     eeglab redraw;
    
    %% reject artifacts (reconstruct signals) and bad channels automatically, i.e., ASR method
    EEG = clean_rawdata(EEG, [],[],[],[],[],[]);
    eeglab redraw
    
    %% ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'options', {'extended', 1, 'maxsteps', 1024, 'stop', 1e-7});
    diary off;
    
    % save1
    EEG = pop_saveset(EEG, 'filename', '_ICA.set', 'filepath', locationsaveica);
    eeglab redraw;
    
    %% ADJUST
%     EEG = interface_ADJ(EEG, [ 'adjust_report.txt']);
    lag=5; %epoch duration (in seconds)
    disp(['Continuous dataset epoched in ' num2str(lag) ' sec long epochs to compute feature statistics over time']);
    % check whether '5sec' events are present
    si=0; 
    for i=1:length(EEG.event) 
        if(strcmp(EEG.event(1,i).type, '5sec')==1) 
            si=1; 
        end
    end
    if si==0 %add events
        ntrials=floor((EEG.xmax-EEG.xmin)/lag);
        nevents=length(EEG.event);
        for index=1:ntrials
            EEG.event(index+nevents).type='5sec';
            EEG.event(index+nevents).latency=1+(index-1)*lag*EEG.srate; %EEG.srate is the sampling frequency
            latency(index)=1+(index-1)*lag*EEG.srate;
        end
        
        EEG=eeg_checkset(EEG,'eventconsistency');
    end
        
    EEGep = pop_epoch( EEG, {  '5sec'  }, [0 lag], 'newname', [EEG.setname '_ep5'] , 'epochinfo', 'yes');
%   % removing baseline
%    EEGep = pop_rmbase( EEGep, [0  0]);
    EEGep = eeg_checkset(EEGep);
        
    % collects ICA data from EEG
    if isempty(EEGep.icaact)
        disp('Warning: EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
        % Next instruction: see eeg_checkset
        EEGep.icaact = reshape(EEGep.icaweights*EEGep.icasphere*reshape(EEGep.data(1:size(EEGep.icaweights,1),:,:),[size(EEGep.icaweights,1) size(EEGep.data,2)*size(EEGep.data,3)]),[size(EEGep.icaweights,1) size(EEGep.data,2) size(EEGep.data,3)]);
    end    
    
    % Now that dataset is epoched, run ADJUST
    [art, horiz, vert, blink, disc,...
        soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin] = ADJUST (EEGep, [EEG.setname '_ADJUST_report.txt']);
    
    % Saving artifacted ICs for further analysis

    nome=['List_' EEG.setname '.mat'];

    save (nome, 'blink', 'horiz', 'vert', 'disc');

    disp(' ')
    disp(['Artifact ICs list saved in ' nome]);

    EEG = pop_subcomp(EEG, art);
    eeglab redraw;  
    
    % save2
    EEG = pop_saveset(EEG, 'filename', '_ICA_ADJUST.set', 'filepath', locationsaveica);
    
    % save3
    EEG = pop_saveset( EEG, 'filename','removeEyesBlink&Move.set', 'filepath', locationsaveica);



    %%  connectivity - Coh/PDC/iCoh
    nfft=50; % number of frequency bins
    fc=250; % sample frequency

    % MVAR coefficients
    % coefficients of extended MVAR model
    %[Bm,B0,Sw]=simuMVARcoeff(simunumber);


    Su=eye(30,30);


    y1=EEG.data;
    gg=size(y1,3);


    for xx=1:gg
        Am=squeeze(y1(:,:,xx));
        % Theoretical spectral functions
        [dc,dtf,pdc,gpdc,coh,pcoh,pcoh2,h,s,pp,f] = fdMVAR_5order(double(Am),Su,nfft,fc);


        coh=abs(pcoh);
        pdc=abs(pdc);
        icoh=abs(gpdc);

        COH(:,:,:,xx)=coh;
        PDC(:,:,:,xx)=pdc;
        iCOH(:,:,:,xx)=icoh;
    end

    COH_OPEN=COH(:,:,:,1:2:end);
    COH_CLOSE=COH(:,:,:,2:2:end);
    PDC_OPEN=PDC(:,:,:,1:2:end);
    PDC_CLOSE=PDC(:,:,:,2:2:end);
    ICOH_OPEN=iCOH(:,:,:,1:2:end);
    ICOH_CLOSE=iCOH(:,:,:,2:2:end);

    band={(1:3),(4:7),(8:12),(13:30),(26:50)};

    for rr=1:5  
        tem=COH_OPEN(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        Coh_open(:,:,rr)=tem3;
    end

    for rr=1:5  
        tem=COH_CLOSE(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        Coh_close(:,:,rr)=tem3;
    end


    for rr=1:5  
        tem=PDC_OPEN(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        Pdc_open(:,:,rr)=tem3;
    end

    for rr=1:5  
        tem=PDC_CLOSE(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        Pdc_close(:,:,rr)=tem3;
    end

    for rr=1:5  
        tem=ICOH_OPEN(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        iCoh_open(:,:,rr)=tem3;
    end

    for rr=1:5  
        tem=ICOH_CLOSE(:,:,band{1,rr},:)
        tem2=mean(tem,4);
        tem3=mean(tem2,3);
        iCoh_close(:,:,rr)=tem3;
    end

    % save
    save([locationsave21],'Coh_open');
    save([locationsave22],'Coh_close');
    save([locationsave23],'Pdc_open');
    save([locationsave24],'Pdc_close');
    save([locationsave25],'iCoh_open');
    save([locationsave26],'iCoh_close');



end

