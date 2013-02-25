function screenProdData(subjID,varargin)
%%
% Overall steps: 
% 1. Manually go through all the trials, label the ones that
% should be througn out, and manually correct identification of landmarks
% such as uTime, jTime, etc. The order of trials is completely randomize in
% this step. The perturbation status of the trials are not marked at this
% tage. 
% 2. Manually go through all the perturbed trials. Mark those ones with bad
% perturbations. The order of the perturbed trials is also randomized.
% However, in one subject, the perturbed trials of the same type are
% blocked. But from subject to subject, the order of the four types (or any
% other number) of perturbations is completely randomized. 

%% Configs
[ret,hostName]=system('hostname');
hostName=deblank(hostName);
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    dacacheDir='E:/speechres/kape/dacache';
    rawDataDir='E:/DATA/KAPE/';
else
    dacacheDir='D:/speechres/kape/dacache';
    rawDataDir='D:/DATA/KAPE/';
end
% elseif isequal(hostName,'glossa')
%     dacacheDir='z:/speechres/apstv2/mcode/dacache';
%     rawDataDir='f:/DATA/APSTV2/';
% else
%     dacacheDir='z:/speechres/apstv2/mcode/dacache';
%     rawDataDir='z:/DATA/APSTV2/';    
% end
mvaWinWidth=21;     % 21 * 1.333 = 28 (ms)
fineParseWin=100e-3;	% sec

ylim=[0,5000];
rmsOneSide=0.01;  % Unit: sec
rmsLBRatio=0.75;
shiraOneSide=0.025; % Unit: sec

%%
expDir=fullfile(rawDataDir,subjID);

if ~isdir(expDir)
    fprintf('Cannot find directory %s. Terminated.\n',expDir);
    return
end
if ~isfile(fullfile(expDir,'expt.mat'))
    fprintf('Cannot find expt.mat in directory %s. Terminated.\n',expDir);
    return
end

%%
load(fullfile(expDir,'expt.mat'));  % gives expt
fprintf('Subject ID: \t%s\n',expt.subject.name);
fprintf('Subject gender: \t%s\n',expt.subject.sex);
% fprintf('Shift direction: %s\n',expt.subject.shiftDirection);
% fprintf('Shift ratio: %f\n',expt.subject.shiftRatio);

idx1=fsic(expt.recPhases,'pre');
idx2=fsic(expt.recPhases,'pract1');
idx3=fsic(expt.recPhases,'pract2');

idx4=fsic(expt.recPhases,'rand');

idx5=fsic(expt.recPhases,'start');
idx6=fsic(expt.recPhases,'ramp');
idx7=fsic(expt.recPhases,'stay');
idx8=fsic(expt.recPhases,'end');

% randPhases = expt.recPhases(setxor(1:numel(expt.recPhases),[idx1, idx2, idx3, idx]));
% nTrainPhases = numel(trainPhases);

% testPhases=expt.recPhases(setxor(1:numel(expt.recPhases),[idx1,idx2,idx3,idx7,idx8,idx9,idx10,idx11]));
% nTestPhases=numel(testPhases);

% prePhases=expt.recPhases(idx12);
% nPrePhases=numel(prePhases);

pdata.subject=expt.subject;
pdata.mvaWinWidth=mvaWinWidth;
pdata.fineParseWin=fineParseWin;

dacacheFN = fullfile(dacacheDir,[pdata.subject.name,'.mat']);
stateFN = fullfile(dacacheDir,[pdata.subject.name,'_state.mat']);

%% Calculate the total number of speech trials
randData.rawDataFNs=cell(1,0);
randData.phases=cell(1,0);
randData.blockNums=[];
randData.trialNums=[];
randData.words=cell(1,0);
randData.datenums=[];
randData.pertType = [];

t_phase = 'rand';

t_dir=fullfile(rawDataDir, subjID, t_phase);
d1 = dir(fullfile(t_dir, 'rep*'));
for i1=1:numel(d1)
    t_repdir = fullfile(t_dir, d1(i1).name);
    
    d2 = dir(fullfile(t_repdir, 'trial-*-1.mat'));
    
    for i2 = 1 : numel(d2)
        t_fn = fullfile(t_repdir, d2(i2).name);
        load(t_fn);     % gives data
    
        randData.rawDataFNs{end+1} = t_fn;
        randData.phases{end+1} = t_phase;
        randData.blockNums(end+1) = str2num(strrep(d1(i1).name,'rep',''));
        randData.trialNums(end+1) = str2num(strrep(strrep(d2(i2).name,'trial-',''),'-1.mat',''));
        randData.words{end+1} = data.params.name;
        randData.datenums(end+1) = datenum(data.timeStamp);
        if data.params.bShift == 0
            randData.pertType(end + 1) = 0;
        else
            if data.params.pertPhi(1) > 0
                randData.pertType(end + 1) = 1; % higher: F1 down, F2 up
            else
                randData.pertType(end + 1) = -1; % lower: F1 up, F2 down
            end
        end
    end
end

% -------------------------------------------
sustPhases = {'start', 'ramp', 'stay', 'end'};
nSustPhases = numel(sustPhases);

sustData.rawDataFNs=cell(1,0);
sustData.phases=cell(1,0);
sustData.blockNums=[];
sustData.trialNums=[];
sustData.words=cell(1,0);
sustData.datenums=[];
for i1 = 1 : nSustPhases
    t_phase = sustPhases{i1};
    if isdir(fullfile(rawDataDir, subjID, t_phase))
        d1 = dir(fullfile(rawDataDir, subjID, t_phase,'rep*'));
        for i2=1:numel(d1)
            t_dir=fullfile(rawDataDir,subjID,t_phase,d1(i2).name);
            d2=dir(fullfile(t_dir,'trial-*-1.mat'));
            for i3=1:numel(d2)
                t_fn=fullfile(t_dir,d2(i3).name);
                load(t_fn);     % gives data
                
                sustData.rawDataFNs{end+1}=t_fn;
                sustData.phases{end+1}=t_phase;
                sustData.blockNums(end+1)=str2num(strrep(d1(i2).name,'rep',''));
                sustData.trialNums(end+1)=str2num(strrep(strrep(d2(i3).name,'trial-',''),'-1.mat',''));
                sustData.words{end+1}=data.params.name;
                sustData.datenums(end+1)=datenum(data.timeStamp);
            end
        end
    end
end

randData.rmsThresh=nan(size(randData.rawDataFNs));
randData.fn1=nan(size(randData.rawDataFNs));
randData.fn2=nan(size(randData.rawDataFNs));
randData.aFact=nan(size(randData.rawDataFNs));
randData.bFact=nan(size(randData.rawDataFNs));
randData.gFact=nan(size(randData.rawDataFNs));
randData.bCepsLift = nan(size(randData.rawDataFNs));

randData.vowelOnset = nan(size(randData.rawDataFNs));
randData.vowelEnd = nan(size(randData.rawDataFNs));
randData.vowelOnsetIdx = nan(size(randData.rawDataFNs));
randData.vowelEndIdx = nan(size(randData.rawDataFNs));
randData.f1Traj = cell(size(randData.rawDataFNs));
randData.f2Traj = cell(size(randData.rawDataFNs));

randData.prodF1=nan(size(randData.rawDataFNs));
randData.prodF2=nan(size(randData.rawDataFNs));
randData.prodF1_LB=nan(size(randData.rawDataFNs));
randData.prodF2_LB=nan(size(randData.rawDataFNs));

randData.prodF1_shira=nan(size(randData.rawDataFNs));
randData.prodF2_shira=nan(size(randData.rawDataFNs));
randData.prodF1_mnlBound=nan(size(randData.rawDataFNs));
randData.prodF2_mnlBound=nan(size(randData.rawDataFNs));

randData.audF1=nan(size(randData.rawDataFNs));
randData.audF2=nan(size(randData.rawDataFNs));
randData.traj_F1=cell(size(randData.rawDataFNs));
randData.traj_F2=cell(size(randData.rawDataFNs));
randData.sigRMS=cell(size(randData.rawDataFNs));
randData.iv1=nan(size(randData.rawDataFNs));
randData.iv2=nan(size(randData.rawDataFNs));
randData.bDiscard=nan(size(randData.rawDataFNs));
randData.rating=nan(size(randData.rawDataFNs));
randData.comments=cell(size(randData.rawDataFNs));

randData.bPertOkay=nan(size(randData.rawDataFNs));


sustData.rmsThresh=nan(size(sustData.rawDataFNs));
sustData.fn1=nan(size(sustData.rawDataFNs));
sustData.fn2=nan(size(sustData.rawDataFNs));
sustData.aFact=nan(size(sustData.rawDataFNs));
sustData.bFact=nan(size(sustData.rawDataFNs));
sustData.gFact=nan(size(sustData.rawDataFNs));
sustData.bCepsLift = nan(size(sustData.rawDataFNs));

sustData.vowelOnset = nan(size(sustData.rawDataFNs));
sustData.vowelEnd = nan(size(sustData.rawDataFNs));
sustData.vowelOnsetIdx = nan(size(sustData.rawDataFNs));
sustData.vowelEndIdx = nan(size(sustData.rawDataFNs));
sustData.f1Traj = cell(size(sustData.rawDataFNs));
sustData.f2Traj = cell(size(sustData.rawDataFNs));

sustData.prodF1=nan(size(sustData.rawDataFNs));
sustData.prodF2=nan(size(sustData.rawDataFNs));
sustData.prodF1_LB=nan(size(sustData.rawDataFNs));
sustData.prodF2_LB=nan(size(sustData.rawDataFNs));

sustData.prodF1_shira=nan(size(sustData.rawDataFNs));
sustData.prodF2_shira=nan(size(sustData.rawDataFNs));
sustData.prodF1_mnlBound=nan(size(sustData.rawDataFNs));
sustData.prodF2_mnlBound=nan(size(sustData.rawDataFNs));

sustData.audF1=nan(size(sustData.rawDataFNs));
sustData.audF2=nan(size(sustData.rawDataFNs));
sustData.traj_F1=cell(size(sustData.rawDataFNs));
sustData.traj_F2=cell(size(sustData.rawDataFNs));
sustData.sigRMS=cell(size(sustData.rawDataFNs));
sustData.iv1=nan(size(sustData.rawDataFNs));
sustData.iv2=nan(size(sustData.rawDataFNs));
sustData.bDiscard=nan(size(sustData.rawDataFNs));
sustData.rating=nan(size(sustData.rawDataFNs));
sustData.comments=cell(size(sustData.rawDataFNs));

sustData.bPertOkay=nan(size(sustData.rawDataFNs));

pdata.randData = randData;
pdata.sustData = sustData;

%% Build a list of all trials
if isfile(stateFN)    
    fprintf('Found state.mat at %s.\n',stateFN);
    a=input('Resume? (0/1): ');
    if a==1
        load(stateFN);  % gives state
        load(dacacheFN);    % gives pdata
        trialList=state.trialList;
        trialListPert=state.trialListPert;
        bNew=0;
    else
        a = input('Are you sure you want to start the screening process over? (0/1): ');
        if a == 1
            bNew=1;
        else
            return
        end
    end
else
    bNew=1;
end
if bNew
    if isfile(dacacheFN)
        delete(dacacheFN);
        fprintf('%s deleted.\n',dacacheFN);
    end
    
    fprintf('\nstate.mat not found.\n');
    fprintf('Getting information about all trials...\n');
    trialList.fn={};
    trialList.phase={};
    trialList.block=[];
    trialList.trialN=[];   
    trialList.word={};
    trialList.allOrderN=[];
    
    fprintf('Randomizing the order of all trials...\n');

    idxRandPerm = randperm(numel(randData.rawDataFNs));
    
    trialListRand.fn=randData.rawDataFNs(idxRandPerm);
    trialListRand.phase=randData.phases(idxRandPerm);
    trialListRand.block=randData.blockNums(idxRandPerm);
    trialListRand.trialN=randData.trialNums(idxRandPerm);
    trialListRand.word=randData.words(idxRandPerm);
    idxOrder = 1:numel(randData.trialNums);
    trialListRand.allOrderN=idxOrder(idxRandPerm);
    
    idxRandPerm2 = randperm(numel(sustData.rawDataFNs));
    
    trialListSust.fn=sustData.rawDataFNs(idxRandPerm2);
    trialListSust.phase=sustData.phases(idxRandPerm2);
    trialListSust.block=sustData.blockNums(idxRandPerm2);
    trialListSust.trialN=sustData.trialNums(idxRandPerm2);
    trialListSust.word=sustData.words(idxRandPerm2);
    idxOrder2 = (1 : numel(sustData.trialNums));
    trialListSust.allOrderN=idxOrder2(idxRandPerm2);
    
    trialListPert=struct;
    idx_ramp = fsic(trialListSust.phase, 'ramp');
    idx_stay = fsic(trialListSust.phase, 'stay');
    idx_pert = [idx_ramp, idx_stay];
    
    trialListPert.fn = trialListSust.fn(idx_pert);
    trialListPert.phase = trialListSust.phase(idx_pert);
    trialListPert.block = trialListSust.block(idx_pert);
    trialListPert.trialN = trialListSust.trialN(idx_pert);
    trialListPert.allOrderN = trialListSust.allOrderN(idx_pert);
    
    % The perturbed trials in the rand phase
    idx_pert = find(randData.pertType ~= 0);
    trialListPert.fn = [trialListPert.fn, trialListRand.fn(idx_pert)];
    trialListPert.phase = [trialListPert.phase, trialListRand.phase(idx_pert)];
    trialListPert.block = [trialListPert.block, trialListRand.block(idx_pert)];
    trialListPert.trialN = [trialListPert.trialN, trialListRand.trialN(idx_pert)];
    trialListPert.allOrderN = [trialListPert.allOrderN, trialListRand.allOrderN(idx_pert)];

    % Combine the training and test lists
    nRandTrials = numel(trialListRand.fn);
    nSustTrials = numel(trialListSust.fn);
    
    flds=fields(trialList);
    for k1=1:numel(flds)
        fld=flds{k1};
        trialList.(fld)=[trialListRand.(fld), trialListSust.(fld)];        
    end
    
    trialList.isRand = [ones(1,nRandTrials), zeros(1, nSustTrials)];
    trialList.isSust = [zeros(1,nRandTrials), ones(1, nSustTrials)];
    
    state.trialList = trialList;
    state.trialListPert = trialListPert;
    
    state.idx=1;
    save(stateFN,'state');
    save(dacacheFN,'pdata');
    
    fprintf('Saved state info to %s\n', stateFN);
    fprintf('Saved pdata to %s\n', dacacheFN);
else
    load(stateFN);
    load(dacacheFN);
end

if ~isempty(fsic(varargin,'phase')) && ~isempty(fsic(varargin,'rep')) && ~isempty(fsic(varargin,'trial'))
    modeSubset=1;
    s_phase=varargin{fsic(varargin,'phase')+1};
    s_repNum=varargin{fsic(varargin,'rep')+1};
    s_trialNum=varargin{fsic(varargin,'trial')+1};
    
    listIdx=NaN;
    for i1=1:numel(trialList.rep)
        if isequal(trialList.phase{i1},s_phase) && trialList.rep{i1}==s_repNum && trialList.trialN(i1)==s_trialNum
            listIdx=i1;
            break;
        end
    end
    
    if isnan(listIdx);
        fprintf('ERROR: rep #%d, trial #%d not found in trialList.\n',s_repNum,s_trialNum);
        return
    end
else
    modeSubset=0;
end

%%
figure('Unit','Normalized','Position',[0.1,0.3,0.8,0.5]);

if modeSubset==0    
    ns=state.idx:numel(state.trialList.fn);
elseif modeSubset==1
    ns=listIdx;
end

persist_rmsThresh = NaN;

bFirstTime=1;
for i1=ns
    this_utter=struct;
    this_utter.phase='main';
    if iscell(state.trialList.block(i1))
        this_utter.blockNum=state.trialList.block{i1};
    else
        this_utter.blockNum=state.trialList.block(i1);
    end
    this_utter.trialNum=state.trialList.trialN(i1);

    if state.trialList.isRand(i1)==1
        dataFld = 'randData';
    else
        dataFld = 'sustData';
    end
    
    load(getRawFN_(rawDataDir,state.trialList.fn{i1}));	% gives data
    dataOrig=data;

%     data=reprocAPSTVData(dataOrig);

    sigIn=data.signalIn;
    fs=data.params.sr;

    this_utter.rmsThresh_orig=data.params.rmsThresh;
    this_utter.rmsThresh=data.params.rmsThresh;
    this_utter.fn1_orig=data.params.fn1;
    this_utter.fn1=data.params.fn1;
    this_utter.fn2_orig=data.params.fn2;
    this_utter.fn2=data.params.fn2;
    this_utter.aFact_orig=data.params.aFact;
    this_utter.aFact=data.params.aFact;
    this_utter.bFact_orig=data.params.bFact;
    this_utter.bFact=data.params.bFact;
    this_utter.gFact_orig=data.params.gFact;
    this_utter.gFact=data.params.gFact;
    this_utter.bCepsLift_orig = data.params.bCepsLift;
    this_utter.bCepsLift = data.params.bCepsLift;
    
    this_utter.nLPC = data.params.nLPC;
    
    if bFirstTime==1
        data=reprocData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
                'aFact',this_utter.aFact,'bFact',this_utter.bFact, 'gFact', this_utter.gFact, 'nLPC', this_utter.nLPC, ...
                'bCepsLift', this_utter.bCepsLift);
        bFirstTime=0;
    end
    
    if ~isnan(persist_rmsThresh)
        this_utter.rmsThresh = persist_rmsThresh;
    end
    
    data=reprocData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
                    'aFact',this_utter.aFact,'bFact',this_utter.bFact,'gFact',this_utter.gFact, 'nLPC', this_utter.nLPC, ...
                    'bCepsLift', this_utter.bCepsLift);
    
    f1v=data.fmts(:,1);
    f2v=data.fmts(:,2);
    sf1v=data.sfmts(:,1);
    sf2v=data.sfmts(:,2);
    sigRMS=data.rms(:,1);
    f1v=mva_nz(f1v,mvaWinWidth,'Hamming');
    f2v=mva_nz(f2v,mvaWinWidth,'Hamming');
    
%     f2v0=f2v;
    frameDur=data.params.frameLen/data.params.sr;
    taxis1=0:(frameDur):(frameDur*(length(f1v)-1));

    [j1,j2,foo1,foo2,iv1,iv2,bPossibleMultiProd]=getFmtPlotBounds(f1v,f2v);
    if bPossibleMultiProd == 1
        fprintf('WARNING: there are possibly multiple productions in this trials\n\tThe way to check: change rmsThresh to a very high value (e.g., 0.1) temporarily.\n')
    end
    if ~isempty(iv1)
        this_utter.iv1=iv1;
    else
        this_utter.iv1=NaN;
    end
    if ~isempty(iv2)
        this_utter.iv2=iv2;
    else
        this_utter.iv2=NaN;
    end
    if (~isempty(iv1) && ~isempty(iv2) && ~isnan(iv1) && ~isnan(iv2))
        this_utter.traj_F1=f1v(iv1:iv2);
        this_utter.traj_F2=f2v(iv1:iv2);
        this_utter.sigRMS=sigRMS(iv1:iv2);
    else
        this_utter.traj_F1=[];
        this_utter.traj_F2=[];
        this_utter.sigRMS=[];
    end

%         if isEMMA
    h1=subplot('Position',[0.1,0.3,0.85,0.65]);
    h2=subplot('Position',[0.1,0.125,0.85,0.175]);
    
%         else
%             h1=subplot('Position',[0.1,0.1,0.85,0.85]);
%         end

    if isequal(state.trialList.phase{i1},'ramp') || isequal(state.trialList.phase{i1},'stay') || isequal(state.trialList.phase{i1},'stay2')
        pertStr='pert';
    else
        pertStr='none';
    end
%     pertStr=trialList.pert{i1};

    if isempty(find(f1v>0))
        this_utter.pertStr=pertStr;
        this_utter.rawDataFN=state.trialList.fn{i1};
        this_utter.bDiscard=1;
        this_utter.rating=0;
        this_utter.comments='';
        this_utter.prodF1=NaN;
        this_utter.prodF2=NaN;
        this_utter.audF1=NaN;
        this_utter.audF2=NaN;

        idx_trial=state.trialList.allOrderN(i1);
%         if isequal(dataFld,'randData')
%             t_allOrderN=trialList.allOrderN(find(trialList.isTrain==1));
%             idx_trial=t_allOrderN(i1);
%         elseif isequal(dataFld,'sustData')
%             t_allOrderN=trialList.allOrderN(find(trialList.isTest==1));
%             idx_trial=t_allOrderN(i1);
%         elseif isequal(dataFld,'preData')
%             t_allOrderN=trialList.allOrderN(find(trialList.isPre==1));
%             idx_trial=t_allOrderN(i1);
%         end
        
        if ~isequal(data.params.name,pdata.(dataFld).words{idx_trial})
            fprintf('WARNING: word mismatch!\n');
        end
        pdata.(dataFld).rmsThresh(idx_trial)=this_utter.rmsThresh;
        pdata.(dataFld).fn1(idx_trial)=this_utter.fn1;
        pdata.(dataFld).fn2(idx_trial)=this_utter.fn2;
        pdata.(dataFld).aFact(idx_trial)=this_utter.aFact;
        pdata.(dataFld).bFact(idx_trial)=this_utter.bFact;
        pdata.(dataFld).gFact(idx_trial)=this_utter.gFact;
        pdata.(dataFld).bCepsLift(idx_trial)=this_utter.bCepsLift;
        
        pdata.(dataFld).prodF1(idx_trial)=this_utter.prodF1;
        pdata.(dataFld).prodF2(idx_trial)=this_utter.prodF2;
        pdata.(dataFld).audF1(idx_trial)=this_utter.audF1;
        pdata.(dataFld).audF2(idx_trial)=this_utter.audF2;
        pdata.(dataFld).traj_F1{idx_trial}=[];
        pdata.(dataFld).traj_F2{idx_trial}=[];
        pdata.(dataFld).sigRMS{idx_trial}=[];
        pdata.(dataFld).iv1(idx_trial)=NaN;
        pdata.(dataFld).iv2(idx_trial)=NaN;
        pdata.(dataFld).bDiscard(idx_trial)=this_utter.bDiscard;
        pdata.(dataFld).rating(idx_trial)=this_utter.rating;
        pdata.(dataFld).comments{idx_trial}=this_utter.comments;

        set(gcf,'CurrentAxes',h1);
        title(sprintf('%s - Rep #%d - Trial #%d: %s: parsing failed --> discarded',strrep(this_utter.rawDataFN,'\','/'),...
            this_utter.blockNum,this_utter.trialNum,this_utter.pertStr),...
            'FontWeight','Bold');
        coord=ginput(1);
        
        state.idx=state.idx+1;
        save(stateFN,'state');
        save(dacacheFN,'pdata');

        continue;
    end

%     if ~isempty(find(data.sentStat==6))
%         xlim=[taxis1(min(find(data.sentStat==1)))-0.2,taxis1(min(find(data.sentStat==6)))+0.8];
%     else
    xlim=[taxis1(j1),taxis1(j2)];
%     end

    this_utter.pertStr=pertStr;
    this_utter.rawDataFN=state.trialList.fn{i1};
    this_utter.bDiscard=0;
    
%     idx_v1=round(iv1+0.4*(iv2-iv1));
%     idx_v2=round(iv1+0.6*(iv2-iv1));
    t_rms=data.rms(iv1:iv2);
    t_rms=mva_nz(t_rms,mvaWinWidth,'Hamming');
    [max_rms,idx_max_rms]=max(t_rms);
    rms_lb=max_rms*rmsLBRatio;
    idx_v1=max([iv1,round(iv1+idx_max_rms-1-rmsOneSide/frameDur)]);
    idx_v2=min([iv2,round(iv1+idx_max_rms-1+rmsOneSide/frameDur)]);
    idx_lb_v1=iv1+min(find(t_rms>rms_lb))-1;
    idx_lb_v2=iv1+max(find(t_rms>rms_lb))-1;
    this_utter.prodF1=nanmean(f1v(idx_v1:idx_v2));
    this_utter.prodF2=nanmean(f2v(idx_v1:idx_v2));
    this_utter.prodF1_LB=nanmean(f1v(idx_lb_v1:idx_lb_v2));
    this_utter.prodF2_LB=nanmean(f2v(idx_lb_v1:idx_lb_v2));
    
    % -- Shira's method --
    idx_max_rms = idx_max_rms + iv1 - 1;
    t_rms = mva_nz(data.rms(:, 1), mvaWinWidth, 'Hamming');
    rms_lb_shira = max_rms * 0.4;
    for k2 = idx_max_rms : -1 : 1
        if t_rms(k2) < rms_lb_shira
            break;
        end
    end
    idx_shira_v1 = k2;
        
    for k2 = idx_max_rms : 1 : length(t_rms)
        if t_rms(k2) < rms_lb_shira
            break;
        end
    end
    idx_shira_v2 = k2;
    
    idx_shira_mid = (idx_shira_v1 + idx_shira_v2) / 2;
    idx_shira_v1 = round(idx_shira_mid - shiraOneSide / frameDur);
    idx_shira_v2 = round(idx_shira_mid + shiraOneSide / frameDur);
    
    t_f1v = f1v(idx_shira_v1 : idx_shira_v2);
    t_f2v = f2v(idx_shira_v1 : idx_shira_v2);    
    this_utter.prodF1_shira = nanmean(t_f1v(t_f1v > 0));
    this_utter.prodF2_shira = nanmean(t_f2v(t_f2v > 0));
    this_utter.rms_lb_shira = rms_lb_shira;
    this_utter.idx_shira_v1 = idx_shira_v1;
    this_utter.idx_shira_v2 = idx_shira_v2;
    % -- ~Shira's method --
    
    this_utter.idx_mnlBound_1 = NaN;
    this_utter.idx_mnlBound_2 = NaN;
    this_utter.prodF1_mnlBound = NaN;
    this_utter.prodF2_mnlBound = NaN;
    
    
    
    this_utter.idx_v1=idx_v1;
    this_utter.idx_v2=idx_v2;
    this_utter.idx_lb_v1=idx_lb_v1;
    this_utter.idx_lb_v2=idx_lb_v2;
    if isequal(pertStr,'none')
        this_utter.audF1=this_utter.prodF1;
        this_utter.audF2=this_utter.prodF2;
    else
        this_utter.audF1=nanmean(sf1v(idx_v1:idx_v2));
        this_utter.audF2=nanmean(sf2v(idx_v1:idx_v2));
    end
    
    this_utter.vowelOnset = NaN;
    this_utter.vowelEnd = NaN;
    this_utter.vowelOnsetIdx = NaN;
    this_utter.vowelEndIdx = NaN;
    this_utter.f1Traj = [];
    this_utter.f2Traj = [];
    
    this_utter.rating=3;
    this_utter.comments='';
    
    t_title = sprintf('Trial #%d / %d - %s (solid: rms LB; dashed: rms-peak-window; dotted: Shira''s method) (randomized order)', ...
                      i1,numel(state.trialList.fn),state.trialList.word{i1});
        
    toRepeat=1;
    while toRepeat
        data=reprocData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
            'aFact',this_utter.aFact,'bFact',this_utter.bFact,'gFact',this_utter.gFact, ...
            'bCepsLift', this_utter.bCepsLift);

        f1v=data.fmts(:,1);
        f2v=data.fmts(:,2);
        f1v=mva_nz(f1v,mvaWinWidth,'Hamming');
        f2v=mva_nz(f2v,mvaWinWidth,'Hamming');
        sf1v=data.sfmts(:,1);
        sf2v=data.sfmts(:,2);
        sigRMS=data.rms(:,1);
        [j1,j2,foo1,foo2,iv1,iv2]=getFmtPlotBounds(f1v,f2v);
        if ~isempty(iv1)
            this_utter.iv1=iv1;
        else
            this_utter.iv1=NaN;
        end
        if ~isempty(iv2)
            this_utter.iv2=iv2;
        else
            this_utter.iv2=NaN;
        end
        if (~isempty(iv1) && ~isempty(iv2) && ~isnan(iv1) && ~isnan(iv2))
            this_utter.traj_F1=f1v(iv1:iv2);
            this_utter.traj_F2=f2v(iv1:iv2);
            this_utter.sigRMS=sigRMS(iv1:iv2);
        else
            this_utter.traj_F1=[];
            this_utter.traj_F2=[];
            this_utter.sigRMS=[];
        end

        xlim=[taxis1(j1),taxis1(j2)];
        
        set(gcf,'CurrentAxes',h1);
        cla;
        [s,f,t]=spectrogram(sigIn,128,96,1024,fs);
        imagesc(t,f,10*log10(abs(s))); hold on;
        axis xy;
        hold on;

        plot(taxis1,f1v,'w-','LineWidth',1);
        hold on;
        plot(taxis1,f2v,'w-','LineWidth',1);

        set(gca,'YLim',ylim,'XLim',xlim);
%         xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        
%         idx_v1=round(iv1+0.4*(iv2-iv1));
%         idx_v2=round(iv1+0.6*(iv2-iv1));
        t_rms=data.rms(iv1:iv2);
        t_rms=mva_nz(t_rms,mvaWinWidth,'Hamming');
        [max_rms,idx_max_rms]=max(t_rms);
        rms_lb=max_rms*rmsLBRatio;
        idx_v1=max([iv1,round(iv1+idx_max_rms-1-rmsOneSide/frameDur)]);
        idx_v2=min([iv2,round(iv1+idx_max_rms-1+rmsOneSide/frameDur)]);
        idx_lb_v1=iv1+min(find(t_rms>rms_lb))-1;
        idx_lb_v2=iv1+max(find(t_rms>rms_lb))-1;
        this_utter.prodF1=nanmean(f1v(idx_v1:idx_v2));
        this_utter.prodF2=nanmean(f2v(idx_v1:idx_v2));
        this_utter.prodF1_LB=nanmean(f1v(idx_lb_v1:idx_lb_v2));
        this_utter.prodF2_LB=nanmean(f2v(idx_lb_v1:idx_lb_v2));
        this_utter.idx_v1=idx_v1;
        this_utter.idx_v2=idx_v2;
        this_utter.idx_lb_v1=idx_lb_v1;
        this_utter.idx_lb_v2=idx_lb_v2;
        
        % -- Shira's method --
        idx_max_rms = idx_max_rms + iv1 - 1;
        t_rms = mva_nz(data.rms(:, 1), mvaWinWidth, 'Hamming');
        rms_lb_shira = max_rms * 0.4;
        for k2 = idx_max_rms : -1 : 1
            if t_rms(k2) < rms_lb_shira
                break;
            end
        end
        idx_shira_v1 = k2;

        for k2 = idx_max_rms : 1 : length(t_rms)
            if t_rms(k2) < rms_lb_shira
                break;
            end
        end
        idx_shira_v2 = k2;
        
        idx_shira_mid = (idx_shira_v1 + idx_shira_v2) / 2;
        idx_shira_v1 = round(idx_shira_mid - shiraOneSide / frameDur);
        idx_shira_v2 = round(idx_shira_mid + shiraOneSide / frameDur);

        t_f1v = f1v(idx_shira_v1 : idx_shira_v2);
        t_f2v = f2v(idx_shira_v1 : idx_shira_v2);    
        this_utter.prodF1_shira = nanmean(t_f1v(t_f1v > 0));
        this_utter.prodF2_shira = nanmean(t_f2v(t_f2v > 0));
        this_utter.rms_lb_shira = rms_lb_shira;
        this_utter.idx_shira_v1 = idx_shira_v1;
        this_utter.idx_shira_v2 = idx_shira_v2;
        % -- ~Shira's method --
        
        if isequal(pertStr,'none')
            this_utter.audF1=this_utter.prodF1;
            this_utter.audF2=this_utter.prodF2;
        else
            this_utter.audF1=nanmean(sf1v(idx_v1:idx_v2));
            this_utter.audF2=nanmean(sf2v(idx_v1:idx_v2));
        end
        
        if ~isempty(idx_v1) && ~isempty(idx_v2)
            plot(repmat(taxis1(idx_v1),1,2),ylim,'k--');
            plot(repmat(taxis1(idx_v2),1,2),ylim,'k--');
        end
        if ~isempty(idx_lb_v1) && ~isempty(idx_lb_v2)
            plot(repmat(taxis1(idx_lb_v1),1,2),ylim,'k-');
            plot(repmat(taxis1(idx_lb_v2),1,2),ylim,'k-');
        end
        if ~isempty(idx_shira_v1) && ~isempty(idx_shira_v2)
            plot(repmat(taxis1(idx_shira_v1),1,50),linspace(ylim(1),ylim(2),50),'k.');
            plot(repmat(taxis1(idx_shira_v2),1,50),linspace(ylim(1),ylim(2),50),'k.');
        end
        if ~isnan(this_utter.idx_mnlBound_1) && ~isnan(this_utter.idx_mnlBound_1)
            plot(repmat(taxis1(this_utter.idx_mnlBound_1),1,50),linspace(ylim(1),ylim(2),50),'g-');
            plot(repmat(taxis1(this_utter.idx_mnlBound_2),1,50),linspace(ylim(1),ylim(2),50),'g-');
        end

        title(t_title);
        
        
        set(gcf,'CurrentAxes',h2);
        cla;
        plot(taxis1,data.rms(:,1),'b-','LineWidth',1);        
        hold on;
        set(gca,'XLim',xlim);
        ys=get(gca,'YLim');
        if ~isempty(idx_v1) && ~isempty(idx_v2)
            plot(repmat(taxis1(idx_v1),1,2),ys,'k--');
            plot(repmat(taxis1(idx_v2),1,2),ys,'k--');
        end
        if ~isempty(idx_lb_v1) && ~isempty(idx_lb_v2)
            plot(repmat(taxis1(idx_lb_v1),1,2),ys,'k-');
            plot(repmat(taxis1(idx_lb_v2),1,2),ys,'k-');
        end
        if ~isempty(idx_shira_v1) && ~isempty(idx_shira_v2)
            plot(repmat(taxis1(idx_shira_v1),1,15),linspace(ys(1),ys(2),15),'k.');
            plot(repmat(taxis1(idx_shira_v2),1,15),linspace(ys(1),ys(2),15),'k.');
        end
        if ~isnan(this_utter.idx_mnlBound_1) && ~isempty(this_utter.idx_mnlBound_2)
            plot(repmat(taxis1(this_utter.idx_mnlBound_1),1,15),linspace(ys(1),ys(2),15),'g-');
            plot(repmat(taxis1(this_utter.idx_mnlBound_2),1,15),linspace(ys(1),ys(2),15),'g-');
        end
        xlabel('Time (s)');
        set(gcf,'CurrentAxes',h1);
        
        
        % Manually set the beginning and end of the vowel
        lblOkay = ~isnan(this_utter.vowelOnset) & ~isnan(this_utter.vowelEnd) & ...
                  this_utter.vowelEnd > this_utter.vowelOnset;
        while ~lblOkay
            title('Set vowel onset...');
            coord = ginput(1);
            this_utter.vowelOnset = coord(1);

            title('Set vowel end...');
            coord = ginput(1);
            this_utter.vowelEnd = coord(1);
            
            lblOkay = ~isnan(this_utter.vowelOnset) & ~isnan(this_utter.vowelEnd) & ...
                  this_utter.vowelEnd > this_utter.vowelOnset;
        end
        [foo, this_utter.vowelOnsetIdx] = min(abs(taxis1 - this_utter.vowelOnset));
        [foo, this_utter.vowelEndIdx] = min(abs(taxis1 - this_utter.vowelEnd));
        
        this_utter.f1Traj = f1v(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx);
        this_utter.f2Traj = f2v(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx);
        
        ylm = get(gca, 'YLim');
        plot(repmat(this_utter.vowelOnset, 1, 2), ylm, 'w--', 'LineWidth', 1.5);
        plot(repmat(this_utter.vowelEnd, 1, 2), ylm, 'w-', 'LineWidth', 1.5);
        % ~Manually set the beginning and end of the vowel
        
        title(t_title);

        % Buttons
        set(gcf,'CurrentAxes',h1);
        xs=get(gca,'XLim'); ys=get(gca,'YLim');        
        x1=xs(1); x2=xs(2); y1=ys(1); y2=ys(2);
        rx=range(xs); ry=range(ys);        

        cmd.items={'Play sigIn', 'Play sigOut', 'rmsThresh',...
                   'fn1', 'fn2', 'aFact', 'bFact', 'gFact', ...
                   'bCepsLift', ...
                   'Manual bounds','Remove mnl. bounds','Relabel','Rating','Discard utter','Save'};
        cmd.x_left=repmat(x1+0.9*rx,1,length(cmd.items));
        cmd.y_bottom=y2-0.05*ry-0.06*(1:length(cmd.items))*ry;
        cmd.x_width=repmat(0.1*rx,1,length(cmd.items));
        cmd.y_height=repmat(0.06*ry,1,length(cmd.items));

        for i3=1:length(cmd.items)            
            rectangle('Position',[cmd.x_left(i3),cmd.y_bottom(i3),cmd.x_width(i3),cmd.y_height(i3)],'FaceColor','w','EdgeColor','k');
            text(cmd.x_left(i3)+0.01*rx,cmd.y_bottom(i3)+0.025*ry,cmd.items{i3},'Color','k');
        end
        % ~Buttons

        coord=ginput(1);
        inBox=0;
        for i3=1:length(cmd.items)
            if (coord(1)>=cmd.x_left(i3) && coord(1)<=cmd.x_left(i3)+cmd.x_width(i3) && coord(2)>=cmd.y_bottom(i3) && coord(2)<=cmd.y_bottom(i3)+cmd.y_height(i3))
                inBox=i3;
                break;
            end
        end

        if inBox==0
            xs=get(gca,'XLim'); ys=get(gca,'YLim');
            if ((coord(1)<xs(1) || coord(1)>xs(2)) || (coord(2)<ys(1) || coord(2)>ys(2)))
                toRepeat=0;
            end
        else
            if isequal(cmd.items{inBox},'Play sigIn')
                soundsc(dataOrig.signalIn,data.params.sr);
            elseif isequal(cmd.items{inBox},'Play sigOut')
                soundsc(dataOrig.signalOut,data.params.sr);
            elseif isequal(cmd.items{inBox},'rmsThresh')
                this_utter.rmsThresh=input(sprintf('[rmsThresh_orig = %.4f; rmsThresh_curr = %.4f] rmsThresh = ',...
                    this_utter.rmsThresh_orig,this_utter.rmsThresh));
                
                persist_rmsThresh = this_utter.rmsThresh;
            elseif isequal(cmd.items{inBox},'fn1')
                this_utter.fn1=input(sprintf('[fn1_orig = %.1f; fn1_curr = %.1f] fn1 = ',...
                    this_utter.fn1_orig,this_utter.fn1));
            elseif isequal(cmd.items{inBox},'fn2')
                this_utter.fn2=input(sprintf('[fn2_orig = %.4f; fn2_curr = %.4f] fn2 = ',...
                    this_utter.fn2_orig,this_utter.fn2));
            elseif isequal(cmd.items{inBox},'aFact')
                this_utter.aFact=input(sprintf('[aFact_orig = %.4f; aFact_curr = %.4f] aFact = ',...
                    this_utter.aFact_orig,this_utter.aFact));
            elseif isequal(cmd.items{inBox},'bFact')
                this_utter.bFact=input(sprintf('[bFact_orig = %.4f; bFact_curr = %.4f] bFact = ',...
                    this_utter.bFact_orig,this_utter.bFact));
            elseif isequal(cmd.items{inBox},'gFact')
                this_utter.gFact=input(sprintf('[gFact_orig = %.4f; gFact_curr = %.4f] gFact = ',...
                    this_utter.gFact_orig,this_utter.gFact));
            elseif isequal(cmd.items{inBox}, 'bCepsLift')
                this_utter.bCepsLift=input(sprintf('[bCepsLift = %d; bCepsLift_curr = %d] bCepsLift = ',...
                    this_utter.bCepsLift_orig, this_utter.bCepsLift));
            elseif isequal(cmd.items{inBox}, 'Manual bounds')
                mnl_t1 = ginput(1); mnl_t1 = mnl_t1(1);
                mnl_t2 = ginput(1); mnl_t2 = mnl_t2(1);
                this_utter.idx_mnlBound_1 = round(mnl_t1 / frameDur);
                this_utter.idx_mnlBound_2 = round(mnl_t2 / frameDur);
                this_utter.prodF1_mnlBound = nanmean(f1v(this_utter.idx_mnlBound_1 : this_utter.idx_mnlBound_2));
                this_utter.prodF2_mnlBound = nanmean(f2v(this_utter.idx_mnlBound_1 : this_utter.idx_mnlBound_2));
            elseif isequal(cmd.items{inBox}, 'Remove mnl. bounds')
                if isnan(this_utter.idx_mnlBound_1) && isnan(this_utter.idx_mnlBound_2) && isnan(this_utter.prodF1_mnlBound) && isnan(this_utter.prodF2_mnlBound)
                    fprintf('ERROR: there are no manual bounds to remove.\n');
                else
                    this_utter.idx_mnlBound_1 = NaN;
                    this_utter.idx_mnlBound_1 = NaN;
                    this_utter.prodF1_mnlBound = NaN;
                    this_utter.prodF2_mnlBound = NaN;
                    fprintf('Maunal bounds removed\n');
                end
            elseif isequal(cmd.items{inBox}, 'Relabel')
                this_utter.vowelOnset = NaN;
                this_utter.vowelEnd = NaN;
                this_utter.vowelOnsetIdx = NaN;
                this_utter.vowelEndIdx = NaN;
                this_utter.f1Traj = [];
                this_utter.f2Traj = [];
            elseif isequal(cmd.items{inBox},'Rating')
                    this_utter.rating=input(['[Old rating = ',num2str(this_utter.rating),'] New rating = ']);
                    this_utter.comments=input(['[e.g., (i1) = problematic i1.] Rating comments = '],'s');
            elseif isequal(cmd.items{inBox},'Discard utter')
                if this_utter.bDiscard==0
                    title('Utter discarded','Color','r','FontSize',13,'FontWeight','Bold');
                    this_utter.bDiscard=1;
                else
                    title('Utter un-discarded','Color','g','FontSize',13,'FontWeight','Bold');
                    this_utter.bDiscard=0;
                end                        
                pause(0.5);
            elseif isequal(cmd.items{inBox},'Save')
                save(dacacheFN,'pdata');
                fprintf('%s saved\n',dacacheFN);                 
            end
            
        end
    end

    bContinue=1;

    if bContinue==1
        this_utter.timeStamp_analysis=clock;
        
        idx_trial=state.trialList.allOrderN(i1);
        
        if ~isequal(state.trialList.fn{i1}, pdata.(dataFld).rawDataFNs{idx_trial})
            fprintf('ERROR: raw data file name mismatch. \n');
            return
        end
        
        if ~isequal(data.params.name,pdata.(dataFld).words{idx_trial})
            fprintf('WARNING: word mismatch!\n');
        end
        pdata.(dataFld).rmsThresh(idx_trial)=this_utter.rmsThresh;
        pdata.(dataFld).fn1(idx_trial)=this_utter.fn1;
        pdata.(dataFld).fn2(idx_trial)=this_utter.fn2;
        pdata.(dataFld).aFact(idx_trial)=this_utter.aFact;
        pdata.(dataFld).bFact(idx_trial)=this_utter.bFact;
        pdata.(dataFld).gFact(idx_trial)=this_utter.gFact;
        pdata.(dataFld).bCepsLift(idx_trial)=this_utter.bCepsLift;
        
        pdata.(dataFld).vowelOnset(idx_trial) = this_utter.vowelOnset;
        pdata.(dataFld).vowelEnd(idx_trial) = this_utter.vowelEnd;
        pdata.(dataFld).vowelOnsetIdx(idx_trial) = this_utter.vowelOnsetIdx;
        pdata.(dataFld).vowelEndIdx(idx_trial) = this_utter.vowelEndIdx;
        pdata.(dataFld).f1Traj{idx_trial} = this_utter.f1Traj;
        pdata.(dataFld).f2Traj{idx_trial} = this_utter.f2Traj;

        pdata.(dataFld).prodF1(idx_trial)=this_utter.prodF1;
        pdata.(dataFld).prodF2(idx_trial)=this_utter.prodF2;
        pdata.(dataFld).prodF1_LB(idx_trial)=this_utter.prodF1_LB;
        pdata.(dataFld).prodF2_LB(idx_trial)=this_utter.prodF2_LB;
        
        pdata.(dataFld).prodF1_shira(idx_trial)=this_utter.prodF1_shira;
        pdata.(dataFld).prodF2_shira(idx_trial)=this_utter.prodF2_shira;
        pdata.(dataFld).prodF1_mnlBound(idx_trial)=this_utter.prodF1_mnlBound;
        pdata.(dataFld).prodF2_mnlBound(idx_trial)=this_utter.prodF2_mnlBound;
        
        pdata.(dataFld).audF1(idx_trial)=this_utter.audF1;
        pdata.(dataFld).audF2(idx_trial)=this_utter.audF2;
        pdata.(dataFld).traj_F1{idx_trial}=this_utter.traj_F1;
        pdata.(dataFld).traj_F2{idx_trial}=this_utter.traj_F2;
        pdata.(dataFld).sigRMS{idx_trial}=this_utter.sigRMS;
        pdata.(dataFld).iv1(idx_trial)=this_utter.iv1;
        pdata.(dataFld).iv2(idx_trial)=this_utter.iv2;
        pdata.(dataFld).bDiscard(idx_trial)=this_utter.bDiscard;
        pdata.(dataFld).rating(idx_trial)=this_utter.rating;
        pdata.(dataFld).comments{idx_trial}=this_utter.comments;
    else
        fprintf('Terminated.\n');
        return
    end
    
    if modeSubset==0
        state.idx=state.idx+1;
        save(stateFN,'state');
    end
    save(dacacheFN,'pdata');
end

if modeSubset==1
    return
end

%% Screening the perturbations
if isfield(state,'pertScrIdx')
    fprintf('Found pertScrIdx in state.\n');
    a=input('Resume perturbation screening? (0/1): ');
    if a==1
        bNew=0;
    else
        bNew=1;
    end
else
    bNew=1;
end
if bNew
    state.pertScrIdx=1;
    
    save(stateFN,'state');
end
% end

%%
ns=state.pertScrIdx:numel(state.trialListPert.fn);
for i1=ns
    load(getRawFN_(rawDataDir,state.trialListPert.fn{i1}));	% gives data
    sigIn=data.signalIn;
    sigOut=data.signalOut;
    fs=data.params.sr;
    f1v=data.fmts(:,1);
    f2v=data.fmts(:,2);
    f1v=mva_nz(f1v,mvaWinWidth,'Hamming');
    f2v=mva_nz(f2v,mvaWinWidth,'Hamming');
    frameDur=data.params.frameLen/data.params.sr;
    taxis1=0:(frameDur):(frameDur*(length(f1v)-1));
    [j1,j2]=getFmtPlotBounds(f1v,f2v);
    h1=subplot('Position',[0.1,0.1,0.85,0.825]);
    set(gcf,'CurrentAxes',h1);
    
    idx1=fsic(state.trialList.fn,state.trialListPert.fn{i1});
    
    f2s=data.sfmts(:,2);
    f1s=data.sfmts(:,1);
    
    cla;
    [s,f,t]=spectrogram(sigOut,128,96,1024,fs);
    imagesc(t,f,10*log10(abs(s))); hold on;
    axis xy;
    hold on;
    plot(taxis1,f1v,'w-','LineWidth',1);
    hold on;
    plot(taxis1,f2v,'w-','LineWidth',1);
    plot(taxis1,f2s,'g-','LineWidth',1); hold on;
    plot(taxis1,f1s,'g-','LineWidth',1); hold on;
    
    xlim=[taxis1(j1),taxis1(j2)];
    
    set(gca,'YLim',ylim,'XLim',xlim);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    a=input(sprintf('Is perturbation okay? (0/1) [Default: 1]: '),'s');
    
    if isempty(a)
        a=1;
    end
    if ~isnumeric(a)
        a=str2num(a);
    end
    
    idx_trial=trialListPert.allOrderN(i1);
    
    if ~isequal(trialListPert.fn{i1},pdata.randData.rawDataFNs{idx_trial})
        fprintf('ERROR: raw data file name mismatch. \n');
        return
    end
    
    pdata.randData.bPertOkay(idx_trial)=a;
    state.pertScrIdx=state.pertScrIdx+1;
    
    save(stateFN,'state');
    save(dacacheFN,'pdata');
end

return

%%
function raw_fn=getRawFN_(expDir,fn)
[path1,fn1]=fileparts(fn);
[path2,fn2]=fileparts(path1);
[path3,fn3]=fileparts(path2);
[path4,fn4]=fileparts(path3);

raw_fn=fullfile(expDir,fn4,fn3,fn2,fn1);
if ~isequal(raw_fn(end-3:end),'.mat')
    raw_fn=[raw_fn,'.mat'];
end