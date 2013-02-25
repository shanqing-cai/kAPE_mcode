function preprocKapeData(subjID,varargin)
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

if ~isempty(fsic(varargin, 'msu'))
    dacacheDir = [dacacheDir, '_msu'];
elseif ~isempty(fsic(varargin, 'boston')) || ~isempty(fsic(varargin, 'bu'))
    dacacheDir = [dacacheDir, '_boston'];
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
otherPhase = {'other'};
otherData = init_data('other', otherPhase, rawDataDir, subjID);

randPhase = {'rand'};
randData = init_data('rand', randPhase, rawDataDir, subjID);

sustPhases = {'start', 'ramp', 'stay', 'end'};
sustData = init_data('sust', sustPhases, rawDataDir, subjID);

pdata.otherData = otherData;
pdata.randData = randData;
pdata.sustData = sustData;

%% Build a list of all trials
if isfile(stateFN)
    fprintf('Found state.mat at %s.\n',stateFN);
    if isempty(fsic(varargin, 'phase')) && isempty(fsic(varargin, 'pdata_fixAutoRMS'))
        a = input('Resume? (0/1): ');
    else
        a = 1;
    end
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
    
    idxRandPerm = randperm(numel(otherData.rawDataFNs));
    
    trialListOther.fn=otherData.rawDataFNs(idxRandPerm);
    trialListOther.phase=otherData.phases(idxRandPerm);
    trialListOther.block=otherData.blockNums(idxRandPerm);
    trialListOther.trialN=otherData.trialNums(idxRandPerm);
    trialListOther.word=otherData.words(idxRandPerm);
    idxOrder_o = 1:numel(otherData.trialNums);
    trialListOther.allOrderN=idxOrder_o(idxRandPerm);

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
    nOtherTrials = numel(trialListOther.fn);
    nRandTrials = numel(trialListRand.fn);
    nSustTrials = numel(trialListSust.fn);
    
    % Combine parts
    flds=fields(trialList);
    for k1=1:numel(flds)
        fld=flds{k1};
        trialList.(fld)=[trialListOther.(fld), trialListRand.(fld), trialListSust.(fld)];        
    end
    
    trialList.isOther = [ones(1, nOtherTrials), zeros(1,nRandTrials), zeros(1, nSustTrials)];
    trialList.isRand = [zeros(1, nOtherTrials), ones(1,nRandTrials), zeros(1, nSustTrials)];
    trialList.isSust = [zeros(1, nOtherTrials), zeros(1,nRandTrials), ones(1, nSustTrials)];
    
    state.trialList = trialList;
    state.trialListPert = trialListPert;
    
    state.rawDataDir = rawDataDir;
    state.dacacheDir = dacacheDir;
    state.expDir = expDir;
    
    state.persist_rmsThresh = NaN;
    state.bFirstTime = 1;
    
    state.stats = zeros(1, numel(state.trialList.fn));
    state.statsPert = zeros(1, numel(state.trialListPert.fn));   
    
    save(stateFN,'state');
    save(dacacheFN,'pdata');
    
    fprintf('Saved state info to %s\n', stateFN);
    fprintf('Saved pdata to %s\n', dacacheFN);
else
    load(stateFN);
    load(dacacheFN);
    
    state.bFirstTime = 1;
end

%%
uihdls = struct;
uihdls.dacacheFN = dacacheFN;

uihdls.hfig = figure('Position', [50, 150, 1500, 600]);
hlist_title = uicontrol('Style', 'text', ...
                  'Unit', 'Normalized', ...
                  'Position', [0.02, 0.93, 0.18, 0.05], ...
                  'String', 'Trial list: (*: Vowel bounds done)', ...
                  'HorizontalAlignment', 'left');
uihdls.hlist_title = hlist_title;

hlist = uicontrol('Style', 'listbox', ...
                  'Unit', 'Normalized', ...
                  'Position', [0.02, 0.15, 0.18, 0.8], ...
                  'BackgroundColor', [1, 1, 1]);
uihdls.hlist = hlist;

uihdls.hmenu_nLPC = uimenu('Parent', uihdls.hfig, 'Label', 'nLPC');
set(uihdls.hfig, 'MenuBar', 'none');
uihdls.hmenu_nLPC_show_overall_best = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Show overall best');
uihdls.hmenu_nLPC_set_overall_best = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Set all trials to overall best');
uihdls.hmenu_nLPC_set_list_1st = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Set all trials to 1st in list', 'Separator', 'on');
uihdls.hmenu_nLPC_restore_user = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Restore user selections', 'Separator', 'on');                                   
                                    
uihdls.hmenu_comments = uimenu('Parent', uihdls.hfig, 'Label', 'Comments');
uihdls.hmenu_comments_recover = uimenu(uihdls.hmenu_comments, ...
                                       'Label', 'Recover from file...');
                                   
uihdls.hmenu_rmsThresh = uimenu('Parent', uihdls.hfig, 'Label', 'rmsThresh');
uihdls.hmenu_rmsThresh_scan = uimenu('Parent', uihdls.hmenu_rmsThresh, ...
                                     'Label', 'Scan for trials with gaps');

hreveal = uicontrol('Style', 'pushbutton', ...
                    'Unit', 'Normalized', ...
                    'Position', [0.02, 0.04, 0.18, 0.04], ...
                    'String', 'Reveal trial details');
uihdls.hreveal = hreveal;

hShowComments = uicontrol('Style', 'pushbutton', ...
                          'Unit', 'Normalized', ...
                          'Position', [0.02, 0.09, 0.18, 0.04], ...
                          'String', 'Show comments');
uihdls.hShowComments = hShowComments;

% htitle = uicontrol('Style', 'text', ...
%                    'Unit', 'Normalized', ...
%                    'Position', [0.24, 0.93, 0.3, 0.04], ...
%                    'String', 'Title', ...
%                    'FontSize', 12);
% uihdls.htitle = htitle;

haxes1 = axes('Unit', 'Normalized', 'Position', [0.24, 0.26, 0.60, 0.68]);
uihdls.haxes1 = haxes1;

haxes2 = axes('Unit', 'Normalized', 'Position', [0.24, 0.13, 0.60, 0.15]);
uihdls.haxes2 = haxes2;

% Zoom buttons
hzo = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.24, 0.025, 0.10, 0.045], ...
                'String', 'Zoom out');
uihdls.hzo = hzo;

hzi = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.34, 0.025, 0.10, 0.045], ...
                'String', 'Zoom in');
uihdls.hzi = hzi;

hpleft = uicontrol('Style', 'pushbutton', ...
                   'Unit', 'Normalized', ...
                   'Position', [0.46, 0.025, 0.10, 0.045], ...
                   'String', 'Pan left');
uihdls.hpleft = hpleft;

hpright = uicontrol('Style', 'pushbutton', ...
                   'Unit', 'Normalized', ...
                   'Position', [0.56, 0.025, 0.10, 0.045], ...
                   'String', 'Pan right');
uihdls.hpright = hpright;

hzd = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.68, 0.025, 0.10, 0.045], ...
                'String', 'Default zoom');
uihdls.hzd = hzd;

bt_playSigIn = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [0.87, 0.88, 0.07, 0.04], ...
                         'String', 'Play sigIn');
uihdls.bt_playSigIn = bt_playSigIn;
bt_playSigOut = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [0.87, 0.83, 0.07, 0.04], ...
                         'String', 'Play sigOut');
uihdls.bt_playSigOut = bt_playSigOut;

rb_alwaysPlaySigIn = uicontrol('Style', 'radiobutton', ...
                               'Unit', 'Normalized', ...
                               'Position', [0.945, 0.88, 0.055, 0.04], ...
                               'String', 'Always play');
uihdls.rb_alwaysPlaySigIn = rb_alwaysPlaySigIn;

lblLeft = 0.845;
lblWidth = 0.04;
editLeft = 0.89;
editWidth = 0.04;
lbl_rmsThresh = uicontrol('Style', 'text', ...
                         'Unit', 'Normalized', ...
                         'Position', [lblLeft, 0.78, lblWidth, 0.04], ...
                         'String', 'rmsThresh: ');
uihdls.lbl_rmsThresh = lbl_rmsThresh;
edit_rmsThresh = uicontrol('Style', 'edit', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft, 0.78, editWidth, 0.04], ...
                         'String', 'rmsThresh', 'HorizontalAlignment', 'left');
uihdls.edit_rmsThresh = edit_rmsThresh;

bt_auto_rmsThresh = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [editLeft + editWidth + 0.004, 0.78, editWidth * 1.65, 0.04], ...
                             'String', 'Auto rmsThresh', 'FontSize', 8);
uihdls.bt_auto_rmsThresh = bt_auto_rmsThresh;

bt_auto_rmsThresh_all = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [editLeft + editWidth + 0.004, 0.73, editWidth * 1.65, 0.04], ...
                             'String', 'Auto rmsThresh all', 'FontSize', 8);
uihdls.bt_auto_rmsThresh_all = bt_auto_rmsThresh_all;

                     
lbl_nLPC = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.72, lblWidth, 0.04], ...
                     'String', 'nLPC: ');
uihdls.lbl_nLPC = lbl_nLPC;
edit_nLPC = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.72, editWidth, 0.04], ...
                      'String', 'nLPC', 'HorizontalAlignment', 'left');
uihdls.edit_nLPC = edit_nLPC;


bt_auto_nLPC = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.65, editWidth * 1.5, 0.04], ...
                         'String', 'Auto nLPC', ...
                         'FontSize', 8);
uihdls.bt_auto_nLPC = bt_auto_nLPC;

lbl_srt_nLPCs = uicontrol('Style', 'text', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.59, editWidth * 1.5, 0.04], ...
                         'String', 'Sorted nLPCs: ', 'HorizontalAlignment', 'left', ...
                         'FontSize', 8);
uihdls.lbl_srt_nLPCs = lbl_srt_nLPCs; 

lst_srt_nLPCs = uicontrol('Style', 'listbox', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.44, editWidth * 1.5, 0.16], ...
                         'String', {}, 'Enable', 'off', ...
                         'FontSize', 8);
uihdls.lst_srt_nLPCs = lst_srt_nLPCs;

bt_auto_nLPC_all = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.38, editWidth * 1.5, 0.04], ...
                         'String', 'Auto nLPC all', ...
                         'FontSize', 8);
uihdls.bt_auto_nLPC_all = bt_auto_nLPC_all;

% bt_best_nLPC = uicontrol('Style', 'pushbutton', ...
%                          'Unit', 'Normalized', ...
%                          'Position', [editLeft + editWidth + 0.004, 0.32, editWidth * 1.5, 0.04], ...
%                          'String', 'Overall best nLPC', ...
%                          'FontSize', 8);
% uihdls.bt_best_nLPC = bt_best_nLPC;

lbl_fn1 = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.67, lblWidth, 0.04], ...
                     'String', 'fn1: ');
uihdls.lbl_fn1 = lbl_fn1;
edit_fn1 = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.67, editWidth, 0.04], ...
                      'String', 'fn1', 'HorizontalAlignment', 'left');
uihdls.edit_fn1 = edit_fn1;
                  
lbl_fn2 = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.62, lblWidth, 0.04], ...
                     'String', 'fn2: ');
uihdls.lbl_fn2 = lbl_fn2;
edit_fn2 = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.62, editWidth, 0.04], ...
                      'String', 'fn2', 'HorizontalAlignment', 'left');
uihdls.edit_fn2 = edit_fn2;
                  
lbl_aFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.56, lblWidth, 0.04], ...
                     'String', 'aFact: ');
uihdls.lbl_aFact = lbl_aFact;
edit_aFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.56, editWidth, 0.04], ...
                      'String', 'aFact', 'HorizontalAlignment', 'left');
uihdls.edit_aFact = edit_aFact;

lbl_bFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.51, lblWidth, 0.04], ...
                     'String', 'bFact: ');
uihdls.lbl_bFact = lbl_bFact;
edit_bFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.51, editWidth, 0.04], ...
                      'String', 'bFact', 'HorizontalAlignment', 'left');
uihdls.edit_bFact = edit_bFact;
                  
lbl_gFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.45, lblWidth, 0.04], ...
                     'String', 'gFact: ');
uihdls.lbl_gFact = lbl_gFact;
edit_gFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.45, editWidth, 0.04], ...
                      'String', 'gFact', 'HorizontalAlignment', 'left');
uihdls.edit_gFact = edit_gFact;

lbl_bCepsLift = uicontrol('Style', 'text', ...
                          'Unit', 'Normalized', ...
                          'Position', [lblLeft, 0.38, lblWidth, 0.04], ...
                          'String', 'bCepsLift: ');
uihdls.lbl_bCepsLift = lbl_bCepsLift;
edit_bCepsLift = uicontrol('Style', 'edit', ...
                           'Unit', 'Normalized', ...
                           'Position', [editLeft, 0.38, editWidth, 0.04], ...
                           'String', 'bCepsLift', 'HorizontalAlignment', 'left');
uihdls.edit_bCepsLift = edit_bCepsLift;

lbl_cepsWinWidth = uicontrol('Style', 'text', ...
                          'Unit', 'Normalized', ...
                          'Position', [lblLeft, 0.32, lblWidth, 0.04], ...
                          'String', 'cepsWinWidth: ');
uihdls.lbl_cepsWinWidth = lbl_cepsWinWidth;
edit_cepsWinWidth = uicontrol('Style', 'edit', ...
                           'Unit', 'Normalized', ...
                           'Position', [editLeft, 0.32, editWidth, 0.04], ...
                           'String', 'cepsWinWidth', 'HorizontalAlignment', 'left');
uihdls.edit_cepsWinWidth = edit_cepsWinWidth;
                  
bt_reproc = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft, 0.27, 0.10, 0.04], ...
                      'String', 'Reprocess');
uihdls.bt_reproc = bt_reproc;

bt_relabel = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft, 0.96, 0.10, 0.035], ...
                      'String', 'ReLabel');
uihdls.bt_relabel = bt_relabel;

bt_relabel_focus = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [lblLeft, 0.925, 0.10, 0.035], ...
                             'String', 'Focus and ReLabel');
uihdls.bt_relabel_focus = bt_relabel_focus;
                  
lbl_rating = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.22, lblWidth, 0.04], ...
                     'String', 'Rating: ');
uihdls.lbl_rating = lbl_rating;
pm_rating = uicontrol('Style', 'popupmenu', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.22, lblWidth, 0.04], ...
                      'String', {'0', '1', '2'}, ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');
uihdls.pm_rating = pm_rating;
                  
lbl_comments = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.18, lblWidth, 0.03], ...
                     'String', 'Comments: ');
uihdls.lbl_comments = lbl_comments;
edit_comments = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft + 0.01, 0.14, 0.14, 0.03], ...
                      'String', 'Comments', ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');                  
uihdls.edit_comments = edit_comments;

lbl_f1mark = uicontrol('Style', 'text', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.86, 0.07, 0.03, 0.03], ...
                      'String', 'F1', ...
                      'HorizontalAlignment', 'Center');

bt_f1_beg = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.89, 0.07, 0.03, 0.03], ...
                      'String', 'beg', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f1_beg = bt_f1_beg;

bt_f1_mid = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.92, 0.07, 0.03, 0.03], ...
                      'String', 'mid', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f1_mid = bt_f1_mid;

bt_f1_end = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.95, 0.07, 0.03, 0.03], ...
                      'String', 'end', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f1_end = bt_f1_end;

lbl_f2mark = uicontrol('Style', 'text', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.86, 0.10, 0.03, 0.03], ...
                      'String', 'F2', ...
                      'HorizontalAlignment', 'Center');

bt_f2_beg = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.89, 0.10, 0.03, 0.03], ...
                      'String', 'beg', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f2_beg = bt_f2_beg;

bt_f2_mid = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.92, 0.10, 0.03, 0.03], ...
                      'String', 'mid', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f2_mid = bt_f2_mid;

bt_f2_end = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [0.95, 0.10, 0.03, 0.03], ...
                      'String', 'end', ...
                      'HorizontalAlignment', 'left');
uihdls.bt_f2_end = bt_f2_end;
                  
bt_next = uicontrol('Style', 'pushbutton', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft + 0.01, 0.01, lblWidth * 2, 0.05], ...
                     'String', 'Next', 'FontSize', 11);
uihdls.bt_next = bt_next;

%% Set upcall back functions
set(uihdls.hlist, 'Callback', {@list_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_playSigIn, 'Callback', {@playSig_cbk, dacacheFN, stateFN, uihdls, 'in'});
set(uihdls.bt_playSigOut, 'Callback', {@playSig_cbk, dacacheFN, stateFN, uihdls, 'out'});
set(uihdls.bt_reproc, 'Callback', {@reproc_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_relabel, 'Callback', {@relabel_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_relabel_focus, 'Callback', {@relabel_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.pm_rating, 'Callback', {@rating_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.edit_comments, 'Callback', {@comments_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_auto_rmsThresh, 'Callback', {@auto_rmsThresh_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_auto_rmsThresh_all, 'Callback', {@auto_rmsThresh_all_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_auto_nLPC, 'Callback', {@auto_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_auto_nLPC_all, 'Callback', {@auto_nLPC_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_nLPC_show_overall_best, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_set_overall_best, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_set_list_1st, 'Callback', {@set_nLPC_list_1st, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_restore_user, 'Callback', {@restore_user_nLPC, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_rmsThresh_scan, 'Callback', {@rmsThresh_scan_cbk, dacacheFN, stateFN, uihdls});
% set(uihdls.bt_best_nLPC, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_comments_recover, 'Callback', {@recover_comments_from_file, dacacheFN, stateFN, uihdls});

set(uihdls.bt_next, 'Callback', {@next_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hreveal, 'Callback', {@reveal_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hShowComments, 'Callback', {@showComments_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_f1_beg, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f1', 'beg'});
set(uihdls.bt_f1_mid, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f1', 'mid'});
set(uihdls.bt_f1_end, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f1', 'end'});
set(uihdls.bt_f2_beg, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f2', 'beg'});
set(uihdls.bt_f2_mid, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f2', 'mid'});
set(uihdls.bt_f2_end, 'Callback', {@mark_cbk, dacacheFN, stateFN, uihdls, 'f2', 'end'});

set(uihdls.lst_srt_nLPCs, 'Callback', {@lst_srt_nLPCs_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_reproc, 'Enable', 'off');
set(uihdls.bt_auto_rmsThresh, 'Enable', 'off');

set(uihdls.hzo, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hzi, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hpleft, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hpright, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hzd, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});

updateTrialList(state, uihdls);

%% Optional: pre-set trial file name.

if ~isempty(fsic(varargin,'phase')) && ~isempty(fsic(varargin,'rep')) && ~isempty(fsic(varargin,'trial'))
    s_phase=varargin{fsic(varargin,'phase')+1};
    s_repNum=varargin{fsic(varargin,'rep')+1};
    s_trialNum=varargin{fsic(varargin,'trial')+1};
    
    updateTrialList(state, uihdls, s_phase, s_repNum, s_trialNum);    
end


%% Optional: pdata_fixAutoRMS
if ~isempty(fsic(varargin, 'pdata_fixAutoRMS'))
    auto_rmsThresh_all_cbk([], [], dacacheFN, stateFN, uihdls, 'noConfirm');
end

%% First show of the spectrogram view
% reproc_cbk(uihdls.bt_reproc, [], dacacheFN, stateFN, uihdls);
list_cbk([], [], dacacheFN, stateFN, uihdls)
return

%%
function rmsThresh_scan_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
lst_str = get(uihdls.hlist, 'String');
load(stateFN); % Gives state
load(dacacheFN); % Gives pdata

a_gapTrialNs = [];
for i1 = 1 : numel(state.trialList.word)
    if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
        dataFld = 'otherData';
    elseif state.trialList.isRand(i1) == 1
        dataFld = 'randData';
    elseif state.trialList.isSust(i1) == 1
        dataFld = 'sustData';
    end
    
    idx_trial = state.trialList.allOrderN(i1);    
    
    if pdata.(dataFld).bDiscard(idx_trial) == 1 || pdata.(dataFld).rating(idx_trial) == 0
        continue;
    end
    
    vowelOnsetIdx = pdata.(dataFld).vowelOnsetIdx(idx_trial);
    vowelEndIdx = pdata.(dataFld).vowelEndIdx(idx_trial);
    if isnan(vowelOnsetIdx) || isnan(vowelEndIdx)
        fprintf(1, 'ERROR: Vowel onset / offset labels have not been set on trial #%d. Aborted.\n', i1);
        return
    end
    
    if ~isempty(find(pdata.(dataFld).f1Traj{idx_trial} == 0));
        a_gapTrialNs(end + 1) = i1;
%         fprintf(1, '\tTrial #%d\n', i1);
    end
end

if length(a_gapTrialNs) == 0
    fprintf(1, 'INFO: Found no trials with gaps in formant trajectories. Good.\n');
else
    fprintf(1, 'INFO: Trials with gaps in formant trajectories (due to incorrect rmsThresh values): \n');
    for i1 = 1 : length(a_gapTrialNs)
        fprintf(1, '\tTrial #%d\n', a_gapTrialNs(i1));
    end
    fprintf(1, 'INFO: Please re-run "auto rmsThresh" and "auto nLPC" on these trials, and then fix the nLPC choices manually.\n');
end
return

%%
function restore_user_nLPC(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % Gives pdata

flds = {'otherData', 'randData', 'sustData'};
bFoundMissing = 0;
a_nLPC_lst = [];
for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if length(pdata.(fld).rawDataFNs) == 0
        continue;
    end
    
    if ~isfield(pdata.(fld), 'srt_nLPCs')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).srt_nLPCs)
        if pdata.(fld).bDiscard(h2)
            continue;
        end

        if pdata.(fld).rating(h2) == 0
            continue;
        end

        if isempty(pdata.(fld).srt_nLPCs{h2})
            bFoundMissing = 1;
            break;
        else
            a_nLPC_lst = [a_nLPC_lst; pdata.(fld).srt_nLPCs{h2}];
        end
    end            
end

if bFoundMissing
    fprintf(1, 'WARNING: This action cannot be taken until auto nLPC has not been run on all trials.\n');
    return
end

if ~isfield(pdata, 'nLPC_status')
    fprintf(1, 'WARNING: User nLPCs have not been overwritten by best nLPC or list-first nLPCs yet.\nCannot perform this restoring at this moment.\n');
    return
else
    if isequal(pdata.nLPC_status, 'user')
        fprintf(1, 'WARNING: User nLPCs have not been overwritten by best nLPC or list-first nLPCs yet.\nCannot perform this restoring at this moment.\n');        
        return
    end
end

fields = {'randData', 'sustData'};
for i1 = 1 : numel(fields)
    fld = fields{i1};

    pdata.(fld).nLPC = pdata.(fld).user_nLPCs; % Restore user selections
end    

pdata.nLPC_status = 'user';
save(dacacheFN, 'pdata');
fprintf(1, 'Restored user nLPCs. \npdata saved to %s\n', dacacheFN);

set(uihdls.hlist, 'Enable', 'off');
lst_str = get(uihdls.hlist, 'String');
for i1 = 1 : numel(lst_str)
    set(uihdls.hlist, 'Value', i1);
    list_cbk(uihdls.hlist, [], dacacheFN, stateFN, uihdls);
    set(uihdls.hlist, 'Enable', 'off');
    drawnow;
end
set(uihdls.hlist, 'Enable', 'on');

fprintf(1, 'Restored user nLPCs. \npdata saved to %s\n', dacacheFN);
return

%%
function set_nLPC_list_1st(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % Gives pdata

flds = {'otherData', 'randData', 'sustData'};
bFoundMissing = 0;
a_nLPC_lst = [];
for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if length(pdata.(fld).rawDataFNs) == 0
        continue;
    end
    
    if ~isfield(pdata.(fld), 'srt_nLPCs')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).srt_nLPCs)
        if pdata.(fld).bDiscard(h2)
            continue;
        end

        if pdata.(fld).rating(h2) == 0
            continue;
        end

        if isempty(pdata.(fld).srt_nLPCs{h2})
            bFoundMissing = 1;
            break;
        else
            a_nLPC_lst = [a_nLPC_lst; pdata.(fld).srt_nLPCs{h2}];
        end
    end            
end

if bFoundMissing
    fprintf(2, 'ERROR: Setting nLPC to 1st in list cannot be done until auto nLPC has not been run on all trials.\n');
    return
end
    
if ~isfield(pdata, 'nLPC_status')
    pdata.nLPC_status = 'user';
else
    if isequal(pdata.nLPC_status, 'list_1st')
        fprintf(1, 'Data are already set to list-first nLPCs. No changes will be made.\n');
        return;
    end
end

fields = {'randData', 'sustData'};
for i1 = 1 : numel(fields)
    fld = fields{i1};

    if isequal(pdata.nLPC_status, 'user')
        pdata.(fld).user_nLPCs = pdata.(fld).nLPC; % Backup user selections
    end

    for i2 = 1 : numel(pdata.(fld).nLPC)
        if ~isnan(pdata.(fld).nLPC)
            if isempty(pdata.(fld).srt_nLPCs{i2})
                continue;
            end
            pdata.(fld).nLPC(i2) = pdata.(fld).srt_nLPCs{i2}(1);
        end
    end
end    

pdata.nLPC_status = 'list_1st';

save(dacacheFN, 'pdata');
fprintf(1, 'Set all nLPCs to list-first. \npdata saved to %s\n', dacacheFN);

set(uihdls.hlist, 'Enable', 'off');
lst_str = get(uihdls.hlist, 'String');
for i1 = 1 : numel(lst_str)
    set(uihdls.hlist, 'Value', i1);
    list_cbk(uihdls.hlist, [], dacacheFN, stateFN, uihdls);
    set(uihdls.hlist, 'Enable', 'off');
    drawnow;
end
set(uihdls.hlist, 'Enable', 'on');
fprintf(1, 'Set all nLPCs to list-first. \npdata saved to %s\n', dacacheFN);

return

%% 
function best_nLPC_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % Gives pdata

flds = {'otherData', 'randData', 'sustData'};
bFoundMissing = 0;
a_nLPC_lst = [];
for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if length(pdata.(fld).rawDataFNs) == 0
        continue;
    end
    
    if ~isfield(pdata.(fld), 'srt_nLPCs')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).srt_nLPCs)
        if pdata.(fld).bDiscard(h2)
            continue;
        end

        if pdata.(fld).rating(h2) == 0
            continue;
        end

        if isempty(pdata.(fld).srt_nLPCs{h2})
            bFoundMissing = 1;
            break;
        else
            a_nLPC_lst = [a_nLPC_lst; pdata.(fld).srt_nLPCs{h2}];
        end
    end            
end

if bFoundMissing
    fprintf(2, 'ERROR: Overall best nLPC cannot be determined until auto nLPC has not been run on all trials.\n');
    return
end

fprintf('\n');
a_nLPCs = unique(a_nLPC_lst(:, 1));
a_nLPCs = sort(a_nLPCs);
for i1 = 1 : numel(a_nLPCs)
    fprintf('nLPC = %d:\tbest in %d of %d trials.\n', ...
            a_nLPCs(i1),  length(find(a_nLPC_lst(:, 1) == a_nLPCs(i1))), ...
            length(a_nLPC_lst(:, 1)));
end

fprintf('======================================\n');
fprintf('Overall best nLPC = %d\n\n', mode(a_nLPC_lst(:, 1)))

if hObject == uihdls.hmenu_nLPC_set_overall_best
    if ~isfield(pdata, 'nLPC_status')
        pdata.nLPC_status = 'user';
    else
        if isequal(pdata.nLPC_status, 'overall_best');
            fprintf(1, 'WARNING: Data are already set to overall best nLPC. No changes will be made.\n');
            return;
        end
    end
    
    fields = {'randData', 'sustData'};
    for i1 = 1 : numel(fields)
        fld = fields{i1};

        if isequal(pdata.nLPC_status, 'user')
            pdata.(fld).user_nLPCs = pdata.(fld).nLPC; % Backup user selections
        end

        for i2 = 1 : numel(pdata.(fld).nLPC)
            if isempty(pdata.(fld).srt_nLPCs{i2})
                continue;
            end
            
            if ~isnan(pdata.(fld).nLPC) 
                pdata.(fld).nLPC(i2) = mode(a_nLPC_lst(:, 1));
            end
        end
    end    
    
    pdata.nLPC_status = 'overall_best';
    save(dacacheFN, 'pdata');
    fprintf(1, 'Set all nLPCs to overall best (%d). \npdata saved to %s\n', ...
            mode(a_nLPC_lst(:, 1)), dacacheFN);
        
    set(uihdls.hlist, 'Enable', 'off');
    lst_str = get(uihdls.hlist, 'String');
    for i1 = 1 : numel(lst_str)
        set(uihdls.hlist, 'Value', i1);
        list_cbk(uihdls.hlist, [], dacacheFN, stateFN, uihdls);
        set(uihdls.hlist, 'Enable', 'off');
        drawnow;
    end
    set(uihdls.hlist, 'Enable', 'on');
    
    fprintf(1, 'Set all nLPCs to overall best (%d). \npdata saved to %s\n', ...
            mode(a_nLPC_lst(:, 1)), dacacheFN);
    
end


return

%%
function zoom_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
zdat = guidata(uihdls.haxes1);

xlim = [NaN, NaN];
if hObject == uihdls.hzo % Zoom out
    xlim(1) = zdat.currXLim(1) - (zdat.currXLim(1) - zdat.tmin) * 0.1;
    xlim(2) = zdat.currXLim(2) + (zdat.tmax - zdat.currXLim(2)) * 0.1;
elseif hObject == uihdls.hzi % Zoom in
    xlim(1) = zdat.currXLim(1) + 0.1 * range(zdat.currXLim);
    xlim(2) = zdat.currXLim(2) - 0.1 * range(zdat.currXLim);
elseif hObject == uihdls.hpleft % Pan left
    xlim(1) = zdat.currXLim(1) - 0.1 * range(zdat.currXLim);    
    if xlim(1) < zdat.tmin
        xlim(1) = zdat.tmin;
    end
    xlim(2) = xlim(1) + range(zdat.currXLim);
elseif hObject == uihdls.hpright
    xlim(2) = zdat.currXLim(2) + 0.1 * range(zdat.currXLim);    
    if xlim(2) > zdat.tmax
        xlim(2) = zdat.tmax;
    end
    xlim(1) = xlim(2) - range(zdat.currXLim);
else % Default zoom
    xlim = zdat.defXLim;    
end

zdat.currXLim = xlim;
guidata(uihdls.haxes1, zdat);

set(gcf, 'CurrentAxes', uihdls.haxes2);
set(gca, 'XLim', xlim);

set(gcf, 'CurrentAxes', uihdls.haxes1);
set(gca, 'XLim', xlim);

return

%%
function lst_srt_nLPCs_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
val = get(uihdls.lst_srt_nLPCs, 'Value');
str = get(uihdls.lst_srt_nLPCs, 'String');
sel_nLPC = str2num(str{val});

set(uihdls.edit_nLPC, 'String', sprintf('%d', sel_nLPC));
reproc_cbk(uihdls.bt_reproc, [], dacacheFN, stateFN, uihdls);

return


function next_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(stateFN);
updateTrialList(state, uihdls, 'next', dacacheFN, stateFN); % gives state

return


function rating_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
    dataFld = 'otherData';
elseif state.trialList.isRand(i1) == 1
    dataFld = 'randData';
elseif state.trialList.isSust(i1) == 1
    dataFld = 'sustData';
end

idx_trial = state.trialList.allOrderN(i1);

val = get(uihdls.pm_rating, 'Value');
items = get(uihdls.pm_rating, 'String');
pdata.(dataFld).rating(idx_trial) = str2num(items{val});
save(dacacheFN, 'pdata');
fprintf('Saved to %s\n', dacacheFN);

fn = state.trialList.fn{i1};
fprintf('INFO: trial: %s: rating -> %d\n', fn, pdata.(dataFld).rating(idx_trial));
return

function comments_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
    dataFld = 'otherData';
elseif state.trialList.isRand(i1) == 1
    dataFld = 'randData';
elseif state.trialList.isSust(i1) == 1
    dataFld = 'sustData';
end


idx_trial = state.trialList.allOrderN(i1);

str = get(uihdls.edit_comments, 'String');
pdata.(dataFld).comments{idx_trial} = str;
save(dacacheFN, 'pdata');

fn = state.trialList.fn{i1};
fprintf('INFO: trial: %s: comments -> [%s]\n', fn, pdata.(dataFld).comments{idx_trial});
return

function playSig_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls, mode)
i1 = get(uihdls.hlist, 'Value');

% load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

hostName=deblank(getHostName);
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    driveLet = 'E:';
else
    driveLet = 'D:';
end

rfn = getRawFN_(state.rawDataDir,state.trialList.fn{i1});
if ~isequal(rfn(1 : 2), driveLet)
    rfn(1 : 2) = driveLet;
end

load(rfn); % gives data

if isequal(mode, 'in')
    soundsc(data.signalIn, data.params.sr);
else
    soundsc(data.signalOut, data.params.sr);
end

return

function reveal_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
revButStr = get(uihdls.hreveal, 'String');
if isequal(revButStr, 'Reveal trial details')
    set(uihdls.hreveal, 'String', 'Hide trial details');
else
    set(uihdls.hreveal, 'String', 'Reveal trial details');
end

load(stateFN);  % gives state;
updateTrialList(state, uihdls)
return

function showComments_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
showCommentsStr = get(uihdls.hShowComments, 'String');
if isequal(showCommentsStr, 'Show comments')
    set(uihdls.hShowComments, 'String', 'Hide comments');
else
    set(uihdls.hShowComments, 'String', 'Show comments');
end

load(stateFN);  % gives state;
updateTrialList(state, uihdls)
return

function mark_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls, fmt, seg)
green = [0, 0.5, 0];
red = [1, 0, 0];

i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
    dataFld = 'otherData';
elseif state.trialList.isRand(i1) == 1
    dataFld = 'randData';
elseif state.trialList.isSust(i1) == 1
    dataFld = 'sustData';
end

idx_trial = state.trialList.allOrderN(i1);

fld = sprintf('mark_%s_%s', fmt, seg);
uifld = sprintf('bt_%s_%s', fmt, seg);

curColor = get(uihdls.(uifld), 'ForegroundColor');

if curColor(1) == 1 % Red
    pdata.(dataFld).(fld)(idx_trial) = 0;
    set(uihdls.(uifld), 'ForegroundColor', green);
else
    pdata.(dataFld).(fld)(idx_trial) = 1;
    set(uihdls.(uifld), 'ForegroundColor', red);
    
    rating_items = get(uihdls.pm_rating, 'String');
    t_val = get(uihdls.pm_rating, 'Value');
    t_val_1 = fsic(rating_items, '1');
    if isequal(rating_items{t_val}, '2')
        set(uihdls.pm_rating, 'Value', t_val_1);
        pdata.(dataFld).rating(idx_trial) = 1;
        fprintf('INFO: Automatically setting rating score to 1\n');
    end
    
end

save(dacacheFN, 'pdata');
fprintf('%s -> %d\n', fld, pdata.(dataFld).(fld)(idx_trial));


return


function auto_rmsThresh_all_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls, varargin)
% Check whether the preproc has finished.
load(stateFN); % gives state
if ~isempty(find(state.stats == 0))
    msgbox('There are some unprocessed trials. Please finish processing all trials before using Auto all', 'Unable to execute Auto ALL');
    return
end

if isempty(fsic(varargin, 'noConfirm'))
    answer = questdlg('Are you sure you want to run auto RMS on all trials of this subject?', 'Confirm Auto RMS all...');
else
    answer = 'Yes';
end

if isequal(answer, 'Yes')
    t_str = get(uihdls.hlist, 'String');    
    for i1 = 1 : numel(t_str)
        set(uihdls.hlist, 'enable', 'off');
        set(uihdls.hlist, 'Value', i1);
        list_cbk([], [], dacacheFN, stateFN, uihdls)
        drawnow;
        
        auto_rmsThresh_cbk([], [], dacacheFN, stateFN, uihdls);
    end
    set(uihdls.hlist, 'enable', 'on');
end
return

%%
function relabel_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
    dataFld = 'otherData';
elseif state.trialList.isRand(i1) == 1
    dataFld = 'randData';
elseif state.trialList.isSust(i1) == 1
    dataFld = 'sustData';
end


idx_trial = state.trialList.allOrderN(i1);

pdata.(dataFld).vowelOnsetIdx(idx_trial) = NaN;
pdata.(dataFld).vowelEndIdx(idx_trial) = NaN;
pdata.(dataFld).vowelOnset(idx_trial) = NaN;
pdata.(dataFld).vowelEnd(idx_trial) = NaN;

save(dacacheFN, 'pdata');

if isequal(uihdls.bt_relabel_focus, hObject) % Focus and relabel
    list_cbk([], [], dacacheFN, stateFN, uihdls, 'focus');
else
    list_cbk([], [], dacacheFN, stateFN, uihdls);
end
return

%%
function recover_comments_from_file(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % gives pdata

[FileName, PathName] = uigetfile('*.mat', 'Select the source .mat file');
if isempty(FileName) || FileName == 0
    fprintf(1, 'recover_comments_from_file: No file selected. No actions taken.\n');
    return
end

pdata_src = load(fullfile(PathName, FileName));

flds = {'otherData', 'randData', 'sustData'};
for i1 = 1 : numel(flds)
    fld = flds{i1};
    nNonEmpty0 = length(pdata.(fld).comments) - length(fsic(pdata.(fld).comments, ''));
    nNonEmpty1 = length(pdata_src.pdata.(fld).comments) - length(fsic(pdata_src.pdata.(fld).comments, ''));
    
    pdata.(fld).comments = pdata_src.pdata.(fld).comments;
    
    fprintf(1, 'INFO: %s: Before: %d non-empty comments; After: %d non-empty comments\n', ...
            fld, nNonEmpty0, nNonEmpty1);
end

save(dacacheFN, 'pdata');

return