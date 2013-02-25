function pdata1 = updateDataUI(uihdls, pdata, dataFld, idx_trial, state, i1, varargin)
%% Data analysis configurations
mvaWinWidth=21;     % 21 * 1.333 = 28 (ms)
fineParseWin=100e-3;	% sec

ylim=[0,5000];
rmsOneSide=0.01;  % Unit: sec
rmsLBRatio=0.75;
shiraOneSide=0.025; % Unit: sec

%% -- Main screening process --
this_utter=struct;
this_utter.phase = state.trialList.phase{i1};
if iscell(state.trialList.block(i1))
    this_utter.blockNum=state.trialList.block{i1};
else
    this_utter.blockNum=state.trialList.block(i1);
end
% this_utter.trialNum=state.trialList.trialN(i1);

[ret, hostName]=system('hostname');
rawfn = getRawFN_(state.rawDataDir,state.trialList.fn{i1});
if ~isequal(lower(deblank(hostName)), 'smcg_w510')
    rawfn = strrep(rawfn, 'E:', 'D:');
end

load(rawfn);	% gives data
dataOrig=data;

sigIn=data.signalIn;
fs=data.params.sr;

if ~isempty(fsic(varargin, 'fromList'))
    if get(uihdls.rb_alwaysPlaySigIn, 'Value') == 1
        soundsc(sigIn, fs);
    end
end

this_utter.rmsThresh_orig = data.params.rmsThresh;
this_utter.fn1_orig = data.params.fn1;
this_utter.fn2_orig = data.params.fn2;
this_utter.aFact_orig = data.params.aFact;
this_utter.bFact_orig = data.params.bFact;
this_utter.gFact_orig = data.params.gFact;
this_utter.bCepsLift_orig = data.params.bCepsLift;
this_utter.cepsWinWidth_orig = data.params.cepsWinWidth;
this_utter.nLPC_orig = data.params.nLPC;

if isnan(pdata.(dataFld).rmsThresh(idx_trial))
    this_utter.rmsThresh = data.params.rmsThresh;
    this_utter.nLPC = data.params.nLPC;
    this_utter.fn1=data.params.fn1;
    this_utter.fn2=data.params.fn2;
    this_utter.aFact=data.params.aFact;
    this_utter.bFact=data.params.bFact;
    this_utter.gFact=data.params.gFact;
    this_utter.bCepsLift = data.params.bCepsLift;
    this_utter.cepsWinWidth = data.params.cepsWinWidth;
    
else
    this_utter.rmsThresh = str2double(get(uihdls.edit_rmsThresh, 'String'));
    this_utter.nLPC = str2num(get(uihdls.edit_nLPC, 'String'));
    this_utter.fn1 = str2double(get(uihdls.edit_fn1, 'String'));
    this_utter.fn2 = str2double(get(uihdls.edit_fn2, 'String'));
    this_utter.aFact = str2double(get(uihdls.edit_aFact, 'String'));
    this_utter.bFact = str2double(get(uihdls.edit_bFact, 'String'));
    this_utter.gFact = str2double(get(uihdls.edit_gFact, 'String'));
    this_utter.bCepsLift = str2num(get(uihdls.edit_bCepsLift, 'String'));
    this_utter.cepsWinWidth = str2num(get(uihdls.edit_cepsWinWidth, 'String'));
end

if state.bFirstTime==1
    data=reprocData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
            'aFact',this_utter.aFact,'bFact',this_utter.bFact, 'gFact', this_utter.gFact, 'nLPC', this_utter.nLPC, ...
            'bCepsLift', this_utter.bCepsLift, 'cepsWinWidth', this_utter.cepsWinWidth);
    state.bFirstTime=0;
end

if ~isnan(state.persist_rmsThresh)
    this_utter.rmsThresh = persist_rmsThresh;
end

data=reprocData(dataOrig,'rmsThresh',this_utter.rmsThresh,'fn1',this_utter.fn1,'fn2',this_utter.fn2,...
                'aFact',this_utter.aFact,'bFact',this_utter.bFact,'gFact',this_utter.gFact, 'nLPC', this_utter.nLPC, ...
                'bCepsLift', this_utter.bCepsLift, 'cepsWinWidth', this_utter.cepsWinWidth);

f1v=data.fmts(:,1);
f2v=data.fmts(:,2);
sf1v=data.sfmts(:,1);
sf2v=data.sfmts(:,2);
sigRMS=data.rms(:,1);
f1v=mva_nz(f1v, mvaWinWidth, 'Hamming');
f2v=mva_nz(f2v, mvaWinWidth, 'Hamming');
            
%% 
frameDur=data.params.frameLen/data.params.sr;
taxis1=0:(frameDur):(frameDur*(length(f1v)-1));

[j1,j2,foo1,foo2,iv1,iv2,bPossibleMultiProd]=getFmtPlotBounds(f1v,f2v);
if bPossibleMultiProd == 1
    fprintf('WARNING: there are possibly multiple productions in this trials\n\tThe way to check: play sigIn or change rmsThresh to a very high value (e.g., 0.1) temporarily.\n')
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
% h1=subplot('Position',[0.1,0.3,0.85,0.65]);
% h2=subplot('Position',[0.1,0.125,0.85,0.175]);
h1 = uihdls.haxes1;
h2 = uihdls.haxes2;


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
    pdata.(dataFld).cepsWinWidth(idx_trial)=this_utter.cepsWinWidth;
    pdata.(dataFld).nLPC(idx_trial)=this_utter.nLPC;

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

    pdata1 = pdata;

    
    return;
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

% ------------ Compute the three F1/F2 averages ------------ %
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
% ------------ ~Compute the three F1/F2 averages ------------ %

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

this_utter.rating = pdata.(dataFld).rating(idx_trial);
this_utter.comments = pdata.(dataFld).comments{idx_trial};

t_title = sprintf('Trial #%d / %d - %s (solid: rms LB; dashed: rms-peak-window; dotted: Shira''s method) (randomized order)', ...
                  i1,numel(state.trialList.fn),state.trialList.word{i1});
              
%%

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


plot(repmat(taxis1(idx_max_rms), 1, 2), ylim, 'k-.');
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
    
    t_xlim = get(gca, 'XLim');
    t_ylim = get(gca, 'YLim');
    plot(t_xlim, repmat(this_utter.prodF1_shira, 1, 2), '--', 'Color', [0.5, 0.5, 0.5]);
    plot(t_xlim, repmat(this_utter.prodF2_shira, 1, 2), '--', 'Color', [0.5, 0.5, 0.5]);
    text(t_xlim(1) + 0.01 * range(t_xlim), this_utter.prodF1_shira + 0.012 * range(t_ylim), ...
         sprintf('F1(shira) = %.1f Hz', this_utter.prodF1_shira), ...
         'Color', 'w');
    text(t_xlim(1) + 0.01 * range(t_xlim), this_utter.prodF2_shira + 0.012 * range(t_ylim), ...
         sprintf('F2(shira) = %.1f Hz', this_utter.prodF2_shira), ...
         'Color', 'w');
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
if isnan(pdata.(dataFld).vowelOnset(idx_trial)) || ...
   isnan(pdata.(dataFld).vowelEnd(idx_trial)) || ...
   isnan(pdata.(dataFld).vowelOnsetIdx(idx_trial)) || ...
   isnan(pdata.(dataFld).vowelEndIdx(idx_trial)) || ...
   pdata.(dataFld).vowelOnset(idx_trial) >= pdata.(dataFld).vowelEnd(idx_trial) || ...
   pdata.(dataFld).vowelOnsetIdx(idx_trial) >= pdata.(dataFld).vowelEndIdx(idx_trial)
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
else
    this_utter.vowelOnsetIdx = pdata.(dataFld).vowelOnsetIdx(idx_trial);
    this_utter.vowelEndIdx = pdata.(dataFld).vowelEndIdx(idx_trial);
    this_utter.vowelOnset = pdata.(dataFld).vowelOnset(idx_trial);
    this_utter.vowelEnd = pdata.(dataFld).vowelEnd(idx_trial);
    
    
end

% Check the validity of Shira bounds
bShiraValid = check_shira_bounds(this_utter);
if bShiraValid == 0
    msgbox('Shira methods bounds are not within the vowel bounds.', 'Vowel bounds error', 'error', 'modal');
end

nzero = numel(find(data.fmts(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx, 1) == 0));
if nzero > 0
    fprintf('WARNING: %d data points between vowelOnsetIdx and vowelEndIdx contain zero formant values.\n', nzero)
end

this_utter.f1Traj = f1v(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx);
this_utter.f2Traj = f2v(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx);

ylm = get(gca, 'YLim');
plot(repmat(this_utter.vowelOnset, 1, 2), ylm, 'w--', 'LineWidth', 1.5);
plot(repmat(this_utter.vowelEnd, 1, 2), ylm, 'w-', 'LineWidth', 1.5);
% ~Manually set the beginning and end of the vowel

%% ------------ Re-compute the three F1/F2 averages ------------ %
if this_utter.vowelOnsetIdx < this_utter.vowelEndIdx
    t_rms = data.rms(this_utter.vowelOnsetIdx : this_utter.vowelEndIdx);
    t_rms = mva_nz(t_rms, mvaWinWidth, 'Hamming');
    [max_rms, idx_max_rms] = max(t_rms);
    rms_lb = max_rms * rmsLBRatio;
    idx_v1 = max([1, round(this_utter.vowelOnsetIdx + idx_max_rms-1-rmsOneSide/frameDur)]);
    idx_v2 = min([length(t_rms), round(this_utter.vowelOnsetIdx+idx_max_rms-1+rmsOneSide/frameDur)]);
    idx_lb_v1 = this_utter.vowelOnsetIdx + min(find(t_rms>rms_lb))-1;
    idx_lb_v2 = this_utter.vowelOnsetIdx + max(find(t_rms>rms_lb))-1;
    this_utter.prodF1 = nanmean(f1v(idx_v1 : idx_v2));
    this_utter.prodF2 = nanmean(f2v(idx_v1 : idx_v2));
    this_utter.prodF1_LB = nanmean(f1v(idx_lb_v1 : idx_lb_v2));
    this_utter.prodF2_LB = nanmean(f2v(idx_lb_v1 : idx_lb_v2));

    % -- Shira's method --
    idx_max_rms = idx_max_rms + this_utter.vowelOnsetIdx - 1;
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
    % ------------ ~Compute the three F1/F2 averages ------------ %
end

title(t_title);

%%
this_utter.timeStamp_analysis=clock;

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
pdata.(dataFld).cepsWinWidth(idx_trial)=this_utter.cepsWinWidth;
pdata.(dataFld).nLPC(idx_trial)=this_utter.nLPC;

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

pdata1 = pdata;
return

%%
function bValid = check_shira_bounds(utter)
bValid = (utter.idx_shira_v1 > utter.vowelOnsetIdx) & ...
         (utter.idx_shira_v2 < utter.vowelEndIdx);

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