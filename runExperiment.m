function varargout = runExperiment(varargin)
DEBUG=0;
[ret,hostName]=system('hostname');

%% ---- Read and parse expt_config.txt ----
configFN = 'expt_config.txt';
expt_config=read_parse_expt_config(configFN);

dlg_fields = {'SUBJECT_ID', 'SUBJECT_GENDER', ...
              'SUBJECT_DOB', 'SUBJECT_GROUP'};
for ii = 1 : numel(dlg_fields)
    dlgf = dlg_fields{ii};
    expt_config.(dlgf) = inputdlg(dlgf, dlgf, 1, {expt_config.(dlgf)});
    if iscell(expt_config.(dlgf))
        expt_config.(dlgf) = expt_config.(dlgf){1};
    end
    expt_config.(dlgf) = deblank(expt_config.(dlgf));
end

%% ---- Subject and experiment information ----
subject.expt_config         =expt_config;
subject.name				=expt_config.SUBJECT_ID; 
subject.sex					=expt_config.SUBJECT_GENDER;  % male / female
% subject.shiftDirectionSust	=expt_config.SHIFT_DIRECTION;  %SC F1Up / F1Down
% subject.shiftRatio			=expt_config.SHIFT_RATIO;    %
subject.age                 =expt_config.SUBJECT_AGE;
subject.group               =expt_config.SUBJECT_GROUP;

subject.mouthMicDist        =expt_config.MOUTH_MIC_DIST;   % cm
subject.closedLoopGain      =0-20*log10(10/subject.mouthMicDist);

subject.dBRange1            =expt_config.SPL_RANGE/0.4;        % is the one-sided dBRange1*0.4
% subject.dBRange2            =12;         % Tightened level range after the initial pract1 training. 

subject.trialLen            =expt_config.TRIAL_LEN;
subject.trialLenMax         =expt_config.TRIAL_LEN_MAX;

subject.hostName            = deblank(hostName);

subject.dataDir             = expt_config.DATA_DIR;
% if isequal(subject.hostName,'smcg-w510') || isequal(subject.hostName,'smcg_w510')
%     subject.dataDir             ='E:\DATA\KIDAPE';
%     subject.percDir             ='E:\DATA\DYS\PERC';
% else
    % subject.dataDir				='C:\CS_2004\PROJECTS\SAP-FMRI\';
%     subject.dataDir				='D:\CS_2004\PROJECTS\DYS';
% end

subject.trigByScanner		=0;
subject.TA					=2.5;
subject.ITI					=6;

subject.vumeterMode         =2;     % 1: 10 ticks; 2: 3 ticks;

if subject.trigByScanner==1
	subject.showProgress		=0;
	subject.showPlayButton      =0;
else
	subject.showProgress		=1;
	subject.showPlayButton      =1;
end

subject.designNum			=2;

subject.lvNoise             =75; % dBA SPL. The level of noise for completely masking speech (mode trialType = 2 or 3).

bAO=~isempty(fsic(varargin,'AO')) | ~isempty(fsic(varargin,'ao'));

bSim = 0;
simDataDir = '';
if ~isempty(fsic(varargin, 'sim'))
    bSim = 1;
    simDataDir = varargin{fsic(varargin, 'sim') + 1};
end

subject.bAO=bAO;
%%
subject.date				=clock;

if (~isempty(findStringInCell(varargin,'subject')))
    clear('subject');
    subject=varargin{findStringInCell(varargin,'subject')+1};
end

% if (~isfield(subject,'pcrKnob'))
%     subject.pcrKnob=input('Phone / Ctrl Room knob = ');
% end
subject.pcrKnob=0;

%%
bNew=true;

dirname=fullfile(subject.dataDir,num2str(subject.name));

if (~isempty(findStringInCell(varargin,'dirname')))
    clear('dirname');
    dirname=varargin{findStringInCell(varargin,'dirname')+1};
end


if isdir(dirname)
    if ~isempty(fsic(varargin, 'forceOverwrite'))
        bNew = true;
        fprintf('Forced to overwrite directory: %s\n', dirname);
    else
        messg={sprintf('The specified directory %s already contains a previously recorded experiment', dirname)
            ''
            'Continue experiment, overwrite  or cancel ?'};
        button1 = questdlg(messg,'DIRECTORY NOT EMPTY','Continue','Overwrite','Cancel','Continue');
        switch button1
            case 'Overwrite'
                button2 = questdlg({sprintf('Are you sure you want to overwrite data in %s?', dirname)} ,'OVERWRITE EXPERIMENT ?');
                switch button2
                    case 'Yes',
                        rmdir(dirname,'s')
                    otherwise,
                        return
                end
            case 'Continue'
                bNew=false;

            otherwise,
                return

        end
    end
end

if bNew % set up new experiment
    if isdir(dirname)
        rmdir(dirname, 's');
    end
    mkdir(dirname)
%     copyfile('expt_config.txt',fullfile(dirname,'expt_config.txt'));
    
	expt.subject=subject;
    
    expt.allPhases={'pre', 'pract1', 'pract2', 'other', 'rand', ...
                    'start', 'ramp', 'stay1', 'noise', 'stay2', 'end'};
%     expt.allPhases={'noise', 'pre', 'pract1', 'pract2', 'other', 'rand', ...
%                     'start', 'ramp', 'stay1', 'stay2', 'end'};
    expt.recPhases={'pre', 'pract1', 'pract2', 'other', 'rand', ...
                    'start', 'ramp', 'stay1', 'noise', 'stay2', 'end'}; % SC The pahses during which the data are recorded
% 	expt.recPhases={'noise', 'pre', 'pract1', 'pract2', 'other', 'rand', ...
%                     'start', 'ramp', 'stay1', 'stay2', 'end'};
    
    if isfield(expt_config,'NATURAL_REPS') && isfield(expt_config,'NATURAL_WORDS') && ...
        ~isempty(expt_config.NATURAL_REPS) && ~isempty(expt_config.NATURAL_WORDS) && ...
        iscell(expt_config.NATURAL_WORDS) && expt_config.NATURAL_REPS>0
        expt.allPhases=[{'natural'},expt.allPhases];
        expt.recPhases=[{'natural'},expt.recPhases];
    end
    
    expt.preWords = expt_config.PRE_WORDS;
    expt.pract1Words = expt.preWords;
    expt.pract2Words = expt.preWords;
    expt.randWords = expt_config.RAND_WORDS;	% 3 3-letters, 5 4-letters
    expt.sustWords = expt_config.SUST_WORDS;
% 	expt.pseudoWords=[];
    
    if isequal(expt.allPhases{1},'natural')
        expt.naturalWords=expt_config.NATURAL_WORDS;
    end

%     if expt_config.TEST1_REPS==0    % No generization tests 
    expt.trialTypes=[1];
    expt.trialOrderRandReps = 1;	%How many reps are randomized together
    expt.script.pre.nReps = expt_config.PRE_REPS;    %SC Numbers of repetitions in the stages   % !!1!!	
    expt.script.pract1.nReps = expt_config.PRACT1_REPS; %SC Default 2   %SC-Mod(09/26/2007)        % !!1!!
    expt.script.pract2.nReps = expt_config.PRACT2_REPS; %SC Default 2   %SC-Mod(09/26/2007)        % !!1!!
    
    expt.script.rand.nBlocks = expt_config.RAND_BLOCKS;  %SC Default 10   %SC-Mod(09/26/2007)       % !!8!!
    expt.script.rand.trialsPerBlock = expt_config.RAND_TRIALS_PER_BLOCK; 
    expt.script.rand.trialsPerBlock_lower = expt_config.RAND_LOWER_TRIALS_PER_BLOCK;
    expt.script.rand.trialsPerBlock_higher = expt_config.RAND_HIGHER_TRIALS_PER_BLOCK;
    
    expt.script.start.nReps = expt_config.SUST_START_REPS;   %SC Default 15   %SC-Mod(09/26/2007)       % !!2!!
    expt.script.ramp.nReps = expt_config.SUST_RAMP_REPS;   %SC Default 15   %SC-Mod(09/26/2007)       % !!2!!
%     expt.script.stay.nReps = expt_config.SUST_STAY_REPS;   %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!
    expt.script.stay1.nReps = expt_config.SUST_STAY1_REPS;
    expt.script.noise.nReps = expt_config.SUST_NOISE_REPS;
    expt.script.stay2.nReps = expt_config.SUST_STAY2_REPS;
    expt.script.end.nReps = expt_config.SUST_END_REPS;    %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!

	expt.trialTypeDesc = cell(1, 5);
	expt.trialTypeDesc{1} = 'Speech with auditory feedback';
	expt.trialTypeDesc{2} = 'Speech with masking noise';
	expt.trialTypeDesc{3} = 'Listen to masking noise, no speech';
	expt.trialTypeDesc{4} = 'Rest (no speech) in silence';
	expt.trialTypeDesc{5} = 'Non-speech bracket task in silence';
	
    expt.script.pre    = genPhaseScript('pre',    ...
                                        expt.script.pre.nReps, expt.preWords);
    expt.script.pract1 = genPhaseScript('pract1', ...
                                        expt.script.pract1.nReps, expt.pract1Words);
    expt.script.pract2 = genPhaseScript('pract2', ...
                                        expt.script.pract2.nReps, expt.pract2Words);
    
    expt.script.other = struct;
    expt.script.other.nReps = 1;
    expt.script.other.nTrials = 8;
    expt.script.other.rep1 = struct;
    expt.script.other.rep1.trialOrder = [1, 1, 1, 1, 1, 1, 1, 1];
    expt.script.other.rep1.word = {'CAT', 'TED', 'PIG', 'CAT', 'TED', 'PIG', 'CAT', 'PIG'};
    
    fprintf('Generating script for the random-perturbation phase... ');
    expt.script.rand = genRandScript(expt.script.rand.nBlocks, expt.script.rand.trialsPerBlock, ...
                                     expt.script.rand.trialsPerBlock_lower, expt.script.rand.trialsPerBlock_higher, ...
                                     expt.randWords);
	fprintf('Done.\n');
                                    
    expt.script.start  = genPhaseScript('start',  ...
                                        expt.script.start.nReps,  expt.sustWords);
%     if isfield(expt.script,'natural')
%         expt.script.natural=   genPhaseScript('natural',    expt.script.natural.nReps,    expt.trialTypes,expt.preWords, expt.trainWords, expt.testWords, expt.pseudoWords, expt.trialOrderRandReps, expt.subject.designNum,expt.naturalWords);
%     end
    
%     if isfield(expt.script,'test1')
%         expt.script.test1= genPhaseScript('test1',  expt.script.test1.nReps,  expt.trialTypes,expt.preWords, expt.trainWords, expt.testWords, expt.pseudoWords, expt.trialOrderRandReps, expt.subject.designNum);
%     end
    expt.script.ramp    = genPhaseScript('ramp', ...
                                         expt.script.ramp.nReps, expt.sustWords);
    expt.script.stay1    = genPhaseScript('stay1', ...
                                          expt.script.stay1.nReps, expt.sustWords);
    expt.script.noise    = genPhaseScript('noise', ...
                                          expt.script.noise.nReps, expt.sustWords);
    expt.script.stay2    = genPhaseScript('stay2', ...
                                          expt.script.stay2.nReps, expt.sustWords);
%     if isfield(expt.script,'test2')
%         expt.script.test2= genPhaseScript('test2',  expt.script.test2.nReps,  expt.trialTypes,expt.preWords, expt.trainWords, expt.testWords, expt.pseudoWords, expt.trialOrderRandReps, expt.subject.designNum);
%     end
    
%     expt.script.stay2=  genPhaseScript('stay2',   expt.script.stay2.nReps,   expt.trialTypes,expt.preWords, expt.trainWords, expt.testWords, expt.pseudoWords, expt.trialOrderRandReps, expt.subject.designNum);
    expt.script.end      = genPhaseScript('end', ...
                                          expt.script.end.nReps, expt.sustWords);
    
%     if isfield(expt.script,'test3')
%         expt.script.test3= genPhaseScript('test3',  expt.script.test2.nReps,  expt.trialTypes,expt.trainWords, expt.testWords, expt.pseudoWords, expt.trialOrderRandReps, expt.subject.designNum);
%     end
    
% 	expt.script.nTotFaces=length(find(expt.trialTypes==5))*(expt.script.pre.nReps+expt.script.pract1.nReps+expt.script.pract2.nReps+...
% 		expt.script.start.nReps+expt.script.ramp.nReps+expt.script.end.nReps);
    p=getTSMDefaultParams(subject.sex,...
        'closedLoopGain',expt.subject.closedLoopGain,...
        'trialLen',expt.subject.trialLen,...
        'trialLenMax', expt.subject.trialLenMax, ...
        'mouthMicDist',expt.subject.mouthMicDist,...
        'sr',expt_config.SAMPLING_RATE/expt_config.DOWNSAMP_FACT,...
        'downFact',expt_config.DOWNSAMP_FACT,...
        'frameLen',expt_config.FRAME_SIZE / expt_config.DOWNSAMP_FACT, ...
        'subject_group', expt_config.SUBJECT_GROUP);
%     p=getTSMDefaultParams(subject.sex,subject.shiftDirection,...
%         'closedLoopGain',expt.subject.closedLoopGain,...
%         'trialLen',expt.subject.trialLen,...
%         'mouthMicDist',expt.subject.mouthMicDist);
%     
    state.phase=1;
    state.rep=1;
    state.params=p;
    rmsPeaks=[];
    
    save(fullfile(dirname,'expt.mat'),'expt');
    save(fullfile(dirname,'state.mat'),'state');
else % load expt
    load(fullfile(dirname,'state.mat'));
    load(fullfile(dirname,'expt.mat'));            
    p=state.params;
%     nPeaks=length(expt.trainWords);
%     if state.phase>1
%     rmsPeaks=ones(length(expt.trainWords),1)*p.rmsMeanPeak;    %SC ***Bug!!***
%     end
    subject=expt.subject;
end


% Make a copy of the configuration file
dconfig = dir(fullfile(dirname, 'config*.txt'));
configBackupFN = fullfile(dirname, sprintf('expt_config_%.2d.txt', length(dconfig) + 1));
copyfile(configFN, configBackupFN);
check_file(configBackupFN);
fprintf(1, 'INFO: Made a backup of the experiment configuration file at: %s\n', configBackupFN);

%% initialize algorithm
MexIO('init',p);      %SC Set the initial (default) parameters

if ((p.frameShift-round(p.frameShift)~=0) || (p.frameShift>p.frameLen))
    uiwait(errordlg(['Frameshift = ' num2str(p.frameShift) ' is a bad value. Set nWin and frameLen appropriately. Frameshift must be an integer & Frameshift <= Framelen'],'!! Error !!'))
    return
else

    fprintf('\n  \n')
    TransShiftMex(0);           %SC Gives input/output device info, and serves as an initialization.

    fprintf('\nSettings : \n')
    fprintf('DMA Buffer    = %i samples \n',p.frameLen) %SC Buffer length after downsampling
    fprintf('Samplerate    = %4.2f kHz \n',p.sr/1000)   %SC sampling rate after downsampling
    fprintf('Analysis win  = %4.2f msec \n',p.bufLen/p.sr*1000)
    fprintf('LPC  window   = %4.2f msec \n',p.anaLen/p.sr*1000)

    fprintf('Process delay = %4.2f msec \n',p.nDelay*p.frameLen/p.sr*1000)
    fprintf('Process/sec   = %4.2f \n',p.sr/p.frameShift)

end

%% Load masking noise
check_file(expt_config.MASKING_NOISE_WAV_FILE_NAME);
[mn_w, mn_fs] = wavread(expt_config.MASKING_NOISE_WAV_FILE_NAME);
if mn_fs ~= expt_config.SAMPLING_RATE
    error('The wav file of masking noise does not have the required sampling rate: %.2f Hz\n', expt_config.SAMPLING_RATE);
end
fprintf(1, 'INFO: Loaded masking-noise waveform from file: %s\n', expt_config.MASKING_NOISE_WAV_FILE_NAME);

mn_w = mn_w * (10 ^ (expt_config.MASKING_NOISE_LEVEL_RELATIVE / 20));
fprintf(1, 'INFO: Scaled masking-noise waveform by %.2f dB\n', expt_config.MASKING_NOISE_LEVEL_RELATIVE);
if ~isempty(find(mn_w > 1, 1)) || ~isempty(find(mn_w < -1, 1))
    fprintf(1, 'WARNING: clipping detected in scaled masking-noise waveform. Consider using lower scaling.');
    mn_w(mn_w > 1) = 1;
    mn_w(mn_w < -1) = -1;
end

TransShiftMex(3, 'datapb', mn_w, 0);

%% Load the multi-talker babble noise
% [x,fs_mtb]=wavread('mtbabble48k.wav');
% lenMTB=round(2.5*fs_mtb);

% gainMTB_fb=dBSPL2WaveAmp(subject.lvNoise,1000)/sqrt(2)/calcMaskNoiseRMS;
% gainMTB_fb3=dBSPL2WaveAmp(subject.lvNoise3,1000,subject.pcrKnob)/sqrt(2)/calcMaskNoiseRMS;
% x_mtb=cell(1,3);
% x_mtb{1}=x(1:lenMTB);               x_mtb{1}=x_mtb{1}-mean(x_mtb{1});
% x_mtb{2}=x(lenMTB+1:lenMTB*2);      x_mtb{2}=x_mtb{2}-mean(x_mtb{2});
% x_mtb{3}=x(lenMTB*2+1:lenMTB*3);    x_mtb{3}=x_mtb{3}-mean(x_mtb{3});
% TransShiftMex(3,'datapb',x_mtb{1});

%% expt
figIdDat=makeFigDataMon;

% wordList=expt.words;

allPhases=expt.allPhases;
recPhases=expt.recPhases;
% nWords=length(wordList);

hgui = UIRecorder('figIdDat', figIdDat, 'dirname', dirname);
set(hgui.UIrecorder, 'Position', [1155, 150, 440, 700]);
% winontop(hgui.UIrecorder, 1);


if ~isempty(fsic(varargin, 'twoScreens'))
    subjFigPos = get(hgui.hkf, 'Position');
    set(hgui.hkf, 'Position', [1800, 100, subjFigPos(3), subjFigPos(4)]);
end

% if (expt.subject.designNum==2)
%     expt.script=addFaceInfo(expt.script,hgui.skin.dFaces);
%     expt.dFaces=hgui.skin.dFaces;
% end

hgui.bSim = bSim;
hgui.simDataDir = simDataDir;
hgui.dirname = dirname;

hgui.pcrKnob=subject.pcrKnob;
hgui.ITI=expt.subject.ITI;
hgui.trigByScanner=expt.subject.trigByScanner;
hgui.TA=expt.subject.TA;
hgui.dBRange=expt.subject.dBRange1;
hgui.trialLen=expt.subject.trialLen;
hgui.trialLenMax = expt.subject.trialLenMax;
% hgui.skin.faceOrder=randperm(length(hgui.skin.dFaces));
hgui.skin.facePnt=1;

hgui.bAO=bAO;
hgui.dScale=p.dScale;

hgui.vumeterMode=expt.subject.vumeterMode;

hgui.rmsTransTarg_spl=getSPLTarg(expt_config.SPL_TARGET,expt.subject.mouthMicDist);
load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));

fprintf('\n');
disp(['Mouth-microphone distance = ',num2str(expt.subject.mouthMicDist),' cm']);
disp(['hgui.rmsTransTarg_spl = ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
fprintf('\n');

hgui.vocaLen=round(expt_config.VOWEL_LEN_TARG*p.sr/(p.frameLen)); % 300 ms, 225 frames
hgui.lenRange=2.5*round(expt_config.VOWEL_LEN_RANGE*p.sr/(p.frameLen));  % single-sided tolerance range: 0.4*250 = 100 ms
disp(['Vowel duration range: [',num2str(300-0.4*250),',',num2str(300+0.4*250),'] ms.']);

hgui.debug=DEBUG;
% hgui.trigKey='equal';

if (isempty(findStringInCell(varargin,'twoScreens')))
% 	set(hgui.UIrecorder,...
% 		'position', [0    5.0000  250.6667   65.8750],...
% 		'toolbar','none');  %SC Set the position of the expt window, partially for the use of multiple monitors.
else
% 	if (expt.subject.trigByScanner==1)
% 		ms=get(0,'MonitorPosition');
% 		set(hgui.UIrecorder,'Position',[ms(2,1),ms(1,4)-ms(2,4),ms(2,3)-ms(2,1)+1,ms(2,4)+20],'toolbar','none','doublebuffer','on','renderer','painters');
% 		pos_win=get(hgui.UIrecorder,'Position');
% 		pos_strh=get(hgui.strh,'Position');
% 		pos_axes_pic=get(hgui.axes_pic,'Position');
% 		pos_rms_axes=get(hgui.rms_axes,'Position');
% 		pos_speed_axes=get(hgui.speed_axes,'Position');
% 		pos_rms_label=get(hgui.rms_label,'Position');
% 		pos_rms_too_soft=get(hgui.rms_too_soft,'Position');
% 		pos_rms_too_loud=get(hgui.rms_too_loud,'Position');
% 		pos_speed_label=get(hgui.speed_label,'Position');
% 		pos_speed_too_slow=get(hgui.speed_too_slow,'Position');
% 		pos_speed_too_fast=get(hgui.speed_too_fast,'Position');
% 		set(hgui.strh,'Position',[(pos_win(3)-pos_strh(3))/2+5,(pos_win(4)-pos_strh(4))/2-15,pos_strh(3),pos_strh(4)*0.9]);
% 		set(hgui.axes_pic,'Position',[(pos_win(3)-pos_axes_pic(3))/2,(pos_win(4)-pos_axes_pic(4))/2,pos_axes_pic(3),pos_axes_pic(4)]);
% 		set(hgui.rms_axes,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_axes(2),pos_rms_axes(3),pos_rms_axes(4)]);
% 		set(hgui.rms_label,'Position',[(pos_win(3)-pos_rms_label(3))/2,pos_rms_label(2),pos_rms_label(3),pos_rms_label(4)]);
% 		set(hgui.rms_too_soft,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_too_soft(2),pos_rms_too_soft(3),pos_rms_too_soft(4)]);
% 		set(hgui.rms_too_loud,'Position',[(pos_win(3)-pos_rms_axes(3))/2+pos_rms_axes(3)-pos_rms_too_loud(3),pos_rms_too_loud(2),pos_rms_too_loud(3),pos_rms_too_loud(4)]);
% 		set(hgui.speed_axes,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_axes(2),pos_speed_axes(3),pos_speed_axes(4)]);		
% 		set(hgui.speed_label,'Position',[(pos_win(3)-pos_speed_label(3))/2,pos_speed_label(2),pos_speed_label(3),pos_speed_label(4)]);
% 		set(hgui.speed_too_slow,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_too_slow(2),pos_speed_too_slow(3),pos_speed_too_slow(4)]);
% 		set(hgui.speed_too_fast,'Position',[(pos_win(3)-pos_speed_axes(3))/2+pos_speed_axes(3)-pos_speed_too_fast(3),pos_speed_too_fast(2),pos_speed_too_fast(3),pos_speed_too_fast(4)]);
%         set(hgui.msgh,'FontSize',17);
% 	else
% 		set(hgui.UIrecorder,'Position',[-1400,180,1254,857],'toolbar','none');
% 	end
	
end

if (subject.showProgress)
	set(hgui.progress_axes,'visible','on');
	set(hgui.progress_imgh,'visible','on');
	progress_meter=0.5*ones(1,100,3);
	progress_mask=zeros(1,100,3);
	set(hgui.progress_imgh,'Cdata',progress_meter.*progress_mask);
    set(hgui.progress_axes, 'YTick', [], 'YColor', [0, 0, 0]);
    set(hgui.progress_axes, 'XTick', [0 : 25 : 100]);
else
	set(hgui.progress_axes,'visible','off');
	set(hgui.progress_imgh,'visible','off');
end

TransShiftMex(2);
if bAO    
    MexIO('reset');
    TransShiftMex(1);
    TransShiftMex(3,'scale',0);
end

rProgress=0;
startPhase=state.phase; %SC For the purpose of resumed experiments
startRep=state.rep;     %SC For the purpose of resumed experiments
for n=startPhase:length(allPhases)
    state.phase=n;
    state.rep=1;
    thisphase=allPhases{1,n};
    subdirname=fullfile(dirname,thisphase);
    mkdir(subdirname);
    
    hgui.phase=thisphase;
    
    % Adjust the number of reps
    if (~isequal(thisphase,'ramp') && ~isequal(thisphase,'stay'))
        disp(['--- Coming up: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
            '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
        nRepsNew=input('(Enter to skip) nRepsNew = ','s');
        nRepsNew=str2num(nRepsNew);
        if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
            expt.script.(thisphase).nReps=nRepsNew;
            expt.script.(thisphase)=genPhaseScript(thisphase, expt.script.(thisphase).nReps,...
                expt.preWords, expt.trialTypes,expt.trainWords,expt.testWords,expt.pseudoWords,...
                expt.trialOrderRandReps,expt.subject.designNum);
            disp(['Changed: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
                '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
            save(fullfile(dirname,'expt.mat'),'expt');
            disp(['Saved ',fullfile(dirname,'expt.mat')]);
        end        
    elseif isequal(thisphase,'ramp')
        disp(['--- Coming up: ','ramp','. nReps = ',num2str(expt.script.ramp.nReps),...
            '; nTrials = ',num2str(expt.script.ramp.nTrials),' ---']);
        disp(['--- Coming up: ','stay','. nReps = ',num2str(expt.script.stay.nReps),...
            '; nTrials = ',num2str(expt.script.stay.nTrials),' ---']);
        disp(['--- Ramp+Stay: nReps = ',num2str(expt.script.ramp.nReps+expt.script.stay.nReps),...
            '; nTrials = ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials)]);
        nRepsNew=input('(Enter to skip) Ramp: nRepsNew = ','s');
        nRepsNew=str2num(nRepsNew);
        if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
            expt.script.ramp.nReps=nRepsNew;
            expt.script.ramp=genPhaseScript('ramp',expt.script.ramp.nReps,...
                expt.trialTypes,expt.preWords, expt.trainWords,expt.testWords,expt.pseudoWords,...
                expt.trialOrderRandReps,expt.subject.designNum);
            disp(['Changed: ramp. ','nReps = ',num2str(expt.script.ramp.nReps),...
                '; nTrials = ',num2str(expt.script.ramp.nTrials),' ---']);
            save(fullfile(dirname,'expt.mat'),'expt');
            disp(['Saved ',fullfile(dirname,'expt.mat')]);
        end
        nRepsNew=input('(Enter to skip) Stay: nRepsNew = ','s');
        nRepsNew=str2num(nRepsNew);
        if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
            expt.script.stay.nReps=nRepsNew;
            expt.script.stay=genPhaseScript('stay',expt.script.stay.nReps,...
                expt.trialTypes,expt.preWords, expt.trainWords,expt.testWords,expt.pseudoWords,...
                expt.trialOrderRandReps,expt.subject.designNum);
            disp(['Changed: stay. ','nReps = ',num2str(expt.script.stay.nReps),...
                '; nTrials = ',num2str(expt.script.stay.nTrials),' ---']);
            save(fullfile(dirname,'expt.mat'),'expt');
            disp(['Saved ',fullfile(dirname,'expt.mat')]);
        end
        disp(['--- Ramp+Stay: nReps = ',num2str(expt.script.ramp.nReps+expt.script.stay.nReps),...
            '; nTrials = ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials)]);
    end
    % Adjust the number of reps
    
    nReps=expt.script.(thisphase).nReps;
%     if ~isequal(thisphase,'stay')
        phaseTrialCnt=1;
%     end

    expt.script.(thisphase).startTime=clock;

% 	set(hgui.rms_axes,'visible','on');
%     set(hgui.rms_imgh,'visible','on');
%     set(hgui.rms_label,'visible','on');
% 	set(hgui.rms_too_soft,'visible','on');
% 	set(hgui.rms_too_loud,'visible','on');	
%     set(hgui.speed_axes,'visible','on');
%     set(hgui.speed_imgh,'visible','on');
%     set(hgui.speed_label,'visible','on');
% 	set(hgui.speed_too_slow,'visible','on');
%     set(hgui.speed_too_fast,'visible','on');
    
    hgui.showSpeedPrompt = 0;
    hgui.showRmsPrompt = 0;
    hgui.bSpeedRepeat=0;
    hgui.bRmsRepeat=0;
	
	if (subject.showPlayButton==0)
		set(hgui.play,'visible','off');
	end
    
    switch(thisphase)
        case 'pre'
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             set(hgui.rms_axes,'visible','off');
% 		    set(hgui.rms_imgh,'visible','off');
% 		    set(hgui.rms_label,'visible','off');
% 			set(hgui.rms_too_loud,'visible','off');
% 		    set(hgui.rms_too_soft,'visible','off');
%             set(hgui.speed_axes,'visible','off');
% 		    set(hgui.speed_imgh,'visible','off');
% 		    set(hgui.speed_label,'visible','off');
% 			set(hgui.speed_too_slow,'visible','off');
% 		    set(hgui.speed_too_fast,'visible','off');
            
            hgui.showRmsPrompt = 0;
            hgui.showSpeedPrompt = 0;
            hgui.bRmsRepeat=0;
            hgui.bSpeedRepeat=0;
            
            p.bDetect=0;
            p.bShift = 0;       %SC No shift in the practice-1 phase
            
        case 'pract1'           
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             if (hgui.vumeterMode==1)
%                 vumeter=hgui.skin.vumeter;
%             elseif (hgui.vumeterMode==2)
%                 vumeter=hgui.skin.vumeter2;
%             end
%             mask=0.5*ones(size(vumeter));
%             % mask(1:50,:,:) = 1;           %SC-Commented(12/11/2007)
%             set(hgui.rms_imgh,'Cdata',vumeter.*mask);
            p.bDetect=0;
            p.bShift = 0;       %SC No shift in the practice-1 phase
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 0;
            hgui.bRmsRepeat=1;
            hgui.bSpeedRepeat=0;
% 		    set(hgui.speed_axes,'visible','off');
% 		    set(hgui.speed_imgh,'visible','off');
% 		    set(hgui.speed_label,'visible','off');
% 			set(hgui.speed_too_slow,'visible','off');
% 		    set(hgui.speed_too_fast,'visible','off');
   
            if exist('rmsPeaks') && ~isempty(rmsPeaks)
                p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
            end

            hgui.showTextCue=1;
            
            subjProdLevel=[];
         case 'pract2'
            if exist('subjProdLevel')
                subjProdLevel=subjProdLevel(find(~isnan(subjProdLevel)));

                if (~isempty(subjProdLevel))
                    hgui.rmsTransTarg_spl=mean(subjProdLevel);
                    load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
                    hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));
                end
            end
            
            fprintf('\n');
            disp(['Target level set as subject mean production level: ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
            fprintf('\n');            
             
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);

            p.bDetect=0;
            p.bShift = 0;
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 1;  %1 
            hgui.bSpeedRepeat = 1;        %SC Make the speed monitor visible %1
            
            if exist('rmsPeaks') && ~isempty(rmsPeaks)
                p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
            end 

            hgui.showTextCue=1;
            
        case 'other'
            if exist('subjProdLevel')
                subjProdLevel=subjProdLevel(find(~isnan(subjProdLevel)));

                if (~isempty(subjProdLevel))
                    hgui.rmsTransTarg_spl=mean(subjProdLevel);
                    load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
                    hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));
                end
            end
            
            fprintf('\n');
            disp(['Target level set as subject mean production level: ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
            fprintf('\n');            
             
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);

            p.bDetect=0;
            p.bShift = 0;
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 1;  %1 
            hgui.bSpeedRepeat = 1;        %SC Make the speed monitor visible %1
            
            if exist('rmsPeaks') && ~isempty(rmsPeaks)
                p.rmsMeanPeak=mean(rmsPeaks);
%                 p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
            end 

            hgui.showTextCue=1;

        case 'rand'
            if bAO
                TransShiftMex(2);
            end
            
            % SC(2008/06/10) Manually determine the optimum tracking params
            % Warning: for consistency, don't change nDelay			
			set(hgui.msgh,'visible','on');
%             set(hgui.msgh_imgh,'CData',CDataMessage.ftparampicking,'visible','on');
			drawnow;
			
			set(hgui.msgh,'string',{'Please stand by...'},'visible','on');
            
            if (hgui.debug==0)
                [vowelF0Mean,vowelF0SD] = getVowelPitches(dirname);
                
                disp(['Vowel meanF0 = ', num2str(vowelF0Mean),' Hz: stdF0 = ',num2str(vowelF0SD),' Hz']);
                disp(['Recommended cepsWinWidth = ',num2str(round(p.sr/vowelF0Mean*0.54))]);
                [vowelF1Mean,vowelF2Mean]=getVowelMeanF1F2(dirname);
                disp(['Vowel meanF1 = ',num2str(vowelF1Mean),' Hz; meanF2 = ',num2str(vowelF2Mean),' Hz']);
                
%                 optimFTParams=compareFormTrackParams(dirname);
%                 p.nLPC=optimFTParams.nLPC;
%                 p.nDelay=optimFTParams.nDelay;
%                 p.bufLen=(2*p.nDelay-1)*(p.frameLen);
%                 p.anaLen=p.frameShift+2*(p.nDelay-1)*p.frameLen;
%                 p.avgLen=optimFTParams.avgLen;
%                 p.bCepsLift=optimFTParams.bCepsLift;
%                 p.cepsWinWidth=optimFTParams.cepsWinWidth;
%                 p.fn1=optimFTParams.fn1;
%                 p.fn2=optimFTParams.fn2;       
%                 p.aFact=optimFTParams.aFact;
%                 p.bFact=optimFTParams.bFact;
%                 p.gFact=optimFTParams.gFact;
            end
            % ~SC(2008/06/10) Manually determine the optimum tracking
            
            if bAO
                TransShiftMex(1);
                TransShiftMex(3,'scale',0);
            end
            
            set(hgui.msgh, 'string', {''}, 'visible', 'on'); 
            if exist('rmsPeaks')
                p.rmsMeanPeak=mean(rmsPeaks);
            end
%             p.rmsThresh=p.rmsMeanPeak/4;    %SC !! Adaptive RMS threshold setting. Always updating
% 			if (p.rmsThresh>0.015)
% 				p.rmsThresh=0.015;
% 				disp('********* Warning: rms too high! Limited at 0.015. *********');
% 			end
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 0;  %1 
            hgui.bSpeedRepeat = 0; 

            p.bDetect=1;
            p.bShift=0;
            
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
            hgui.showTextCue=1;
        case 'start'
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 0;  %1 
            hgui.bSpeedRepeat = 0; 
            
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
            hgui.showTextCue=1;
        case 'ramp'             %SC !! Notice that adaptive RMS threshold updating is no longer done here.           			
            %====================================================
            % Perform the Vowel Identification task
%             if bAO
%                 TransShiftMex(2);
%             end
            
%             [meanF1_EH,meanF2_EH,meanF1_AE,meanF2_AE,medianF0,bKeep]=extract_EH_AE_formants(dirname,{'pract2','rand','start'});
%             ehaeInfo=struct;
%             ehaeInfo.expDir=''; ehaeInfo.stage=''; ehaeInfo.rep=NaN; ehaeInfo.word=''; 
%             ehaeInfo.meanF1_EH=meanF1_EH;
%             ehaeInfo.meanF2_EH=meanF2_EH;
%             ehaeInfo.meanF1_AE=meanF1_AE;
%             ehaeInfo.meanF2_AE=meanF2_AE;
%             ehaeInfo.F0=medianF0;    % TODO
%             ehaeInfo.arbitraryF0=0;
%             ehaeInfo.bAO=bAO;
%             ehaeInfo.bKeep=bKeep;
%             if isempty(ehaeInfo.F0) || isnan(ehaeInfo.F0)
%                 if isequal(expt.subject.sex,'male')
%                     ehaeInfo.F0=120;
%                     
%                 else
%                     ehaeInfo.F0=240;
%                 end
%                 ehaeInfo.arbitraryF0=1;
%             end
%             expt.ehaeInfo=ehaeInfo;
%             save(fullfile(dirname,'expt.mat'),'expt');
%             if ~isnan(ehaeInfo.meanF1_EH) && ~isempty(ehaeInfo.meanF1_EH) &&...
%                ~isnan(ehaeInfo.meanF2_EH) && ~isempty(ehaeInfo.meanF2_EH) &&...
%                ~isnan(ehaeInfo.meanF1_AE) && ~isempty(ehaeInfo.meanF1_AE) &&...
%                ~isnan(ehaeInfo.meanF2_AE) && ~isempty(ehaeInfo.meanF2_AE) 
%                 expt.bDoVID1=input('Proceed with vowel identification #1? (0/1) (Default=1): ','s');
%                 if isempty(expt.bDoVID1)
%                     expt.bDoVID1=1;
%                 else
%                     if ischar(expt.bDoVID1)
%                         expt.bDoVID1=str2num(expt.bDoVID1);
%                     end
%                 end
%                 save(fullfile(dirname,'expt.mat'),'expt');
% %                 if expt.bDoVID1
% %                     eh_discrim([expt.subject.name,'_VID1'],'eh',expt.subject.percDir,11,ehaeInfo,[],'twoScreens');
% %                 end
%             end
%             input('Press enter to proceed to "ramp" phase: ');
            % ~Perform the Vowel Identification task
            %====================================================
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 0;  %1 
            hgui.bSpeedRepeat = 0; 
            
            p.bDetect = 1;
			p.bShift = 1;
            
            if bAO
                TransShiftMex(1);
                TransShiftMex(3,'scale',0);
            end
            
            set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
			hgui.showTextCue=1;
%             if doPlot
%                 uiwait(gcf,10);
% 			end
        case 'stay'         	
			set(hgui.msgh,'visible','on');
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 0;  %1 
            hgui.bSpeedRepeat = 0; 
            
            p.bDetect = 1;
            p.bShift = 1;
            hgui.showTextCue=1;
%         case 'test2'
%             set(hgui.msgh,'visible','on');
%             
%             hgui.showRmsPrompt = 1;
%             hgui.showSpeedPrompt = 1;
%             hgui.bRmsRepeat = 0;  %1 
%             hgui.bSpeedRepeat = 0; 
%             
%             p.bDetect = 0;
%             p.bShift = 0;
%             hgui.showTextCue=1;
%         case 'stay2';
%             if bAO
%                 TransShiftMex(2);
%             end
%             
%             %====================================================
%             % Perform the Vowel Identification task: repetition after adaptation
%             expt.bDoVID2=input('Proceed with vowel identification #2? (0/1) (Default=1): ','s');
%             if isempty(expt.bDoVID2)
%                 expt.bDoVID2=1;
%             else
%                 if ischar(expt.bDoVID2)
%                     expt.bDoVID2=str2num(expt.bDoVID2);
%                 end
%             end
%             save(fullfile(dirname,'expt.mat'),'expt');
% %             if expt.bDoVID2
% %                 eh_discrim([expt.subject.name,'_VID2'],'eh',expt.subject.percDir,11,ehaeInfo,[],'twoScreens');
% %             end
%             input('Press enter to proceed to "stay2" phase: ');
%             % ~Perform the Vowel Identification task
%             %====================================================
%             
%             set(hgui.msgh,'visible','on');
%             p.bDetect = 1;
%             p.bShift = 1;
%             
%             if bAO
%                 TransShiftMex(1);
%                 TransShiftMex(3,'scale',0);
%             end
%             
%             set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
%             hgui.showTextCue=1;
        case 'end'       			
			set(hgui.msgh,'visible','on');
            
            hgui.showRmsPrompt = 1;
            hgui.showSpeedPrompt = 1;
            hgui.bRmsRepeat = 0;  %1 
            hgui.bSpeedRepeat = 0; 
            
            p.bDetect = 0;
            p.bShift = 0;
            hgui.showTextCue=1;
        case 'test3'
            
    end
    drawnow    
    

    set(hgui.msgh,'string',getMsgStr(thisphase),'visible','on');        

    if ~bAO
        if isfile(fullfile(dirname, 'p.mat'))
            load(fullfile(dirname, 'p.mat'))    % gives p;            
        end
        MexIO('init',p);  %SC Inject p to TransShiftMex
    else
        p0=p; 
        p0.dScale=0;
        MexIO('init',p0);
    end

    bSkip = 0;
    for i0=startRep:nReps    %SC Loop for the reps in the phase
        if isequal(thisphase, 'pre') || isequal(thisphase, 'pract1') || isequal(thisphase, 'pract2')
            if bSkip
                fprintf('Skipping the rest of the %s phase...\n', thisphase);
                bSkip = 0;
                break;
            end
        end
        
        repString=['rep',num2str(i0)];
        state.rep=i0;
        state.params=p;
        save(fullfile(dirname,'state.mat'),'state');
        
        nTrials=length(expt.script.(thisphase).(repString).trialOrder);

        subsubdirname=fullfile(subdirname,repString);
        mkdir(subsubdirname);
		
		% --- Perturbation field ---
		p.pertF2=linspace(p.F2Min, p.F2Max, p.pertFieldN);
        t_amp = norm([subject.expt_config.SHIFT_RATIO_SUST_F1, ...
                      subject.expt_config.SHIFT_RATIO_SUST_F2]);
        t_angle = atan2(subject.expt_config.SHIFT_RATIO_SUST_F2, ...
                         subject.expt_config.SHIFT_RATIO_SUST_F1);
        switch (thisphase)
            case 'start'
                p.pertAmp = zeros(1, p.pertFieldN);
                p.pertPhi = zeros(1, p.pertFieldN);
                p.bShift = 0;
            case 'ramp'
				p.pertAmp = i0 / (expt.script.ramp.nReps+1) * t_amp * ones(1, p.pertFieldN);
                p.pertPhi = t_angle * ones(1, p.pertFieldN);
                p.bShift = 1;
            case 'stay'
                p.pertAmp = t_amp * ones(1, p.pertFieldN);
                p.pertPhi = t_angle * ones(1, p.pertFieldN);
                p.bShift = 1;
%             case 'stay2'
%                 p.pertAmp=subject.shiftRatio*ones(1,p.pertFieldN);
            otherwise,
				p.pertAmp = zeros(1, p.pertFieldN);
                p.pertPhi = zeros(1, p.pertFieldN);
                p.bShift = 0;	
		end
% 		if isequal(subject.shiftDirection,'F1Up')
% 			p.pertPhi=0*ones(1,p.pertFieldN);
% 		elseif isequal(subject.shiftDirection,'F1Down')
% 			p.pertPhi=pi*ones(1,p.pertFieldN);
%         end
        
		if ~bAO
            MexIO('init',p);  %SC Inject p to TransShiftMex
        else
            p0=p; 
            p0.dScale=0;
            MexIO('init',p0);
        end
		% --- ~Perturbation field ---

        for k=1:nTrials           
            thisTrial=expt.script.(thisphase).(repString).trialOrder(k); % 0: silent; 1: no noise; 2: noise only; 			
            thisWord=expt.script.(thisphase).(repString).word{k};     %SC Retrieve the word from the randomly shuffled list

            if isequal(thisphase, 'rand')   % Configure random perturbation
                thispert = expt.script.(thisphase).(repString).pertType(k);
                if thispert == 1 % Upward perturbation                    
                    p.pertAmp = norm([subject.expt_config.SHIFT_RATIO_RAND_HIGHER_F1, ...
                                      subject.expt_config.SHIFT_RATIO_RAND_HIGHER_F2]) * ...
                                ones(1, p.pertFieldN);
                    t_angle =  atan2(subject.expt_config.SHIFT_RATIO_RAND_HIGHER_F2, ...
                                     subject.expt_config.SHIFT_RATIO_RAND_HIGHER_F1);
                    p.pertPhi = t_angle * ones(1, p.pertFieldN);
                    p.bShift = 1;
                elseif thispert == -1 % Downward perturbation
                    p.pertAmp = norm([subject.expt_config.SHIFT_RATIO_RAND_LOWER_F1, ...
                                      subject.expt_config.SHIFT_RATIO_RAND_LOWER_F2]) * ...
                                ones(1, p.pertFieldN);
                    t_angle =  atan2(subject.expt_config.SHIFT_RATIO_RAND_LOWER_F2, ...
                                     subject.expt_config.SHIFT_RATIO_RAND_LOWER_F1);
                    p.pertPhi = t_angle * ones(1, p.pertFieldN);
                    p.bShift = 1;
                else % No perturbation
                    p.pertAmp = zeros(1, p.pertFieldN);
                    p.pertPhi = zeros(1, p.pertFieldN);
                    p.bShift = 0;
                end
                
                MexIO('init',p);
            end

			hgui.trialType=thisTrial;
			hgui.word=thisWord;
            hgui.phase = thisphase;
            hgui.repNum = i0;
            hgui.trialNum = k;

%             if (hgui.trialType==2 || hgui.trialType==3)	% Speech with masking noise or passively listening to masking noise
%                 TransShiftMex(3,'datapb',gainMTB_fb*x_mtb{3-mod(k,3)},0);
% 			end
               
			disp('');
            if (ischar(thisWord))
    			disp([thisphase,' - ',repString,', k = ',num2str(k),': trialType = ',num2str(hgui.trialType),' - ',thisWord]);
            else
                disp([thisphase,' - ',repString,', k = ',num2str(k),': trialType = ',num2str(hgui.trialType),' - Pseudoword-',num2str(thisWord)]);
            end
            
            % Count down    
            if ~(isequal(thisphase,'start') || isequal(thisphase,'ramp') || isequal(thisphase,'stay') || isequal(thisphase,'end'))
                disp(['Left: ',num2str(expt.script.(thisphase).nTrials-phaseTrialCnt+1),'/',num2str(expt.script.(thisphase).nTrials)]);
            else 
                if ~(isequal(thisphase,'ramp') || isequal(thisphase,'stay'))
                    disp(['Left: ',num2str(expt.script.(thisphase).nTrials-phaseTrialCnt+1),'/',num2str(expt.script.(thisphase).nTrials),...
                        ', ',num2str((expt.script.(thisphase).nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
                else
                    disp(['Left: ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1),'/',...
                        num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials),...
                        ', ',num2str((expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
                end
            end
            % ~Count down
            
            if (hgui.trialType>=2)   %SC The distinction between train and test words                             
                TransShiftMex(3, 'bdetect', 0, 1);
                TransShiftMex(3, 'bshift', 0, 1);
			else
                TransShiftMex(3, 'bdetect', p.bDetect, 1);
                TransShiftMex(3, 'bshift', p.bShift, 1);                
            end
            
%             if (thisTrial==5)
% 				hgui.skin.facePnt=expt.script.(thisphase).(repString).face(k);
%             end
            
            updateParamDisp(p, hgui);
            set(hgui.button_reproc, 'enable', 'off');
            set(hgui.button_markBad, 'enable', 'off', 'String', 'Mark as bad');
            
            UIRecorder('singleTrial', hgui.play, 1, hgui);
            data=get(hgui.UIrecorder, 'UserData');           %SC Retrieve the data
            % -- Attach uiConfig info --
            load(hgui.uiConfigFN);
            data.uiConfig = uiConfig;
            clear('uiConfig');
            
            % -- Update the parameter settings --
            checkParams = {'rmsThresh', 'nLPC', 'fn1', 'fn2', 'aFact', 'bFact', 'gFact', 'bCepsLift', 'cepsWinWidth'};
            for k0 = 1 : numel(checkParams)
                t_param = checkParams{k0};
                if ~isequal(p.(t_param), data.params.(t_param))
                    fprintf('Updating parameter %s: %f --> %f\n', t_param, p.(t_param), data.params.(t_param));
                    p.(t_param) = data.params.(t_param);
                    
                end
            end
            
            data.timeStamp=clock;
            data.subject=expt.subject;
            data.params.name=thisWord;
            data.params.trialType=thisTrial;
            
            if (thisTrial==1)
                if ~exist('rmsPeaks')
                    rmsPeaks = [];
                end
                
                if ~isempty(data.rms)
                    switch (thisphase)  %SC Record the RMS peaks in the bout
                        case 'pre'
                            rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];
                        case 'pract1',
                            rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
                        case 'pract2',
                            rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
                        otherwise,
                    end
                end
            end
            
            if (isequal(thisphase,'pract1'))
                if (thisTrial==1 || thisTrial==2)
                    if (isfield(data,'vowelLevel') && ~isempty(data.vowelLevel) && ~isnan(data.vowelLevel) && ~isinf(data.vowelLevel))
                        subjProdLevel=[subjProdLevel,data.vowelLevel];
                    end
                end
            end
			
            save(fullfile(subsubdirname,['trial-',num2str(k),'-',num2str(thisTrial)]),'data');
            disp(['Saved ',fullfile(subsubdirname,['trial-',num2str(k),'-',num2str(thisTrial)]),'.']);
            disp(' ');
            
            phaseTrialCnt=phaseTrialCnt+1;

            % Calculate and show progress
            if (subject.showProgress)
                [rProgress, nDoneTrials, nTotTrials] = calcExpProgress(expt, thisphase, i0, k, rProgress);
				if (~isnan(rProgress))
	                progress_mask=zeros(size(progress_meter));
		            progress_mask(:,1:round(rProgress*100),:)=1;            
			        set(hgui.progress_imgh, 'Cdata', progress_meter .* progress_mask);
                    set(hgui.txt_prog, 'String', sprintf('Completed: %d / %d', nDoneTrials, nTotTrials));
				end
            end
            
            if isequal(thisphase, 'pre') || isequal(thisphase,'pract1') || isequal(thisphase,'pract2') && ...
                ~isempty(data) && ~isempty(data.fmts)
                items={'rmsThresh', 'nLPC', 'fn1', 'fn2'};
                toRepeat=1;
                while toRepeat
                    statMsg = '[';
                    for k1  = 1 : numel(items)
                        statMsg = [statMsg, sprintf('%s=%.4f', items{k1}, p.(items{k1}))];
                        if k1 ~= numel(items)
                            statMsg = [statMsg, '; '];
                        else
                            statMsg = [statMsg, ']'];
                        end
                    end
                    
                    fprintf('%s\n', statMsg);
                    cmdAdj = input(sprintf('cmd (or skip): '), 's');
                    cmdAdj = strrep(cmdAdj, ' ', '');
                    
                    if isequal(lower(deblank(cmdAdj)), 'skip')
                        bSkip = 1; % Skip the rest                        
                    else
                        bSkip = 0;
                       
                    end
                    
                    for j1=1:length(items)
                        cmdAdj=strrep(cmdAdj,';','');
                        idx=strfind(cmdAdj,[items{j1},'=']);
                        if ~isempty(idx)
                            tVal=str2num(cmdAdj(idx+length([items{j1},'=']) : end));
                            if ~isfield(p, items{j1})
                                fprintf('WARNING: p does not have the field %s\n', items{j1})
                                cmdAdj = '';
                            else
                                p.(items{j1})=tVal;
                            end
                            
%                             data=reprocAPSTVData(data,'iF2LB',p.iF2LB,'uF2UB',p.uF2UB,'rmsThresh',p.rmsThresh);
%                             taxis1=0 : (data.params.frameLen/data.params.sr) : (data.params.frameLen/data.params.sr)*(length(data.sentStat)-1);
%                             plot(taxis1,data.sentStat*250,'w','LineWidth',2); hold on;
%                             set(gca,'XLim',xs);
%                             save(fullfile(subsubdirname,['trial-',num2str(k),'-',num2str(thisTrial)]),'data');
%                             disp(['Saved NEW ',fullfile(subsubdirname,['trial-',num2str(k),'-',num2str(thisTrial)]),'.']);
                        end
                    end
                    if isempty(cmdAdj) || isequal(lower(cmdAdj),'q') || bSkip == 1
                        toRepeat=0;
                    end
                end
                
            else
                bSkip = 0;
            end
            
            if bSkip
                break;
            end
            
        end
    end
    startRep=1;
end
set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
set(hgui.msgh,'string',...
	{'Congratulations!';...
	'You have finished the expt.'},'visible','on');
% set(hgui.msgh_imgh,'CData',CDataMessage.finish,'visible','on');
% pause(1);
close(hgui.hkf);
close(hgui.figIdDat(1));
% close(hgui.UIrecorder);
% pause(1);
% saveExperiment(dirname);

if bAO
    TransShiftMex(2);
end

save(fullfile(dirname,'expt.mat'),'expt');
save(fullfile(dirname,'state.mat'),'state');

%% F1 JND (perceptual acuity) test
% percTokenInfo=chooseEHPercToken(dirname,'HEAD');
% percTokenInfo.prodF0=expt.ehaeInfo.F0;
% percTokenInfo.x0=1.0;
% percTokenInfo.bAO=bAO;

% expt.percTokenInfo=percTokenInfo;
save(fullfile(dirname,'expt.mat'),'expt');

% eh_discrim([expt.subject.name,'_updown1'],'eh',expt.subject.percDir,6,expt.percTokenInfo,expt.percTokenInfo.x0,'twoScreens');
% eh_discrim([expt.subject.name,'_updown2'],'eh',expt.subject.percDir,6,expt.percTokenInfo,expt.percTokenInfo.x0,'twoScreens');
% eh_discrim([expt.subject.name,'_updown3'],'eh',expt.subject.percDir,6,percTokenInfo,x0,'twoScreens');
% eh_discrim([expt.subject.name,'_updown4'],'eh',expt.subject.percDir,6,percTokenInfo,x0,'twoScreens');
% eh_discrim([expt.subject.name,'_updown5'],'eh',expt.subject.percDir,6,percTokenInfo,x0,'twoScreens');
% eh_discrim([expt.subject.name,'_updown6'],'eh',expt.subject.percDir,6,percTokenInfo,x0,'twoScreens');

return

%% 
% function phaseScript = genPhaseScript(stage, nReps, preWords, trainWords, testWords, varargin)
function phaseScript = genPhaseScript(stage, nReps, words, varargin)
phaseScript=struct();
phaseScript.nReps=nReps;
phaseScript.nTrials=0;

for n=1:nReps
%     bt=[zeros(1,length(trainWords)),ones(1,round(length(trainWords)/2))];
%     if isequal(stage,'natural')
%         wordsUsed=varargin{1};
%         bt=[zeros(1,length(preWords))];
%     if isequal(stage,'pre')
%         wordsUsed=preWords(randperm(length(preWords)));
%         bt=[zeros(1,length(preWords))];
%     elseif isequal(stage,'test1') || isequal(stage,'test2')
%         wordsUsed=testWords(randperm(length(testWords)));
%         bt=[zeros(1,length(testWords))];
%     else
%         wordsUsed=trainWords(randperm(length(trainWords)));
%         bt=[zeros(1,length(trainWords))];
%     end
    wordsUsed = words(randperm(length(words)));
    bt = zeros(1, length(words));
    
%     pseudoWordsUsed=pseudoWords(randperm(length(pseudoWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));            
    twCnt=1;
    bt=bt(randperm(length(bt)));
    oneRep=struct;
    oneRep.trialOrder=[];
    oneRep.word=cell(1,0);
    cntTW=1;
    for m=1:length(bt)
%         if (bt(m)==0)
        if isequal(stage, 'noise')
            oneRep.trialOrder=[oneRep.trialOrder, 2];
        else
            oneRep.trialOrder=[oneRep.trialOrder, 1];
        end
            
        oneRep.word{length(oneRep.word)+1}=wordsUsed{twCnt};
        twCnt=twCnt+1;
%         elseif (bt(m)==1)
%             oneRep.trialOrder=[oneRep.trialOrder,[5,4,4]];
%             oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
%             oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW);
%             oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
%             cntTW=cntTW+2;
%         end
    end

    if (isequal(stage,'pract1') || isequal(stage,'pract2'))
        idx=find(oneRep.trialOrder<4);
        oneRep.trialOrder=oneRep.trialOrder(idx);
        oneRep.word=oneRep.word(idx);
    end

%     if n==nReps
%         oneRep.trialOrder=[oneRep.trialOrder,4];    % Dummy trial at the end
%         oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(1);
%     end
    phaseScript.(['rep',num2str(n)])=oneRep;
    phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);
end
return