function p = getTSMDefaultParams(sex,varargin)
switch sex
    case 'male'
        p.nLPC          = 11; 
        p.fn1           = 591;
        p.fn2           = 1314;
    case 'female'
        p.nLPC          = 9;
        p.fn1           = 700;
        p.fn2           = 1800;
    otherwise,
        error('specify sex (male / female');
end

if ~isempty(fsic(varargin, 'subject_group'))
    subjectGroup = varargin{fsic(varargin, 'subject_group') + 1};
    if isequal(subjectGroup, 'CNS') || isequal(subjectGroup, 'CWS')
        p.nLPC = 9;
        p.fn1 = 1000;
        p.fn2 = 2000;
    end
end

p.aFact         = 1;
p.bFact         = 0.8;
p.gFact         = 1;

p.downfact      = 4;
if ~isempty(fsic(varargin,'downFact'))
    p.downfact=varargin{fsic(varargin,'downFact')+1};
end

% p.closedLoopGain= 15;    % dB
p.closedLoopGain= 0;    % dB

if ~isempty(findStringInCell(varargin,'closedLoopGain'))
    p.closedLoopGain=varargin{findStringInCell(varargin,'closedLoopGain')+1};
end

p.dScale        = 10^((p.closedLoopGain-calcClosedLoopGain)/20);
p.preempFact    = 0.98;% preemp factor

p.sr            = 48000/p.downfact;
if ~isempty(fsic(varargin,'sr'))
    p.sr=varargin{fsic(varargin,'sr')+1};
end

% Frame structure
p.nWin          = 1;% 1 2 4  8 16 32 64 (max=p.framLen) Number of windows per frame  !!

p.frameLen      = 64 / p.downfact;% must be a valid DMA Buffer size (64 128 256 ..)
if ~isempty(fsic(varargin,'frameLen'))
    p.frameLen=varargin{fsic(varargin,'frameLen')+1};
end

p.nDelay        = 5;% the total process delay is: p.frameLen*p.nDelay/p.sr
p.frameShift    = p.frameLen/p.nWin;% 
p.bufLen        = (2*p.nDelay-1)*p.frameLen;
p.anaLen        = p.frameShift+2*(p.nDelay-1)*p.frameLen;
p.avgLen      = 8;    %ceil(p.sr/(f0*p.frameShift));

p.bCepsLift     = 0;

p.minVowelLen   = 60;

if (isequal(sex,'male'))
    p.cepsWinWidth  = 50;
elseif (isequal(sex,'female'))
    p.cepsWinWidth  = 30;    
end

% Formant tracking 
p.nFmts         = 2;
p.nTracks       = 4;
p.bTrack        = 1;
p.bWeight       = 1; % weigthing (short time rms) of moving average formant estimate o

% RMS calculation
if ~isempty(findStringInCell(varargin,'mouthMicDist'))
    p.mouthMicDist=varargin{findStringInCell(varargin,'mouthMicDist')+1};
else
    p.mouthMicDist=10;
end

p.rmsThresh     = 0.06 * 10^((getSPLTarg(p.mouthMicDist)-85)/20); % Before: 0.04*10^((getSPLTarg('prod')-85)/20); % 2009/11/27, changed from 0.04* to 0.032*
p.rmsRatioThresh= 1.3;% threshold for sibilant / vowel detection
p.rmsMeanPeak   = 6*p.rmsThresh;
p.rmsForgFact   = 0.95;% forgetting factor for rms computation

% Vowel detection
p.bDetect       = 1;
%p.fmtsForgFact  = 0.97;% formants forgetting factor [1 2 .. nFmts]
p.dFmtsForgFact = 0.93;% formants forgetting factor for derivate calculation

% Shifting
p.bShift        = 0;
p.bRatioShift   = 1;
p.bMelShift     = 0;    % Use mel as the unit

p.gainAdapt     =0;

% Trial-related
p.fb=1; % Voice only;

p.trialLen=3; %SC(2008/06/22) % Maximum allowable: 

if ~isempty(findStringInCell(varargin,'trialLen'))
    p.trialLen=varargin{findStringInCell(varargin,'trialLen')+1};
end

p.rampLen=0.25; %SC(2008/06/22)

% SC(2009/02/06) RMS clipping protection
p.bRMSClip=1;
% p.rmsClipThresh=1.0;

load('../../signals/leveltest/micRMS_100dBA.mat');
p.rmsClipThresh=micRMS_100dBA / (10^((100-100)/20));	% 100 dB maximum input level

%% Perturbation-related variables: these are for the mel (bMelShift=1) or Hz (bMelShift=0) frequency space
p.F2Min=0;
p.F2Max=5000;
p.F1Min=0;
p.F1Max=5000;
p.LBb=0;
p.LBk=0;

p.pertFieldN=257;
p.pertF2=zeros(1,p.pertFieldN);
p.pertAmp=zeros(1,p.pertFieldN);
p.pertPhi=zeros(1,p.pertFieldN);

return
