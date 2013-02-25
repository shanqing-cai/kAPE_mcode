function pdata_fixAutoRMS(subjID, varargin)
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

%%
dacacheFN = fullfile(dacacheDir, [subjID, '.mat']);
stateFN = strrep(dacacheFN, '.mat', '_state.mat');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

%%
ptFlds = {'otherData', 'randData', 'sustData'};
for i1 = 1 : numel(ptFlds)
    ptFld = ptFlds{i1};
    if ~isfield(pdata, ptFld)
        fprintf('WARNING: pdata has no field %s. Skipped this part\n', ptFld);
    else
        nTrials = numel(pdata.(ptFld).rawDataFNs);
        for i2 = 1 : nTrials
            t_rawDataFN = pdata.(ptFld).rawDataFNs{i2};
            t_phase = pdata.(ptFld).phases{i2};
            t_blockNum = pdata.(ptFld).blockNums(i2);
            t_trialNum = pdata.(ptFld).trialNums(i2);
            
            idx_trial = i2;
            auto_rmsThresh_cbk([], [], dacacheFN, stateFN, ...
                               {t_rawDataFN, ptFld, idx_trial, dacacheFN, stateFN, t_phase, t_blockNum, t_trialNum});
        end
    end
end
return