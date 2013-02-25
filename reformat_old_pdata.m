function reformat_old_pdata(pdataFN)
load(pdataFN); % gives pdata

datFlds = {'otherData', 'randData', 'sustData'};
nTotTrials = 0;
for i1 = 1 : numel(datFlds)
    fld = datFlds{i1};
    if ~isfield(pdata, fld)
        fprintf('WARNING: not field %s in pdata. Skipped.\n', fld);
        continue;
    end
    
    nTrials = numel(pdata.(fld).rawDataFNs);
    nTotTrials = nTotTrials + nTrials;
    
    pdata.(fld).bCepsLift = nan(1, nTrials);
    pdata.(fld).cepsWinWidth = nan(1, nTrials);
    pdata.(fld).mark_f1_beg = nan(1, nTrials);
    pdata.(fld).mark_f1_mid = nan(1, nTrials);
    pdata.(fld).mark_f1_end = nan(1, nTrials);
    pdata.(fld).mark_f2_beg = nan(1, nTrials);
    pdata.(fld).mark_f2_mid = nan(1, nTrials);
    pdata.(fld).mark_f2_end = nan(1, nTrials);
    pdata.(fld).nLPC = nan(1, nTrials);
    
    for i2 = 1 : nTrials
        rawfn = pdata.(fld).rawDataFNs{i2};
        if ~isempty(strfind(rawfn, 'D:'))
            if isequal(getHostName, 'smcg_w510')
                rawfn = strrep(rawfn, 'D:', 'E:');
            end
        end
        load(rawfn); % gives data
        
        pdata.(fld).bCepsLift(i2) = data.params.bCepsLift;
        pdata.(fld).cepsWinWidth(i2) = data.params.cepsWinWidth;
        pdata.(fld).mark_f1_beg(i2) = 0;
        pdata.(fld).mark_f1_mid(i2) = 0;
        pdata.(fld).mark_f1_end(i2) = 0;
        pdata.(fld).mark_f2_beg(i2) = 0;
        pdata.(fld).mark_f2_mid(i2) = 0;
        pdata.(fld).mark_f2_end(i2) = 0;
        pdata.(fld).nLPC(i2) = data.params.nLPC;
    end
end

%% Re-format state
stateFN = strrep(pdataFN, '.mat', '_state.mat');
load(stateFN); % gives state

if isfield(state, 'idx')
    rmfield(state, 'idx');
end
state.rawDataDir = fileparts(fileparts(fileparts(fileparts(rawfn))));
state.dacacheDir = 'E:/speechres/kape/dacache';
state.expDir = fileparts(fileparts(fileparts(rawfn)));
state.persist_rmsThresh = NaN;
state.bFirstTime = 1;
state.stats = ones(1, length(state.trialList.fn));
state.statsPert = zeros(1, length(state.trialListPert.fn));


%% Write to file
save(pdataFN, 'pdata');
fprintf('New pdata written to file: %s\n', pdataFN);

save(stateFN, 'state');
fprintf('New state written to file: %s\n', stateFN);

return