function rData = init_data(type, phases, rawDataDir, subjID)
rData.rawDataFNs=cell(1,0);
rData.phases=cell(1,0);
rData.blockNums=[];
rData.trialNums=[];
rData.words=cell(1,0);
rData.datenums=[];

if isequal(type, 'rand')
    rData.pertType = [];
end

% sustData.rawDataFNs=cell(1,0);
% sustData.phases=cell(1,0);
% sustData.blockNums=[];
% sustData.trialNums=[];
% sustData.words=cell(1,0);
% sustData.datenums=[];
for i1 = 1 : numel(phases)
    t_phase = phases{i1};
    if isdir(fullfile(rawDataDir, subjID, t_phase))
        d1 = dir(fullfile(rawDataDir, subjID, t_phase,'rep*'));
        for i2=1:numel(d1)
            t_dir=fullfile(rawDataDir,subjID,t_phase,d1(i2).name);
            d2=dir(fullfile(t_dir,'trial-*-1.mat'));
            for i3=1:numel(d2)
                t_fn=fullfile(t_dir,d2(i3).name);
                load(t_fn);     % gives data
                
                rData.rawDataFNs{end+1}=t_fn;
                rData.phases{end+1}=t_phase;
                rData.blockNums(end+1)=str2num(strrep(d1(i2).name,'rep',''));
                rData.trialNums(end+1)=str2num(strrep(strrep(d2(i3).name,'trial-',''),'-1.mat',''));
                rData.words{end+1}=data.params.name;
                rData.datenums(end+1)=datenum(data.timeStamp);
                
                if isequal(type, 'rand')
                    t_fn = fullfile(t_dir, d2(i3).name);
                    load(t_fn);
                    
                    if data.params.bShift == 0
                        rData.pertType(end + 1) = 0;
                    else
                        if data.params.pertPhi(1) > 0
                            rData.pertType(end + 1) = 1; % higher: F1 down, F2 up
                        else
                            rData.pertType(end + 1) = -1; % lower: F1 up, F2 down
                        end
                    end
                end
                
            end
        end
    end
end

%% Data fields
dataFlds = {'rmsThresh',  0;
          'fn1',        0;
          'fn2',        0;
          'aFact',      0;
          'bFact',      0;
          'gFact',      0;
          'bCepsLift',  0;
          'cepsWinWidth', 0;
          'nLPC',       0;
          'vowelOnset',     0;
          'vowelEnd',       0;
          'vowelOnsetIdx',  0;
          'vowelEndIdx',    0;
          'f1Traj',         1;
          'f2Traj',         1;
          'prodF1',         0;
          'prodF2',         0;
          'prodF1_LB',      0;
          'prodF2_LB',      0;
          'prodF1_shira',   0;
          'prodF2_shira',   0;
          'prodF1_mnlBound',    0;
          'prodF2_mnlBound',    0;
          'audF1',      0;
          'audF2',      0;
          'traj_F1',    1;
          'traj_F2',    1;
          'sigRMS',     1;
          'iv1',        0;
          'iv2',        0;
          'bDiscard',   0;
          'bPertOkay',  0};

for i1 = 1 : size(dataFlds, 1)
    fld = dataFlds{i1, 1};
    isc = dataFlds{i1, 2};
    
    if isc == 1
        rData.(fld) = cell(size(rData.rawDataFNs));
    else
        rData.(fld) = nan(size(rData.rawDataFNs));
    end
end

rData.mark_f1_beg = 0 * ones(size(rData.rawDataFNs));
rData.mark_f1_mid = 0 * ones(size(rData.rawDataFNs));
rData.mark_f1_end = 0 * ones(size(rData.rawDataFNs));
rData.mark_f2_beg = 0 * ones(size(rData.rawDataFNs));
rData.mark_f2_mid = 0 * ones(size(rData.rawDataFNs));
rData.mark_f2_end = 0 * ones(size(rData.rawDataFNs));
 
rData.rating=2 * ones(size(rData.rawDataFNs));
rData.comments=cell(size(rData.rawDataFNs));

rData.bPertOkay=nan(size(rData.rawDataFNs));

return