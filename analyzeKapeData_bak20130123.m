function out = analyzeKapeData(subjID, varargin)
%% Options: data analysis
RAND_MAX_T = 0.250; % Unit: sec
FRAME_DUR = 24 / 12000;

%% Options: colors
colors.higher = [1, 0, 0];
colors.lower = [0, 0, 1];
colors.noPert = [0, 0, 0];

%%

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
elseif ~isempty(fsic(varargin, 'bu')) || ~isempty(fsic(varargin, 'boston'))
    dacacheDir = [dacacheDir, '_boston'];
end

bPlot = 1;
if ~isempty(fsic(varargin, 'noPlot'))
    bPlot = 0;
end

%%
pdataFN = fullfile(dacacheDir, [subjID, '.mat']);

load(pdataFN);


%% Analyze and visualize rand data
TRAJ_N = RAND_MAX_T / FRAME_DUR;
f1Mat = struct;
f2Mat = struct;

f1Mat.noPert = nan(TRAJ_N, 0);
f1Mat.higher = nan(TRAJ_N, 0);
f1Mat.lower = nan(TRAJ_N, 0);

f2Mat.noPert = nan(TRAJ_N, 0);
f2Mat.higher = nan(TRAJ_N, 0);
f2Mat.lower = nan(TRAJ_N, 0);

nRndNoPert = numel(find(pdata.randData.pertType == 0));
rndNoPert.prodF1_shira = nan(1, nRndNoPert);
rndNoPert.prodF2_shira = nan(1, nRndNoPert);
rndNoPert.rawDataFNs = cell(1, nRndNoPert);
% f1Mat. = nan(0, TRAJ_N);

nrnp_cnt = 0;
for i1 = 1 : numel(pdata.randData.rawDataFNs)
    if pdata.randData.pertType(i1) == 0
        nrnp_cnt = nrnp_cnt + 1;
    end
    
    if pdata.randData.bDiscard(i1) == 1
        continue;
    end
    if pdata.randData.rating(i1) == 0
        continue;
    end
    
    if pdata.randData.pertType(i1) == 0
        pert = 'noPert';
    elseif pdata.randData.pertType(i1) == 1
        pert = 'higher';
    else
        pert = 'lower';
    end
    
    if length(pdata.randData.f1Traj{i1}) >= TRAJ_N
        t_f1Traj = pdata.randData.f1Traj{i1}(1 : TRAJ_N);
        t_f2Traj = pdata.randData.f2Traj{i1}(1 : TRAJ_N);
        f1Mat.(pert) = [f1Mat.(pert), t_f1Traj];
        f2Mat.(pert) = [f2Mat.(pert), t_f2Traj];   
    end
    
    if pdata.randData.pertType(i1) == 0
        if pdata.randData.rating(i1) > 0
            rndNoPert.prodF1_shira(nrnp_cnt) = pdata.randData.prodF1_shira(i1);
            rndNoPert.prodF2_shira(nrnp_cnt) = pdata.randData.prodF2_shira(i1);
        else
            rndNoPert.prodF1_shira(nrnp_cnt) = NaN;
            rndNoPert.prodF2_shira(nrnp_cnt) = NaN;
        end
        rndNoPert.rawDataFNs{nrnp_cnt} = pdata.randData.rawDataFNs{i1};
    end
end

avg_f1Mat = struct;
avg_f2Mat = struct;
pertFlds = fields(f1Mat);
for i1 = 1 : numel(pertFlds)
    fld = pertFlds{i1};
    
    avg_f1Traj.(fld) = mean(transpose(f1Mat.(fld)));
    avg_f2Traj.(fld) = mean(transpose(f2Mat.(fld)));
end

if bPlot
    figure('Name', 'Rand: F1 changes');
    hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (TRAJ_N - 1);
    plot(tAxis, avg_f1Traj.higher - avg_f1Traj.noPert, ...
         '-', 'Color', colors.higher);
    plot(tAxis, avg_f1Traj.lower - avg_f1Traj.noPert, ...
         '-', 'Color', colors.lower);
    plot([tAxis(1), tAxis(end)], [0, 0], '-', 'color', [0.5, 0.5, 0.5]);
    xlabel('Time (s)');
    ylabel('F1 change from noPert (Hz)');
    legend({'Higher', 'Lower'});
    
    figure('Name', 'Rand: F2 changes');
    hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (TRAJ_N - 1);
    plot(tAxis, avg_f2Traj.higher - avg_f2Traj.noPert, ...
         '-', 'Color', colors.higher);
    plot(tAxis, avg_f2Traj.lower - avg_f2Traj.noPert, ...
         '-', 'Color', colors.lower);
    plot([tAxis(1), tAxis(end)], [0, 0], '-', 'color', [0.5, 0.5, 0.5]);
    xlabel('Time (s)');
    ylabel('F2 change from noPert (Hz)');
    legend({'Higher', 'Lower'});
end

%% Analyze and visualize sust data
idxStart = strmatch('start', pdata.sustData.phases, 'exact');
idxRamp = strmatch('ramp', pdata.sustData.phases, 'exact');
idxStay = strmatch('stay', pdata.sustData.phases, 'exact');
idxEnd = strmatch('end', pdata.sustData.phases, 'exact');
ntStart = numel(idxStart);
ntRamp = numel(idxRamp);
ntStay = numel(idxStay);
ntEnd = numel(idxEnd);
ntTot = ntStart + ntRamp + ntStay + ntEnd;

sust_prodF1_shira = [];
sust_prodF2_shira = [];

hft = nan(1, 2);
for i1 = 1 : 2
    if i1 == 1
        meas = 'prodF1_shira';
        measName = 'F1 (Shira)';
    else
        meas = 'prodF2_shira';
        measName = 'F2 (Shira)';
    end
    unit = 'Hz';
    
    pdata.sustData.(meas)(pdata.sustData.rating == 0) = NaN;
    
    idx_start = fsic(pdata.sustData.phases, 'start');
    idx_ramp = fsic(pdata.sustData.phases, 'ramp');
    idx_stay = fsic(pdata.sustData.phases, 'stay');
    idx_end = fsic(pdata.sustData.phases, 'end');
    
    if isequal(meas, 'prodF1_shira')
        sust_prodF1_shira = [nanmean(pdata.sustData.(meas)(idx_start)), ...
                             nanmean(pdata.sustData.(meas)(idx_ramp)), ...
                             nanmean(pdata.sustData.(meas)(idx_stay)), ...
                             nanmean(pdata.sustData.(meas)(idx_end))];
    elseif isequal(meas, 'prodF2_shira')
        sust_prodF2_shira = [nanmean(pdata.sustData.(meas)(idx_start)), ...
                             nanmean(pdata.sustData.(meas)(idx_ramp)), ...
                             nanmean(pdata.sustData.(meas)(idx_stay)), ...
                             nanmean(pdata.sustData.(meas)(idx_end))];
    end
    
    if bPlot
        hft(i1) = figure('Name', measName);
        subplot(2, 1, 1);
        prodFmt_baseLine = nanmean([rndNoPert.(meas), pdata.sustData.(meas)(idxStart)]);
        
        plot([rndNoPert.(meas), pdata.sustData.(meas)], 'o-');
        set(gca, 'XLim', [0, nRndNoPert + ntTot + 1]);
        hold on;
        plot([0, nRndNoPert + ntTot + 1], repmat(prodFmt_baseLine, 1, 2), '-', ...
             'Color', [0.5, 0.5, 0.5]);
        ys = get(gca, 'YLim');
        plot(repmat(nRndNoPert + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
        plot(repmat(nRndNoPert + ntStart + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
        plot(repmat(nRndNoPert + ntStart + ntRamp + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
        plot(repmat(nRndNoPert + ntStart + ntRamp + ntStay + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);

        xlabel('Trial #');
        ylabel(sprintf('%s (%s)', measName, unit));
    end
    
    %% Use ginput to look at the outliers.
    if ~isempty(fsic(varargin, 'pick'))
        xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
        bMoveOn = 0;
        while (bMoveOn == 0)
            coord = ginput(1);
            if coord(1) > xs(2) || coord(1) < xs(1) || coord(2) > ys(2) || coord(2) < ys(1)
                bMoveOn = 1;
            else
                trn = round(coord(1));
                if trn <= length(rndNoPert.(meas))
                    trn1 = trn;
                    c_rawDataFN = rndNoPert.rawDataFNs{trn1};
                    c_F1 = rndNoPert.prodF1_shira(trn1);
                    c_F2 = rndNoPert.prodF2_shira(trn1);
                else
                    trn1 = trn - length(rndNoPert.(meas));
                    c_rawDataFN = pdata.sustData.rawDataFNs{trn1};
                    c_F1 = pdata.sustData.prodF1_shira(trn1);
                    c_F2 = pdata.sustData.prodF2_shira(trn1);
                end

                c_rawDataFN = strrep(c_rawDataFN, '.mat', '');
                [fp0, fp1] = fileparts(c_rawDataFN);
                c_trialN = str2num(strrep(strrep(fp1, 'trial-', ''), '-1', ''));
                [fp0, fp1] = fileparts(fp0);
                c_repN = str2num(strrep(fp1, 'rep', ''));
                [fp0, fp1] = fileparts(fp0);
                c_phase = fp1;

                cmd = sprintf('preprocKapeData(''%s'', ''phase'', ''%s'', ''rep'', %d, ''trial'', %d);', ...
                              subjID, c_phase, c_repN, c_trialN);
                fprintf('F1 = %f, F2 = %f, cmd = \n\t%s\n', c_F1, c_F2, cmd);
            end
        end
    end
end

sustF1s_shira = pdata.sustData.prodF1_shira;
sustF2s_shira = pdata.sustData.prodF2_shira;

extF1s_shira = [rndNoPert.prodF1_shira, sustF1s_shira];
extF2s_shira = [rndNoPert.prodF2_shira, sustF2s_shira];

%% Compute by-phase means and SEs in the sust data
for i1 = 1 : 2
    if i1 == 1
        bph_mn = [nanmean(extF1s_shira(1 : nRndNoPert + ntStart)), ...
                  nanmean(extF1s_shira(nRndNoPert + ntStart + 1 : nRndNoPert + ntStart + ntRamp)), ...
                  nanmean(extF1s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay)), ...
                  nanmean(extF1s_shira(nRndNoPert + ntStart + ntRamp + ntStay + 1 : nRndNoPert + ntStart + ntRamp + ntStay + ntEnd))];
        bph_se = [nanste(extF1s_shira(1 : nRndNoPert + ntStart)), ...
                  nanste(extF1s_shira(nRndNoPert + ntStart + 1 : nRndNoPert + ntStart + ntRamp)), ...
                  nanste(extF1s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay)), ...
                  nanste(extF1s_shira(nRndNoPert + ntStart + ntRamp + ntStay + 1 : nRndNoPert + ntStart + ntRamp + ntStay + ntEnd))];
              
        [h_t, p_t, ci, stats_t] = ttest2(extF1s_shira(1 : nRndNoPert + ntStart), ...
                                         extF1s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay));
        ylab = 'F1 (Shira) (Hz)';
    else
        bph_mn = [nanmean(extF2s_shira(1 : nRndNoPert + ntStart)), ...
                  nanmean(extF2s_shira(nRndNoPert + ntStart + 1 : nRndNoPert + ntStart + ntRamp)), ...
                  nanmean(extF2s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay)), ...
                  nanmean(extF2s_shira(nRndNoPert + ntStart + ntRamp + ntStay + 1 : nRndNoPert + ntStart + ntRamp + ntStay + ntEnd))];
        bph_se = [nanste(extF2s_shira(1 : nRndNoPert + ntStart)), ...
                  nanste(extF2s_shira(nRndNoPert + ntStart + 1 : nRndNoPert + ntStart + ntRamp)), ...
                  nanste(extF2s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay)), ...
                  nanste(extF2s_shira(nRndNoPert + ntStart + ntRamp + ntStay + 1 : nRndNoPert + ntStart + ntRamp + ntStay + ntEnd))];
              
        [h_t, p_t, ci, stats_t] = ttest2(extF2s_shira(1 : nRndNoPert + ntStart), ...
                                         extF2s_shira(nRndNoPert + ntStart + ntRamp + 1 : nRndNoPert + ntStart + ntRamp + ntStay));
        ylab = 'F2 (Shira) (Hz)';
    end
         
    set(0, 'CurrentFigure', hft(i1));
    subplot(2, 1, 2);
    errorbar(bph_mn, bph_se, 'b-');
    hold on;
    plot([0.5, length(bph_mn) + 0.5], repmat(bph_mn(1), 1, 2), '-', 'Color', [0.5, 0.5, 0.5]);
    ylabel(ylab);
    xlabel('Phase');
    
    set(gca, 'XLim', [0.5, length(bph_mn) + 0.5]); 
    set(gca, 'XTick', [1 : length(bph_mn)], 'XTickLabel', {'Base', 'Ramp', 'Stay', 'End'});
    
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.10 * range(ys), ...
         sprintf('Base v. Stay: t = %.3f, p = %.3f', stats_t.tstat, p_t));
end

%% Output
out = struct;
out.rndNoPert = rndNoPert;
out.avg_f1Traj = avg_f1Traj;
out.avg_f2Traj = avg_f2Traj;
out.extF1s_shira = extF1s_shira;
out.extF2s_shira = extF2s_shira;
out.sust_prodF1_shira = sust_prodF1_shira;
out.sust_prodF2_shira = sust_prodF2_shira;

return