function g_kapeAnalysis(subjListFN, varargin)
%% Constants:
SHIFT_RATIO_SUST_F1 = +0.25;
SHIFT_RATIO_SUST_F2 = -0.125;

infoXlsFN = 'E:/DATA/KAPE/kAPE_Info_&_Data_6-19-13_LW.xls';

fontSize = 12;
CORP_THRESH = 0.05;

%% Options: subjects
% SID.ANS = {'PILOT_M01', 'PILOT_F03', 'PILOT_M02', ...
%            'ANS_F05', 'ANS_F06', 'ANS_F07', 'ANS_F08', 'ANS_F09', 'ANS_F10','ANS_F11', ...
%            'ANS_M03', 'ANS_M04', 'ANS_M05', 'ANS_M06'};  
%          % Noise interference, 'PILOT_F04'; old paradigm: 'PILOT_F01', 'PILOT_F02'
% 
% SID.CNS = {'CNS_F01', 'CNS_F02', 'CNS_F03', 'CNS_F04', ...
%            'CNS_F05', 'CNS_F06', 'CNS_F07', ...
%            'CNS_M02', 'CNS_M03', 'CNS_M04', 'CNS_M05', ...
%            'CNS_M06', 'CNS_M07', 'CNS_M08'}; 
%        
%          % Bad formant tracking and old paradigm: 'CNS_M01'   
% 
% SID.CWS = {'CWS_F01', 'CWS_F02', 'CWS_F03', 'CWS_F04', ...
%            'CWS_F05', 'CWS_F06', 'CWS_F07', ...
%            'CWS_M01', 'CWS_M02', 'CWS_M03', 'CWS_M04', ...
%            'CWS_M05', 'CWS_M06', 'CWS_M06'};

% --- Read subject lists --- %
sList = read_struct_from_text(subjListFN);
sgrps = fields(sList);
SID = struct;
for i1 = 1 : numel(sgrps)
    sgrp = sgrps{i1};
    
    SID.(sgrp) = splitstring(strrep(sList.(sgrp), ',', ' '));
    
    fprintf(1, 'INFO: Subject list file contains %d subjects in group %s\n', ...
            numel(SID.(sgrp)), sgrp);
end      

%% Load additional information from xls file
check_file(infoXlsFN);
[N, T] = xlsread(infoXlsFN);

age_mo = struct;
for i1 = 1 : numel(sgrps)
    grp = sgrps{i1};
    age_mo.(grp) = nan(1, length(SID.(grp)));
end
SSI_tot.CWS = nan(1, length(SID.CWS));

for i0 = 1 : numel(sgrps)
    grp = sgrps{i0};
    
    for i1 = 1 : length(SID.(grp))
        if isequal(grp, 'CNS') || isequal(grp, 'CWS')
            age_mo.(grp)(i1) = get_kape_subj_info(N, T, SID.(grp){i1}, 'Age: kAPE (mo)');
        end
        
        if isequal(grp, 'CWS')
            SSI_tot.(grp)(i1) = get_kape_subj_info(N, T, SID.(grp){i1}, 'SSI - Total');
        end
    end
end

%% Options: colors and other visualization options
colors.higher = [1, 0, 0];
colors.lower = [0, 0, 1];
colors.noPert = [0, 0, 0];

colors.ANS = [0, 0, 0];
colors.CNS = [0, 0, 1];
colors.CWS = [1, 0, 0];

FRAME_DUR = 24 / 12000;
EXT_LEN = 132;
EXT_BASE_LEN = 49 + 3;
LEN_RAND_NOPERT = 49;
LEN_START = 3;
LEN_RAMP = 18;
LEN_STAY = 36;
LEN_END = 18;
randShiftDirs = {'higher', 'lower'};

FRAME_DUR = 24 / 12000;

barW = 0.5;

%% Process input arguments
gend = 'both';
if ~isempty(fsic(varargin, 'maleOnly'))
    gend = 'male';
    fprintf(1, 'INFO: including data from only %s subjects\n', gend);
elseif ~isempty(fsic(varargin, 'femaleOnly'))
    gend = 'female';
    fprintf(1, 'INFO: including data from only %s subjects\n', gend);
end

grps = fields(SID);
if isequal(gend, 'male') || isequal(gend, 'female')
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        bKeep = zeros(1, numel(SID.(grp)));
        for i2 = 1 : numel(SID.(grp))
            idxh = strfind(SID.(grp){i2}, '_');
            
            gchar = lower(SID.(grp){i2}(idxh + 1));
            if isequal(gchar, lower(gend(1)))
                bKeep(i2) = 1;
            end
        end
        
        SID.(grp) = SID.(grp)(find(bKeep));
    end
end

bNoXLS = ~isempty(fsic(varargin, '--noXLS'));

nPerm = 0;
if ~isempty(fsic(varargin, '--perm'))
    nPerm = varargin{fsic(varargin, '--perm') + 1};
end

testType = 'rs';     % {'t', 'rs'} (t-test or rank-sum / signed-rank tests)
if ~isempty(fsic(varargin, '--test'))
    testType = varargin{fsic(varargin, '--test') + 1};
end

%%
g_rndF1TrajChg = struct;

g_extF1_shira = struct;
g_extF2_shira = struct;
g_extFp_shira = struct;

g_extF1Chg_shira = struct;
g_extF2Chg_shira = struct;

gp_extF1Chg_shira = struct;
gp_extF2Chg_shira = struct;
gp_extF12Chg_shira = struct;

meanVDur_rand = struct;
meanVDur_sust = struct;

meanF1s_ext = struct;
meanF2s_ext = struct;


% --- Composite F1-F2 compensation under higher and lower perturbations at
% the last time point --- %
g_rndFp_comp = struct;
g_rndF12_comp = struct;

age_yo = struct;

ntp = NaN;
pertVec = [SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2];
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    age_yo.(grp) = [];
    
    meanF1s_ext.(grp) = [];
    meanF2s_ext.(grp) = [];
    
    meanVDur_rand.(grp) = nan(0, 3); % [noPert, higher, lower]    
    meanVDur_sust.(grp) = nan(0, 4); % [Start, Ramp, Stay, End]
    
    g_rndF1TrajChg.(grp).higher = [];
    g_rndF1TrajChg.(grp).lower = [];
    g_rndF2TrajChg.(grp).higher = [];
    g_rndF2TrajChg.(grp).lower = [];
    
    g_rndFpTrajChg.(grp).higher = [];
    g_rndFpTrajChg.(grp).lower = [];
    
    g_rndFpTrajChgContra.(grp) = [];
    
    g_rndFp_comp.(grp) = [];
    g_rndF12_comp.(grp) = [];
    
    g_extF1_shira.(grp) = [];
    g_extF2_shira.(grp) = [];
    g_extFp_shira.(grp) = [];
    
    g_extF1Chg_shira.(grp) = [];
    g_extF2Chg_shira.(grp) = [];
    
    gp_extF1Chg_shira.(grp) = [];
    gp_extF2Chg_shira.(grp) = [];
    
    for i2 = 1 : numel(SID.(grp))
        sID = SID.(grp){i2};
        
        out = analyzeKapeData(sID, 'noPlot');
        age_yo.(grp)(end + 1) = out.age_yo;
        if ~isequal(grp, out.group)
            error('Mismatch between out.group = %ds and grp = %s\n', out.group, grp);
        end
        meanF1s_ext.(grp)(end + 1) = nanmean(out.extF1s_shira);
        meanF2s_ext.(grp)(end + 1) = nanmean(out.extF2s_shira);
        
        meanVDur_rand.(grp)(end + 1, :) = [nanmean(out.vDur_rand.noPert), nanmean(out.vDur_rand.higher), nanmean(out.vDur_rand.lower)];
        meanVDur_sust.(grp)(end + 1, :) = [nanmean(out.vDur_sust.start), nanmean(out.vDur_sust.ramp),  nanmean(out.vDur_sust.stay), nanmean(out.vDur_sust.end)];
        
%         if i1 == 2 && i2 == 14 % DEBUG
%             pause(0);
%         end            
        if length(out.avg_f1Traj.noPert) ~= length(out.avg_f1Traj.higher) || ...
           length(out.avg_f1Traj.noPert) ~= length(out.avg_f1Traj.lower) % DEBUG
            fprintf(2, 'WARNING: Mismatch in vector lengths in out.avg_f1Traj in group %s, subject %s\n', grp, sID);
        end
        
        if ~isempty(find(isnan(out.avg_f1Traj.noPert)))
            fprintf(2, 'WARNING: NaN found in out.avg_f1Traj in group %s, subject %s\n', grp, sID);
        end
        
        % --- Projection of the rand data --- %
        flds = fields(out.avg_f1Traj);
        for i1 = 1 : numel(flds)
            fld = flds{i1};
            out.avg_fpTraj.(fld) = proj2PertLine(SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2, ...
                                       out.avg_f1Traj.(fld), out.avg_f2Traj.(fld));
        end
        % --- ~Projection of the rand data --- %
        
        
        for i3 = 1 : numel(randShiftDirs)
            sDir = randShiftDirs{i3};
            g_rndF1TrajChg.(grp).(sDir) = [g_rndF1TrajChg.(grp).(sDir); ...
                   (out.avg_f1Traj.(sDir) - out.avg_f1Traj.noPert) ./ out.avg_f1Traj.noPert];
            g_rndF2TrajChg.(grp).(sDir) = [g_rndF2TrajChg.(grp).lower; ...
                   (out.avg_f2Traj.(sDir) - out.avg_f2Traj.noPert) ./ out.avg_f2Traj.noPert];
               
            g_rndFpTrajChg.(grp).(sDir) = [g_rndFpTrajChg.(grp).(sDir); ...
                   out.avg_fpTraj.(sDir) - out.avg_fpTraj.noPert];            
        end
        
        g_rndFpTrajChgContra.(grp) = [g_rndFpTrajChgContra.(grp); ...
                   out.avg_fpTraj.higher - out.avg_fpTraj.lower];
        g_rndFp_comp.(grp)(end + 1) = g_rndFpTrajChgContra.(grp)(end, end);
        g_rndF12_comp.(grp)(end + 1, :) = pertVec * [(nanmean(out.avg_f1Traj.higher(end - 0 : end)) - nanmean(out.avg_f1Traj.lower(end - 0 : end))) / nanmean(out.avg_f1Traj.noPert), ...
                                                  (nanmean(out.avg_f2Traj.higher(end - 0 : end)) - nanmean(out.avg_f2Traj.lower(end - 0 : end))) / nanmean(out.avg_f2Traj.noPert)]' ...
                                       / norm(pertVec);
        
        g_extF1_shira.(grp) = [g_extF1_shira.(grp); out.sust_prodF1_shira];
        g_extF2_shira.(grp) = [g_extF2_shira.(grp); out.sust_prodF2_shira];
        
        approx_pert_mag = norm([mean(out.sust_prodF1_shira(2 : 3)) * SHIFT_RATIO_SUST_F1, ...
                                mean(out.sust_prodF2_shira(2 : 3)) * SHIFT_RATIO_SUST_F2]);
%         g_extFp_shira.(grp) = [g_extFp_shira.(grp); ...
%                                proj2PertLine(SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2, ...
%                                              out.sust_prodF1_shira, out.sust_prodF2_shira) / approx_pert_mag];
        g_extFp_shira.(grp) = [g_extFp_shira.(grp); ...
                               proj2PertLine(SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2, ...
                                             out.sust_prodF1_shira, out.sust_prodF2_shira)];
        
        if length(out.extF1s_shira) == EXT_LEN
            t_f1s = out.extF1s_shira;
            t_f1c = (t_f1s - nanmean(t_f1s(1 : EXT_BASE_LEN))) / nanmean(t_f1s(1 : EXT_BASE_LEN));
            t_f2s = out.extF2s_shira;
            t_f2c = (t_f2s - nanmean(t_f2s(1 : EXT_BASE_LEN))) / nanmean(t_f2s(1 : EXT_BASE_LEN));
        else
            t_f1s = [nan(1, EXT_LEN - length(out.extF1s_shira)), out.extF1s_shira];
            t_f1c = (t_f1s - nanmean(t_f1s(1 : EXT_BASE_LEN))) / nanmean(t_f1s(1 : EXT_BASE_LEN));
            t_f2s = [nan(1, EXT_LEN - length(out.extF2s_shira)), out.extF2s_shira];
            t_f2c = (t_f2s - nanmean(t_f2s(1 : EXT_BASE_LEN))) / nanmean(t_f2s(1 : EXT_BASE_LEN));
        end
        g_extF1Chg_shira.(grp) = [g_extF1Chg_shira.(grp); t_f1c];
        g_extF2Chg_shira.(grp) = [g_extF2Chg_shira.(grp); t_f2c];
        
        gp_extF1Chg_shira.(grp) = [gp_extF1Chg_shira.(grp); ...
            (out.sust_prodF1_shira - out.sust_prodF1_shira(1)) / out.sust_prodF1_shira(1)];
        gp_extF2Chg_shira.(grp) = [gp_extF2Chg_shira.(grp); ...
            (out.sust_prodF2_shira - out.sust_prodF2_shira(1)) / out.sust_prodF2_shira(1)];

        if isnan(ntp)
            ntp = length(out.avg_f1Traj.noPert);
        end
    end
    
    % --- Plot Rand F1, F2, and Fp change trajectories --- %
    figure('Name', sprintf('rand - F1, F2 & Fp - %s', grp), ...
           'Position', [100, 100, 1350, 500]);
       
	for j1 = 1 : 3 
        subplot(1, 3, j1);

        if j1 == 1
            t_dat = g_rndF1TrajChg.(grp);
            dat_name = 'F1';
        elseif j1 == 2
            t_dat = g_rndF2TrajChg.(grp);
            dat_name = 'F2';
        else
            t_dat = g_rndFpTrajChg.(grp);
            dat_name = 'Fp';
        end
        
        hold on;
        tAxis = 0 : FRAME_DUR : FRAME_DUR * (ntp - 1);
        plot(tAxis, zeros(size(tAxis)), '-', 'color', [0.5, 0.5, 0.5]);
        for i2 = 1 : numel(randShiftDirs)
            sDir = randShiftDirs{i2};
            plot(tAxis, nanmean(t_dat.(sDir)), ...
                 'Color', colors.(sDir), 'LineWidth', 1.5);
            plot(tAxis, nanmean(t_dat.(sDir)) + nanste1(t_dat.(sDir)), ...
                 '--', 'Color', colors.(sDir), 'LineWidth', 1);
            plot(tAxis, nanmean(t_dat.(sDir)) - nanste1(t_dat.(sDir)), ...
                 '--', 'Color', colors.(sDir), 'LineWidth', 1);
        end
        xlabel('Time (s)');
        ylabel(sprintf('%s change (ratio)', dat_name));
        title(sprintf('rand - %s - %s', dat_name, grp));

        ezlegend([NaN, 0.6, 0.75, 0.3, 0.2], 0.45, ...
                 randShiftDirs, {colors.(randShiftDirs{1}), colors.(randShiftDirs{2})}, ...
                 [12, 12], {'-', '-', '.-'}, {colors.(randShiftDirs{1}), colors.(randShiftDirs{2})}, ...
                 [1.5, 1.5], [1, 1]);
    end
    
    % --- Plot Rand F2 change trajectories --- %
%     figure('Name', sprintf('rand - F2 - %s', grp));
%     subplot(1, 2, 2);
% 
%     hold on;
%     tAxis = 0 : FRAME_DUR : FRAME_DUR * (ntp - 1);
%     plot(tAxis, zeros(size(tAxis)), '-', 'color', [0.5, 0.5, 0.5]);
%     for i2 = 1 : numel(randShiftDirs)
%         sDir = randShiftDirs{i2};
%         plot(tAxis, nanmean(g_rndF2TrajChg.(grp).(sDir)), ...
%              'Color', colors.(sDir), 'LineWidth', 1.5);
%         plot(tAxis, nanmean(g_rndF2TrajChg.(grp).(sDir)) + nanste1(g_rndF2TrajChg.(grp).(sDir)), ...
%              '--', 'Color', colors.(sDir), 'LineWidth', 1);
%         plot(tAxis, nanmean(g_rndF2TrajChg.(grp).(sDir)) - nanste1(g_rndF2TrajChg.(grp).(sDir)), ...
%              '--', 'Color', colors.(sDir), 'LineWidth', 1);
%     end
%     xlabel('Time (s)');
%     ylabel('F2 change (ratio)');
%     title(sprintf('rand - F2 - %s', grp));
    
%     ezlegend([NaN, 0.6, 0.75, 0.3, 0.2], 0.45, ...
%              randShiftDirs, {colors.(randShiftDirs{1}), colors.(randShiftDirs{2})}, ...
% 	         [12, 12], {'-', '-', '.-'}, {colors.(randShiftDirs{1}), colors.(randShiftDirs{2})}, ...
%              [1.5, 1.5], [1, 1]);
    
    % --- Plot Sust F1 changes --- %
    figure('Name', sprintf('sust - F1 - %s', grp), ...
           'Position', [100, 100, 1200, 500]);
    subplot(1, 2, 1);
    hold on;
    
    plot([0, EXT_LEN + 1], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    plot(nanmean(g_extF1Chg_shira.(grp)), '-');
    plot(nanmean(g_extF1Chg_shira.(grp)) - nanste1(g_extF1Chg_shira.(grp)), ...
         '--');
    plot(nanmean(g_extF1Chg_shira.(grp)) + nanste1(g_extF1Chg_shira.(grp)), ...
         '--');
    set(gca, 'XLim', [0, EXT_LEN + 1]);
    ys = get(gca, 'YLim');
    plot(repmat(LEN_RAND_NOPERT + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    xlabel('Trial #');
    ylabel('Sust. F1 change (fraction)');
%     plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + LEN_END + 1, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);

    % --- Plot Sust F2 changes --- %
%     figure('Name', sprintf('sust - F2 - %s', grp));
    subplot(1, 2, 2);    
    hold on;
    
    plot([0, EXT_LEN + 1], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    plot(nanmean(g_extF2Chg_shira.(grp)), '-');
    plot(nanmean(g_extF2Chg_shira.(grp)) - nanste1(g_extF2Chg_shira.(grp)), ...
         '--');
    plot(nanmean(g_extF2Chg_shira.(grp)) + nanste1(g_extF2Chg_shira.(grp)), ...
         '--');
    set(gca, 'XLim', [0, EXT_LEN + 1]);
    ys = get(gca, 'YLim');
    plot(repmat(LEN_RAND_NOPERT + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    xlabel('Trial #');
    ylabel('Sust. F2 change (fraction)');
    
    % Combined F1/F2 amount of compensation   
    gp_extF12Chg_shira.(grp) = nan(size(gp_extF1Chg_shira.(grp)));
    
    for k0 = 1 : size(gp_extF12Chg_shira.(grp), 1)
        for k1 = 1 : size(gp_extF12Chg_shira.(grp), 2)
            gp_extF12Chg_shira.(grp)(k0, k1) = ...
                (-pertVec * [gp_extF1Chg_shira.(grp)(k0, k1), gp_extF2Chg_shira.(grp)(k0, k1)]') / norm(pertVec);
        end
    end
     
    % --- In-phase averages --- %
    figure('Name', sprintf('sust - F1 & F2 chg - %s', grp), ...
           'Position', [150, 150, 600, 600]);
    subplot(4, 1, 1); hold on;
    errorbar([1 : 4], mean(gp_extF1Chg_shira.(grp)), ste(gp_extF1Chg_shira.(grp)), 'bo-');
    ylabel('F1 fraction change');
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
    
    subplot(4, 1, 2); hold on;
    errorbar([1 : 4], mean(gp_extF2Chg_shira.(grp)), ste(gp_extF2Chg_shira.(grp)), 'bo-');
    ylabel('F2 fraction change');
    xlabel('Phase'); 
    set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
    
    subplot(4, 1, 3); hold on;
    errorbar([1 : 4], mean(gp_extF12Chg_shira.(grp)), ste(gp_extF12Chg_shira.(grp)), 'bo-');
    ylabel('Combined F1/F2 compens.');
    xlabel('Phase'); 
    set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
    
    subplot(4, 1, 4); hold on;
    errorbar([1 : 4], mean(g_extFp_shira.(grp)), ste(g_extFp_shira.(grp)), 'bo-');
end

%% Subject-related info: age
for i1 = 1 : numel(grps)
    grp = grps{i1};
    fprintf(1, 'Group %s: age (Y.O.) range = [%.1f, %.1f], mean = %.1f, SD = %.1f, median = %.1f, IQR = %.1f\n', ...
            grp, min(age_yo.(grp)), max(age_yo.(grp)), mean(age_yo.(grp)), std(age_yo.(grp)), median(age_yo.(grp)), iqr(age_yo.(grp)));
end

%% Descriptive statistics: Duration - rand.
figure('Name', 'Vowel duration - rand.');
hold on;
for i1= 1 : numel(grps)
    grp = grps{i1};
    errorbar(1 : 3, ...
             1e3 * mean(meanVDur_rand.(grp)), ...
             1e3 * ste(meanVDur_rand.(grp)), ...
             'Color', colors.(grp));    
end
ylabel('Vowel duration (ms) (mean\pm1 SEM)');
set(gca, 'XLim', [0, 4], 'XTick', [1 : 3], 'XTickLabel', {'noPert', 'higher', 'lower'});
legend(grps);

%% Descriptive statistics: Duration - sust.
figure('Name', 'Vowel duration - sust.');
hold on;
for i1= 1 : numel(grps)
    grp = grps{i1};
    errorbar(1 : 4, ...
             1e3 * mean(meanVDur_sust.(grp)), ...
             1e3 * ste(meanVDur_sust.(grp)), ...
             'Color', colors.(grp));    
end
ylabel('Vowel duration (ms) (mean\pm1 SEM)');
set(gca, 'XLim', [0, 5], 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
legend(grps);

%% Group comparison: rand. projected F1/F2 compens.
figure('Name', 'Groups: rand.: projected F1/F2 (Fp) compens.');
hold on;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    len = size(g_rndFpTrajChgContra.(grp), 2);
    tAxis = (0 : FRAME_DUR : FRAME_DUR * (len - 1)) * 1e3;
   
    plot(tAxis, nanmean(g_rndFpTrajChgContra.(grp)), '-', 'Color', colors.(grp), 'LineWidth', 1.5);
    plot(tAxis, nanmean(g_rndFpTrajChgContra.(grp)) - nanste1(g_rndFpTrajChgContra.(grp)), ... 
          '--', 'Color', colors.(grp), 'LineWidth', 0.5);
	plot(tAxis, nanmean(g_rndFpTrajChgContra.(grp)) + nanste1(g_rndFpTrajChgContra.(grp)), ... 
          '--', 'Color', colors.(grp), 'LineWidth', 0.5);
end
xs = get(gca, 'XLim');
plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
xlabel('Time from vowel onset (ms)');
ylabel('Fp compensation');

ezlegend([NaN, 0.6, 0.05, 0.3, 0.2], 0.45, ...
         grps, {colors.(grps{1}), colors.(grps{2}), colors.(grps{3})}, ...
         [12, 12, 12], {'-', '-', '-'}, {colors.(grps{1}), colors.(grps{2}), colors.(grps{3})}, ...
         [1.5, 1.5, 1.5], [1, 1, 1]);
     
%% Rand. projected F1 / F2 at the last analysis time point
figure('Name', sprintf('Groups: rand.: Normalized composite F12 compens. (%.1f ms)', tAxis(end)));
hold on;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    bar(i1, nanmean(g_rndF12_comp.(grp)), barW, 'EdgeColor', colors.(grp), 'FaceColor', 'none');
    plot(repmat(i1, 1, 2), nanmean(g_rndF12_comp.(grp)) + [-1, 1] * nanste(g_rndF12_comp.(grp)), ...
         'Color', 'k');
end
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', grps);
ylabel('Composite F12 compensation to rand. pert. (normalized)');

%% Group comparison: sust. combined F1/F2 compens.
% Perform statistical comparisons

% check_dir('perm_files', '-create');
% permMatFN = sprintf('%s_%s_perm%d_gend%s.mat', mfilename, testType, nPerm, upper(gend(1)));
% permMatFN = fullfile('perm_files', permMatFN);
% 
% if isfile(permMatFN)
%     load(permMatFN); 
% else
%     ps_wg = nan(numel(grps), size(gp_extF12Chg_shira.ANS, 2) - 1);
%     rp_ps_wg = nan(nPerm, numel(grps), size(gp_extF12Chg_shira.ANS, 2) - 1);
%     for i0 = 1 : 1 + nPerm    
%         for i1 = 1 : numel(grps)
%             grp = grps{i1};
%             if i0 > 1
%                 signPerm = (rand(size(gp_extF12Chg_shira.(grp), 1), 1) > 0.5) * 2 - 1;
%                 dat = gp_extF12Chg_shira.(grp)(:, 2 : end) .* repmat(signPerm, 1, size(gp_extF12Chg_shira.(grp), 2) - 1);
%             else
%                 dat = gp_extF12Chg_shira.(grp)(:, 2 : end);
%             end
% 
%             for i2 = 1 : size(gp_extF12Chg_shira.ANS, 2) - 1
%                 if isequal(testType, 't')
%                     [~, t_p] = ttest(dat(:, i2));
%                 elseif isequal(testType, 'rs')
%                     t_p = signrank(dat(:, i2));
%                 end
% 
%                 if i0 == 1
%                     ps_wg(i1, i2) = t_p;
%                 else
%                     rp_ps_wg(i0 - 1, i1, i2) = t_p;
%                 end
%             end
%         end
%     end
% 
%     save(permMatFN, 'rp_ps_wg', 'ps_wg');
%     fprintf(1, 'INFO: Saved permutation data to file: %s\n', permMatFN);
% end
% 
% assert(size(rp_ps_wg, 1) == nPerm);
% min_rp_ps_wg = min(rp_ps_wg, [], 3);   
%     
% corr_ps_wg = nan(size(ps_wg));
% for i1 = 1 : numel(grps)
%     for i2 = 1 : size(gp_extF12Chg_shira.ANS, 2) - 1
%         corr_ps_wg(i1, i2) = length(find(rp_ps_wg(:, i2) < ps_wg(i1, i2))) / nPerm;
%     end
% end

%% Perform random permutation 
if nPerm > 0
    check_dir('perm_files', '-create');
    permMatFN = fullfile('perm_files', ...
                         sprintf('%s_%s_perm%d_gend%s.mat', mfilename, testType, nPerm, upper(gend(1))));
end

if isfile(permMatFN)
    load(permMatFN);    % gives gp_res
    
    assert(exist('gp_res', 'var') == 1);
    fprintf(1, 'INFO: Loaded permutation data and results from file: %s\n', permMatFN);
else
    gp_res = gp_stats(gp_extF12Chg_shira, nPerm, '--test', testType);
    save(permMatFN, 'gp_res');
    check_file(permMatFN);
    fprintf(1, 'INFO: Saved permutation data and results to file: %s\n', permMatFN);
end


%% Visualization of gp_extF12Chg_shira 
figure('Name', 'Groups: sust.: combined F1/F2 compens.');
hold on;

horizPad = 0.05;
NP = 4;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    mns = mean(gp_extF12Chg_shira.(grp));
    stes = ste(gp_extF12Chg_shira.(grp));
    errorbar([1 : NP], mns, stes, ...
             'o-', 'Color', colors.(grp));
         
	for i2 = 2 : NP
        if gp_res.corp_wg(i1, i2) < CORP_THRESH
            fontWeight = 'bold';
        else
            fontWeight = 'normal';
        end
        text(i2 + horizPad, mns(i2), sprintf('%.3f', gp_res.corp_wg(i1, i2)), ...
             'FontSize', fontSize * 0.7, 'Color', colors.(grp), 'FontWeight', fontWeight);
    end       
end

gcPad = 0.2;
for i1 = 1 : numel(gp_res.gc)
    g1 = gp_res.gc{i1}{1};
    g2 = gp_res.gc{i1}{2};
    
    ig1 = fsic(grps, g1);
    ig2 = fsic(grps, g2);    
    assert(length(ig1) == 1 & length(ig2) == 1);
    
    mn_g1 = mean(gp_extF12Chg_shira.(g1));
    mn_g2 = mean(gp_extF12Chg_shira.(g2));
    clr = (colors.(g1) + colors.(g2)) / 2;
    for i2 = 2 : NP
        plot(repmat(i2 - gcPad * i1, 1, 2), [mn_g1(i2), mn_g2(i2)], ...
             '-', 'Color', clr);
        plot(i2 - gcPad * [i1, i1 - 0.5], repmat(mn_g1(i2), 1, 2), ...
             '-', 'Color', clr);
        plot(i2 - gcPad * [i1, i1 - 0.5], repmat(mn_g2(i2), 1, 2), ...
             '-', 'Color', clr);
        
        if gp_res.corp_bg(i1, i2) < CORP_THRESH
            fontWeight = 'bold';
        else
            fontWeight = 'normal';
        end
        
        ht = text(i2 - gcPad * (i1 + 0.3), mean([mn_g1(i2), mn_g2(i2)]) - 0.2 * (mn_g2(i2) - mn_g1(i2)), ...
                  sprintf('%.3f', gp_res.corp_bg(i1, i2)), ...
                  'Color', clr, 'FontSize', fontSize * 0.7, 'FontWeight', fontWeight);
        set(ht, 'rotation', 90);
    end
end
plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
set(gca, 'XLim', [0.5, 4.5]);
xlabel('Phase'); 
set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
ylabel('Composite F1/F2 adaptation');
legend(grps, 'Location', 'Southeast');


%% Group comparison: combined F1/F2 compens.
figure('Name', 'Groups: sust.: combined F1/F2 compens.');
hold on;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    errorbar([1 : 4], mean(g_extFp_shira.(grp)), ste(g_extFp_shira.(grp)), ...
             'o-', 'Color', colors.(grp))
end
plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
set(gca, 'XLim', [0.5, 4.5]);
xlabel('Phase'); 
set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
ylabel('Fp adaptation');
legend(grps, 'Location', 'Southeast');

%% Correlation between SSI total score and compensation 
meta_corr(SSI_tot, gp_extF12Chg_shira, {'CWS'}, 3, ...
          'SSI total score', 'F12 change (Stay phase)', ...
          colors);

%% Correlation betwen age and compensation 
meta_corr(age_mo, gp_extF12Chg_shira, {'CNS', 'CWS'}, 3, ...
          'Age (m.o.)', 'F12 change (Stay phase)', ...
          colors);

figure;
hold on;
plot(age_mo.CNS, gp_extF12Chg_shira.CNS(:, 3), 'o', 'Color', colors.CNS);
plot(age_mo.CWS, gp_extF12Chg_shira.CWS(:, 3), 'o', 'Color', colors.CWS);
grid on; box on;
legend({'CNS', 'CWS'});

xlabel('Age (mo)');
ylabel('F12 change (Stay phase)');


%% Write data to an xls file, which can be used by other stats software, e.g., SPSS, SYSTAT
if ~bNoXLS
    xlsFNs.sustF1 = sprintf('../data_sets/sust_%s_F1_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(g_extF1Chg_shira.ANS, 1), ...
                            size(g_extF1Chg_shira.CNS, 1), ...
                            size(g_extF1Chg_shira.CWS, 1));
    xlsFNs.sustF2 = sprintf('../data_sets/sust_%s_F2_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(g_extF1Chg_shira.ANS, 1), ...
                            size(g_extF1Chg_shira.CNS, 1), ...
                            size(g_extF1Chg_shira.CWS, 1));
    xlsFNs.sustF12 = sprintf('../data_sets/sust_%s_F12_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(gp_extF12Chg_shira.ANS, 1), ...
                            size(gp_extF12Chg_shira.CNS, 1), ...
                            size(gp_extF12Chg_shira.CWS, 1));
    xlsFNs.sustFp = sprintf('../data_sets/sust_%s_Fp_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(g_extFp_shira.ANS, 1), ...
                            size(g_extFp_shira.CNS, 1), ...
                            size(g_extFp_shira.CWS, 1));
    xlsFNs.rndF12_comp = sprintf('../data_sets/randF12_comp_%s_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(g_rndF12_comp.ANS, 1), ...
                            size(g_rndF12_comp.CNS, 1), ...
                            size(g_rndF12_comp.CWS, 1));
    xlsFNs.meanVDur_rand = sprintf('../data_sets/meanVDur_rand_%s_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(meanVDur_rand.ANS, 1), ...
                            size(meanVDur_rand.CNS, 1), ...
                            size(meanVDur_rand.CWS, 1));
    xlsFNs.meanVDur_sust = sprintf('../data_sets/meanVDur_sust_%s_%dANS_%dCNS_%dCWS.xls', ...
                            gend(1), ...
                            size(meanVDur_sust.ANS, 1), ...
                            size(meanVDur_sust.CNS, 1), ...
                            size(meanVDur_sust.CWS, 1));
    if ~isdir('../data_sets')
        mkdir('../data_sets');
    end


    flds = fields(xlsFNs);

    for i1 = 1 : numel(flds)
        fld = flds{i1};

        if isequal(fld, 'sustF1')
            meas = g_extF1_shira;
        elseif isequal(fld, 'sustF2')
            meas = g_extF2_shira;
        elseif isequal(fld, 'sustF12')
            meas = gp_extF12Chg_shira;
        elseif isequal(fld, 'sustFp')
            meas = g_extFp_shira;
        elseif isequal(fld, 'rndF12_comp')
            meas = g_rndF12_comp;
        elseif isequal(fld, 'meanVDur_rand')
            meas = meanVDur_rand;
        elseif isequal(fld, 'meanVDur_sust')
            meas = meanVDur_sust;
        else
            error('Unexpected field name: %s', fld);
        end

        if ~isempty(strfind(fld, 'sust'));
            title_row = {'NUMBER', 'GRP$', 'SID$', 'START', 'RAMP', 'STAY', 'END'};
        elseif isequal(fld, 'meanVDur_rand')
            title_row = {'NUMBER', 'GRP$', 'SID$', 'NOPERT', 'HIGHER', 'LOWER'};
        elseif isequal(fld, 'rndF12_comp')
            title_row = {'NUMBER', 'GRP$', 'SID$', 'RNDF12COMP'};
        end
        nums = num2cell([meas.ANS; meas.CNS; meas.CWS]);
        sNums = num2cell([1 : size(nums, 1)]');
        t_SIDs = [SID.ANS'; SID.CNS'; SID.CWS'];
        grps = [repmat({'ANS'},  size(g_extF1Chg_shira.ANS, 1), 1); 
                repmat({'CNS'},  size(g_extF1Chg_shira.CNS, 1), 1); 
                repmat({'CWS'},  size(g_extF1Chg_shira.CWS, 1), 1)];

        dat = [title_row; [sNums, grps, t_SIDs, nums]];
        xlswrite(xlsFNs.(fld), dat, 1);

        if isfile(xlsFNs.(fld))
            fprintf(1, 'INFO: Wrote %s data to file %s\n', fld, xlsFNs.(fld));
        end
    end

end
return
