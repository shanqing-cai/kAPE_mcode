function g_kapeAnalysis()
%% Options: subjects
SID.ANS = {'PILOT_M01', 'PILOT_F03', ...
           'PILOT_M02', 'PILOT_F04', 'ANS_F05', 'ANS_F06', ...
           'ANS_M03', 'ANS_M04', 'ANS_F07', 'ANS_F08'};  % Noise interference, 'PILOT_F04'; old paradigm: 'PILOT_F01',

SID.CNS = {'CNS_F01', 'CNS_F02', 'CNS_F03', ...
           'CNS_M02', 'CNS_M03', 'CNS_M04'}; 
       
         % Bad formant tracking and old paradigm: 'CNS_M01'
         % preproc not done yet: 'CNS_M06', 'CNS_M07';       

SID.CWS = {'CWS_F01', 'CWS_F02', 'CWS_F03', ...
           'CWS_M01', 'CWS_M02', 'CWS_M03'};
          % Where is CWS_M04?
          
       
SHIFT_RATIO_SUST_F1 = +0.25;
SHIFT_RATIO_SUST_F2 = -0.125;

%% Options: colors
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

grps = fields(SID);
g_rndF1TrajChg = struct;

g_extF1_shira = struct;
g_extF2_shira = struct;
g_extFp_shira = struct;

g_extF1Chg_shira = struct;
g_extF2Chg_shira = struct;

gp_extF1Chg_shira = struct;
gp_extF2Chg_shira = struct;
gp_extF12Chg_shira = struct;

ntp = NaN;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    g_rndF1TrajChg.(grp).higher = [];
    g_rndF1TrajChg.(grp).lower = [];
    g_rndF2TrajChg.(grp).higher = [];
    g_rndF2TrajChg.(grp).lower = [];
    
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
        
        for i3 = 1 : numel(randShiftDirs)
            sDir = randShiftDirs{i3};
            g_rndF1TrajChg.(grp).(sDir) = [g_rndF1TrajChg.(grp).(sDir); ...
                   (out.avg_f1Traj.(sDir) - out.avg_f1Traj.noPert) ./ out.avg_f1Traj.noPert];
            g_rndF2TrajChg.(grp).(sDir) = [g_rndF2TrajChg.(grp).lower; ...
                   (out.avg_f1Traj.(sDir) - out.avg_f1Traj.noPert) ./ out.avg_f1Traj.noPert];
        end
        
        g_extF1_shira.(grp) = [g_extF1_shira.(grp); out.sust_prodF1_shira];
        g_extF2_shira.(grp) = [g_extF2_shira.(grp); out.sust_prodF2_shira];
        
        approx_pert_mag = norm([mean(out.sust_prodF1_shira(2 : 3)) * SHIFT_RATIO_SUST_F1, ...
                                mean(out.sust_prodF2_shira(2 : 3)) * SHIFT_RATIO_SUST_F2]);
        g_extFp_shira.(grp) = [g_extFp_shira.(grp); ...
                               proj2PertLine(SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2, ...
                                             out.sust_prodF1_shira, out.sust_prodF2_shira) / approx_pert_mag];
        
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
    
    figure('Name', sprintf('rand - F1 - %s', grp));
    hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (ntp - 1);
    plot(tAxis, zeros(size(tAxis)), '-', 'color', [0.5, 0.5, 0.5]);
    for i2 = 1 : numel(randShiftDirs)
        sDir = randShiftDirs{i2};
        plot(tAxis, mean(g_rndF1TrajChg.(grp).(sDir)), ...
             'Color', colors.(sDir), 'LineWidth', 1.5);
        plot(tAxis, mean(g_rndF1TrajChg.(grp).(sDir)) + ste(g_rndF1TrajChg.(grp).(sDir)), ...
             '--', 'Color', colors.(sDir), 'LineWidth', 1);
        plot(tAxis, mean(g_rndF1TrajChg.(grp).(sDir)) - ste(g_rndF1TrajChg.(grp).(sDir)), ...
             '--', 'Color', colors.(sDir), 'LineWidth', 1);
    end
    xlabel('Time (s)');
    ylabel('F1 change (ratio)');
    
    figure('Name', sprintf('rand - F2 - %s', grp));
    hold on;
    tAxis = 0 : FRAME_DUR : FRAME_DUR * (ntp - 1);
    plot(tAxis, zeros(size(tAxis)), '-', 'color', [0.5, 0.5, 0.5]);
    for i2 = 1 : numel(randShiftDirs)
        sDir = randShiftDirs{i2};
        plot(tAxis, mean(g_rndF2TrajChg.(grp).(sDir)), ...
             'Color', colors.(sDir), 'LineWidth', 1.5);
        plot(tAxis, mean(g_rndF2TrajChg.(grp).(sDir)) + ste(g_rndF2TrajChg.(grp).(sDir)), ...
             '--', 'Color', colors.(sDir), 'LineWidth', 1);
        plot(tAxis, mean(g_rndF2TrajChg.(grp).(sDir)) - ste(g_rndF2TrajChg.(grp).(sDir)), ...
             '--', 'Color', colors.(sDir), 'LineWidth', 1);
    end
    xlabel('Time (s)');
    ylabel('F2 change (ratio)');
    
    figure('Name', sprintf('sust - F1 - %s', grp));
    hold on;
    plot([0, EXT_LEN + 1], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    plot(nanmean(g_extF1Chg_shira.(grp)), '-');
    plot(nanmean(g_extF1Chg_shira.(grp)) - nanste(g_extF1Chg_shira.(grp)), ...
         '--');
    plot(nanmean(g_extF1Chg_shira.(grp)) + nanste(g_extF1Chg_shira.(grp)), ...
         '--');
    set(gca, 'XLim', [0, EXT_LEN + 1]);
    ys = get(gca, 'YLim');
    plot(repmat(LEN_RAND_NOPERT + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    xlabel('Trial #');
    ylabel('F1 change (fraction)');
%     plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + LEN_END + 1, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);

    figure('Name', sprintf('sust - F2 - %s', grp));
    hold on;
    plot([0, EXT_LEN + 1], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    plot(nanmean(g_extF2Chg_shira.(grp)), '-');
    plot(nanmean(g_extF2Chg_shira.(grp)) - nanste(g_extF2Chg_shira.(grp)), ...
         '--');
    plot(nanmean(g_extF2Chg_shira.(grp)) + nanste(g_extF2Chg_shira.(grp)), ...
         '--');
    set(gca, 'XLim', [0, EXT_LEN + 1]);
    ys = get(gca, 'YLim');
    plot(repmat(LEN_RAND_NOPERT + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    plot(repmat(LEN_RAND_NOPERT + LEN_START + LEN_RAMP + LEN_STAY + 0.5, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    xlabel('Trial #');
    ylabel('F2 change (fraction)');
    
    % Combined F1/F2 amount of compensation   
    gp_extF12Chg_shira.(grp) = nan(size(gp_extF1Chg_shira.(grp)));
    pertVec = [SHIFT_RATIO_SUST_F1, SHIFT_RATIO_SUST_F2];
    for k0 = 1 : size(gp_extF12Chg_shira.(grp), 1)
        for k1 = 1 : size(gp_extF12Chg_shira.(grp), 2)
            gp_extF12Chg_shira.(grp)(k0, k1) = ...
                (-pertVec * [gp_extF1Chg_shira.(grp)(k0, k1), gp_extF2Chg_shira.(grp)(k0, k1)]') / norm(pertVec);
        end
    end
     
    % --- In-phase averages --- %
    figure('Name', sprintf('sust - F1 & F2 chg - %s', grp));
    subplot(3, 1, 1); hold on;
    errorbar([1 : 4], mean(gp_extF1Chg_shira.(grp)), ste(gp_extF1Chg_shira.(grp)), 'bo-');
    ylabel('F1 fraction change');
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
    
    subplot(3, 1, 2); hold on;
    errorbar([1 : 4], mean(gp_extF2Chg_shira.(grp)), ste(gp_extF2Chg_shira.(grp)), 'bo-');
    ylabel('F2 fraction change');
    xlabel('Phase'); 
    set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
    
    subplot(3, 1, 3); hold on;
    errorbar([1 : 4], mean(gp_extF12Chg_shira.(grp)), ste(gp_extF12Chg_shira.(grp)), 'bo-');
    ylabel('Combined F1/F2 compens.');
    xlabel('Phase'); 
    set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
    plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', [0.5, 4.5]);
end

%% Group comparison: combined F1/F2 compens.
figure('Name', 'Groups: combined F1/F2 compens.');
hold on;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    errorbar([1 : 4], mean(gp_extF12Chg_shira.(grp)), ste(gp_extF12Chg_shira.(grp)), ...
             'o-', 'Color', colors.(grp))
end
plot([0.5, 4.5], [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
set(gca, 'XLim', [0.5, 4.5]);
xlabel('Phase'); 
set(gca, 'XTick', [1 : 4], 'XTickLabel', {'Start', 'Ramp', 'Stay', 'End'});
ylabel('Composite F1/F2 adaptation');
legend(grps, 'Location', 'Southeast');

%% Write data to an xls file, which can be used by other stats software, e.g., SPSS, SYSTAT
xlsFNs.sustF1 = sprintf('sust_F1_%dANS_%dCNS_%dCWS.xls', ...
                        size(g_extF1Chg_shira.ANS, 1), ...
                        size(g_extF1Chg_shira.CNS, 1), ...
                        size(g_extF1Chg_shira.CWS, 1));
xlsFNs.sustF2 = sprintf('sust_F2_%dANS_%dCNS_%dCWS.xls', ...
                        size(g_extF1Chg_shira.ANS, 1), ...
                        size(g_extF1Chg_shira.CNS, 1), ...
                        size(g_extF1Chg_shira.CWS, 1));
xlsFNs.sustFp = sprintf('sust_Fp_%dANS_%dCNS_%dCWS.xls', ...
                        size(g_extF1Chg_shira.ANS, 1), ...
                        size(g_extF1Chg_shira.CNS, 1), ...
                        size(g_extF1Chg_shira.CWS, 1));
                    
flds = fields(xlsFNs);

for i1 = 1 : numel(flds)
    fld = flds{i1};
    
    if isequal(fld, 'sustF1')
        meas = g_extF1_shira;
    elseif isequal(fld, 'sustF2')
        meas = g_extF2_shira;
    elseif isequal(fld, 'sustFp')
        meas = g_extFp_shira;
    else
        error('Unexpected field name: %s', fld);
    end
    
    title_row = {'NUMBER', 'GRP$', 'SID$', 'START', 'RAMP', 'STAY', 'END'};
    nums = num2cell([meas.ANS; meas.CNS; meas.CWS]);
    sNums = num2cell([1 : size(nums, 1)]');
    t_SIDs = [SID.ANS'; SID.CNS'; SID.CWS'];
    grps = [repmat({'ANS'},  size(g_extF1Chg_shira.ANS, 1), 1); 
            repmat({'CNS'},  size(g_extF1Chg_shira.CNS, 1), 1); 
            repmat({'CWS'},  size(g_extF1Chg_shira.CWS, 1), 1)];
        
	dat = [title_row; [sNums, grps, t_SIDs, nums]];
    xlswrite(xlsFNs.(fld), dat, 1);
    
    if isfile(xlsFNs.(fld))
        fprintf(1, 'INFO: Wrote %d data to file %s\n', fld, xlsFNs.(fld));
    end
end
return
