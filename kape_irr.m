function kape_irr(subjID, varargin)
[ret,hostName]=system('hostname');
hostName=deblank(hostName);
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    dacacheDir='E:/speechres/kape/dacache';
    rawDataDir='E:/DATA/KAPE/';
else
    dacacheDir='D:/speechres/kape/dacache';
    rawDataDir='D:/DATA/KAPE/';
end


dacacheDir_msu = [dacacheDir, '_msu'];
dacacheDir_bu = [dacacheDir, '_boston'];

pdataFN_msu = fullfile(dacacheDir_msu, [subjID, '.mat']);
pdataFN_bu = fullfile(dacacheDir_bu, [subjID, '.mat']);

%% Load pdata 
if ~isfile(pdataFN_msu)
    error('Cannot find MSU pdata file: %s', pdataFN_msu);
end
if ~isfile(pdataFN_bu)
    error('Cannot find BU pdata file: %s', pdataFN, bu);
end

load(pdataFN_msu);
pdatam = pdata;
clear pdata;

load(pdataFN_bu);
pdatab = pdata;
clear pdata;

%% QA: Compare phase, block and trial number entries
ptflds = {'otherData', 'randData', 'sustData'};
numflds = {'phases', 'blockNums', 'trialNums'};
for i1 = 1 : numel(ptflds)
    ptfld = ptflds{i1};
    for i2 = 1 : numel(numflds)
        numfld = numflds{i2};
        if isequal(pdatab.(ptfld).(numfld), pdatam.(ptfld).(numfld));
            fprintf('%s.%s match. Good.\n', ptfld, numfld);
        else
            error('Mismatch in %s.%s', ptfld, numfld);
        end
    end
    
end

%% 
for i1 = 1 : 2
    figure; hold on;
    if i1 == 1
        x_b = pdatab.sustData.prodF1_shira;
        x_m = pdatam.sustData.prodF1_shira;        
        measName = 'F1';
    else
        x_b = pdatab.sustData.prodF2_shira;
        x_m = pdatam.sustData.prodF2_shira;
        measName = 'F2';
    end
    nLPC_b = pdatab.sustData.nLPC;
    nLPC_m = pdatam.sustData.nLPC;
    
    plot(x_b, x_m, 'o');
    xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
    xl = [min([xs(1), ys(1)]), max([xs(2), ys(2)])];
    yl = xl;
    plot(xl, yl, '-', 'Color', [0.5, 0.5, 0.5]);
    set(gca, 'XLim', xl, 'YLim', yl);
    axis square;
    
    if i1 == 1
        xlabel('F1 (BU)'); ylabel('F1 (MSU)');
    else
        xlabel('F2 (BU)'); ylabel('F2 (MSU)');
    end
    
    % ginput: look up trial number
    if ~isempty(fsic(varargin, 'ginput'))
        bContinue = 0;
        while ~bContinue
            coord = ginput(1);
            
            if coord(1) < xl(1) || coord(1) > xl(2) || coord(2) < yl(1) || coord(2) > yl(2)
                bContinue = 1;
            else
                dist = sqrt((coord(1) - x_b) .^ 2 + (coord(2) - x_m) .^ 2);
                [foo, idxMinDist] = min(dist);

                t_phase = pdatab.sustData.phases{idxMinDist};
                t_blockNum = pdatab.sustData.blockNums(idxMinDist);
                t_trialNum = pdatab.sustData.trialNums(idxMinDist);

                fprintf('Phase %s - block %d - trial %d:\n\t%s (BU) = %.2f Hz (nLPC = %d);\n\t%s (MSU) = %.2f Hz (nLPC = %d)\n', ...
                        t_phase, t_blockNum, t_trialNum, ...
                        measName, x_b(idxMinDist), nLPC_b(idxMinDist), ...
                        measName, x_m(idxMinDist), nLPC_m(idxMinDist));
            end
        end
    end
end

return