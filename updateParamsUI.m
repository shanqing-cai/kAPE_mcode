function updateParamsUI(uihdls, varargin)
%% 
DEFAULT_RATING = 2;

%%
red = [1, 0, 0];
green = [0, 0.5, 0];

if nargin == 2 % data
    data = varargin{1};
    
    rmsThresh = data.params.rmsThresh;
    nLPC = round(data.params.nLPC);
    fn1 = data.params.fn1;
    fn2 = data.params.fn2;
    aFact = data.params.aFact;
    bFact = data.params.bFact;
    gFact = data.params.gFact;
    bCepsLift = data.params.bCepsLift;
    cepsWinWidth = data.params.cepsWinWidth;
    
    srt_nLPCs = [];
    str_srt_nLPCs = {};
    idx_nLPC = NaN;
    
    mark_f1_beg = 0;
    mark_f1_mid = 0;
    mark_f1_end = 0;
    mark_f2_beg = 0;
    mark_f2_mid = 0;
    mark_f2_end = 0;
    
    rating = DEFAULT_RATING;
    comments = '';
    
elseif nargin == 4 % pdata
    pdata = varargin{1};
    dataFld = varargin{2};
    idx = varargin{3};
    
    rmsThresh = pdata.(dataFld).rmsThresh(idx);
    nLPC = round(pdata.(dataFld).nLPC(idx));
    
    if isfield(pdata.(dataFld), 'srt_nLPCs') ...
            && ~isempty(pdata.(dataFld).srt_nLPCs{idx})
        srt_nLPCs = pdata.(dataFld).srt_nLPCs{idx};
        str_srt_nLPCs = cell(size(srt_nLPCs));
        
        for k1 = 1 : length(srt_nLPCs)
            str_srt_nLPCs{k1} = sprintf('%d', srt_nLPCs(k1));
        end
        
        idx_nLPC = find(srt_nLPCs == nLPC);        
    else
        srt_nLPCs = [];
        str_srt_nLPCs = {};
        idx_nLPC = NaN;
    end
    
    fn1 = pdata.(dataFld).fn1(idx);
    fn2 = pdata.(dataFld).fn2(idx);
    aFact = pdata.(dataFld).aFact(idx);
    bFact = pdata.(dataFld).bFact(idx);
    gFact = pdata.(dataFld).gFact(idx);
    bCepsLift = pdata.(dataFld).bCepsLift(idx);
    cepsWinWidth = pdata.(dataFld).cepsWinWidth(idx);
    
    mark_f1_beg = pdata.(dataFld).mark_f1_beg(idx);
    mark_f1_mid = pdata.(dataFld).mark_f1_mid(idx);
    mark_f1_end = pdata.(dataFld).mark_f1_end(idx);
    mark_f2_beg = pdata.(dataFld).mark_f2_beg(idx);
    mark_f2_mid = pdata.(dataFld).mark_f2_mid(idx);
    mark_f2_end = pdata.(dataFld).mark_f2_end(idx);
    
    rating = pdata.(dataFld).rating(idx);
    comments = pdata.(dataFld).comments{idx};
else
    return;
end

% set(uihdls.edit_rmsThresh, 'String', sprintf('%.5f', rmsThresh));
set(uihdls.edit_rmsThresh, 'String', sprintf('%.8f', rmsThresh));
set(uihdls.edit_nLPC, 'String', sprintf('%d', nLPC));

if ~isempty(srt_nLPCs)
    set(uihdls.lst_srt_nLPCs, 'String', str_srt_nLPCs, 'Enable', 'on', ...
        'Value', idx_nLPC);
else
    set(uihdls.lst_srt_nLPCs, 'String', str_srt_nLPCs, 'Enable', 'off');
end

set(uihdls.edit_fn1, 'String', sprintf('%.1f', fn1));
set(uihdls.edit_fn2, 'String', sprintf('%.1f', fn2));
set(uihdls.edit_aFact, 'String', sprintf('%.1f', aFact));
set(uihdls.edit_bFact, 'String', sprintf('%.1f', bFact));
set(uihdls.edit_gFact, 'String', sprintf('%.1f', gFact));
set(uihdls.edit_bCepsLift, 'String', sprintf('%d', bCepsLift));
set(uihdls.edit_cepsWinWidth, 'String', sprintf('%d', cepsWinWidth));

if mark_f1_beg == 1
    set(uihdls.bt_f1_beg, 'ForegroundColor', red);
else
    set(uihdls.bt_f1_beg, 'ForegroundColor', green);
end

if mark_f1_mid == 1
    set(uihdls.bt_f1_mid, 'ForegroundColor', red);
else
    set(uihdls.bt_f1_mid, 'ForegroundColor', green);
end

if mark_f1_end == 1
    set(uihdls.bt_f1_end, 'ForegroundColor', red);
else
    set(uihdls.bt_f1_end, 'ForegroundColor', green);
end

if mark_f2_beg == 1
    set(uihdls.bt_f2_beg, 'ForegroundColor', red);
else
    set(uihdls.bt_f2_beg, 'ForegroundColor', green);
end

if mark_f2_mid == 1
    set(uihdls.bt_f2_mid, 'ForegroundColor', red);
else
    set(uihdls.bt_f2_mid, 'ForegroundColor', green);
end

if mark_f2_end == 1
    set(uihdls.bt_f2_end, 'ForegroundColor', red);
else
    set(uihdls.bt_f2_end, 'ForegroundColor', green);
end

items = get(uihdls.pm_rating, 'String');
for i1 = 1 : numel(items)
    if isequal(str2num(items{i1}), rating)
        break;
    end
end
set(uihdls.pm_rating, 'Value', i1);

set(uihdls.edit_comments, 'String', comments);

return