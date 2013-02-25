function list_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls, varargin)
set(uihdls.hlist, 'Enable', 'off');
set(uihdls.bt_reproc, 'Enable', 'off');
set(uihdls.bt_auto_rmsThresh, 'Enable', 'off');

i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

%%
if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
    dataFld = 'otherData';
elseif state.trialList.isRand(i1) == 1
    dataFld = 'randData';
elseif state.trialList.isSust(i1) == 1
    dataFld = 'sustData';
end
% fprintf('dataFld = %s\n', dataFld);

idx_trial = state.trialList.allOrderN(i1);
if isnan(pdata.(dataFld).rmsThresh(idx_trial)) % New
    load(getRawFN_(state.rawDataDir,state.trialList.fn{i1}));	% gives data
    updateParamsUI(uihdls, data);
else
    updateParamsUI(uihdls, pdata, dataFld, idx_trial);
end

if ~isempty(fsic(varargin, 'focus'))
    pdata1 = updateDataUI(uihdls, pdata, dataFld, idx_trial, state, i1, ...
                          'fromList', 'focus');
else
    pdata1 = updateDataUI(uihdls, pdata, dataFld, idx_trial, state, i1, ...
                          'fromList');
end
pdata = pdata1;
state.stats(i1) = 1;

save(stateFN, 'state');
save(dacacheFN, 'pdata');
fprintf('Saved to %s\n', dacacheFN);

updateTrialList(state, uihdls);
set(uihdls.hlist, 'Enable', 'on');
set(uihdls.bt_reproc, 'Enable', 'on');
set(uihdls.bt_auto_rmsThresh, 'Enable', 'on');
return
