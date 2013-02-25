function reproc_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

if ~iscell(uihdls)
    i1 = get(uihdls.hlist, 'Value');
    
    if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
        dataFld = 'otherData';
    elseif state.trialList.isRand(i1) == 1
        dataFld = 'randData';
    elseif state.trialList.isSust(i1) == 1
        dataFld = 'sustData';
    end
    % fprintf('dataFld = %s\n', dataFld);

    idx_trial = state.trialList.allOrderN(i1);
    
    pdata1 = updateDataUI(uihdls, pdata, dataFld, idx_trial, state, i1);
else
    dataFld = uihdls{2};
    idx_trial = uihdls{3};
    
    pdata1 = updateDataUI(uihdls, pdata, dataFld, idx_trial, state, i1);
end


pdata = pdata1;

if ~iscell(uihdls)
    state.stats(i1) = 1;
end

save(stateFN, 'state');
save(dacacheFN, 'pdata');

if ~iscell(uihdls)
    updateTrialList(state, uihdls);
end
return

