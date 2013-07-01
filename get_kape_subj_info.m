function meas = get_kape_subj_info(N, T, sID, measName)
%%
assert(length(strfind(sID, '_')) == 1);

strItems = splitstring(sID, '_');
grpID = strItems{1};
kapeID = strItems{2};

%%
cGrpID = fsic(T(1, :), 'Group');
cKapeID = fsic(T(1, :), 'kAPE ID');
cMeas = fsic(T(1, :), measName);
assert(length(cMeas) == 1);


bFoundSubj = 0;
for i1 = 1 : size(T, 1)
    if isequal(T{i1, cGrpID}, grpID) && isequal(T{i1, cKapeID}, kapeID)
        bFoundSubj = 1;
        break;
    end
end

if bFoundSubj == 0
    meas = NaN;
    fprintf(2, 'WARNING: Failed to find entry [] for subject %s\n', ...
            measName, sID);
else
    meas = N(i1 - 1, cMeas - 4);
end

return