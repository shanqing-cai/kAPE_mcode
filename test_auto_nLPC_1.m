function test_auto_nLPC_1(varargin)
tfn = 'E:\DATA\KAPE\ANS_M03\rand\rep1\trial-3-1.mat';
a_nLPCs = [11, 13];

%%

load(tfn); % gives data
dataOrig = data;

fprintf('Original nLPC = %d\n', dataOrig.params.nLPC);

for i1 = 1 : numel(a_nLPCs)    
    for k = 1 : 2
        t_data = reprocData(dataOrig, 'nLPC', a_nLPCs(i1));
    end
    a_data(i1) = t_data;
end

return