function phaseScript = genRandScript(nBlocks, trialsPerBlock, trialsPerBlock_lower, trialsPerBlock_higher, words)
%% 
MIN_PERT_DIST = 2;  % 0 corresponds to adjacent, 1 corresponds to one intervening noPert trials, ...

%%
phaseScript = struct();
phaseScript.nReps = nBlocks;
phaseScript.nTrials = 0;

%% What word to perturb in each block
nWords = numel(words);
pertVec = [];
while (length(pertVec) < nBlocks)
    nToGo = nBlocks - length(pertVec);
    if (nToGo >= length(words))
        pertVec = [pertVec, [1 : length(words)]];
    elseif (nToGo > 0)
        idx = randperm(length(words));
        pertVec = [pertVec, idx(1 : nToGo)];
    end
end
pertVec = pertVec(randperm(numel(pertVec)));

%% Order: headPad, pertTrial1, interPad, pertTrial2, tailPad
prevDist = Inf;
for n = 1 : nBlocks
    oneRep = struct;
    oneRep.trialOrder = ones(1, trialsPerBlock);
    oneRep.word = cell(1, trialsPerBlock);
    oneRep.pertType = [];
    
    % Determine the minimum and maximum length of the headPad
    minhp = max([MIN_PERT_DIST - prevDist, 0]);
    maxhp = trialsPerBlock - 2 - MIN_PERT_DIST;
    
    hpLen = floor(rand * (maxhp - minhp + 1)) + minhp;
    headPad = zeros(1, hpLen);
    
    minip = MIN_PERT_DIST;
    maxip = trialsPerBlock - 2 - hpLen;
    
    ipLen = floor(rand * (maxip - minip + 1)) + minip;    
    interPad = zeros(1, ipLen);
    
    tpLen = trialsPerBlock - 2 - hpLen - ipLen;
    tailPad = zeros(1, tpLen);
    
    if rand > 0.5   % [1, -1]
       oneRep.pertType = [headPad, 1, interPad, -1, tailPad];
    else
       oneRep.pertType = [headPad, -1, interPad, 1, tailPad];
    end
    
    prevDist = tpLen;
    
    % Supply the words
    useCnt = [];
    for i1 = 1 : nWords
        if i1 < nWords
            useCnt(end + 1) = ceil(trialsPerBlock / nWords);
        else
            useCnt(end + 1) = trialsPerBlock - sum(useCnt);
        end
    end
    
    idxAll = 1 : numel(oneRep.pertType);
    idxPert = find(oneRep.pertType ~= 0);
    idxNoPert = setxor(idxAll, idxPert);
    for m = 1 : numel(idxPert)
        oneRep.word{idxPert(m)} = words{pertVec(n)};
        useCnt(pertVec(n)) =  useCnt(pertVec(n)) - 1;
    end
    
    t_idx = [];
    for i1 = 1 : nWords
        t_idx = [t_idx, i1 * ones(1, useCnt(i1))];
    end
    t_idx = t_idx(randperm(length(t_idx)));
    
    for i1 = 1 : numel(t_idx)
        oneRep.word{idxNoPert(i1)} = words{t_idx(i1)};
    end
    
    phaseScript.(['rep',num2str(n)]) = oneRep;
    phaseScript.nTrials = phaseScript.nTrials + length(oneRep.trialOrder);
end
return