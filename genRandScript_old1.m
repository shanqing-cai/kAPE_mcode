function phaseScript = genRandScript_old1(nBlocks, trialsPerBlock, trialsPerBlock_lower, trialsPerBlock_higher, words)
phaseScript = struct();
phaseScript.nReps = nBlocks;
phaseScript.nTrials = 0;

for n=1 : nBlocks
    wordsUsed = {};
    while (length(wordsUsed) < trialsPerBlock)
        nToGo = trialsPerBlock - length(wordsUsed);
        if (nToGo >= length(words))
            wordsUsed = [wordsUsed, words{randperm(length(words))}];
        elseif (nToGo > 0)
            idx = randperm(length(words));            
            wordsUsed = [wordsUsed, words(idx(1 : nToGo))];
        end
    end
    
%     wordsUsed = wordsUsed(randperm(length(wordsUsed)));
    bt = zeros(1, length(wordsUsed));
    
%     pseudoWordsUsed=pseudoWords(randperm(length(pseudoWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));            
    twCnt=1;
    bt=bt(randperm(length(bt)));
    oneRep=struct;
    oneRep.trialOrder=[];
    oneRep.word=cell(1,0);
    cntTW=1;
    for m=1:length(bt)
        if (bt(m)==0)					
            oneRep.trialOrder=[oneRep.trialOrder,1];
            oneRep.word{length(oneRep.word)+1}=wordsUsed{twCnt};
            twCnt=twCnt+1;
        elseif (bt(m)==1)
            oneRep.trialOrder=[oneRep.trialOrder,[5,4,4]];
            oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
            oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW);
            oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
            cntTW=cntTW+2;
        end
    end

    oneRep.pertType = zeros(1, numel(oneRep.trialOrder));
    phaseScript.(['rep',num2str(n)])=oneRep;
    phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);
end

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

pertLastRepEnd = 0;
for i1 = 1 : nBlocks
    t_word = words(pertVec(i1));
    idx_word = fsic(phaseScript.(['rep',num2str(i1)]).word, t_word);
    if length(idx_word) < 2
        error('length(idx_word) should be at least 2.')
    end
    
    iOkay = 0;
    while iOkay == 0
        idx_perts = idx_word(randperm(length(idx_word)));
        if abs(idx_perts(1) - idx_perts(2)) > 2 && ~(pertLastRepEnd + min([idx_perts(1), idx_perts(2)]) <= 2)
            iOkay = 1;
        else
            if abs(idx_perts(2) - idx_perts(3)) > 2 && ~(pertLastRepEnd + min([idx_perts(2), idx_perts(3)]) <= 2)
                iOkay = 2;
            else
                if abs(idx_perts(1) - idx_perts(3)) > 2 && ~(pertLastRepEnd + min([idx_perts(1), idx_perts(3)]) <= 2)
                    iOkay = 3;
                end
            end
        end
    end    
    
    if iOkay == 1
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(1)) = -1;   % Lower
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(2)) = 1;   % Higher
        
        pertLastRepEnd = length(phaseScript.(['rep',num2str(i1)]).pertTyp) - max([idx_perts(1), idx_perts(2)]);
    elseif iOkay == 2
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(2)) = -1;   % Lower
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(3)) = 1;   % Higher
        
        pertLastRepEnd = length(phaseScript.(['rep',num2str(i1)]).pertTyp) - max([idx_perts(2), idx_perts(3)]);
    else
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(1)) = -1;   % Lower
        phaseScript.(['rep',num2str(i1)]).pertType(idx_perts(3)) = 1;   % Higher
        
        pertLastRepEnd = length(phaseScript.(['rep',num2str(i1)]).pertTyp) - max([idx_perts(1), idx_perts(3)]);
    end
    
    if pertLastRepEnd
        pause(0);
    end
    
%     phaseScript.(['rep',num2str(i1)]).
end
return