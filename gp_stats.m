function res = gp_stats(dat, nPerm, varargin)
%% Constants
NP = 4; % Number of phases
ALL_TEST_TYPES = {'t', 'rs'}; % {t-test | signed-rank; rank-sum}

%% Additional input arguments
testType = 'rs';
if ~isempty(fsic(varargin, '--test'))
    testType = varargin{fsic(varargin, '--test') + 1};
end
assert(length(fsic(ALL_TEST_TYPES, testType)) == 1);

%% Check input data
grps = fields(dat);
for i1 = 1 : numel(grps)
    grp = grps{i1};
    assert(size(dat.(grp), 2) == 4);
end

assert(nPerm >= 0);
assert(mod(nPerm, 1.0) == 0.0);

%% Within group difference from zero
uncp_wg = nan(numel(grps), NP); % Uncorrected p-values
corp_wg = nan(numel(grps), NP); % Corrected p-values

rp_p_wg = nan(numel(grps), NP, nPerm);
rp_sig_wg = nan(numel(grps), NP, nPerm);

% --- Perform test and random permutations --- %
if nPerm > 0
    nc = print_progress_bar(0, nPerm, sprintf('Performing permutation test: Within-group'));
end

for i0 = 0 : nPerm
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        if i0 == 0
            sgn = ones(size(dat.(grp), 1), 1);
        else
            sgn = (rand(size(dat.(grp), 1), 1) > 0.5) * 2 - 1;
        end
        t_dat = dat.(grp) .* repmat(sgn, 1, NP);
        
        for i2 = 2 : NP
            if isequal(testType, 'rs')
                [p, ~, stats] = signrank(t_dat(:, i2));
                sig = -log10(p) * sign(median(t_dat(:, i2)));
            elseif isequal(testType, 't')
                [~, p, ~, stats] = ttest(t_dat(:, i2));
                sig = -log10(p) * sign(mean(t_dat(:, i2)));
            end
            
            if i0 == 0
                uncp_wg(i1, i2) = p;
            else
                rp_p_wg(i1, i2, i0) = p;
                rp_sig_wg(i1, i2, i0) = sig;
            end
        end
    end
    
    if nPerm > 0 && mod(i0, round(nPerm / 10)) == 0
        for k1 = 1 : nc; fprintf(1, '\b'); end 
        print_progress_bar(i0, nPerm, sprintf('Performing permutation test: Within-group'));
    end
end

if nPerm > 0
    fprintf(1, '\n');
end

% --- Calculate corrected p-values --- %
for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    % -- Across phases -- %
    t_ps = squeeze(rp_p_wg(i1, :, :));
    min_ps = min(t_ps, [], 1);
    
    % -- Across groups and phases -- %
%     t_ps = rp_p_wg;
%     min_ps = squeeze(min(t_ps, [], 1));
%     min_ps = min(min_ps, [], 1);
    
    for i2 = 2 : NP
        corp_wg(i1, i2) = numel(find(min_ps <= uncp_wg(i1, i2))) / nPerm;
    end
end

%% Between-group comparisons
gc = {};    % Group crossings
for i1 = 1 : numel(grps)
    grp1 = grps{i1};
    for i2 = i1 + 1 : numel(grps)
        grp2 = grps{i2};
        
        gc{end + 1} = {grp1, grp2};
    end
end

ngc = numel(gc);    % Number of group crossings

uncp_bg = nan(ngc, NP); % Uncorrected p-values
corp_bg = nan(ngc, NP); % Corrected p-values

rp_p_bg = nan(ngc, NP, nPerm);
rp_sig_bg = nan(ngc, NP, nPerm);

if nPerm > 0
    nc = print_progress_bar(0, nPerm, sprintf('Performing permutation test: Between-group'));
end

for i0 = 0 : nPerm    
    for i1 = 1 : ngc
        grp1 = gc{i1}{1};
        grp2 = gc{i1}{2};
        ns1 = size(dat.(grp1), 1);
        ns2 = size(dat.(grp2), 1);
        
        a_dat = [dat.(grp1); dat.(grp2)];
        
        if i0 > 0 % Random permutation 
            rpidx = randperm(size(a_dat, 1));
            
            a_dat = a_dat(rpidx, :);
        end
        
        g1_dat = a_dat(1 : ns1, :);
        g2_dat = a_dat(ns1 + 1 : end, :);
        
        for i2 = 2 : NP
            if isequal(testType, 'rs')
                [p, ~, stats] = ranksum(g1_dat(:, i2), g2_dat(:, i2));
                sig = -log10(p) * sign(median(g1_dat(:, i2)) - median(g2_dat(:, i2)));
            elseif isequal(testType, 't')
                [~, p, ~, stats] = ttest2(g1_dat(:, i2), g2_dat(:, i2));
                sig = -log10(p) * sign(mean(g1_dat(:, i2)) - mean(g2_dat(:, i2)));
            end
            
            if i0 == 0
                uncp_bg(i1, i2) = p;
            else
                rp_p_bg(i1, i2, i0) = p;
                rp_sig_bg(i1, i2, i0) = sig;
            end
        end
    end
    
    if nPerm > 0 && mod(i0, round(nPerm / 10)) == 0
        for k1 = 1 : nc; fprintf(1, '\b'); end 
        print_progress_bar(i0, nPerm, sprintf('Performing permutation test: Between-group'));
    end
end

if nPerm > 0
    fprintf(1, '\n');
end

% --- Calculate corrected p-values --- %
for i1 = 1 : ngc
    
    % -- Across phases -- %
    t_ps = squeeze(rp_p_bg(i1, :, :));
    min_ps = min(t_ps, [], 1);
    
    % -- Across groups and phases -- %
%     t_ps = rp_p_bg;
%     min_ps = squeeze(min(t_ps, [], 1));
%     min_ps = min(min_ps, [], 1);
    
    for i2 = 2 : NP
        corp_bg(i1, i2) = numel(find(min_ps <= uncp_bg(i1, i2))) / nPerm;
    end
end

%% Furnish output
res = struct('uncp_wg', uncp_wg, ...
             'corp_wg', corp_wg, ...
             'rp_p_wg', rp_p_wg, ...
             'rp_sig_wg', rp_sig_wg, ...             
             'uncp_bg', uncp_bg, ...
             'corp_bg', corp_bg, ...
             'rp_p_bg', rp_p_bg, ...
             'rp_sig_bg', rp_sig_bg);
res.gc = gc;
return