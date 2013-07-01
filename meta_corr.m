function meta_corr(meas1, meas2, grps, col2, xlab, ylab, colors)
figure;
hold on;

a_x = [];
a_y = [];

lc_ks = nan(0, 2);  % Linear correlation
lc_rs = nan(0, 1);
lc_ps = nan(0, 1);

sp_ps = nan(0, 1);  % Spearman's rank-order correlation

for i1 = 1 : numel(grps)
    grp = grps{i1};
    
    x = meas1.(grp);
    y = meas2.(grp)(:, col2);
    x = x(:);
    y = y(:);
    a_x = [a_x; x];
    a_y = [a_y; y];
    
    plot(x, y, 'o', 'Color', colors.(grp));
    
    [k, r2, p] = lincorr(x, y);
    r = sqrt(r2) * sign(k(2));
    lc_ks = [lc_ks; k(:)'];
    lc_rs = [lc_rs; r];
    lc_ps = [lc_ps; p];
    
    [rho, t, p] = spear(x, y);
    sp_ps = [sp_ps; p];
end


xs = get(gca, 'XLim');
ys = get(gca, 'YLim');
for i1 = 1 : numel(grps)
    grp = grps{i1};
    plot(xs, lc_ks(i1, 1) + lc_ks(i1, 2) * xs, '--', 'Color', colors.(grp));    
end

set(gca, 'XLim', xs, 'YLim', ys);


xlabel(xlab);
ylabel(ylab);
grid on; box on;
legend(grps);

txtHoriPad = 0.05;
txtVertPad = 0.05;
for i1 = 1 : numel(grps)
    grp = grps{i1};
    text(xs(1) + txtHoriPad * range(xs), ys(2) - txtVertPad * (i1 * 2 - 1) * range(ys), ...
         sprintf('Lin. corr. (%s): R=%.3f, p=%.3f', grp, lc_rs(i1), lc_ps(i1)), ...
         'Color', colors.(grp))
    text(xs(1) + txtHoriPad * range(xs), ys(2) - txtVertPad * (i1 * 2) * range(ys), ...
         sprintf('Spearman (%s): p=%.3f', grp, sp_ps(i1)), ...
         'Color', colors.(grp))
end

% xlabel('Age (mo)');
% ylabel('F12 change (Stay phase)');

return