function fig = plot_pc_spike_structure(t, pclus, pidx, nclus, nidx, coeff, score, options)
%PLOT_PC_SPIKE_STRUCTURE  Plots the spike structure for PC eigenvectors.
%
% Syntax:
%   fig = plot_pc_spike_structure(t, pclus, pidx, nclus, nidx, coeff, score, 'Name', value, ...);
%
% Inputs:
%   t - Time vector. Should have same number of elements as `score`
%   pclus - Cluster indices for each positive score peak
%   pidx  - Time index for each positive score peak
%   nclus - Cluster indices for each negative score peak
%   nidx  - Time index for each negative score peak
%   coeff - Coefficients for this eigenvector (or whatever basis)
%   score - The projection of data onto the basis defined by coeff.
%   options:
%     'FigureName' {mustBeTextScalar} = 'NN DRG Depth-Wise PCA';
%     'Position' (1,4) double = [680   193   560   786];
%     'CoefficientsName' {mustBeTextScalar} = 'Eigenvector Coefficients'
%     'CLim' (1,2) double = [-0.75 0.75];
%     'PClusSyntax' {mustBeTextScalar} = 'PClus %d';
%     'NClusSyntax' {mustBeTextScalar} = 'NClus %d';
%
% Output:
%   fig - Figure handle for the structure plot
%
% See also: Contents, example_generate_session_spike_templates

arguments
    t     (1, :) double
    pclus (:, 1) double
    pidx  (:, 1) double
    nclus (:, 1) double
    nidx  (:, 1) double
    coeff (32,1) double
    score (:, 1) double
    options.FigureName {mustBeTextScalar} = 'NN DRG Depth-Wise PCA';
    options.Position (1,4) double = [680   193   560   786];
    options.CoefficientsName {mustBeTextScalar} = 'Eigenvector Coefficients'
    options.CLim (1,2) double = [-0.75 0.75];
    options.PClusSyntax {mustBeTextScalar} = 'PClus %d';
    options.NClusSyntax {mustBeTextScalar} = 'NClus %d';
end


fig = figure(...
    'Color','w', ...
    'Name',options.FigureName, ...
    'Position',options.Position);
% Determine column widths.
if isempty(pclus)
    if isempty(nclus)
        w_depth = 3;
        w_pclus = 0;
        w_nclus = 0;
        n_pclus = 0;
        n_nclus = 0;
    else
        n_pclus = 0;
        n_nclus = max(nclus);
        w_depth = 1;
        w_pclus = 0;
        w_nclus = 2;
    end
elseif isempty(nclus)
    w_depth = 1;
    w_pclus = 2;
    w_nclus = 0;
    n_pclus = max(pclus);
    n_nclus = 0;
else
    n_pclus = max(pclus);
    n_nclus = max(nclus);
    w_depth = 1;
    w_pclus = 1;
    w_nclus = 1;
end
n_rows = max(max(n_pclus, n_nclus)+1,2);
L = tiledlayout(fig, n_rows, 3);

% Plot the coefficients
ax = nexttile(L, 1, [n_rows-1 w_depth]);
set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',turbo(32), ...
    'YDir','reverse', ...
    'YTick',[250, 500, 750], ...
    'YLim',[180, 840], ...
    'XLim',[0 1], ...
    'CLim', options.CLim);
imagesc(ax, [0 1], [180 840], coeff);
ylabel(ax,'Depth (\mum)');
colorbar(ax, 'Location', 'southoutside');
title(ax, options.CoefficientsName, 'FontName','Tahoma','Color','k');

% Plot the scores (projections) time-series from the 10-second epoch.
ax_ts = nexttile(L, n_rows*3 - 2, [1 3]);
set(ax_ts,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
plot(ax_ts, t, score, 'LineWidth', 1.25, 'Color','k');
ylabel(ax_ts, 'Score', 'FontName','Tahoma','Color','k');
xlabel(ax_ts, 'Time (s)', 'FontName','Tahoma','Color','k');

% Plot the positive-clusters
if w_pclus > 0

    [snips,t_snip,i_snip] = signal_2_snippets(score,pidx);
    cdata = cm.umap([0.2 0.2 0.8]);
    cmobj = cm.cmap([1,numel(pidx)], cdata(73:200, :));
    i_col_start = 0;
    for ik = 1:n_pclus

        idx = (pclus(i_snip)==ik);
        if sum(idx) == 0
            continue;
        end
        ax = nexttile(L, 2 + (ik-1)*3, [1 w_pclus]);
        set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',double(cmobj(i_col_start + (1:sum(idx))))/255.0);
        i_col_start = i_col_start + sum(idx);
        plot(ax, t_snip, snips(idx,:), 'LineWidth', 0.75);
        ylabel(ax, 'Score', 'FontName','Tahoma','Color','k');
        xlabel(ax, 'Time (ms)', 'FontName','Tahoma','Color','k');
        title(ax, sprintf(options.PClusSyntax, ik), 'FontName','Tahoma','Color','k');
    end
end

% Plot the negative-clusters
if w_nclus > 0
    [snips,t_snip,i_snip] = signal_2_snippets(score,nidx);
    cdata = cm.umap([0.8 0.2 0.2]);
    cmobj = cm.cmap([1,numel(nidx)], cdata(73:200, :));
    i_col_start = 0;
    for ik = 1:n_nclus
        idx = (nclus(i_snip)==ik);
        if sum(idx) == 0
            continue;
        end
        ax = nexttile(L, 4 - w_nclus + (ik-1)*3, [1 w_nclus]);
        set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',double(cmobj(i_col_start + (1:sum(idx))))/255.0);
        i_col_start = i_col_start + sum(idx);
        plot(ax, t_snip, snips(idx,:), 'LineWidth', 0.75);
        ylabel(ax, 'Score', 'FontName','Tahoma','Color','k');
        xlabel(ax, 'Time (ms)', 'FontName','Tahoma','Color','k');
        title(ax, sprintf(options.NClusSyntax, ik), 'FontName','Tahoma','Color','k');
    end
end

end