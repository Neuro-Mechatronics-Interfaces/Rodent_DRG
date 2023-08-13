function fig = plot_correlogram(tau, rho, options)
%PLOT_CORRELOGRAM  Plots correlation function and optional templates for correlated events.
%
% Syntax:
%   fig = plot_correlogram(tau, rho, 'Name', value, ...);
%
% Inputs:
%   tau - Correlation function time lages (milliseconds)
%   rho - Correlation function value at each time lag in tau
%   options:
%     'Name' {mustBeTextScalar} = 'Correlogram';
%     'Position' (1,4) double = [371   558   870   420];
%     'TemplateI' (1,1) struct = struct('mean', [], 'var', []);
%     'TemplateIYLabel' {mustBeTextScalar} = 'Score';
%     'TemplateIName' {mustBeTextScalar} = 'Template_i';
%     'TemplateJ' (1,1) struct = struct('mean', [], 'var', []);
%     'TemplateJYLabel' {mustBeTextScalar} = 'Score';
%     'TemplateJName' {mustBeTextScalar} = 'Template_j';
%     'TemplateTime' (1,:) double = linspace(-0.8,0.4,37); % (milliseconds)
%     'YLabel' {mustBeTextScalar} = 'Counts / Spike';
%
% Output:
%   fig - Figure handle to correlogram figure.
%
% See also: Contents, example_generate_session_spike_templates

arguments
    tau (1,:) double % Correlation function time lags (milliseconds)
    rho (1,:) double % Correlation function value at each time lag
    options.Name {mustBeTextScalar} = 'Correlogram';
    options.Position (1,4) double = [371   558   870   420];
    options.TemplateI (1,1) struct = struct('mean', [], 'var', []);
    options.TemplateILabel {mustBeTextScalar} = 'Score';
    options.TemplateIName {mustBeTextScalar} = 'Template_i';
    options.TemplateJ (1,1) struct = struct('mean', [], 'var', []);
    options.TemplateJLabel {mustBeTextScalar} = 'Score';
    options.TemplateJName {mustBeTextScalar} = 'Template_j';
    options.TemplateTime (1,:) double = linspace(-0.8,0.4,37); % (milliseconds)
    options.YLabel {mustBeTextScalar} = 'Counts / Spike';
end


fig = figure('Color','w', ...
    'Name',options.Name, ...
    'Position',options.Position);
if ~isempty(options.TemplateJ.mean)
    nRow = 2;
    nCol = 2;
elseif ~isempty(options.TemplateI.mean)
    nRow = 1;
    nCol = 2;
else
    nRow = 1;
    nCol = 1;
end
L = tiledlayout(fig, nRow, nCol);

title(L, options.Name, ...
    'FontName','Tahoma','Color','k');
ax = nexttile(L,1,[nRow 1]);
set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
bar(ax, tau, rho, 1, 'EdgeColor','none','FaceColor','b','DisplayName','Observation');
yline(ax, mean(rho), 'Label', 'Mean', 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
if (mean(rho)-std(rho)) > 0
    yline(ax, [mean(rho)-std(rho), mean(rho)+std(rho)], 'Label', "1 SD", 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1.25);
else
    yline(ax, mean(rho)+std(rho), 'Label', "+1 SD", 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1.25);
end
xlabel(ax,'Time Lag (ms)', 'FontName','Tahoma','Color','k');
ylabel(ax, options.YLabel, 'FontName','Tahoma','Color','k');
title(ax, 'Correlogram', 'FontName','Tahoma','Color','k');

if ~isempty(options.TemplateI.mean)
    ax = nexttile(L,2,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    errorbar(ax,  options.TemplateTime, options.TemplateI.mean, options.TemplateI.var);
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax, options.TemplateILabel, 'FontName','Tahoma','Color','k');
    title(ax, options.TemplateIName,'FontName','Tahoma','Color','k');
end
if ~isempty(options.TemplateJ.mean)
    ax = nexttile(L,4,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    errorbar(ax,  options.TemplateTime, options.TemplateJ.mean, options.TemplateJ.var);
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax, options.TemplateJLabel, 'FontName','Tahoma','Color','k');
    title(ax, options.TemplateJName,'FontName','Tahoma','Color','k');
end
end