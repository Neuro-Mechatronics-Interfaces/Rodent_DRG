%EXAMPLE_GENERATE_SESSION_SPIKE_TEMPLATES Example of generating spike templates for each channel, for a session (tank), which can then be used for each separate recording to correspond clusters (for different stimulus conditions etc.) 
close all force;
clear;
clc;

SUBJ = 'DW100';
RAW_DATA_FOLDER = parameters('raw_data_folder');
GEN_DATA_FOLDER = parameters('generated_data_folder');
YYYY = 2023;
MM = 8;
DD = 9;
CHANNELS = 1:32; % Amplifier channels to retain

TEMPLATE_BLOCK = 2; % Block index to use for generating templates
EPOCH = [10, 70]; % Seconds to use for detection/template generation

K_LIST = 2:8;                % K-means to try for Calinski-Harabasz-criterion cluster evaluation
SILHOUETTE_THRESHOLD = 0.25; % Minimum silhouette score to include in cluster

P_COLS = [validatecolor("#F0A202"); ...
          validatecolor("#F19504"); ...
          validatecolor("#F18805")];

N_COLS = [validatecolor("#0E1428"); ...
          validatecolor("#2A3741"); ...
          validatecolor("#455959")];

%% Load data
% Load data (see comment below first!)
% Need to have raptor datashare (\\raptor.ni.cmu.edu\NML) mapped -- see the
% lab wiki for instructions on how to do that.
x = io.load_data(SUBJ, YYYY, MM, DD, "*", TEMPLATE_BLOCK, '.rhd', RAW_DATA_FOLDER);
tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD);


%% Extract relevant endpoints and clear data struct from memory
fs = x.frequency_parameters.amplifier_sample_rate;
t = x.t_amplifier((x.t_amplifier >= EPOCH(1)) & (x.t_amplifier < EPOCH(2)));
data = x.amplifier_data(CHANNELS, (x.t_amplifier >= EPOCH(1)) & (x.t_amplifier < EPOCH(2)));

spike_band = struct;
[spike_band.b,spike_band.a] = butter(4, [300 3000]./(fs/2), 'bandpass');
data_f = filtfilt(spike_band.b, spike_band.a, data')';

%% Get the channel indexing and a masking vector for plotting data snapshot
[i_intan, i_NN, depth_NN] = intan_2_neuronexus();
i_plot = (t >= 55) & (t < 65);
ch_exclude = [10, 26]; % 1-indexed versions of the Intan channels (e.g. A-000 becomes "1")
i_exclude = nan(size(ch_exclude));
for ii = 1:numel(i_exclude)
    i_exclude(ii) = find(i_intan == ch_exclude(ii), 1, 'first');
end

%% Make a figure with the data streams shown as they would be on probe
labs = strings(32,1);
for iL = 1:numel(labs)
    labs(iL) = string(sprintf('A-%03d: %d \\mum', i_intan(iL)-1, depth_NN(iL)));
end
fig = figure('Color','w','Name','NN DRG Depth-wise Layout','Position',[362   356   878   623]);
ax = axes(fig,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',turbo(32), ...
    'YDir','reverse','YTick',depth_NN,'YLim',[depth_NN(1)-10, depth_NN(end)+10],'YTickLabel',labs);
xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');
data_depth = -data_f(i_intan,i_plot); 
h = plot(ax, t(i_plot), data_depth./6 + depth_NN); % Depths in microns
for ii = 1:numel(i_exclude)
    h(i_exclude(ii)).Color = [0.85 0.85 0.85];
end
out_folder = fullfile(GEN_DATA_FOLDER, tank);
if exist(out_folder,'dir')==0
    mkdir(out_folder);
end
default.savefig(fig, fullfile(out_folder, sprintf('%s_%d_Depths_Snippet', tank, TEMPLATE_BLOCK)));

%% Estimate principal components of the signals shown in this plot.
data_blanked = data_depth';
data_blanked(:,i_exclude) = 0; % Blank the noisy channels
spike_d = struct;

warning('off','stats:pca:ColRankDefX'); % We expect this warning due to the blanking step.
[spike_d.coeff, spike_d.score, spike_d.latent, spike_d.tsquared, spike_d.explained, spike_d.mu] = pca(data_blanked,'Algorithm','eig');
warning('on','stats:pca:ColRankDefX');

%% Apply thresholding on the PCA scores. Use peak-finding to get timestamps.
spike = struct;
spike.idx = cell(1,12);
spike.pidx = cell(1,12);
spike.nidx = cell(1,12);
spike.clus = cell(1,12);
spike.pclus = cell(1,12);
spike.nclus = cell(1,12);
warning('off','signal:findpeaks:largeMinPeakHeight');
fprintf(1,'Please wait, finding spikes and clustering using Calinski-Harabasz-optimal means...000%%\n');
for ii = 1:12
    [~, spike.idx{ii}] = findpeaks(abs(spike_d.score(:,ii)), 'MinPeakHeight', 25, 'MinPeakDistance', 60);
    spike.clus{ii} = sign(spike_d.score(spike.idx{ii}, ii));
    spike.pidx{ii} = spike.idx{ii}(spike.clus{ii}==1);
    if numel(spike.pidx{ii}) < 8
        spike.pidx{ii} = [];
    else
        snips = signal_2_snippets(spike_d.score(:,ii), spike.pidx{ii});
        eva = evalclusters(snips,'kmeans','CalinskiHarabasz','KList',K_LIST);
        spike.pclus{ii} = kmeans(snips,eva.OptimalK);
        s = silhouette(snips, spike.pclus{ii});
        i_remove = s < SILHOUETTE_THRESHOLD;
        spike.pclus{ii}(i_remove) = [];
        spike.pidx{ii}(i_remove) = [];
    end
    spike.nidx{ii} = spike.idx{ii}(spike.clus{ii}==-1);
    if numel(spike.nidx{ii}) < 8
        spike.nidx{ii} = [];
    else
        snips = signal_2_snippets(spike_d.score(:,ii), spike.nidx{ii});
        eva = evalclusters(snips,'kmeans','CalinskiHarabasz','KList',K_LIST);
        spike.nclus{ii} = kmeans(snips,eva.OptimalK);
        s = silhouette(snips, spike.nclus{ii});
        i_remove = s < SILHOUETTE_THRESHOLD;
        spike.nclus{ii}(i_remove) = [];
        spike.nidx{ii}(i_remove) = [];
    end
    fprintf(1,'\b\b\b\b\b%03d%%\n', round(100*ii/12));
end
warning('on','signal:findpeaks:largeMinPeakHeight');

%% Plot the coefficients and corresponding eigenvector projections.
for ii = 1:12
    fig = plot_pc_spike_structure( t(i_plot), ...
        spike.pclus{ii}, spike.pidx{ii}, ...
        spike.nclus{ii}, spike.nidx{ii}, ...
        spike_d.coeff(:,ii), spike_d.score(:,ii), ...
        'CoefficientsName', sprintf('PC-%02d: Coefficients',ii), ...
        'PClusSyntax', sprintf('PC-%02d: PClus %%d', ii), ...
        'NClusSyntax', sprintf('PC-%02d: NClus %%d', ii));

    out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'Structure');
    if exist(out_folder,'dir')==0
        mkdir(out_folder);
    end
    default.savefig(fig, fullfile(out_folder, sprintf('%s_%d_PC-%02d_Structure', tank, TEMPLATE_BLOCK, ii)));
end

%% Plot the spikes
c_order = turbo(32);
c_order(i_exclude,:) = repmat([0.85 0.85 0.85],numel(i_exclude),1);
for ip = 1:7 
    [snips,t_snip,i_snip] = signal_2_snippets(spike_d.score(:,ip),spike.pidx{ip},'PreMilliseconds',3,'PostMilliseconds',3);
    pidx = spike.pidx{ip}(i_snip);
    for iclus = 1:max(spike.pclus{ip})
        tmp = pidx(spike.pclus{ip}(i_snip)==iclus);
        if size(tmp, 1) == 0
            continue;
        end
        i_sample = sort(randsample(numel(tmp),min(10,numel(tmp)),false),'ascend');
        for ik = 1:numel(i_sample)
            waves = nan(32,181);
            t_waves = ((-90):90)/30; % ms
            for ii = 1:32
                waves(ii,:) = signal_2_snippets(data_blanked(:,ii),tmp(i_sample(ik)),'PreMilliseconds',3,'PostMilliseconds',3);
            end
            proj = snips(i_sample(ik),:)' * spike_d.coeff(:,ip)';
            fig = figure('Position',[680   509   217   470], ...
                'Color','w'); 
            L = tiledlayout(fig, 1, 2);
            title(L, sprintf('PC-%02d P%02d', ip, iclus), 'FontName','Tahoma','Color','k');
            subtitle(L, sprintf('Snippet-%03d', i_sample(ik)), 'FontName','Tahoma','Color',[0.65 0.65 0.65]);
            % Plot the original waveforms
            ax = nexttile(L, 1, [1 1]);
            set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',c_order, ...
                'YDir','reverse','YTick',[250, 500, 750],'YLim',[180, 840],'XLim',[-3 3]);
            plot(ax, t_snip, waves./6 + ((0:20:620)+200)', 'LineWidth', 1.25);
            title(ax, 'Original', 'FontName','Tahoma','Color','k');
            xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
            ylabel(ax, 'Depth (\mum)','FontName','Tahoma','Color','k');

            % Plot the reconstruction, using just this eigenvector and its projection
            ax = nexttile(L, 2, [1 1]);
            set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',c_order, ...
                'YDir','reverse','YTick',[250, 500, 750],'YLim',[180, 840],'XLim',[-3 3]);
            plot(ax, t_snip, proj./3 + ((0:20:620)+200), 'LineWidth', 1.25);
            title(ax, 'Reconstructed', 'FontName','Tahoma','Color','k');
            xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
            ylabel(ax, 'Depth (\mum)','FontName','Tahoma','Color','k');

            out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'Snippets', sprintf('PC-%02d_PClus-%02d', ip, iclus));
            if exist(out_folder,'dir')==0
                mkdir(out_folder);
            end
            default.savefig(fig, fullfile(out_folder, sprintf('%s_%d_PC-%02d_PClus-%02d_Snippet-%03d', tank, TEMPLATE_BLOCK, ip, iclus, i_sample(ik))));
        end
    end

    [snips,t_snip,i_snip] = signal_2_snippets(spike_d.score(:,ip),spike.nidx{ip},'PreMilliseconds',3,'PostMilliseconds',3);
    nidx = spike.nidx{ip}(i_snip);
    for iclus = 1:max(spike.nclus{ip})
        tmp = nidx(spike.nclus{ip}(i_snip)==iclus);
        if size(tmp, 1) == 0
            continue;
        end
        i_sample = sort(randsample(numel(tmp),min(10,numel(tmp)),false),'ascend');
        for ik = 1:numel(i_sample)
            waves = nan(32,181);
            for ii = 1:32
                waves(ii,:) = signal_2_snippets(data_blanked(:,ii),tmp(i_sample(ik)),'PreMilliseconds',3,'PostMilliseconds',3);
            end
            proj = snips(i_sample(ik),:)' * spike_d.coeff(:,ip)';
            fig = figure('Position',[380   309   217   470], ...
                'Color','w'); 
            L = tiledlayout(fig, 1, 2);
            title(L, sprintf('PC-%02d N%02d', ip, iclus), 'FontName','Tahoma','Color','k');
            subtitle(L, sprintf('Snippet-%03d', i_sample(ik)), 'FontName','Tahoma','Color',[0.65 0.65 0.65]);
            % Plot the original waveforms
            ax = nexttile(L, 1, [1 1]);
            set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',c_order, ...
                'YDir','reverse','YTick',[250, 500, 750],'YLim',[180, 840],'XLim',[-3 3]);
            plot(ax, t_snip, waves./6 + ((0:20:620)+200)', 'LineWidth', 1.25);
            title(ax, 'Original', 'FontName','Tahoma','Color','k');
            xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
            ylabel(ax, 'Depth (\mum)','FontName','Tahoma','Color','k');

            % Plot the reconstruction, using just this eigenvector and its projection
            ax = nexttile(L, 2, [1 1]);
            set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',c_order, ...
                'YDir','reverse','YTick',[250, 500, 750],'YLim',[180, 840],'XLim',[-3 3]);
            plot(ax, t_snip, proj./3 + ((0:20:620)+200), 'LineWidth', 1.25);
            title(ax, 'Reconstructed', 'FontName','Tahoma','Color','k');
            xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
            ylabel(ax, 'Depth (\mum)','FontName','Tahoma','Color','k');

            out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'Snippets', sprintf('PC-%02d_NClus-%02d', ip, iclus));
            if exist(out_folder,'dir')==0
                mkdir(out_folder);
            end
            default.savefig(fig, fullfile(out_folder, sprintf('%s_%d_PC-%02d_NClus-%02d_Snippet-%03d', tank, TEMPLATE_BLOCK, ip, iclus, i_sample(ik))));
        end
    end
end

%% Save the spike data struct along with the filtered data.
% Beyond top-7 they start not looking so great.
clu = struct('id',{},'idx', {}, 'pc', {}, 'mean', {}, 'var', {});
id = 0;
for ii = 1:7
    pclus = spike.pclus{ii};
    pidx = spike.pidx{ii};
    [snips,~,i_snip] = signal_2_snippets(spike_d.score(:,ii),pidx);
    pclus = pclus(i_snip);
    pidx = pidx(i_snip);
    for ik = 1:min(max(pclus),2)
        idx = pclus == ik;
        clu(end+1) = struct('id',{id},'idx',{pidx(idx)},'pc',{ii},'mean',{mean(snips(idx,:),1)},'var',{var(snips(idx,:),0,1)}); %#ok<SAGROW> 
        id = id + 1;
    end
    nclus = spike.nclus{ii};
    nidx = spike.nidx{ii};
    [snips,~,i_snip] = signal_2_snippets(spike_d.score(:,ii),nidx);
    nclus = nclus(i_snip);
    nidx = nidx(i_snip);
    for ik = 1:min(max(nclus),2)
        idx = nclus == ik;
        clu(end+1) = struct('id',{id},'idx',{nidx(idx)},'pc',{ii},'mean',{mean(snips(idx,:),1)},'var',{var(snips(idx,:),0,1)}); %#ok<SAGROW> 
        id = id + 1;
    end
end
out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Data');
if exist(out_folder,'dir')==0
    mkdir(out_folder);
end
save(fullfile(out_folder, sprintf('%s_clu.mat', tank)), 'data_blanked', 'clu', '-v7.3');

%% Look at the spike time correlograms for all the units.
XC_VEC = -100:100; % ms
tau_xc = XC_VEC(2:end) - 0.5.*mean(diff(XC_VEC));
comparison = nchoosek(1:numel(clu),2);

nc = size(comparison,1);
xc = struct('comparison',mat2cell(comparison,ones(nc,1),2),'rho',cell(size(comparison,1),1));
isi = zeros(numel(clu),numel(tau_xc));

ix = 0;
fprintf(1,'Please waiting, computing isi & cross-correlations...000%%\n');
for ii = 1:(numel(clu)-1)
    ts_ii = 1e3 .* clu(ii).idx ./ fs; % ms
    for ik = 1:numel(ts_ii)
        isi(ii,:) = isi(ii,:) + histcounts(setdiff(ts_ii,ts_ii(ik)) - ts_ii(ik), XC_VEC);
    end
    isi(ii,:) = isi(ii,:) ./ numel(ts_ii);

    for ij = (ii+1):numel(clu)
        ix = ix + 1;
        ts_ij = 1e3 .* clu(ij).idx ./ fs; % ms
        xc(ix).rho = zeros(1,numel(tau_xc));
        for ik = 1:numel(ts_ii)
            xc(ix).rho = xc(ix).rho + histcounts(ts_ij - ts_ii(ik), XC_VEC);
        end
        xc(ix).rho = xc(ix).rho ./ numel(ts_ii);
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(100 * ix/nc));
    end
end
ts_ii = 1e3 .* clu(end).idx ./ fs; % ms
for ik = 1:numel(ts_ii)
    isi(end,:) = isi(end,:) + histcounts(setdiff(ts_ii,ts_ii(ik)) - ts_ii(ik), XC_VEC);
end
isi(end,:) = isi(end,:) ./ numel(ts_ii);

%% Generate ISI Plots
out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'ISI',sprintf('%dms-to-%dms',XC_VEC(1),XC_VEC(end)));
if exist(out_folder,'dir')==0
    mkdir(out_folder);
end

t_template = linspace(-0.4,0.8,37);
for ii = 1:size(isi,1)
    fig = figure('Color','w','Name','ISI','Position',[371   558   870   420]);
    L = tiledlayout(fig, 1, 2);
    name = sprintf('Unit-%02d--PC-%02d',clu(ii).id,clu(ii).pc);
    title(L, name,'FontName','Tahoma','Color','k');
    ax = nexttile(L,1,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    bar(ax, tau_xc, isi(ii,:), 1, 'EdgeColor','none','FaceColor','b');
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax,'Counts / Spike', 'FontName','Tahoma','Color','k');
    title(ax,'ISI', 'FontName','Tahoma','Color','k');

    ax = nexttile(L,2,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    errorbar(ax,  t_template, clu(ii).mean, clu(ii).var);
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax,sprintf('PC-%02d: Score', clu(ii).pc), 'FontName','Tahoma','Color','k');
    title(ax, 'Template','FontName','Tahoma','Color','k');
    default.savefig(fig,fullfile(out_folder, sprintf('%s_ISI',name)));
end

%% Generate Correlogram Plots
out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'Correlograms',sprintf('%dms-to-%dms',XC_VEC(1),XC_VEC(end)));
if exist(out_folder,'dir')==0
    mkdir(out_folder);
end

t_template = linspace(-0.4,0.8,37);
for ii = 1:numel(xc)
    fig = figure('Color','w','Name','Correlogram','Position',[371   558   870   420]);
    L = tiledlayout(fig, 2, 2);
    name = sprintf('Unit-%02d--PC-%02d-Unit-%02d--PC-%02d',clu(xc(ii).comparison(1)).id,clu(xc(ii).comparison(1)).pc,clu(xc(ii).comparison(2)).id,clu(xc(ii).comparison(2)).pc);
    title(L, name,'FontName','Tahoma','Color','k');
    ax = nexttile(L,1,[2 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    bar(ax, tau_xc, xc(ii).rho, 1, 'EdgeColor','none','FaceColor','b');
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax,'Counts / Spike', 'FontName','Tahoma','Color','k');
    title(ax,'Correlogram', 'FontName','Tahoma','Color','k');

    ax = nexttile(L,2,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    errorbar(ax,  t_template, clu(xc(ii).comparison(1)).mean, clu(xc(ii).comparison(1)).var);
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax,sprintf('PC-%02d: Score', clu(xc(ii).comparison(1)).pc), 'FontName','Tahoma','Color','k');
    title(ax, sprintf('Template_i (%d)', xc(ii).comparison(1)),'FontName','Tahoma','Color','k');

    ax = nexttile(L,4,[1 1]);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    errorbar(ax,  t_template, clu(xc(ii).comparison(2)).mean, clu(xc(ii).comparison(2)).var);
    xlabel(ax,'Time (ms)', 'FontName','Tahoma','Color','k');
    ylabel(ax,sprintf('PC-%02d: Score', clu(xc(ii).comparison(2)).pc), 'FontName','Tahoma','Color','k');
    title(ax, sprintf('Template_j (%d)', xc(ii).comparison(2)),'FontName','Tahoma','Color','k');
    default.savefig(fig,fullfile(out_folder, sprintf('%s_Correlogram',name)));
end

%% After checking correlograms, look to see how PC-1 and PC-4 combined checks out.
pidx = union(spike.pidx{1}(spike.pclus{1}==1), spike.pidx{4}(spike.pclus{4}==1));
pclus = ones(size(pidx));
fig = plot_pc_spike_structure( t(i_plot), ...
    pclus, pidx, ...
    [], [], ...
    sum(spike_d.coeff(:,[1,2,4]),2), sum(spike_d.score(:,[1,2,4]),2), ...
    'CoefficientsName', '\SigmaPC_{1,2,4}: Coefficients');

out_folder = fullfile(GEN_DATA_FOLDER, tank, 'Spike-Band PCA', 'Structure');
if exist(out_folder,'dir')==0
    mkdir(out_folder);
end
default.savefig(fig, fullfile(out_folder, sprintf('%s_%d_PC-01_PC-02_PC-04_combined_Structure', tank, TEMPLATE_BLOCK)));