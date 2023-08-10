%EDA__DW100_2023_08_09 Exploratory data analysis (EDA) for DW100 experiment conducted on 2023-08-09. Pertains specifically to "flappy/tappy" pilot thing.
%
% Note: running this script will produce a folder named `figures` in the
% current workspace folder, with a subfolder corresponding to a specific
% "block name" which is just the subject name, date in ISO-8601 order, and
% block number which is a key indexing which experiment it was (i.e. "zero"
% "one" "two" and "three" in the original filenames). 

close all force;
clear;
clc;

%% Set parameters
SUBJ = 'DW100';
BLOCK = 2;

YYYY = 2023;
MM = 8;
DD = 9;

BLOCK_NAME = sprintf("%s_%04d_%02d_%02d_%d", SUBJ, YYYY, MM, DD, BLOCK);
OUT_FOLDER = fullfile(pwd, 'figures', BLOCK_NAME);

N_DECIMAL_ROUNDING = 1; % For block-1 use 0, all others use 2 (for categorical grouping of sync pulses based on estimated frequency of stimuli)
SYNC_EDGE = 'rising'; % 'rising' | 'falling' | 'both'

%% Load data (see comment below first!)
% Need to have raptor datashare (\\raptor.ni.cmu.edu\NML) mapped -- see the
% lab wiki for instructions on how to do that.
x = io.load_data(SUBJ, YYYY, MM, DD, "*", BLOCK, '.rhd', 'R:/NMLShare/raw_data/rodent');

%% Parse data used in each plot
% Get timing vector
fs = x.frequency_parameters.amplifier_sample_rate;
t = x.t_amplifier;
acc = (x.aux_input_data - mean(x.aux_input_data,2))/0.34; % Acceleration approximately (g's)
p = cumtrapz(x.t_aux_input, acc, 2);

if exist(OUT_FOLDER, 'dir')==0
    mkdir(OUT_FOLDER);
end

% Get sync array
[i_sync, g_sync, id_sync] = parse_edges(x.board_dig_in_data(4,:), SYNC_EDGE, ...
    'NDecimalPointsRound', N_DECIMAL_ROUNDING);
cdata = turbo(max(g_sync));
vec = (-60:450)'; % -2 : +15 ms
t_snips = 1e3 .* vec ./ fs;
mask = i_sync + vec;
mask(:, any(mask < 1, 1) | any(mask > numel(t), 1)) = []; % Remove out-of-bounds samples
[b,a] = butter(4, [300 3000]./(fs/2), 'bandpass');

%% Plot figures for each channel
for iCh = 1:32
    fig = figure('Color','w', ...
                 'Name',sprintf('Amplifier Channel A%03d EDA', iCh-1), ...
                 'Position', [200 100 500 600]);
    L = tiledlayout(fig, 2, 1);
    title(L, sprintf('Channel-A%03d EDA', iCh), 'FontName','Tahoma','Color','k');
    xlabel(L, 'Time (s)', 'FontName','Tahoma','Color','k');
    ax = nexttile(L);
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
    plot(ax, t, x.amplifier_data(iCh,:), 'b-', ...           % Plot channel data
             t, x.board_dig_in_data(4,:)*50 + 300, 'k-');   % Offset DIG-SYNC
    title(ax, 'Amplifier and Sync Data', 'FontName','Tahoma','Color','k');
    ylabel(ax, 'DRG Amplitude (\muV)', 'FontName','Tahoma','Color','k');
    
    ax = nexttile(L);
    yyaxis(ax, 'left');
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k', 'ColorOrder', [1 0 0; 0 1 0; 0 0 1], 'LineStyleOrder', {'-'});
    plot(ax, x.t_aux_input, acc);
    ylabel(ax, 'Headstage Acceleration (g)', 'FontName','Tahoma','Color','k');
    
    yyaxis(ax, 'right');
    set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k', 'ColorOrder', [1 0 0; 0 1 0; 0 0 1], 'LineStyleOrder', {'--'});
    plot(ax, x.t_aux_input, p);
    title(ax, 'Headstage Movement', 'FontName','Tahoma','Color','k');
    ylabel(ax, '(\intAcceleration)', 'FontName','Tahoma','Color','k');
    legend(ax, ["a_x", "a_y", "a_z", "\inta_x", "\inta_y", "\inta_z"], 'FontName','Tahoma');
    
    default.savefig(fullfile(OUT_FOLDER, sprintf('Channel-A-%03d_All', iCh-1)));
    
    % Identify sync events
    y = x.amplifier_data(iCh,:);
    yf = filtfilt(b,a,y);
    snips = yf(mask);
    n = size(snips, 2);
    k = min(100, size(snips,2));
    
    i_sample = randsample(n, k, false);
    i_sample = sort(i_sample, 'ascend');
    
    fig = figure('Color', 'w', 'Name', 'Rising Edge Snippets');
    L = tiledlayout(fig, 2, 1);
    ax = nexttile(L, 1, [2 1]);
    set(ax,'NextPlot','add','FontName','Tahoma', 'ColorOrder', cdata(g_sync(i_sample), :), ...
        'XLim', [t_snips(1), t_snips(end)], 'YLim',[0, k]);
    plot(ax, t_snips, snips(:,i_sample)./120+(1:k));
    ylabel(ax, 'Rising Edge Number', 'FontName','Tahoma',"Color",'k');
    xlabel(ax, 'Time From Rising Edge (ms)', 'FontName','Tahoma',"Color",'k');
    
    ax = nexttile(L, 'east');
    set(ax, 'NextPlot','add','XLim',[0 1], 'YLim', [0.5 size(cdata,1)+0.5],...
        'YTick',1:size(cdata,1),'YTickLabel', id_sync);
    ylabel(ax, 'Frequency (Hz)','FontName','Tahoma','Color','k');
    imagesc(ax, [0 1], [1 size(cdata,1)], reshape(cdata,[],1,3));
    default.savefig(fullfile(pwd, 'figures', BLOCK_NAME, sprintf('Channel-A-%03d_Snippets', iCh-1)));

end