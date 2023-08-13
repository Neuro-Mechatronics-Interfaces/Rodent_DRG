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

% Set parameters
N_MS_PRE_SYNC = 6;
N_MS_POST_SYNC = 6;

SUBJ = 'DW100';
RAW_DATA_FOLDER = parameters('raw_data_folder');
YYYY = 2023;
MM = 8;
DD = 9;
BLOCK = 1:3; % options for block
SYNC_EDGE = ["rising", "falling"];
STIM_MODE = [0,1,3];

TRIAL_AMPLITUDE_SCALING = 120; % microvolts, this should roughly correspond to the largest expected deviations so that trials fit nicely with most values scaled between ~-0.5 to 0.5, with outliers maybe going as large as -1.0 to 1.0.

PLOT_HEADSTAGE_MOVEMENT = true;
N_EXAMPLE_TRIALS = 10;   % How many trials to plot, per frequency block?
        
for sync_edge = SYNC_EDGE
    for block = BLOCK
        block_name = sprintf("%s_%04d_%02d_%02d_%d", SUBJ, YYYY, MM, DD, block);
        out_folder = fullfile(pwd, 'figures', block_name, sync_edge);

        % Load data (see comment below first!)
        % Need to have raptor datashare (\\raptor.ni.cmu.edu\NML) mapped -- see the
        % lab wiki for instructions on how to do that.
        x = io.load_data(SUBJ, YYYY, MM, DD, "*", block, '.rhd', RAW_DATA_FOLDER);
        
        % Parse data used in each plot
        % Get timing vector
        fs = x.frequency_parameters.amplifier_sample_rate;
        t = x.t_amplifier;
        acc = (x.aux_input_data - mean(x.aux_input_data,2))/0.34; % Acceleration approximately (g's)
        p = cumtrapz(x.t_aux_input, acc, 2);
        
        if exist(out_folder, 'dir')==0
            mkdir(out_folder);
        end
        
        % Get sync array
        [i_sync, g_sync, id_sync] = parse_edges(x.board_dig_in_data(4,:), sync_edge, ...
            'StimMode', STIM_MODE(block));
        n_groups = max(g_sync);
        cdata = turbo(n_groups);
        vec = (-round(N_MS_PRE_SYNC*0.001*fs):round(N_MS_POST_SYNC*0.001*fs))';
        t_snips = 1e3 .* vec ./ fs;
        mask = i_sync + vec;
        mask(:, any(mask < 1, 1) | any(mask > numel(t), 1)) = []; % Remove out-of-bounds samples
        spike_band = struct;
        [spike_band.b,spike_band.a] = butter(4, [300 3000]./(fs/2), 'bandpass');
        
        i_sample = cell(n_groups,1);
        for ii = 1:n_groups
            tmp = find(g_sync == ii);
            i_sample{ii} = tmp(1:min(numel(tmp), N_EXAMPLE_TRIALS));
        end
        i_sample = horzcat(i_sample{:});
        k = numel(i_sample);
        mask = mask(:, i_sample);
        
        % Plot the one example figure for the full block
        if PLOT_HEADSTAGE_MOVEMENT
            fig = figure('Color','w', ...
                         'Name',sprintf('Amplifier Channel A%03d EDA', 23), ...
                         'Position', [200 100 500 600]); 
            L = tiledlayout(fig, 2, 1);
            title(L, sprintf('Channel-A%03d EDA', 23), 'FontName','Tahoma','Color','k');
            xlabel(L, 'Time (s)', 'FontName','Tahoma','Color','k');
            ax = nexttile(L);
            set(ax,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
            plot(ax, t, x.amplifier_data(24,:), 'b-', ...           % Plot channel data
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
            
            default.savefig(fullfile(out_folder, sprintf('Channel-A-%03d_All', 23)));
        end

        % Plot figures for each channel
        for iCh = 1:numel(x.amplifier_channels)  
            % Identify sync events
            y = x.amplifier_data(iCh,:);
            if iCh < 33
                yf = filtfilt(spike_band.b,spike_band.a,y);
            else
                yf = y;
            end
            snips = yf(mask);      
            if iCh > 32
                snips = snips - median(snips,1);
            end
        
            fig = figure('Color', 'w', 'Name', 'Rising Edge Snippets','Position',[680   224   560   751]);
            L = tiledlayout(fig, 2, 1);
            ax = nexttile(L, 1, [2 1]);
            set(ax,'NextPlot','add','FontName','Tahoma', 'ColorOrder', cdata(g_sync(i_sample), :), ...
                'XLim', [t_snips(1), t_snips(end)], 'YLim',[0, k+1]);
            plot(ax, t_snips, snips./TRIAL_AMPLITUDE_SCALING+(1:k));
            ylabel(ax, 'Rising Edge Number', 'FontName','Tahoma',"Color",'k');
            xlabel(ax, 'Time From Rising Edge (ms)', 'FontName','Tahoma',"Color",'k');
            
            ax = nexttile(L, 'east');
            set(ax, 'NextPlot','add','XLim',[0 1], 'YLim', [0.5 size(cdata,1)+0.5],...
                'YTick',1:size(cdata,1),'YTickLabel', id_sync);
            ylabel(ax, 'Frequency (Hz)','FontName','Tahoma','Color','k');
            imagesc(ax, [0 1], [1 size(cdata,1)], reshape(cdata,[],1,3));
            title(L, strrep(block_name,'_','\_'), 'FontName','Tahoma','Color','k');
            subtitle(L, sprintf('Channel-%s Snippets', x.amplifier_channels(iCh).native_channel_name), 'FontName','Tahoma','Color',[0.65 0.65 0.65]);
            default.savefig(fullfile(out_folder, sprintf('Channel-%s_Snippets', x.amplifier_channels(iCh).native_channel_name)));
        
        end
    end
end