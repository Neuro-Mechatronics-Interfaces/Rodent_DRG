%EXAMPLE_GENERATE_SESSION_SPIKE_TEMPLATES Example of generating spike templates for each channel, for a session (tank), which can then be used for each separate recording to correspond clusters (for different stimulus conditions etc.) 
close all force;
clear;
clc;

SUBJ = 'DW100';
RAW_DATA_FOLDER = parameters('raw_data_folder');
YYYY = 2023;
MM = 8;
DD = 9;
CHANNELS = 1:32; % Amplifier channels to retain

TEMPLATE_BLOCK = 2; % Block index to use for generating templates
EPOCH = [10, 70]; % Seconds to use for detection/template generation

%% Load data
% Load data (see comment below first!)
% Need to have raptor datashare (\\raptor.ni.cmu.edu\NML) mapped -- see the
% lab wiki for instructions on how to do that.
x = io.load_data(SUBJ, YYYY, MM, DD, "*", TEMPLATE_BLOCK, '.rhd', RAW_DATA_FOLDER);

%% Extract relevant endpoints and clear data struct from memory
fs = x.frequency_parameters.amplifier_sample_rate;
t = x.t_amplifier((x.t_amplifier >= EPOCH(1)) & (x.t_amplifier < EPOCH(2)));
data = x.amplifier_data(CHANNELS, (x.t_amplifier >= EPOCH(1)) & (x.t_amplifier < EPOCH(2)));
clear x;

spike_band = struct;
[spike_band.b,spike_band.a] = butter(4, [300 3000]./(fs/2), 'bandpass');
data_f = filtfilt(spike_band.b, spike_band.a, data')';

%% Get the channel indexing and a masking vector for plotting data snapshot
[i_intan, i_NN, depth_NN] = intan_2_neuronexus();
i_plot = (t >= 55) & (t < 65);

%% Make a figure with the data streams shown as they would be on probe
fig = figure('Color','w','Name','NN DRG Depth-wise Layout');
ax = axes(fig,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k','ColorOrder',turbo(32), ...
    'YDir','reverse','YTick',[250, 500, 750],'YLim',[180, 840]);
ylabel(ax, 'Depth (\mum)', 'FontName','Tahoma','Color','k');
xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');
data_depth = -data_f(i_intan,i_plot)./6 + ((0:20:620) + 200); % Depths in microns
plot(ax, t(i_plot), data_depth);
