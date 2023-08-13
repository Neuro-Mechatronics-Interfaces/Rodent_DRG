function [i_sync, g_sync, id_sync, freq] = parse_edges(signal, edge, options)
%PARSE_EDGES Detects edge transition indices (first HIGH for rising, last HIGH for falling).
%
% Syntax:
%   [i_sync, g_sync, id_sync, freq] = parse_edges(signal, edge);
%
% Inputs:
%   sync - Digital logic vector
%   edge - 'rising' (default) | 'falling'
%
% Output:
%   i_sync - Sample indices for edge transition
%   g_sync - "Frequency grouping" for each sync pulse
%   id_sync - Frequency-ID corresponding to elements of sync groupings
%   freq    - The actual frequency estimates (should match up with id_sync at least somewhat closely).
%
% See also: Contents

arguments
    signal (1,:) double {mustBeNumeric, mustBeVector}
    edge {mustBeTextScalar, mustBeMember(edge, {'rising', 'falling'})} = 'rising'
    options.SampleRate (1,1) double = 30000;
    options.StimMode (1,1) double {mustBeMember(options.StimMode, 0:3)} = 1;
    options.MaxIterations (1,1) double = 10;
    options.NDecimalPointsRound (1,1) double {mustBeNumeric, mustBeInteger} = 0;
end

switch options.StimMode
    case 0
        stim_freqs = [5, 10, 20];
    case {1, 2}
        stim_freqs = [1, 2, 4, 8, 16]; 
    case 3
        stim_freqs = [2, 5, 10, 15, 20, 25, 30, 50, 80, 200];
end

d = [0, diff(signal)];

switch edge
    case 'rising'
        i_sync = find(d > 0) + 1;
    case 'falling'
        i_sync = find(d < 0);
end

t_diff = diff(i_sync);
t_diff = [t_diff(1), t_diff];

has_overfrequencies = true;
i_count = 0;
min_freq = min(stim_freqs);
max_freq = max(stim_freqs);

while has_overfrequencies && (i_count < options.MaxIterations)
    i_slow = find(t_diff > (options.SampleRate / min_freq));
    for ii = 1:numel(i_slow)
        if (i_slow(ii)+1) <= numel(t_diff)
            t_diff(i_slow(ii)) = t_diff(i_slow(ii)+1);
        end
    end
    freq = round(1 ./ (t_diff ./ options.SampleRate),options.NDecimalPointsRound);
    i_bad = freq > max_freq;
    t_diff(i_bad) = [];
    i_sync(i_bad) = [];
    freq(i_bad) = [];
    has_overfrequencies = sum(i_bad) > 0;
    i_count = i_count + 1;
end
if i_count == options.MaxIterations
    warning('Iteration limit reached, frequencies may still contain over-frequencies.');
end
freq_prog = freq;
for ii = 1:numel(freq)
    [~,i_prog] = min(abs(stim_freqs - freq(ii)));
    freq_prog(ii) = stim_freqs(i_prog);
end
[g_sync, id_sync] = findgroups(freq_prog);

end