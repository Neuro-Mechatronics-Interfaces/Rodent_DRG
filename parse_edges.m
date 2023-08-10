function [i_sync, g_sync, id_sync] = parse_edges(signal, edge, options)
%PARSE_EDGES Detects edge transition indices (first HIGH for rising, last HIGH for falling).
%
% Syntax:
%   [i_sync, g_sync] = parse_edges(signal, edge);
%
% Inputs:
%   sync - Digital logic vector
%   edge - 'rising' (default) | 'falling'
%
% Output:
%   i_sync - Sample indices for edge transition
%   g_sync - "Frequency grouping" for each sync pulse
%   id_sync - Frequency-ID corresponding to elements of sync groupings
%
% See also: Contents

arguments
    signal (1,:) double {mustBeNumeric, mustBeVector}
    edge {mustBeTextScalar, mustBeMember(edge, {'rising', 'falling'})} = 'rising'
    options.SampleRate (1,1) double = 30000;
    options.MinFrequency (1,1) double = 0.25;
    options.MaxFrequency (1,1) double = 300;
    options.MaxIterations (1,1) double = 10;
    options.NDecimalPointsRound (1,1) double {mustBeNumeric, mustBeInteger} = 1;
end

d = [0, diff(signal)];

switch edge
    case 'rising'
        i_sync = find(d < 0);
    case 'falling'
        i_sync = find(d > 0) - 1;
end

t_diff = diff(i_sync);
t_diff = [t_diff(1), t_diff];

has_overfrequencies = true;
i_count = 0;
while has_overfrequencies && (i_count < options.MaxIterations)
    i_slow = find(t_diff > (options.SampleRate / options.MinFrequency));
    for ii = 1:numel(i_slow)
        if (i_slow(ii)+1) <= numel(t_diff)
            t_diff(i_slow(ii)) = t_diff(i_slow(ii)+1);
        end
    end
    freq = round(1 ./ (t_diff ./ options.SampleRate),options.NDecimalPointsRound);
    i_bad = freq > options.MaxFrequency;
    t_diff(i_bad) = [];
    has_overfrequencies = sum(i_bad) > 0;
    i_count = i_count + 1;
end
if i_count == options.MaxIterations
    warning('Iteration limit reached, frequencies may still contain over-frequencies.');
end
[g_sync, id_sync] = findgroups(freq);

end