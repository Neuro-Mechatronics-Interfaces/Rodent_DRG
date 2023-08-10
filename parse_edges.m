function [i_sync, g_sync, id_sync] = parse_edges(signal, edge, options)
%PARSE_EDGES Detects edge transition indices (first HIGH for rising, last HIGH for falling).
%
% Syntax:
%   [i_sync, g_sync] = parse_edges(signal, edge);
%
% Inputs:
%   sync - Digital logic vector
%   edge - 'rising' (default) | 'falling' | 'both'
%
% Output:
%   i_sync - Sample indices for edge transition
%   g_sync - "Frequency grouping" for each sync pulse
%   id_sync - Frequency-ID corresponding to elements of sync groupings
%
% See also: Contents

arguments
    signal (1,:) double {mustBeNumeric, mustBeVector}
    edge {mustBeTextScalar, mustBeMember(edge, {'rising', 'falling', 'both'})} = 'rising'
    options.SampleRate (1,1) double = 30000;
    options.MinFrequency (1,1) double = 0.25;
    options.NDecimalPointsRound (1,1) double {mustBeNumeric, mustBeInteger} = 1;
end

d = [0, diff(signal)];

switch edge
    case 'rising'
        i_sync = find(d < 0);
    case 'falling'
        i_sync = find(d > 0) - 1;
    case 'both'
        i_sync = union(find(d < 0), find(d > 0)-1);
end

t_diff = diff(i_sync);
t_diff = [t_diff(1), t_diff];

i_slow = find(t_diff > (options.SampleRate / options.MinFrequency));

for ii = 1:numel(i_slow)
    if (i_slow(ii)+1) <= numel(t_diff)
        t_diff(i_slow(ii)) = t_diff(i_slow(ii)+1);
    end
end
freq = round(1 ./ (t_diff ./ options.SampleRate),options.NDecimalPointsRound);
[g_sync, id_sync] = findgroups(freq);

end