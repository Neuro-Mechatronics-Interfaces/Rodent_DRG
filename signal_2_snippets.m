function [snips, t_snip, i_snip] = signal_2_snippets(signal, i_align, options)
%SIGNAL_2_SNIPPETS  Chop up signal into snippets and return array.
%
% Syntax:
%   [snips, t_snip, i_snip] = signal_2_snippets(signal, i_align, 'Name', value, ...);
%
% Inputs:
%   signal (:,1) double - The signal to chop up
%   i_align (:,1) double - Samples to use for alignment
%
%   Options:
%       SampleRate (1,1) double = 30000
%       PreMilliseconds (1,1) double = 0.4
%       PostMilliseconds (1,1) double = 0.8
%
% Output:
%   snips - Array where rows are observations, columns are time samples
%   t_snip - Times (ms) of each column
%
% See also: Contents

arguments
    signal (:,1) double
    i_align (:,1) double
    options.SampleRate (1,1) double = 30000;
    options.PreMilliseconds (1,1) double = 0.4;
    options.PostMilliseconds (1,1) double = 0.8;
end

vec = round(-options.PreMilliseconds*1e-3*options.SampleRate):round(options.PostMilliseconds*1e-3*options.SampleRate);

t_snip = 1e3 .* vec ./ options.SampleRate; % ms

mask = i_align + vec;
i_snip = (1:numel(i_align))';
if isempty(i_align)
    snips = zeros(size(t_snip));
    return;
end
i_remove = any(mask < 1, 2) | any(mask > numel(signal), 2);
mask(i_remove,:) = [];
i_snip(i_remove) = [];

snips = signal(mask);

end