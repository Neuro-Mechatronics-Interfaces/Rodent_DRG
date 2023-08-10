function [i_intan, i_NN, depth_NN, shank_NN] = intan_2_neuronexus(options)
%INTAN_2_NEURONEXUS  Convert NeuroNexus channels to Intan channels
%
% Syntax:
%   [i_intan, i_NN, depth_NN, shank_NN] = intan_2_neuronexus('Name', value, ...);
%
% Output:
%   i_intan     - Indices for corresponding elements of i_NN. Arrangement
%                   corresponds with shank layout, but generally this
%                   should be a 1-indexed vector that goes in order from
%                   1:number of electrodes.
%   i_NN        - Indices corresponding to the Neuro-Nexus electrode. Use
%                   this like amplifier_data(i_NN(1),:) to index the top
%                   electrode (most-superficial), on the first shank. The
%                   depth info will correspond with the value returned in
%                   `depth_NN`. 
%   depth_NN    - Depths corresponding to each element in `i_NN`
%   shank_NN    - Shank indicating which shank the electrode is on. 
%
% See also: Contents, 
%           https://www.neuronexus.com/files/probemapping/32-channel/CM32-Maps.pdf
%           https://www.neuronexus.com/files/wiringconfiguration/Wiring_CM32.pdf

arguments
    options.ConnectorType {mustBeTextScalar, mustBeMember(options.ConnectorType, {'CM32', 'CM16', 'HC16', 'HC32'})} = 'CM32';
    options.ProbeType {mustBeTextScalar, mustBeMember(options.ProbeType, {'A1x32', 'A1x16', 'A4x8', 'A4x4', 'A1x32-Poly2', 'A2x16', 'A2x8'})} = 'A1x32';
    options.TipDepth    (1,1) double = 800; % Microns at tip (900 in, minus 100)
    options.ShankLength (1,1) double = 620; % Total distance from top electrode to bottom electrode (center-to-center), microns.
end

switch options.ProbeType % "Shallowest" probe is listed first in the array. Multi-shank arrays will have 1 column per shank.
    case 'A1x32'
        i_NN = [17,16,18,15,19,14,20,13,21,12,22,11,23,10,24,9,25,8,26,7,27,6,28,5,29,4,30,3,31,2,32,1]';
    case 'A1x32-Poly2'
        i_NN = (32:-1:1)';
    otherwise
        error("ProbeType: '%s' has not yet been implemented, sorry. See: https://www.neuronexus.com/files/probemapping/32-channel/CM32-Maps.pdf to fill in the corresponding 'case' statement above.", options.ProbeType);
end

switch options.ConnectorType 
    case 'CM32'
        N_Per_Shank = 32;
        mapping = [16, 15, 1, 2, 4, 6, 8 10, 12, 11, 9, 7, 5, 3, 13, 14, 19, 20, 30, 28, 26, 24, 22, 21, 23, 25, 27, 29, 31, 32, 18, 17]';
    otherwise
        error("ConnectorType: '%s' has not yet been implemented, sorry. See: https://www.neuronexus.com/files/wiringconfiguration/Wiring_CM32.pdf to fill in the corresponding 'case' statement above.", options.ConnectorType);
end

depths = linspace(options.TipDepth-options.ShankLength, options.TipDepth, N_Per_Shank);
depth_NN = nan(size(i_NN));
shank_NN = nan(size(i_NN));
i_intan = nan(size(i_NN));
for ii = 1:numel(i_NN)
    i_intan(ii) = find(mapping == i_NN(ii), 1, 'first');
    [row, col] = ind2sub(size(i_NN), i_intan(ii));
    depth_NN(ii) = depths(row);
    shank_NN(ii) = col;
end

end