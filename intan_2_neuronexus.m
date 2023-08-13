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
    options.ConnectorType {mustBeTextScalar, mustBeMember(options.ConnectorType, {'CM32', 'CM16', 'HC16', 'HC32', 'A32-OM32', 'OM32'})} = 'A32-OM32';
    options.ProbeType {mustBeTextScalar, mustBeMember(options.ProbeType, {'A1x32', 'A1x32-Edge', 'A1x16', 'A4x8', 'A4x4', 'A1x32-Poly2', 'A2x16', 'A2x8'})} = 'A1x32-Edge';
    options.TipDepth    (1,1) double = 800; % Microns at tip (900 in, minus 100)
    options.ShankLength (1,1) double = 620; % Total distance from top electrode to bottom electrode (center-to-center), microns.
end

switch options.ProbeType % "Shallowest" probe is listed first in the array. Multi-shank arrays will have 1 column per shank.
    case 'A1x32'
        N_Per_Shank = 32;
        i_NN = [17,16,18,15,19,14,20,13,21,12,22,11,23,10,24,9,25,8,26,7,27,6,28,5,29,4,30,3,31,2,32,1]';
    case 'A1x32-Edge'
        N_Per_Shank = 32;
        i_NN = (32:-1:1)';
    otherwise
        error("ProbeType: '%s' has not yet been implemented, sorry. See: https://www.neuronexus.com/files/probemapping/32-channel/CM32-Maps.pdf to fill in the corresponding 'case' statement above.", options.ProbeType);
end

% Array element 1 == IN0, Array element 2 == IN1, etc.
switch options.ConnectorType 
    case 'CM32' 
        headstage_conn_2_electrode_conn_mapping = [17, 18, 32, 31, 29, 27, 25, 23, 21, 22, 24, 26, 28, 30, 20, 19, 14, 13, 3, 5, 7, 9, 11, 12, 10, 8, 6, 4, 2, 1, 15, 16]';
    case 'A32-OM32' % A32-OM32 adapter
        headstage_conn_2_adaptor_mapping = [22, 18, 20, 32, 30, 28, 26, 24, 23, 25, 27, 29, 31, 19, 17, 21, 11, 15, 13, 1, 3, 5, 7, 9, 10, 8, 6, 4, 2, 14, 16, 12]';
        electrode_2_adaptor_conn_mapping = [16, 6, 5, 15, 4, 7, 3, 8, 2, 9, 1, 10, 14, 13, 12, 11, 22, 21, 20, 19, 23, 25, 24, 18, 26, 17, 27, 29, 28, 31, 30, 32]'; % Index-1 is NN-1, mapping gives which adaptor pin connects with NN-1 contact.
        headstage_conn_2_electrode_conn_mapping = nan(size(headstage_conn_2_adaptor_mapping));
        for ii = 1:numel(headstage_conn_2_adaptor_mapping)
            headstage_conn_2_electrode_conn_mapping(ii) = find(electrode_2_adaptor_conn_mapping == headstage_conn_2_adaptor_mapping(ii), 1, 'first');
        end
    case 'OM32' % Jordyn's pinout from ACUTE_OM32_v2
        headstage_conn_2_electrode_conn_mapping = [22, 18, 20, 32, 30, 28, 26, 24, 23, 25, 27, 29, 31, 19, 17, 21, 11, 15, 13, 1, 3, 5, 7, 9, 10, 8, 6, 4, 2, 14, 16, 12]';
    otherwise
        error("ConnectorType: '%s' has not yet been implemented, sorry. See: https://www.neuronexus.com/files/wiringconfiguration/Wiring_CM32.pdf to fill in the corresponding 'case' statement above.", options.ConnectorType);
end

depths = linspace(options.TipDepth-options.ShankLength, options.TipDepth, N_Per_Shank);
depth_NN = nan(size(i_NN));
shank_NN = nan(size(i_NN));
i_intan = nan(size(i_NN));
for ii = 1:numel(i_NN)
    i_intan(ii) = find(headstage_conn_2_electrode_conn_mapping == i_NN(ii), 1, 'first');
    [row, col] = ind2sub(size(i_NN), ii);
    depth_NN(ii) = depths(row);
    shank_NN(ii) = col;
end

end