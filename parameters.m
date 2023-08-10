function varargout = parameters(name)
%PARAMETERS Return parameters struct, which sets default values for things like epoch durations etc.
%
% Example 1
%   f = parameters('raw_data_folder'); % Returns raw data folder value
%
% See also: Contents, io.load_data

arguments (Repeating)
    name {mustBeTextScalar, mustBeMember(name, {'raw_data_folder', 'generated_data_folder'})}
end

pars = struct;
% Add fields to `pars` struct here, and make sure to add to the validator
% in repeating `arguments` block above as well. The arguments block allows
% you to use "tab-complete" which is super-helpful when the list of
% possible parameters is long!
pars.raw_data_folder = 'R:/NMLShare/raw_data/rodent';
pars.generated_data_folder = 'R:/NMLShare/generated_data/rodent/Acute_DRG';

N = numel(name);
if nargout == 1
    if rem(N, 2) == 1
        varargout = {pars.(name{end})};
        return;
    else
        f = fieldnames(pars);
        for iV = 1:2:N
            idx = strcmpi(f, name{iV});
            if sum(idx) == 1
               pars.(f{idx}) = name{iV+1}; 
            end
        end
        varargout = {pars};
        return;
    end
else
    f = fieldnames(pars);
    varargout = cell(1, nargout);
    for iV = 1:numel(varargout)
        idx = strcmpi(f, name{iV});
        if sum(idx) == 1
            varargout{iV} = pars.(f{idx}); 
        else
            error('Could not find parameter: %s', name{iV}); 
        end
    end
end

end
