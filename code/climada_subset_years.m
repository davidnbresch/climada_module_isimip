function obj_sub=climada_subset_years(obj, obj_type, years_in, silent_mode)
% climada template
% MODULE:
%   isimip
% NAME:
%   climada_subset_years
% PURPOSE:
%   Subset an object (entity, hazard, em_dat, ...) to include only specific
%   years
%
%   previous call: -
%   next call: -
% CALLING SEQUENCE:
%   em_data=emdat_read('','USA',['FL';'F1';'F2'],2005,0);
%   
% EXAMPLE:
%   em_data=emdat_read('','USA',['FL';'F1';'F2'],2005,0);
%   em_data_subset=climada_subset_years(em_data, 'obs', 1995:2005)
% INPUTS:
%   obj: object to be subset
%   obj_type: type of object, can be 'hazard', 'entity', 'obs'
%   years_in: vector of years to be retained in the output
% OPTIONAL INPUT PARAMETERS:
%   silent_mode: =1 prints a confirmation, prevented if =0 (default), 
% OUTPUTS:
%   obj_sub: as input 'obj' but with only the years specified by years_in
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180918, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180920, bug fixes
%-

damage_data=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('obj','var'),error('Error: obj is not given');return;end
if ~exist('obj_type','var'),error('Error: obj_type is not given');return;end
if ~exist('years_in','var'),error('Error: years_in is not given');return;end
if ~exist('silent_mode','var'),silent_mode=1;end

% by default or if anything fails, returns the original object
obj_sub = obj;

% check that 'obj_type' has a valid value
if ~ismember(obj_type, {'hazard', 'entity', 'obs'})
    error('***** Error: obj_type had unexpected value *****')
    return
end

errorStruct.message = ['*** Error: not all requested years are included in the ' obj_type ' ***'];
errorStruct.identifier = 'climada_subset_years:missingYear';


if strcmp(obj_type, 'entity')
    obj_yyyy_pos=find(ismember(obj.assets.Values_yyyy, years_in));
    if length(obj_yyyy_pos) < length(years_in)
        error(errorStruct);
    end
    obj_sub.assets.Values=obj_sub.assets.Values(obj_yyyy_pos,:);
    obj_sub.assets.Values_yyyy=obj_sub.assets.Values_yyyy(obj_yyyy_pos,:);
elseif strcmp(obj_type, 'hazard')
    obj_yyyy_pos=find(ismember(obj_sub.yyyy, years_in));
    if length(obj_yyyy_pos) < length(years_in)
        error(errorStruct);
    end
    obj_sub.intensity       =obj_sub.intensity(obj_yyyy_pos,:);
    obj_sub.fraction        =obj_sub.fraction(obj_yyyy_pos,:);
    obj_sub.frequency       =ones(1,length(obj_yyyy_pos))/length(obj_yyyy_pos);
    obj_sub.event_count     =length(obj_yyyy_pos);
    obj_sub.event_ID        =obj_sub.event_ID(obj_yyyy_pos);
    obj_sub.yyyy            =obj_sub.yyyy(obj_yyyy_pos);
    obj_sub.datenum         =obj_sub.datenum(obj_yyyy_pos);
    obj_sub.mm              =obj_sub.mm(obj_yyyy_pos);
    obj_sub.dd              =obj_sub.dd(obj_yyyy_pos);
    obj_sub.name            =obj_sub.name(obj_yyyy_pos);
elseif strcmp(obj_type, 'obs')
    obj_yyyy_pos=find(ismember(obj.year, years_in));
    if length(obj_yyyy_pos) < length(years_in)
        error(errorStruct);
    end
    obj_sub.values=obj_sub.values(obj_yyyy_pos);
    obj_sub.year==obj_sub.year(obj_yyyy_pos);
else
    error('***** Error: obj_type had unexpected value *****')
end

if ~silent_mode
    fprintf('\nYears successfully subset for object');
end

return

end % climada_subset_years