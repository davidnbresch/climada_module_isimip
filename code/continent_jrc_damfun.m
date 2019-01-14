function [continent,jrc_file]=continent_jrc_damfun(country, params)
% climada isimip find continent and JRC flood damage function file
% MODULE:
%   isimip
% NAME:
%   continent_jrc_damfun
% PURPOSE:
%   Identify the continent to which the country belongs to and the file
%   name of the damage function
%
%   next call: isimip...
% CALLING SEQUENCE:
%   [continent,jrc_file]=continent_jrc_damfun(country, params)
% EXAMPLE:
%   clear params;params.filename_suffix='PAA_1';
%   [continent,damfun_file]=continent_jrc_damfun('DEU',params);
% INPUTS:
%   country: country name (full name or ISO3)
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%    filename_suffix: a suffix to the JRC filename (e.g. for different
%    PAA).
%    filepath: the path where the damage function file is located.
% OUTPUTS:
%   continent: the name of the continent to which the country belongs
%   jrc_file: the name of the jrc file containing the respective JRC damage
%       function.
%
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180316, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180316, minor bug fixes
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('country','var'),   country   = '';end
if ~exist('params','var'), params = struct;end

if ischar(params) % special case, where we pass only the name suffix in params
    filename_suffix=params;clear params;params = struct;
    params.filename_suffix=filename_suffix;clear filename_suffix;
end

% check for some parameter fields we need
if ~isfield(params,'filename_suffix'),      params.filename_suffix='';end
if ~isempty(params.filename_suffix)
    if ~strcmp(params.filename_suffix(1),'_'),params.filename_suffix=['_' params.filename_suffix];end
end
if ~isfield(params,'filepath'),             params.filepath='';end


% get country and ISO3 names
[country, country_iso3] =  climada_country_name(country);

shapes = climada_shaperead(climada_global.map_border_file);
shapes_country = string({shapes(:).NAME});
ish = find(shapes_country==country);
if isempty(ish)
    fprintf('*** ERROR: no country match, continent cannot be identified %s\n\n',country_iso3);
    return
elseif length(ish)>1
    fprintf('*** ERROR: more than one country match, continent cannot be identified %s\n\n',country_iso3);
    return
end
continent = shapes(ish).CONTINENT;
% Hack continent:
% 1) Only US and Canada belong to North America, central
% america belong to South America in JRC damage functions
if strcmp(continent, 'North America')
    if ~(strcmp(country, 'United States') || strcmp(country, 'Canada'))
        continent = 'South America';
    end
elseif strcmp(continent, 'Antarctica')
    continent = 'Oceania';
elseif strcmp(continent, 'Seven seas (open ocean)')
    if ismember(string(country), ['French Southern and Antarctic Lands', 'Mauritius', 'Saint Helena', 'Seychelles'])
        continent = 'Africa';
    elseif ismember(country, ['Clipperton Island', 'British Indian Ocean Territory', 'South Georgia and South Sandwich Islands'])
        continent = 'Europe';
    elseif ismember(country, ['Heard Island and McDonald Islands'])
        continent = 'Oceania';
    elseif ismember(country, ['Maldives'])
        continent = 'Asia';
    else
        fprintf('*** ERROR: continent not found %s\n\n',country);
    end
end

continent(isspace(continent)) = [];
jrc_file=[continent '_FL_JRCdamagefunc_residential' params.filename_suffix '.xls'];
if ~strcmp(params.filepath,'')
    jrc_file=[params.filepath filesep jrc_file];
end

end