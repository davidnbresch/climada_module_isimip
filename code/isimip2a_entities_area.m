function country_area=isimip2a_entities_area(country, params)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip2a_entities_area
% PURPOSE:
%   for a given country, load isimip entity and determine country area
%
% CALLING SEQUENCE:
%   country_area=isimip2a_FL_countrydata_for_PIK(country,params)
% EXAMPLE:
%   country='Peru';
%   clear params;
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   country_area=isimip2a_FL_countrydata_for_PIK(country,params)
% INPUTS:
%   country: country name (full name or ISO3)
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%    entity_folder: the folder where the entities are located (default:
%       [climada_global.data_dir filesep 'isimip/entities'] ).
%    entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%       entity_prefix='Try1' will result in Try1_DEU_0150as.mat
% OUTPUTS:
%   country_area: country area in km**2
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180827, initial, based
% on isimip2a_FL_countrydata_for_PIK
%   
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('country','var'),                 country=                 '';end
if ~exist('params','var'),                  params=              struct;end

% check for some parameter fields we need
if ~isfield(params,'entity_folder'),    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];end
if ~isfield(params,'entity_prefix'),     params.entity_prefix='';end

if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end

% define variables and paths
[country country_iso3] =  climada_country_name(country);

% define entity files for asset source (at resolution 'as0150')
entity_file_isimip=[params.entity_folder filesep params.entity_prefix strtrim(country_iso3) '_0150as_entity'];

% -----------------
% ASSET LOADING
% -----------------
% load (or, if necessary, construct) the asset base
% ISIMIP assets
entity_isimip=climada_entity_load(entity_file_isimip,1); % try to load, flag to 1 to avoir overwrite
if isempty(entity_isimip)
    fprintf('*** ERROR: entity file not found %s\n\n',entity_file_isimip);
    return
%     clear params;params.val_filename='0150as';
%     [entity_isimip,params]=isimip_gdp_entity(country_iso3,params);
%     % fix the grid rounding errors
%     entity_isimip=fix_coords_isimip(entity_isimip, '0150as');
%     entity_isimip.assets.centroid_index = 1:length(entity_isimip.assets.centroid_index);
%     entity = entity_isimip;
%     save(entity_file_isimip,'entity');
%     clear entity;
end
entity_isimip.assets.centroid_index = 1:length(entity_isimip.assets.centroid_index);
    
% compute area of centroids
centroids_area = zeros(1,length(entity_isimip.assets.lon));
ellipsoid = referenceEllipsoid('wgs84','kilometers');
dlon = 360/8640;
dlat   = 180/4320;
for i=1:length(entity_isimip.assets.lon)
    centroids_area(i) = areaquad(entity_isimip.assets.lat(i)-dlon/2,...
        entity_isimip.assets.lon(i)-dlon/2,entity_isimip.assets.lat(i)+dlon/2,...
        entity_isimip.assets.lon(i)+dlon/2,ellipsoid);
end
country_area = sum(centroids_area);



end % isimip2a_FL_countrydata_for_PIK
