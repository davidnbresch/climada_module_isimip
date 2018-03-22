function output=isimip2a_FL_countrydata_for_PIK(country, ghm, forcing, params, params_damfun)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip2a_FL_countrydata_for_PIK
% PURPOSE:
%   for a given country and ISIMIP-a2 model run, create a csv file
%   containing information on flood, asset (total value and exposed value),
%   and damage using the JRC impact functions but with PAA=1. Files for all
%   countries can later on be merged together.
%
% CALLING SEQUENCE:
%   [status,output_filename]=isimip2a_FL_countrydata_for_PIK(country,ghm,forcing)
% EXAMPLE:
%   country='Peru';
%   ghm='CLM';
%   forcing='gswp3';
%   clear params;
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   params_damfun.filename_suffix='PAA1';
%   params_damfun.filepath=[climada_global.data_dir filesep 'isimip/entities/damfun'];
%   [status,output_filename]=isimip2a_FL_countrydata_for_PIK(country,ghm,forcing,params,params_damfun)
% INPUTS:
%   country: country name (full name or ISO3)
%   ghm: Global hydroloical model name
%   forcing: observational forcing name
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%    entity_folder: the folder where the entities are located (default:
%       [climada_global.data_dir filesep 'isimip/entities'] ).
%    entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%       entity_prefix='Try1' will result in Try1_DEU_0150as.mat
%   params_damfun: parameters passed to 'continent_jrc_damfun': a structure
%       with fields:
%    filename_suffix: the suffix to append to the damage function
%       file name (see 'continent_jrc_damfun'). Default is 'PAA1'.
%    filepath: the path which contains the damage function file.
%       Default is [climada_global.data_dir filesep
%       'isimip/entities/damfun'.
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .csv file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180308, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180314, every
%   computation in, still need to allow loading the JRC damage function and
%   save the file
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180316, loading the JRC
%   damage function, still needs to save the file
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('country','var'),                 country=                 '';end
if ~exist('ghm','var'),                     ghm=                     '';end
if ~exist('forcing','var'),                 forcing=                 '';end
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_damfun','var'),           params_damfun=       struct;end

% check for some parameter fields we need
if ~isfield(params,'entity_folder'),    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];end
if ~isfield(params,'entity_prefix'),     params.entity_prefix='';end
if ~isfield(params_damfun,'filename_suffix'),params_damfun.filename_suffix='PAA1';end
if ~isfield(params_damfun,'filepath'),     params_damfun.filepath='';end

if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end

% define variables and paths
isimip_simround = '2a';
isimip_data_dir = [climada_global.data_dir filesep 'isimip'];
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
    clear params;params.val_filename='0150as';
    [entity_isimip,params]=isimip_gdp_entity(country_iso3,params);
    % fix the grid rounding errors
    entity_isimip=fix_coords_isimip(entity_isimip, '0150as');
    entity_isimip.assets.centroid_index = 1:length(entity_isimip.assets.centroid_index);
    entity = entity_isimip;
    save(entity_file_isimip,'entity');
    clear entity;
end
entity_isimip.assets.centroid_index = 1:length(entity_isimip.assets.centroid_index);

% Replace damage function with JRC
[continent,damfun_file]=continent_jrc_damfun(country, params_damfun);
[~,entity_isimip]=climada_damagefunctions_read(damfun_file,entity_isimip);
fprintf('* damage function from continent %s is used\n\n',continent);
% fprintf('*** ERROR: DAMAGE FUNCTION NOT REPLACED, IMPLEMENT THIS ***');
% damfun_file = sprintf('%s_%s',cont,'entity_residential.xls');               % Get damagefunction
%   % to first read an entity, then replace its damagefunctions
%   entity=climada_entity_read('entity_template','TCNA_today_small',1);
%   [~,entity]=climada_damagefunctions_read('entity_template',entity);
% Read the JRC damage function
% entity.damagefunctions = climada_damagefunctions_read(damfun_file,entity);
% Adjust damage function parameters
% entity.damagefunctions.PAA = entity.damagefunctions.PAA*sc;

    
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


% -----------------
% HAZARD PARAMETERS
% isimip data files
% -----------------

% generic parameters
protection_levels = {'0', '100', 'flopros'};
time_period = 'historical'; % not needed for ISIMIP-2a
years_range = [0 0]; % all years available for each simulation

% loop on protection levels
for i=1:length(protection_levels)
    protection = protection_levels{i};
    
    % define the climada hazard set files (to be generated below)
    [flddph_filename,fldfrc_filename,fld_path] = isimip_get_flood_filename(isimip_simround, ghm, forcing, protection, time_period);
    flood_filename=[fld_path filesep flddph_filename];
    hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity_isimip,isimip_simround,years_range);

    % switch to FULL RESOLUTION OUTPUT, i.e. event damage at each centroid
    initial_damage_at_centroid=climada_global.damage_at_centroid; % the one used globally
    climada_global.damage_at_centroid=1; % pass on, reset at the end
    
    % -----------------
    % HAZARD LOADING AND CONVERTING TO ASSET CENTROIDS
    % -----------------
    hazard_FL=climada_hazard_load(hazard_FL_file);
    if isempty(hazard_FL)
        fprintf('*** NOTE: generating FL hazard from %s\n\n',flood_filename);
        figure % new figure for the check_plot of isimip_flood_load
        hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity_isimip,0,isimip_simround,years_range,'nearest');
    end
    hazard_FL.yyyy = double(string(hazard_FL.yyyy));

    % -----------------
    % COMPUTATIONS
    % -----------------
    
    % -----------------
    
    % to do for each year available in hazard
    all_years = hazard_FL.yyyy;
    % which indice of the entity corresponds to that year
    ind_year_entity = find(entity_isimip.assets.Values_yyyy == all_years(1)):find(entity_isimip.assets.Values_yyyy == all_years(end));
    ind_year_2005 = find(entity_isimip.assets.Values_yyyy==2005);
    
    
    % First, using entity 2005 (is the default, but replace to be sure)
    entity_isimip.assets.Value = entity_isimip.assets.Values(ind_year_2005,:);
    
    % prepare matrices if first in protection loop
    if i == 1
        % arrays that do not depend on protection level
        asset_value_dim = [length(all_years) 1];
        total_asset_value_2005 = sum(entity_isimip.assets.Value)*ones(asset_value_dim)*entity_isimip.assets.currency_unit;
        total_asset_value = zeros(size(total_asset_value_2005));
        % arrays that depend on protection level
        per_protection_dim = [length(all_years) length(protection_levels)];
        exposed_asset_value_2005 = zeros(per_protection_dim);
        exposed_asset_value = zeros(per_protection_dim);
        mean_flddph = zeros(per_protection_dim);
        affected_area = zeros(per_protection_dim);
        damage_2005 = zeros(per_protection_dim);
        damage = zeros(per_protection_dim);
    end
    
    
    % 1) severity indices (hazard only)
    % total affected area
    affected_area(:,i) = hazard_FL.fraction*centroids_area';
    % mean flood depth where flooded: done per year in the loop below -
    % ignore the next few lines (kept just in case)
    %     % area affected per grid cell (also sums up to affected_area)
    %     centroids_area_relative = centroids_area/sum(centroids_area);
    %     hazard_fraction_area_per_gridcell = times(hazard_FL.fraction, repmat(centroids_area, [40 1]));
    %     % flood height average weighted by area affected
    %     mean_flddph = times(hazard_area_per_gridcell, hazard_FL.intensity)/mean(hazard_area_per_gridcell, 1);
    %     % compare to simple (no lat correction) average
    %     temp = times(hazard_FL.intensity,hazard_FL.fraction);
    %     mean_flddph2=mean(temp(temp>0),1);
    %     frac_intens_prod = times(hazard_FL.fraction,hazard_FL.intensity);
    %     mean_flddph = (times(hazard_FL.fraction,hazard_FL.intensity)*repmat(centroids_area', []))/;

    
    % total and exposed asset value
%     total_asset_value_2005 = sum(entity_isimip.assets.Value)*ones(length(all_years),1)*entity_isimip.assets.currency_unit;
%     total_asset_value = zeros(size(total_asset_value_2005));
%     exposed_asset_value_2005 = zeros(size(total_asset_value_2005)); exposed_asset_value = zeros(size(total_asset_value_2005));
%     mean_flddph = zeros(size(total_asset_value_2005));
    for j=1:length(all_years)
        % compute average flood height where it is flooded
        frac_non_0 = find(hazard_FL.fraction(j,:)~=0);
        mean_flddph(j,i) = (times(hazard_FL.fraction(j,frac_non_0), ...
            centroids_area(frac_non_0))*hazard_FL.intensity(j,frac_non_0)')...
            /(hazard_FL.fraction(j,frac_non_0)*centroids_area(frac_non_0)');
        exposed_asset_value_2005(j,i) = sum(entity_isimip.assets.Value*hazard_FL.fraction(j,:)')*entity_isimip.assets.currency_unit;
        exposed_asset_value(j,i) = sum(entity_isimip.assets.Values(ind_year_entity(j),:)*hazard_FL.fraction(j,:)')*entity_isimip.assets.currency_unit;
        if i==1
            total_asset_value(j,:) = sum(entity_isimip.assets.Values(ind_year_entity(j),:))*entity_isimip.assets.currency_unit;
        end
    end
    
    EDS_FL_2005=climada_EDS_calc(entity_isimip,hazard_FL); % the damage calculation
    damage_2005(:,i) = EDS_FL_2005.damage*EDS_FL_2005.currency_unit;
    % for damage, look into isimip_YDS_calc?
    YDS_FL=isimip_YDS_calc(entity_isimip,hazard_FL);
    damage(:,i) = YDS_FL.damage*EDS_FL_2005.currency_unit;
    
end

% create final matrix

output = cat(2, all_years, repmat(string(country_iso3), [length(all_years) 1]),...
    repmat(string(continent), [length(all_years) 1]),...
    affected_area,  mean_flddph,...
    total_asset_value, total_asset_value_2005, ...
    exposed_asset_value, exposed_asset_value_2005, ...
    damage, damage_2005);

output_names = [string('year') string('country') string('continent')...
    string(['affected_area_',protection_levels{1}])...
    string(['affected_area_',protection_levels{2}])...
    string(['affected_area_',protection_levels{3}])...
    string(['mean_flddph_',protection_levels{1}])...
    string(['mean_flddph_',protection_levels{2}])...
    string(['mean_flddph_',protection_levels{3}])...
    string('total_asset_value')...
    string('total_asset_value_2005')...
    string(['exposed_asset_value',protection_levels{1}])...
    string(['exposed_asset_value',protection_levels{2}])...
    string(['exposed_asset_value',protection_levels{3}])...
    string(['exposed_asset_value_2005',protection_levels{1}])...
    string(['exposed_asset_value_2005',protection_levels{2}])...
    string(['exposed_asset_value_2005',protection_levels{3}])...
    string(['damage',protection_levels{1}])...
    string(['damage',protection_levels{2}])...
    string(['damage',protection_levels{3}])...`
    string(['damage_2005',protection_levels{1}])...
    string(['damage_2005',protection_levels{2}])...
    string(['damage_2005',protection_levels{3}])];

output = cat(1, output_names, output);
    
fprintf('Note that we set climada_global.damage_at_centroid (back to) %i\n',initial_damage_at_centroid);
climada_global.damage_at_centroid=initial_damage_at_centroid; % reset


end % isimip2a_FL_countrydata_for_PIK
