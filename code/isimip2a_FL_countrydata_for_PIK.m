function [status,output_filename]=isimip2a_FL_countrydata_for_PIK(country, ghm, forcing)
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
%   country='Switzerland';
%   ghm='CLM';
%   forcing='gswp3';
%   [status,output_filename]=isimip2a_FL_countrydata_for_PIK(country,ghm,forcing)
% INPUTS:
%   country: country name
%   ghm: Global hydroloical model name
%   forcing: observational forcing name
% OPTIONAL INPUT PARAMETERS:
%   none
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .csv file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180308, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180314, every
%   computation in, still need to allow loading the JRC damage function and
%   save the file
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('country','var'),                country=                 '';end
if ~exist('ghm','var'),                    ghm=                     '';end
if ~exist('forcing','var'),                forcing=                 '';end

% define variables and paths
isimip_simround = '2a';
isimip_data_dir = [climada_global.data_dir filesep 'isimip'];
[country country_iso3] =  climada_country_name(country);

% define entity files for asset source (at resolution 'as0150')
entity_folder = [climada_global.data_dir filesep 'isimip/entities'];
entity_file_isimip=[entity_folder filesep country_iso3 '_0150as_entity_isimip_v0'];

% -----------------
% ASSET LOADING
% -----------------
% load (or, if necessary, construct) the asset base
% ISIMIP assets
entity_isimip=climada_entity_load(entity_file_isimip); % try to load
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

% TO DO: replace damage function with JRC
fprintf('*** ERROR: DAMAGE FUNCTION NOT REPLACED, IMPLEMENT THIS ***');

    
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
        hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity_isimip,1,isimip_simround,years_range,'nearest');
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
        affected_area,  mean_flddph,...
        total_asset_value, total_asset_value_2005, ...
        exposed_asset_value, exposed_asset_value_2005, ...
        damage, damage_2005);
    
fprintf('Note that we set climada_global.damage_at_centroid (back to) %i\n',initial_damage_at_centroid);
climada_global.damage_at_centroid=initial_damage_at_centroid; % reset


end % isimip2a_FL_countrydata_for_PIK



















% -----------------
% ASSET PARAMETERS
% -----------------
% define the assets file (constructed if not existing)
% either create from night-light as follows:
% country_ISO3='USA';country_name='UnitedStates';admin1_name='Florida'; % admin1 default='', whole country, but useful for e.g US states
% %country_ISO3='JPN';country_name='Japan';admin1_name=''; % default='', whole country, but useful for e.g US states
% %entity_file           ='USA_UnitedStates_Florida'; % one can also specify direcly
% if isempty(admin1_name)
%     entity_file           =[country_ISO3 '_' strrep(country_name,' ','')];
% else
%     entity_file           =[country_ISO3 '_' strrep(country_name,' ','') '_' admin1_name];
% end
% %
% nightlight_params.resolution_km=10; % 10x10km resolution
% or specify file directly

% Set of countries, set icountry to choose which one to do
countries = {'Germany', 'France', 'United Kingdom', 'Vietnam'};
icountry = 1;
country = countries{icountry};
[country country_iso3] =  climada_country_name(country);

% define entity files for various asset sources (at resolution 'as0150')
entity_folder = [climada_global.data_dir filesep 'isimip/entities'];
entity_file_isimip=[entity_folder filesep country_iso3 '_0150as_entity_isimip'];
entity_file_nl=[entity_folder filesep country_iso3 '_0150as_entity_nightlight'];
entity_file_bm=[entity_folder filesep country_iso3 '_0150as_entity_blackmarble'];


% -----------------
% ASSET LOADING
% -----------------
% load (or, if necessary, construct) the asset base
% ISIMIP assets
entity_isimip=climada_entity_load(entity_file_isimip); % try to load
if isempty(entity_isimip)
    clear params;params.val_filename='0150as';
    [entity_isimip,params]=isimip_gdp_entity(country_iso3,params);
    % fix the grid rounding errors
    entity_isimip=fix_coords_isimip(entity_isimip, '0150as');
    entity_isimip.assets.centroid_index = 1:length(entity_isimip.assets.centroid_index);
    entity = entity_isimip;
    save(entity_file_isimip,'entity');
    clear entity;
end
    
% compute area of centroids
centroids_area = zeros(1,length(entity_isimip.assets.lon));
ellipsoid = referenceEllipsoid('wgs84','kilometers');
for i=1:length(entity_isimip.assets.lon)
    centroids_area(i) = areaquad(entity_isimip.assets.lat(i)-dlon/2,...
        entity_isimip.assets.lon(i)-dlon/2,entity_isimip.assets.lat(i)+dlon/2,...
        entity_isimip.assets.lon(i)+dlon/2,ellipsoid);
end
country_area = sum(centroids_area);

% figure;climada_entity_plot(entity) % the assets plot
% figure;climada_damagefunctions_plot(entity,'FL 001'); % the damage function plot


%
% -----------------
% HAZARD PARAMETERS
% isimip data files
% -----------------
% for flood (FL): define isimip simulations round (in folder ../climada_data/isimip)
isimip_simround = '2a';
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
ghm = ghms{1};
forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};
forcing = forcings{1};
protections = {'0', '100', 'flopros'};
protection = protections{3};
time_period = 'historical'; % not needed for ISIMIP-2a
years_range = [1971 2010]; % if you need to do only for a subset of years - isimip-2a data contain 1971-2010
% flood_filename        ='USA_UnitedStates_Florida_FL.nc';
% isimip_simround = '2b';
% flood_filename        ='merged_LPJmL_miroc5_historical_flopros_gev_0.1.nc';
% years_range = [0 0];
%

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
% Damage computation for floods
hazard_FL=climada_hazard_load(hazard_FL_file);
if isempty(hazard_FL)
    fprintf('*** NOTE: generating FL hazard from %s\n\n',flood_filename);
    figure % new figure for the check_plot of isimip_flood_load
    hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity_isimip,1,isimip_simround,years_range,'nearest');
end
% 
% centroids = hazard_FL;
% country_names={'Vietnam'};
% 
% entity2=climada_centroids_generate_blackmarble_entity(centroids, country_names)
% calculate the from ground up damage for each event at each centroid

EDS_FL=climada_EDS_calc(entity_isimip,hazard_FL); % the damage calculation

EDS_FL2=climada_EDS_calc(entity_bm,hazard_FL); % the damage calculation

EDS_all=[EDS_FL EDS_FL2];

figure;climada_EDS_DFC(EDS_all);

figure;climada_EDS_DFC({EDS_FL, EDS_FL2}); % plot flood DFC
    hold on;
climada_EDS_DFC(EDS_FL2);

fprintf('Note that we set climada_global.damage_at_centroid (back to) %i\n',initial_damage_at_centroid);
climada_global.damage_at_centroid=initial_damage_at_centroid; % reset
