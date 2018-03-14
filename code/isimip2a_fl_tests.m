% isimip_fl_tests
% climada isimip-2a flood tests
% MODULE:
%   isimip
% NAME:
%   isimip2a_fl_tests
% PURPOSE:
%   batch script to run isimip input into climada and core climada
%   calculations.
%
%   Mainly demonstrates the use of isimip_flood_load
%
%   this batch job reads already existing (intermediate step) .mat files,
%   otherwise creates them (faster on subsequent calls).
%
%   needs modules:
%   climada         https://github.com/davidnbresch/climada
%   country_risk    https://github.com/davidnbresch/climada_module_country_risk
%   isimip          https://github.com/davidnbresch/climada_module_isimip
%
% CALLING SEQUENCE:
%   isimip_fl_tests
% EXAMPLE:
%   isimip_fl_tests
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   some figures and all data loaded into MATLAB session (accessible from
%       comand line, obviously)
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180118, initial whereby
%       the basis of the code was isimip_fl_tests.m, modified to load
%       isimip-2a hydrological model output
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%


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
    
% % nightlight assets (old version) 
% entity_nl=climada_entity_load(entity_file_nl); % try to load
% if isempty(entity)
%     fprintf('*** NOTE: generating asset base %s\n\n',entity_file);
%     nightlight_params.entity_filename=entity_file; % pass on
%     entity=climada_nightlight_entity(country_ISO3,admin1_name,nightlight_params); % create assets from nightlight
% %     entity.assets.Value=entity.assets.Value*893189e6*5; % scale to GDP*income_group_factor
%     entity.assets.Cover=entity.assets.Value; % technical step, ignore
%     save(entity.assets.filename,'entity');
% end

% Black Marble night light assets
entity_bm=climada_entity_load(entity_file_bm); % try to load
if isempty(entity_bm)
    entity_bm=climada_centroids_generate_blackmarble_entity(entity_isimip.assets, {country_iso3});
    entity = entity_bm;
    save(entity_file_bm, 'entity');
    clear entity;
end


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
