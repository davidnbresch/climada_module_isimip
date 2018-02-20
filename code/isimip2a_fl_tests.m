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
entity_file='VNM_0150as_entity'

% -----------------
% ASSET LOADING
% -----------------
% load (or, if necessary, construct) the asset base
entity=climada_entity_load(entity_file); % try to load
if isempty(entity)
    fprintf('*** NOTE: generating asset base %s\n\n',entity_file);
    nightlight_params.entity_filename=entity_file; % pass on
    entity=climada_nightlight_entity(country_ISO3,admin1_name,nightlight_params); % create assets from nightlight
%     entity.assets.Value=entity.assets.Value*893189e6*5; % scale to GDP*income_group_factor
    entity.assets.Cover=entity.assets.Value; % technical step, ignore
    save(entity.assets.filename,'entity');
end
figure;climada_entity_plot(entity) % the assets plot
figure;climada_damagefunctions_plot(entity,'FL 001'); % the damage function plot

% based on the entity, decide whether we prefer the range -180..180 or
% 0..360 for the tracks, currently always -180..180
maxlon=180;
% decide on hemisphere
hemisphere='N';if max(entity.assets.lat)<0,hemisphere='S';end


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
hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity,isimip_simround,years_range);


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
    hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity,1,isimip_simround,years_range,'nearest');
end

% calculate the from ground up damage for each event at each centroid
EDS_FL=climada_EDS_calc(entity,hazard_FL); % the damage calculation
figure;climada_EDS_DFC(EDS_FL); % plot flood DFC

fprintf('Note that we set climada_global.damage_at_centroid (back to) %i\n',initial_damage_at_centroid);
climada_global.damage_at_centroid=initial_damage_at_centroid; % reset
