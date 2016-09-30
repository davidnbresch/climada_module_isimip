% isimip_step_by_step
% climada isimip step by step
% MODULE:
%   isimip
% NAME:
%   isimip_step_by_step
% PURPOSE:
%   batch script to run isimip input into climada and core climada
%   calculations. 
%
%   Mainly demonstrates the use of isimip_tc_track_load and isimip_flood_load
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
%   isimip_step_by_step
% EXAMPLE:
%   isimip_step_by_step
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   some figures and all data loaded into MATLAB session (accessible from
%       comand line, obviously)
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160929, initial
% David N. Bresch, david.bresch@gmail.com, 20160929, hazard_dummy_std_file added
% David N. Bresch, david.bresch@gmail.com, 20160930, flood added
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
% define the assets file (constructed if not existing)
entity_file           ='USA_UnitedStates_Florida';
nightlight_params.resolution_km=10; % set 10 km asset resolution (10x10km), other option=1
%
% define default climada tropical cyclone (TC) hazard set
hazard_std_file       ='USA_UnitedStates_atl_TC';
hazard_dummy_std_file ='TCNA_today_small'; % if hazard_std_file does not exist
%
% for tropical cyclones (TC): define isimip input data (in folder ../climada_data/isimip)
tc_track_20th_file    ='temp_mpi20thcal'; % the file with Kerry's tracks
tc_track_rcp85_file   ='temp_mpircp85cal_full'; % the file with Kerry's tracks
%
% for flood (FL): define isimip input data (in folder ../climada_data/isimip)
flood_fraction_filename='fldfrc_max.nc';
flood_depth_filename   ='flddph_max.nc';
%
% switch to FULL RESOLUTION OUTPUT, i.e. event damage at each centroid
climada_global.EDS_at_centroid=1; % climada default=0, set =1 specifically

% define the climada hazard set files (to be generated below)
hazard_20th_file       =[entity_file '_' tc_track_20th_file];
hazard_rcp85_file      =[entity_file '_' tc_track_rcp85_file];
hazard_FL_file         =[entity_file '_FL'];

% load (or, if necessary, construct) the asset base
entity=climada_entity_load(entity_file); % try to load
if isempty(entity)
    fprintf('*** NOTE: generating asset base %s\n\n',entity_file);
    nightlight_params.entity_filename=entity_file; % pass on
    entity=climada_nightlight_entity('USA','Florida',nightlight_params); % create assets from nightlight
    entity.assets.Value=entity.assets.Value*893189e6*5; % scale to GDP*income_group_factor
    entity.assets.Cover=entity.assets.Value; % technical step, ignore
    save(entity.assets.filename,'entity');
end
figure;climada_entity_plot(entity) % the assets plot
figure;climada_damagefunctions_plot(entity,'TC 001'); % the damage function plot

% load (or, if necessary, construct) the hazard set for 20th century
hazard_20th=climada_hazard_load(hazard_20th_file);
if isempty(hazard_20th)
    fprintf('*** NOTE: imparting TC tracks from %s\n\n',tc_track_20th_file);
    % load Kerry's TC tracks and generate the hazard set
    tc_track_20th=isimip_tc_track_load(tc_track_20th_file,0); % load the tracks (the only new code)
    hazard_20th=climada_tc_hazard_set(tc_track_20th,hazard_20th_file,entity); % generate the new hazard set (4 min)
end

% load the climada default hazard set for comparison
hazard_std=climada_hazard_load(hazard_std_file); % load default hazard
if isempty(hazard_std)
    fprintf('Warning: using DUMMY standard hazard\n');
    hazard_std=climada_hazard_load(hazard_dummy_std_file); % load default hazard
end

% calculate the from ground up damage for each event at each centroid
EDS(1)=climada_EDS_calc(entity,hazard_std); % the default damage calculation for comparison (0.5 sec)
EDS(2)=climada_EDS_calc(entity,hazard_20th); % the damage calculation with Kerry?s tracks (0.5 sec)
figure;climada_EDS_DFC(EDS); % plot climada std and 20th century DFCs

% load (or, if necessary, construct) the hazard set for rcp85
hazard_rcp85=climada_hazard_load(hazard_rcp85_file);
if isempty(hazard_rcp85)
    fprintf('*** NOTE: imparting TC tracks from %s\n\n',tc_track_20th_file);
    % load Kerry's TC tracks and generate the hazard set
    tc_track_rcp85=isimip_tc_track_load(tc_track_rcp85_file,0); % load rcp85 tracks, convert
    hazard_rcp85=climada_tc_hazard_set(tc_track_rcp85,hazard_rcp85_file,entity); % generate the new hazard set (4 min)
end
if ~isempty(hazard_rcp85) % as we do not provide this hazard in the TEST suite, hence one might not run this step and just proceed to flood
    % calculate the from ground up damage for each event at each centroid
    EDS(3)=climada_EDS_calc(entity,hazard_rcp85); % the damage calculation with Kerry?s tracks (0.5 sec)
    figure;climada_EDS_DFC(EDS); % plot climada std, 20th century and rcp85 DFCs
end

fprintf('\n*** NOTE: all basic tests for TC done ***\n');

hazard_FL=climada_hazard_load(hazard_FL_file);
if isempty(hazard_FL)
    fprintf('*** NOTE: generating FL hazard from %s and %s\n\n',flood_fraction_filename,flood_depth_filename);
    figure % new figure for the check_plot of isimip_flood_load
    hazard=isimip_flood_load(flood_fraction_filename,flood_depth_filename,'auto',entity,1);
end

% calculate the from ground up damage for each event at each centroid
EDS_FL=climada_EDS_calc(entity,hazard_FL); % the damage calculation with Kerry?s tracks (0.5 sec)
figure;climada_EDS_DFC(EDS_FL); % plot flood DFC