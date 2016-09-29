% isimip_step_by_step
% climada isimip step by step
% MODULE:
%   isimip
% NAME:
%   isimip_step_by_step
% PURPOSE:
%   batch script to run isimip input into climada and core climada
%   calculations
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
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
% define the assets file (constructed if not existing)
entity_file      ='USA_UnitedStates_Florida';
%
% define default climada tropical cyclone (TC) hazard set
hazard_std_file  ='USA_UnitedStates_atl_TC';
%
hazard_20th_file ='USA_UnitedStates_Florida_temp_mpi20thcal';
hazard_rcp85_file='USA_UnitedStates_Florida_temp_mpircp85cal';
%
% switch to FULL RESOLUTION OUTPUT, i.e. event damage at each centroid
climada_global.EDS_at_centroid=1;


entity=climada_entity_load(entity_file); % try to load
if isempty(entity)
    p.resolution=10; % set 10 km resolution
    entity=climada_nightlight_entity('USA','Florida',p); % create assets from nightlight
    entity.assets.Value=entity.assets.Value*893189e6*5; % scale to GDP*income_group_factor
    entity.assets.Cover=entity.assets.Value; % technical step, ignore
    save(entity.assets.filename,'entity');
end
figure;climada_entity_plot(entity) % the assets plot
figure;climada_damagefunctions_plot(entity,'TC 001'); % the damage function plot

hazard_20th=climada_hazard_load(hazard_20th_file);
if isempty(entity)
    % load Kerry's TC tracks and generate the hazard set
    tc_track_20th=isimip_load_tc_tracks('temp_mpi20thcal',0); % load the tracks (the only new code)
    hazard_20th=climada_tc_hazard_set(tc_track_20th,hazard_20th_file,entity); % generate the new hazard set (4 min)
end

hazard_std=climada_hazard_load(hazard_std_file); % load default hazard
EDS(1)=climada_EDS_calc(entity,hazard_std); % the default damage calculation for comparison (0.5 sec)
EDS(2)=climada_EDS_calc(entity,hazard_20th); % the damage calculation with Kerry?s tracks (0.5 sec)
figure;climada_EDS_DFC(EDS); % plot climada std and 20th century DFCs

hazard_rcp85=climada_hazard_load(hazard_rcp85_file);
if isempty(entity)
    % load Kerry's TC tracks and generate the hazard set
    tc_track_rcp85=isimip_load_tc_tracks('temp_mpircp85cal_full',0); % load rcp85 tracks, convert
    hazard_rcp85=climada_tc_hazard_set(tc_track_rcp85,hazard_rcp85_file,entity); % generate the new hazard set (4 min)
end

EDS(3)=climada_EDS_calc(entity,hazard_rcp85); % calculate rcp85 century
figure;climada_EDS_DFC(EDS); % plot climada std, 20th century and rcp85 DFCs

fprintf('\n*** NOTE: all basic tests for TC done ***\n');


