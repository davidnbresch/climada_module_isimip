%function BGD_TS_hist_kerry_batch
% BGD_TS_hist_kerry_batch
% MODULE:
%   isimip
% NAME:
%   BGD_TS_hist_kerry_batch
% PURPOSE:
%   generate entities and TC as well as TS hazard sets for Bangladesh (or
%   in fact any other country), both in 0360as and 0150as (arcsec) resolution. 
%
%   Highly recommended to set climada_global.parfor=1 for speedup
%
%   First, TEST with 0360as, as already this needs to map  53'913'600 points
%   to 2'303 points for SRTM in Bangladesh, which takes abut 6 minutes. Please
%   also conside to set FAST_TEST (see PARAMETERS) in order to run only a
%   subset of all tracks for check first.
%
%   We use the abbreviation ETOP for ETOPO1, such that it has same length
%   as SRTM (nicer for filenames)
%
%   previous call: none
%   next call: climada_EDS_calc
% CALLING SEQUENCE:
%   BGD_TS_hist_kerry_batch
% EXAMPLE:
%   BGD_TS_hist_kerry_batch
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   batch job, stdout
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20180104, initial
% David N. Bresch, david.bresch@gmail.com, 20180105, figures also inivisbly produced
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(fileparts(mfilename('fullpath')))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% define the country (see climada_country_name)
country_ISO3='BGD';
country_ISO3='BRB';
%
% define entity resolution and spacing of regular ocean grid in degrees (and markersize to plot):
% in a first got, set markersize=[] for it to be automatically determined, you might later adjust
entity_params.val_filename='0150as';dlonlat=0.1;markersize=[];
%entity_params.val_filename='0360as';dlonlat=0.2;markersize=3;
%
% FAST_TEST to only process a FAST_TEST number of tracks
% =0: default, all tracks (and preserving event_ID as track_i etc. (recommended)
% =100: only first 100 tracks hitting at least one centroid
% =1e9: to select ALL tracks hitting at least one centroid
FAST_TEST           =  0; % =100 for only first 100 tracks to check, default=0
%
% define the check plots
centroids_check_plot=  1;
ETOP_check_plots    = 10; % =0 to suppress plots, =10 for plots (to show ETOPO1 height in range 0..10m)
SRTM_check_plots    = 10; % =0 to suppress plots, =10 for plots (to show SRTM height in range 0..10m)
Kerry_check_plots   =  1; % =0 to suppress plots,  =1 for plots
DFC_check_plots     =  1; % =0 to suppress plots,  =1 for plots
fig_save            =  1; % whether we save the figures (=1, default) or not (=0)
%
% whether we show the figures on screen (='on') or not (='off'). In the
% case of 'off', some figures (ETOP and SRTM check plots) are not generated
fig_visible         ='off'; % default='on', ='off' for processing w/o figures
%
% more technical PARAMETERS follow:
% to improve temporal resolution, consider climada_global.tc.default_min_TimeStep
%climada_global.tc.default_min_TimeStep=.1;
%
% whether we shift the regular grid relative to integer lat/lons
dlonlat_offset=0; % default=0 (in degrees)
%
% the file with all IBTrACS alread read (see isimip_ibtracs_read, automatically invoked below)
ibtracs_save_file=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat'];
%
% the list of TC Kerry's track files to be processed
track_files={
    'Trial1_GB_dkgfdl_rcp60cal'
    'Trial1_GB_dkhad_rcp60cal'
    'Trial1_GB_dkipsl_rcp60cal'
    'Trial1_GB_dkmiroc_rcp60cal'
    };
%
fig_dir =[climada_global.results_dir filesep 'isimip_' country_ISO3];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';

if fig_save,fprintf('> saving all figures to %s\n\n',fig_dir);end

% define a unique name (used to name hazard event sets
unique_name=[country_ISO3 '_' entity_params.val_filename];

% suppress check_plots (since generated in functions) if figures not visible
if strcmpi(fig_visible,'off'),ETOP_check_plots=0;SRTM_check_plots=0;end
    
clear centroids EDS % aviod troubles with preceding calls (since a batch job)

t0=clock;

% prepare the centroids
% =====================

entity_file=[climada_global.entities_dir filesep unique_name '_entity.mat'];
if exist(entity_file,'file')
    % load isimip entity
    fprintf('< loading entity %s\n',entity_file);
    entity=climada_entity_load(entity_file);
else
    % generate isimip entity
    entity=isimip_gdp_entity(country_ISO3,entity_params);
end

% add coarse grid
entity.BoundingBox(1,1)=floor(min(entity.assets.lon)); % shape-file convention
entity.BoundingBox(2,1)= ceil(max(entity.assets.lon));
entity.BoundingBox(1,2)=floor(min(entity.assets.lat));
entity.BoundingBox(2,2)= ceil(max(entity.assets.lat));

centroids.lon=entity.assets.lon;centroids.lat=entity.assets.lat; % copy
n_centroids_entity=length(centroids.lon);
fprintf('> adding coarse (1x1deg) ocean grid\n')
for lon_i=entity.BoundingBox(1,1):dlonlat:entity.BoundingBox(2,1) % add regular coarse grid
    for lat_i=entity.BoundingBox(1,2):dlonlat:entity.BoundingBox(2,2)
        centroids.lon(end+1)=lon_i+dlonlat_offset; % to keep this points unique
        centroids.lat(end+1)=lat_i+dlonlat_offset;
    end
end
centroids.centroid_ID=1:length(centroids.lon);
% set the ocean-point centroids apart:
centroids.centroid_ID(n_centroids_entity+1:end)=centroids.centroid_ID(n_centroids_entity+1:end)+(3e6-n_centroids_entity);
centroids.on_land=centroids.centroid_ID*0+1;centroids.on_land(n_centroids_entity+1:end)=0; % mask ocean points
centroids.on_land=logical(centroids.on_land);
centroids.distance2coast_km = climada_distance2coast_km(centroids.lon,centroids.lat); % add distance to coast in km
centroids.elevation_m       = etopo_elevation_m(        centroids.lon,centroids.lat); % add ETOPO1 elevation
open_water_points           = centroids.elevation_m<0 & ~centroids.on_land; % zero or below elevation and not on land
centroids.distance2coast_km(open_water_points) = -centroids.distance2coast_km(open_water_points); % open water points neg. distance

if centroids_check_plot
    fig_centroids=figure('Name','centroids','Visible',fig_visible);
    climada_entity_plot(entity,7);
    hold on
    plot(centroids.lon,centroids.lat,'xr','MarkerSize',2*3); % larger marker, since overlays
    plot(centroids.lon(open_water_points),centroids.lat(open_water_points),'ob','MarkerSize',3);
    plot(centroids.lon(centroids.on_land),centroids.lat(centroids.on_land),'.g','MarkerSize',1);
    if fig_save,saveas(fig_centroids,[fig_dir filesep unique_name '_centroids'],fig_ext);delete(fig_centroids);end % save and delete figure
end


% calculate tropical cyclone wind and surge hazard set
% ====================================================

% prepare IBTrACS
% ---------------
if exist(ibtracs_save_file,'file')
    load(ibtracs_save_file) % contains tc_track
else
    tc_track=isimip_ibtracs_read('all','',1,1);
end
if FAST_TEST
    % select tracks hitting centroids
    tc_track=climada_tc_track_info(tc_track,0,[],centroids);
    if length(tc_track)>FAST_TEST,tc_track=tc_track(1:FAST_TEST);end
end

% calculations for IBTrACS
% ------------------------

% generate the TC windfield
hazard_TC      = climada_tc_hazard_set(tc_track,[unique_name '_TC_hist'],centroids);
fig_TC=figure('Name','TC','Visible',fig_visible);
climada_hazard_plot(hazard_TC,0,markersize);title([strrep(unique_name,'_',' ') ' TC maxintens']);
if fig_save,saveas(fig_TC,[fig_dir filesep unique_name '_TC'],fig_ext);delete(fig_TC);end % save and delete figure

    
% generate the simply coarse-resolution (ETOPO)  TS surge field
[hazard_TS_ETOP,ETOP_elevation] = climada_ts_hazard_set(hazard_TC,[unique_name '_TS_ETOP'],'ETOPO',ETOP_check_plots,0);
if ETOP_check_plots
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_hist'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_maxintens'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_elevation'],fig_ext);delete(gcf);end % save and delete figure
    if strcmpi(get(gcf,'Name'),'etopo')
        if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_elevation_all'],fig_ext);delete(gcf);end % save and delete figure
    else
        close all
    end
    climada_hazard_plot(hazard_TS_ETOP,0,markersize);
    if fig_save,saveas(gcf,     [fig_dir filesep unique_name '_TS_ETOP_maxintens'],fig_ext);delete(gcf);end % save and delete figure
else
    fig_ETOP=figure('Name','ETOP','Visible',fig_visible);
    climada_hazard_plot(hazard_TS_ETOP,0,markersize);title([strrep(unique_name,'_',' ') ' TS ETOP maxintens']);
    if fig_save,saveas(fig_ETOP,[fig_dir filesep unique_name '_TS_ETOP_maxintens'],fig_ext);delete(fig_ETOP);end % save and delete figure
end

% generate the high-resolution (SRTM) TS surge field
% to show the SRTM mapping, set 2nd last parameter to -SRTM_check_plots (plot takes a long time)
[hazard_TS_SRTM,SRTM_elevation] = climada_ts_hazard_set(hazard_TC,[unique_name '_TS_SRTM'],'SRTM',SRTM_check_plots,1);
if SRTM_check_plots
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_hist'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_maxintens'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_elevation'],fig_ext);delete(gcf);end % save and delete figure
    climada_hazard_plot(hazard_TS_SRTM,0,markersize);
    if fig_save,saveas(gcf,     [fig_dir filesep unique_name '_TS_SRTM_maxintens'],fig_ext);delete(gcf);end % save and delete figure
else
    fig_SRTM=figure('Name','SRTM','Visible',fig_visible);
    climada_hazard_plot(hazard_TS_SRTM,0,markersize);title([strrep(unique_name,'_',' ') ' TS SRTM maxintens']);
    if fig_save,saveas(fig_SRTM,[fig_dir filesep unique_name '_TS_SRTM_maxintens'],fig_ext);delete(fig_SRTM);end % save and delete figure
end

entity    =climada_assets_encode(entity,hazard_TC); % encode once
EDS       =climada_EDS_calc(entity,hazard_TC,     [unique_name '_TC_hist']);
EDS(end+1)=climada_EDS_calc(entity,hazard_TS_ETOP,[unique_name '_TS_ETOP_hist']);
EDS(end+1)=climada_EDS_calc(entity,hazard_TS_SRTM,[unique_name '_TS_SRTM_hist']);

% calculations for Kerry's simulations
% ------------------------------------

for file_i=1:length(track_files)
    
    tc_track=isimip_tc_track_load(track_files{file_i},'N',180,-1); % Northern hemisphere
    if FAST_TEST
        tc_track=climada_tc_track_info(tc_track,0,[],centroids);
        if length(tc_track)>FAST_TEST,tc_track=tc_track(1:FAST_TEST);end
    end
    hazard_name_TC=[unique_name '_' track_files{file_i} '_TC'];
    hazard_TC_kerry=isimip_tc_hazard_set(tc_track,hazard_name_TC,centroids,1,hazard_name_TC);
    
    hazard_name_ETOP=[unique_name '_' track_files{file_i} '_TS_ETOP'];
    hazard_name_SRTM=[unique_name '_' track_files{file_i} '_TS_SRTM'];
    hazard_TS_ETOP_kerry = climada_ts_hazard_set(hazard_TC_kerry,hazard_name_ETOP,'ETOPO',0,1,ETOP_elevation);
    hazard_TS_SRTM_kerry = climada_ts_hazard_set(hazard_TC_kerry,hazard_name_SRTM,'SRTM', 0,1,SRTM_elevation);
    
    if Kerry_check_plots
        fig_Kerry_ETOP=figure('Name','Kerry ETOP','Visible',fig_visible);
        climada_hazard_plot(hazard_TS_ETOP_kerry,0,markersize);title([strrep(hazard_name_ETOP,'_',' ') ' maxintens']);
        if fig_save,saveas(fig_Kerry_ETOP,[fig_dir filesep hazard_name_ETOP],fig_ext);delete(fig_Kerry_ETOP);end % save and delete figure
        fig_Kerry_SRTM=figure('Name','Kerry SRTM','Visible',fig_visible);
        climada_hazard_plot(hazard_TS_SRTM_kerry,0,markersize);title([strrep(hazard_name_SRTM,'_',' ') ' maxintens']);
        if fig_save,saveas(fig_Kerry_SRTM,[fig_dir filesep hazard_name_SRTM],fig_ext);delete(fig_Kerry_SRTM);end % save and delete figure
    end
    
    EDS(end+1)=climada_EDS_calc(entity,hazard_TC_kerry,     hazard_name_TC);
    EDS(end+1)=climada_EDS_calc(entity,hazard_TS_ETOP_kerry,hazard_name_ETOP);
    EDS(end+1)=climada_EDS_calc(entity,hazard_TS_SRTM_kerry,hazard_name_SRTM);

    % create probabilistic tracks
    %[~,p_rel] = climada_tc_track_wind_decay_calculate(tc_track,0); % wind speed decay at track nodes after landfall
    %tc_track  = climada_tc_random_walk(tc_track); % overwrites tc_track to save memory
    %tc_track  = climada_tc_track_wind_decay(tc_track, p_rel,0); % add the inland decay correction to all probabilistic nodes
    
end % file_i

if DFC_check_plots
    % for a show all EDSs
    fig_DFC=figure('Name','DFC','Visible',fig_visible);
    climada_EDS_DFC(EDS);title(strrep(unique_name,'_',' '))
    if fig_save,saveas(fig_DFC,[fig_dir filesep unique_name '_DFC'],fig_ext);delete(fig_DFC);end % save and delete figure
end

fprintf('all calculations took %f sec.\n',etime(clock,t0));

%end % BGD_TS_hist_kerry_batch