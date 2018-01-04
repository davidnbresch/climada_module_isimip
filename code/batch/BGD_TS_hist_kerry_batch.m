%function BGD_TS_hist_kerry_batch
% BGD_TS_hist_kerry_batch
% MODULE:
%   isimip
% NAME:
%   BGD_TS_hist_kerry_batch
% PURPOSE:
%   generate entities and TC as well as TS hazard sets for Bangladesh
%
%   highly recommended to set climada_global.parfor=1
%
%   first, TEST with 0360as, as already this needs to map  53'913'600 points
%   to 2303 points for SRTM, which takes abut 6 min.
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
country_ISO3='BGD';
%
% define entity resolution and spacing of regular ocean grid (and markersize to plot):
%entity_params.val_filename='0150as';dlonlat=0.1;markersize=4;
entity_params.val_filename='0360as';dlonlat=0.2;markersize=7;
%
% to improve temporal resolution, consider climada_global.tc.default_min_TimeStep
%%climada_global.tc.default_min_TimeStep=.1;
%
% whether we shift the regular grid relative to integer lat/lons
dlonlat_offset=0;
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
% define the check plots
centroids_check_plot=1;
ETOP_check_plots=10; % =0 to suppress plots, =10 for plots (to show ETOPO1 height in range 0..10m)
SRTM_check_plots=10; % =0 to suppress plots, =10 for plots (to show SRTM height in range 0..10m)
%
fig_dir =[climada_global.results_dir filesep 'isimip_BGD'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
fig_save=1;


% define a unique name (used to name hazard event sets
unique_name=[country_ISO3 '_' entity_params.val_filename];

% 50km coastal cut-off

% prepare the centroids
% =====================

entity_file=[climada_global.entities_dir filesep unique_name '_entity.mat'];
if exist(entity_file,'file')
    % load isimip entity
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
centroids.distance2coast_km=climada_distance2coast_km(centroids.lon,centroids.lat);
n_centroids=length(centroids.lon);
fprintf('> adding coarse (1x1deg) ocean grid\n')
for lon_i=entity.BoundingBox(1,1):dlonlat:entity.BoundingBox(2,1) % add regular coarse grid
    for lat_i=entity.BoundingBox(1,2):dlonlat:entity.BoundingBox(2,2)
        centroids.lon(end+1)=lon_i+dlonlat_offset; % to keep this points unique
        centroids.lat(end+1)=lat_i+dlonlat_offset;
    end
end
centroids.centroid_ID=1:length(centroids.lon);
centroids.centroid_ID(n_centroids+1:end)=centroids.centroid_ID(n_centroids+1:end)+(3e6-n_centroids); % to set the ocean-point centroids apart
centroids.on_land=centroids.centroid_ID*0+1;centroids.on_land(n_centroids+1:end)=0; % mask ocean points
centroids.on_land=logical(centroids.on_land);
centroids.distance2coast_km(n_centroids+1:end)=0; % add

if centroids_check_plot
    figure('Name','centroids')
    climada_entity_plot(entity,7);
    hold on
    plot(centroids.lon,centroids.lat,'xr','MarkerSize',0.1);
    %plot(centroids.lon(centroids.on_land),centroids.lat(centroids.on_land),'.x','MarkerSize',0.1);
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_centroids'],fig_ext);delete(gcf);end % save and delete figure
end

% calculate tropical cyclone hazard set
% =====================================

% prepare IBTrACS
% ---------------
if exist(ibtracs_save_file,'file')
    load(ibtracs_save_file) % contains tc_track
else
    tc_track=isimip_ibtracs_read('all','',1,1);
end

% generate the TC windfield
hazard_TC      = climada_tc_hazard_set(tc_track,[unique_name '_TC_hist'],centroids);

% generate the simply coarse-resolution (ETOPO)  TS surge field
hazard_TS_ETOP = climada_ts_hazard_set(hazard_TC,[unique_name '_TS_ETOP'],'ETOPO',ETOP_check_plots);
if ETOP_check_plots
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_hist'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_maxintens'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_elevation'],fig_ext);delete(gcf);end % save and delete figure
    if strcmpi(get(gcf,'Name'),'etopo')
        if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_ETOP_elevation_all'],fig_ext);delete(gcf);end; % save and delete figure
    else
        close all
    end
end

% generate the high-resolution (SRTM) TS surge field
% to show the SRTM mapping, set 2nd last parameter to -SRTM_check_plots (plot takes a long time)
hazard_TS_SRTM = climada_ts_hazard_set(hazard_TC,[unique_name '_TS_SRTM'],'SRTM',SRTM_check_plots,1);
if SRTM_check_plots
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_hist'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_maxintens'],fig_ext);delete(gcf);end % save and delete figure
    if fig_save,saveas(gcf,[fig_dir filesep unique_name '_TS_SRTM_elevation'],fig_ext);delete(gcf);end % save and delete figure
end

%end % BGD_TS_hist_kerry_batch