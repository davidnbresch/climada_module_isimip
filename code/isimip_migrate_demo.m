% isimip_migrate_demo
% isimip_migrate_demo
% MODULE:
%   isimip
% NAME:
%   isimip_migrate_demo
% PURPOSE:
%   batch job to check migration in Bangladesh
%
%   FIRST, set PARAMETERS in code below
%
% CALLING SEQUENCE:
%   isimip_migrate_demo
% EXAMPLE:
%   isimip_migrate_demo
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   entity (assets) and hazard sets to respective folders
%   to stdout and figures stored as .png files
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171013, initial
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% PARAMETERS
%
% define the country (ISO3 code)
country_ISO3='BGD'; % e.g. 'BGD', see climada_country_name
%
% define single events to provide detail for
% Sidr: 2007314N10093, damage 3.775 USD bn
% Aila: 2009143N17089, damage 0.270 USD bn
single_event_ID_no={2007314010093,2009143017089}; % N-->0, S-->1
%
prob_str='prob'; % default ='' for historic only
%
% LATER, start from global ismip entity, see isimip_gdp_entity
%GLB_entity_file=[climada_global.hazards_dir filesep 'GLB_0360as_entity.mat'];
% in this case, revise also hazard generation (might not be needed etc.)
%
% define the country asset base
country_entity_file=[climada_global.entities_dir filesep country_ISO3 '_0360as_entity.mat'];
%
% define the tc tracks
%tc_track_file=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat'];
tc_track_file=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs_' prob_str '.mat'];
hazard_set_file_TC=[climada_global.hazards_dir filesep country_ISO3 '_isimip_TC' prob_str '.mat'];
hazard_set_file_TS=[climada_global.hazards_dir filesep country_ISO3 '_isimip_TS' prob_str '.mat'];
%
% local folder to write the figures
fig_dir =[climada_global.results_dir filesep 'isimip_migrate'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
%
% heavy switch to force recaulation of hazard sets
recalc_all=0;


% get the assets
% --------------
if ~exist(country_entity_file,'file')
    entity=isimip_gdp_entity(country_ISO3); % single country entity
    % SPECIAL to reset centroid_index (which points to global isimip centroids
    entity.assets.centroid_index=1:length(entity.assets.lon);
    save(country_entity_file,'entity',climada_global.save_file_version);
else
    entity=climada_entity_load(country_entity_file);
end


% get the TC (tropical cyclone) tracks
% ------------------------------------
if ~exist(tc_track_file,'file')
    [tc_track,tc_track_file]=isimip_ibtracs_read('all','',1,1); % all tracks globally, quality checked
else
    tc_track=climada_tc_track_load(tc_track_file);
end


% generate the historic TC (tropical cyclone wind) event set
% ----------------------------------------------------------
if ~exist(hazard_set_file_TC,'file') || recalc_all
    hazard_TC=isimip_tc_hazard_set(tc_track,hazard_set_file_TC,entity);
else
    hazard_TC=climada_hazard_load(hazard_set_file_TC);
end


% generate the historic TS (tropical cyclone surge) event set
% -----------------------------------------------------------
if ~exist(hazard_set_file_TS,'file') || recalc_all
    hazard_TS=climada_ts_hazard_set(hazard_TC,hazard_set_file_TS);
else
    hazard_TS=climada_hazard_load(hazard_set_file_TS);
end


% plot maximum hazard intensity
figure('Name','Maximum intensities');
subplot(2,1,1);climada_hazard_plot(hazard_TC,0);
subplot(2,1,2);climada_hazard_plot(hazard_TS,0);
        

% calculate damages
% -----------------
EDS   =climada_EDS_calc(entity,hazard_TC,'TC');
EDS(2)=climada_EDS_calc(entity,hazard_TS,'TS');
EDS(3)=climada_EDS_combine(EDS);

    
% damage exceedance plot
% ----------------------
figure;[~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS);

% add EM-DAT for comparison
em_data=emdat_read('',country_ISO3,'-TC',1,1); % last parameter =1 for verbose
if isfield(entity.assets,'currency_unit')
    em_data.damage      = em_data.damage/entity.assets.currency_unit;
    em_data.damage_orig = em_data.damage_orig/entity.assets.currency_unit;
    em_data.DFC.damage  = em_data.DFC.damage/entity.assets.currency_unit;
    em_data.DFC.ED      = em_data.DFC.ED/entity.assets.currency_unit;
    em_data.DFC_orig.damage = em_data.DFC_orig.damage/entity.assets.currency_unit;
    em_data.DFC_orig.ED     = em_data.DFC_orig.ED/entity.assets.currency_unit;
    em_data.YDS.damage  = em_data.YDS.damage/entity.assets.currency_unit;
end
[legend_str,legend_handle]=emdat_barplot(em_data,'dm','om','EM-DAT',legend_str,legend_handle,'SouthEast');

saveas(gcf,[fig_dir filesep country_ISO3 '_DFC' prob_str '.' fig_ext],fig_ext);


% single events
% -------------
for event_i=1:length(single_event_ID_no)
    ID_no=single_event_ID_no{event_i};
    event_index=find(hazard_TC.ID_no==ID_no);
    if ~isempty(event_index)
        event_index=event_index(1); % take first one (original)
        fprintf('event ID_no %i, index %i, TC: %2.2g, TS: %2.2g, combined: %2.2g\n',...
            ID_no,event_index,...
            EDS(1).damage(event_index),EDS(2).damage(event_index),EDS(3).damage(event_index))
        figure('Name',num2str(ID_no));
        subplot(2,1,1);climada_hazard_plot(hazard_TC,event_index);
        subplot(2,1,2);climada_hazard_plot(hazard_TS,event_index);
        saveas(gcf,[fig_dir filesep country_ISO3 '_event_' num2str(ID_no) fig_ext],fig_ext);
    end % found
end % event_i
