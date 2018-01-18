%function BGD_TS_batch(IBTrACS_ID,entity,interest_area,marker_size)
% climada BGD TS storm surge batch
% MODULE:
%   isimip
% NAME:
%   BGD_TS_batch
% PURPOSE:
%   run 2007 SIDR (2007314N10093), create TS surge field etc.
%
%   Set climada_global.parfor=1; % for substabtial speedup (in
%   climada_ts_hazard_set, also in climada_regrid)
%
%   previous call: isimip_ibtracs_read and isimip_ibtracs_load
%   next call: climada_EDS_calc
% CALLING SEQUENCE:
%   BGD_TS_batch(IBTrACS_ID)
% EXAMPLE:
%   BGD_TS_batch('2007314N10093')
% INPUTS:
%   IBTrACS_ID: either an original IBTrACS ID, such as '2007314N10093' or a
%       climada-compatible one as (long) integer,such as 2007314010093 (just
%       S->1, N->0)
%       Usually obtained from
%       http://www.atms.unca.edu/ibtracs/ibtracs_v03r07/browse-ibtracs/browseIbtracs.php 
%       Aila: http://www.atms.unca.edu/ibtracs/ibtracs_v03r07/browse-ibtracs/index.php?name=v03r07-2009143N17089
%   entity: a climada entity, used for centroids, loaded using
%       climada_entity_load, hence either a filename (even w/o path) or a
%       full entity structure.
% OPTIONAL INPUT PARAMETERS:
%   interest_area: [minlon maxlon minlat maxlat] the area to show on plots
%       (to zoom in by default).
%   marker_size: default =2
% OUTPUTS:
%   command window and figures
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20171231, intial
% David N. Bresch, david.bresch@gmail.com, 20180102, save_figures
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('IBTrACS_ID','var'),IBTrACS_ID='2007314N10093';end
if ~exist('entity','var'),entity='BGD_Bangladesh_Barisal_01x01';end
if ~exist('interest_area','var'),interest_area=[];end
if ~exist('marker_size','var'),marker_size=2;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% the file with all IBTrACS alread read (see isimip_ibtracs_read, automatically invoked below)
ibtracs_save_file=[climada_global.data_dir filesep 'tc_tracks' filesep 'ibtracs' filesep 'ibtracs.mat'];
%
%if isempty(interest_area),interest_area=[89.8 90.8 21.8 23];end % all Barisal
if isempty(interest_area),interest_area=[90 90.2 22 22.2];end % SouthWest Barisal
%
fig_dir =[climada_global.results_dir filesep 'isimip_BGD'];
if ~isdir(fig_dir),[fP,fN]=fileparts(fig_dir);mkdir(fP,fN);end % create it
fig_ext ='png';
save_figures=1;
%
% TEST
IBTrACS_ID='2007314N10093';fig_name ='Sidr_BGD_Barisal'; % Sidr, 2007
%IBTrACS_ID='2009143N17089';fig_name ='Aila_BGD_Barisal'; % Aila, 2009, track does not hit Bangladesh...
entity='BGD_Bangladesh_Barisal_01x01';
%interest_area=[90 90.2 22 22.2]; % SouthWest Barisal
interest_area=[90.30 90.60 22.65 22.95]; % with Barisal city in lower left
marker_size=10; % for such a small area
%
% BIG
% IBTrACS_ID='2007314N10093';fig_name ='Sidr_BGD_Bangladesh'; % Sidr, 2007
% entity='BGD_Bangladesh_01x01.mat';
% IBTrACS_ID='2007314N10093';
% interest_area=[89 92 21 23.5];
% marker_size=1;

if exist(ibtracs_save_file,'file')
    load(ibtracs_save_file) % contains tc_track
else
    tc_track=isimip_ibtracs_read('all','',1,1);
end

% find the requested track
if ischar(IBTrACS_ID) % convert to integer for faster search
    IBTrACS_ID_no=str2double(strrep(strrep(IBTrACS_ID,'S','1'),'N','0'));
    fprintf('%s --> %i\n',IBTrACS_ID,IBTrACS_ID_no) % S->1, N->0
else
    IBTrACS_ID_no=IBTrACS_ID;
end

n_tracks=length(tc_track);IBTrACS_index=[];
for track_i=1:n_tracks
    if ~isempty(find(tc_track(track_i).ID_no==IBTrACS_ID_no, 1))
        IBTrACS_index=track_i;
    end
end % track_i

if ~isempty(IBTrACS_index)
    fprintf('processing %s (ID_no %i):\n',tc_track(IBTrACS_index).ID_str,IBTrACS_ID_no);
    
    % load the high-resolution entity
    entity=climada_entity_load(entity);if isempty(entity),return;end
    if ~isempty(findstr(entity.assets.filename,'BGD_Bangladesh_Barisal_01x01'))
        % scale Value to inhabintants of Barisal division (admin1) https://en.wikipedia.org/wiki/Barisal_Division
        f_pop=8325666; % scale with inhabintants of Barisal division (admin1) https://en.wikipedia.org/wiki/Barisal_Division
        f_gdp=1358.8;  % scale with GDP per capita https://data.worldbank.org/indicator/NY.GDP.PCAP.CD
        f_inc=2;  % scale with Wordl bank income group (=1) +1
        fprintf('> scaling asset values with population %i, GDP/capita %5.0f and factor %i\n',f_pop,f_gdp,f_inc)
        entity.assets.Value=entity.assets.Value/sum(entity.assets.Value)*f_pop*f_gdp*f_inc;
        entity.assets.Cover=entity.assets.Value;
    end
    
    % restrict to area of interest (especially for time/memory resons)
    if ~isempty(interest_area)
        fprintf('> restricting to interest_area\n');
        interest_area_pos=(entity.assets.lon>interest_area(1) & entity.assets.lon<interest_area(2)) ...
            & (entity.assets.lat>interest_area(3) & entity.assets.lat<interest_area(4));
        entity.assets=climada_subarray(entity.assets,interest_area_pos);
    end
    
    % generate the TC windfield
    % to improve temporal resolution, consider climada_global.tc.default_min_TimeStep
    climada_global.tc.default_min_TimeStep=.1;
    hazard_TC = climada_tc_hazard_set(tc_track(IBTrACS_index),['_' fig_name '_TC'],entity);
    
    % generate the simply coarse-resolution (ETOPO)  TS surge field
    hazard_TS_ETOP = climada_ts_hazard_set(hazard_TC,['_' fig_name '_TS_ETOP'],'ETOPO',10);
    
    % generate the high-resolution (SRTM) TS surge field
    % to show the SRTM mapping, set 2nd last parameter to -10 (plot takes a long time)
    hazard_TS_SRTM = climada_ts_hazard_set(hazard_TC,['_' fig_name '_TS_SRTM'],'SRTM', 10,1);
    
    % check plots
    % -----------
    figure('Name','BGD_TS_batch','Color',[1 1 1]);
    if ~save_figures,subplot(2,2,1);end
    entity_params.unit_scale=1e9;entity_params.blue_ocean=1;entity_params.title_str='assets';
    climada_entity_plot(entity,marker_size,entity_params); % to check assets
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep 'BGD_assets'],fig_ext);end
    
    if ~save_figures
        subplot(2,2,2);
    else
        figure('Name','BGD_TS_batch','Color',[1 1 1]);
    end
    climada_hazard_plot(hazard_TC,1,marker_size); % to check TC
    title('wind field (TC)');
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep fig_name '_TC'],fig_ext);end
    
    if ~save_figures
        subplot(2,2,3);
    else
        figure('Name','BGD_TS_batch','Color',[1 1 1]);
    end
    climada_hazard_plot(hazard_TS_ETOP,1,marker_size); % to check TS ETOPO1
    title('surge field (TS) - ETOPO1');
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep fig_name '_TS_ETOP'],fig_ext);end
    
    if ~save_figures
        subplot(2,2,4);
    else
        figure('Name','BGD_TS_batch','Color',[1 1 1]);
    end
    climada_hazard_plot(hazard_TS_SRTM,1,marker_size); % to check TS SRTM
    title('surge field (TS) - SRTM');
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep fig_name '_TS_SRTM'],fig_ext);end
    
    figure('Name','BGD_TS_batch difference','Color',[1 1 1]);
    d_hazard_TS_SRTM=hazard_TS_SRTM;d_hazard_TS_SRTM.intensity=d_hazard_TS_SRTM.intensity-hazard_TS_ETOP.intensity;
    params.difference=1;params.title_str='surge field (TS) difference SRTM minus ETOPO1';
    climada_hazard_plot(d_hazard_TS_SRTM,1,marker_size,params); % to check TS SRTM
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep fig_name '_TS_SRTM-ETOP'],fig_ext);end
    
    figure('Name','BGD_TS_batch fraction','Color',[1 1 1]);
    hazard_TS_SRTM.peril_ID='fraction';climada_hazard_plot(hazard_TS_SRTM,1,marker_size); % to check TS SRTM
    title('surge field (TS) SRTM fraction');hazard_TS_SRTM.peril_ID='TS'; % reset!
    if ~isempty(interest_area),xlim(interest_area(1:2)),ylim(interest_area(3:4));end
    if save_figures,saveas(gcf,[fig_dir filesep fig_name '_TS_SRTM_fraction'],fig_ext);end
    
    % eoncode assets
    entity=climada_assets_encode(entity,hazard_TC);
    EDS   =climada_EDS_calc(entity,hazard_TC);
    EDS(2)=climada_EDS_calc(entity,hazard_TS_ETOP);
    EDS(3)=climada_EDS_calc(entity,hazard_TS_SRTM);
    fprintf('modeled damage: TC %f TS ETOPO %f SRTM %f USD mio\n',EDS(1).ED/1e6,EDS(2).ED/1e6,EDS(3).ED/1e6)
    
else
    fprintf('ID_no %i not found in %s\n',IBTrACS_ID_no,ibtracs_save_file);
end

%end % BGD_TS_batch