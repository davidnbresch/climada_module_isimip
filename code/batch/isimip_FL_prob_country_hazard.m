function hazard=isimip_FL_prob_country_hazard(ISO3,twoab)
% climada isimip_FL_prob_country_hazard
% MODULE:
%   _LOCAL
% NAME:
%   isimip_FL_prob_country_hazard
% PURPOSE:
%   batch job to generate a pragmatically probabilistic hazard set for
%   country. Construct one hazard set with all isimip simulations
%
%   see PARAMETERS
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%
%   copy data from cluster: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/hazards/???_FL.mat /Users/bresch/Documents/_GIT/climada_data/hazards/.
% CALLING SEQUENCE:
%   hazard=isimip_FL_prob_country_hazard(ISO3,twoab)
% EXAMPLE:
%   hazard=isimip_FL_prob_country_hazard('GBR')
% INPUTS:
%   ISO3: the ISO3 code of the country. If ='all', all countries are
%       processed recursively
% OPTIONAL INPUT PARAMETERS:
%   twoab: whether we use the isimip 2a (='2a', default) simulations or the 2b (='2b')
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   a hazard set, also stored as {ISO3}_FL.mat in hazards
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20190312, initial from isimip_FL_prob_TEST
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('ISO3','var'),     ISO3 = 'GBR';   end % default for TEST
if ~exist('twoab','var'),    twoab = '2a';   end % default today (2a)

if strcmp(ISO3,'all') % to process all automatically
    [country_name,country_ISO3] = climada_country_name('all');
    for country_i=1:length(country_name)
        isimip_FL_prob_country_hazard(country_ISO3{country_i});
        return
    end % country_i
end

% PARAMETERS
%
FAST_TEST=1;
%
% the folder with isimip hazard sets on the cluster
%cluster_data_folder='/cluster/work/climate/dbresch/climada_data/isimip/';
cluster_data_folder='/cluster/work/climate/bguillod/climada_data/hazards/isimip';
%
models={
    'FLOOD_CLM_gswp3'
    'FLOOD_CLM_princeton'
    'FLOOD_CLM_watch'
    'FLOOD_CLM_wfdei'
    'FLOOD_DBH_gswp3'
    'FLOOD_DBH_princeton'
    'FLOOD_DBH_watch'
    'FLOOD_DBH_wfdei'
    'FLOOD_H08_gswp3'
    'FLOOD_H08_princeton'
    'FLOOD_H08_watch'
    'FLOOD_H08_wfdei'
    'FLOOD_JULES-TUC_gswp3'
    'FLOOD_JULES-TUC_princeton'
    'FLOOD_JULES-TUC_wfdei'
    'FLOOD_JULES-UoE_gswp3'
    'FLOOD_JULES-UoE_princeton'
    'FLOOD_JULES-UoE_watch'
    'FLOOD_LPJmL_gswp3'
    'FLOOD_LPJmL_princeton'
    'FLOOD_LPJmL_watch'
    'FLOOD_LPJmL_wfdei'
    'FLOOD_MATSIRO_gswp3'
    'FLOOD_MATSIRO_princeton'
    'FLOOD_MATSIRO_watch'
    'FLOOD_MPI-HM_gswp3'
    'FLOOD_MPI-HM_princeton'
    'FLOOD_MPI-HM_watch'
    'FLOOD_MPI-HM_wfdei'
    'FLOOD_ORCHIDEE_gswp3'
    'FLOOD_ORCHIDEE_princeton'
    'FLOOD_ORCHIDEE_watch'
    'FLOOD_ORCHIDEE_wfdei'
    'FLOOD_PCR-GLOBWB_gswp3'
    'FLOOD_PCR-GLOBWB_princeton'
    'FLOOD_PCR-GLOBWB_watch'
    'FLOOD_PCR-GLOBWB_wfdei'
    'FLOOD_VEGAS_gswp3'
    'FLOOD_VIC_gswp3'
    'FLOOD_VIC_princeton'
    'FLOOD_VIC_watch'
    'FLOOD_VIC_wfdei'
    'FLOOD_WaterGAP_gswp3'
    'FLOOD_WaterGAP_princeton'
    'FLOOD_WaterGAP_watch'
    'FLOOD_WaterGAP_wfdei'
    };
%
if FAST_TEST,models=models(1:3);end % only first three models
%
%ext_nomg=[      '_0_gev_0.1_FL1950_' ISO3 '_0150as_mFRCmatsiro_FL.mat'];
ext_mngd=['_flopros_gev_0.1_FL1950_' ISO3 '_0150as_mFRCmatsiro_FL.mat'];


n_models=length(models);
info.n_models=n_models;
info.models=models;
info.ext_mngd=ext_mngd;
info.model_i=[];
info.event_ID=[];
info.frequency=[];
hazard=[];
% first, test for country being available
model_i=1;
hazard_file=[cluster_data_folder filesep twoab filesep models{model_i} ext_mngd];
if exist(hazard_file,'file')
    fprintf('construct %s hazard set with all %i simulations:\n',ISO3,n_models);
    for model_i=1:n_models
        fprintf('- dealing with model %s (%i of %i)\n',models{model_i},model_i,n_models);
        hazard_file=[cluster_data_folder filesep twoab filesep models{model_i} ext_mngd];
        hazard_mngd=climada_hazard_load(hazard_file,1); % no save
        info.model_i = [info.model_i   ones(1,hazard_mngd.event_count)*model_i];
        info.event_ID= [info.event_ID  hazard_mngd.event_ID];
        info.frequency=[info.frequency hazard_mngd.frequency];
        hazard_mngd.yyyy   =str2double(hazard_mngd.yyyy)';
        hazard_mngd.mm     =str2double(hazard_mngd.mm)';
        hazard_mngd.dd     =str2double(hazard_mngd.dd)';
        hazard_mngd.datenum=hazard_mngd.datenum';
        if isempty(hazard)
            hazard=hazard_mngd;
        else
            hazard=climada_hazard_merge(hazard,hazard_mngd,'events');
        end
    end % model_i
    % re-define frequency
    hazard.orig_years=hazard.orig_event_count; % we know each event is a year
    hazard.frequency=hazard.frequency*0+1/hazard.orig_years;
    hazard.comment='a collection of isimip FL hazars sets, see hazard.info';
    fprintf('  NOTE: event frequency redefined as 1/%i years (%i models combined)\n',hazard.orig_years,n_models);
    hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity); % update
    hazard.info=info; % add info
    hazard_filename=[climada_global.hazards_dir filesep ISO3 '_FL'];
    fprintf('saving combined hazard as %s ..',hazard_filename);
    save(hazard_filename,'hazard');
    fprintf(' done\n');
else
    fprintf('NOTE: for %s no FL hazard exists\n',ISO3);
end

end % isimip_FL_prob_country_hazard

% % and to assess a ptf locally (after copying {ISO3}_FL_test_merged.mat to local)
% entity=climada_entity_load('GBR_UnitedKingdom_10x10');
% hazard=climada_hazard_load([climada_global.results_dir  filesep 'GBR_FL_test' filesep 'GBR_FL_test_merged']);
% entity=climada_assets_encode(entity,hazard);
% EDS=climada_EDS_calc(entity,hazard);