%function res=isimip_FL_prob_country_EDS
% climada isimip_FL_prob_country_EDS
% MODULE:
%   _LOCAL
% NAME:
%   isimip_FL_prob_country_EDS
% PURPOSE:
%   batch job to calculate the EDS for each country based on the pragmatic
%   probabilistic hazard (see isimip_FL_prob_country_hazard) and compare
%   with EM-DAT (where available)
%
%   see PARAMETERS
%   for speedup, consider climada_global.parfor=1
%   to run this on the cluster as batch, see batch_job_template
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%
%   copy job to cluster: scp -r Documents/_GIT/climada_modules/isimip/code/batch/isimip_FL_prob_country_EDS.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   copy data from cluster: scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/results/???_FL_test_DFC_combined*.png /Users/bresch/Documents/_GIT/climada_data/results/FL_test/.
%
%   for more, see http://www.clusterwiki.ethz.ch/brutus/Parallel_MATLAB_and_Brutus
% CALLING SEQUENCE:
%   isimip_FL_prob_country_EDS
% EXAMPLE:
%   bsub -W 2:00 -R "rusage[mem=4000]" -n 24 matlab -nodisplay -singleCompThread -r isimip_FL_prob_country_EDS
%   -W: time, here 2 hours
%   mem: memory, for large jobs, request e.g. 9000
%   -n: number of cluster workers, here 2
%
%   run_on_desktop=1; % to test the job on a desktop
%   isimip_FL_prob_country_EDS
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   run_on_desktop: if you set =1 before calling GBR_FL_test (all in
%       MATLAB command window), it sets the number of parallel pool workers
%       to two, does not delete the pool after execution and does not quit
%       MATLAB and hence allows to TEST a job on a local desktop. This
%       allows to TEST a job without editing it.
%       Default=0 for cluster.
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to stdout and figures
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20190312, initial from isimip_FL_prob_TEST
%-


% CLUSTER JOB PARAMETERS and setup - see PARAMETERS below for specifics
%
FAST_TEST=0; % default=0, if =1, set -R "rusage[mem=500]"
%
cluster_climada_root_dir='/cluster/home/dbresch/climada'; % to make sure the cluster finds climada
cluster_N_pool_workers= 24; % number of parpool workers on pool (same as argument in bsub -n ..)
desktop_N_pool_workers= 2; % number of parpool workers on desktop

% aaa: some admin to start with (up to % eee standard code)
if ~exist('run_on_desktop','var'),run_on_desktop=[];end
if isempty(run_on_desktop),run_on_desktop=0;end % default=0, =1 to run job on mac
if run_on_desktop % for parpool on desktop
    N_pool_workers=desktop_N_pool_workers;
    pool_active=gcp('nocreate');
    if isempty(pool_active),pool=parpool(N_pool_workers);end
else
    cd(cluster_climada_root_dir)
    N_pool_workers=cluster_N_pool_workers;
    pool=parpool(N_pool_workers);
end
startup % climada_global exists afterwards
if exist('startupp','file'),startupp;end
fprintf('executing in %s\n',pwd) % just to check where the job is running from
climada_global.parfor=0; % for parpool, see e.g. climada_tc_hazard_set
t0=clock;
% eee: end of admin (do not edit until here)


% PARAMETERS
%
% whether we use the isimip 2a (='2a', default) simulations or the 2b (='2b')
twoab='2a';
%
% the folder with isimip data on the cluster
cluster_data_folder='/cluster/work/climate/dbresch/climada_data/isimip/';


[country_name,country_ISO3] = climada_country_name('all');
% keep only coutries with a hazard (see isimip_FL_prob_country_hazard)
n_countries=length(country_name);
valid_pos=zeros(1,n_countries,'logical');
entity_file=country_name; % init
hazard_file=country_name; % init
for country_i=1:length(country_name)
    entity_file{country_i}=[cluster_data_folder filesep 'entities' filesep 'FL1950_' country_ISO3{country_i} '_0150as_entity.mat'];
    hazard_file{country_i}=[climada_global.hazards_dir filesep country_ISO3{country_i} '_' twoab '_FL.mat'];
    if exist(entity_file{country_i},'file') && exist(hazard_file{country_i},'file'),valid_pos(country_i)=1;end
end % country_i
country_name=country_name(valid_pos);
country_ISO3=country_ISO3(valid_pos);
entity_file=entity_file(valid_pos);
hazard_file=hazard_file(valid_pos);

n_countries=length(country_name);

entity_temp=climada_entity_load('entity_template'); % load template damage functions
temp_damagefunctions=entity_temp.damagefunctions; % convert to easy broadcast variables
climada_global_results_dir=climada_global.results_dir; % dito

fprintf('processing %i countries with FL hazard:\n',n_countries)

parfor country_i=1:n_countries
    
    saveas_filename=[climada_global_results_dir filesep country_ISO3{country_i} '_FL_test_DFC_combined.png'];
    saveas_filename2=[climada_global_results_dir filesep country_ISO3{country_i} '_FL_test_DFC_combined_emdat.png'];
    if ~exist(saveas_filename,'file') && ~exist(saveas_filename2,'file')
        
        fprintf('\n*** processing %s %s (%i/%i) ***\n',country_ISO3{country_i},country_name{country_i},country_i,n_countries)
        
        entity=climada_entity_load(entity_file{country_i},1); % no-save
        entity.assets.Cover=entity.assets.Value;
        entity.damagefunctions=temp_damagefunctions; % replace to include FL damage function
        hazard=climada_hazard_load(hazard_file{country_i},1); % % load hazard, no-save
        entity=climada_assets_encode(entity,hazard);
        EDS=climada_EDS_calc(entity,hazard);
        EDS(1).annotation_name=[country_ISO3{country_i} ' ' country_name{country_i}];
        
        fig1=figure('Visible','off');
        [~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS);
        xlim([0 250]);
        
        % get EM-DAT for ISO3 FL (and scale to USD bn, as standard in isimip)
        em_data=emdat_read('',country_ISO3{country_i},'FL',1,0);
        if isempty(em_data)
            saveas(fig1,[climada_global_results_dir filesep country_ISO3{country_i} '_FL_test_DFC_combined'],'png'); % in /cluster/work/climate/dbresch/climada_data/results
        else
            em_data.damage          = em_data.damage          /entity.assets.currency_unit;
            em_data.damage_orig     = em_data.damage_orig     /entity.assets.currency_unit;
            em_data.DFC.damage      = em_data.DFC.damage      /entity.assets.currency_unit;
            em_data.DFC_orig.damage = em_data.DFC_orig.damage /entity.assets.currency_unit;
            
            fig1=figure('Visible','off');
            [~,~,legend_str,legend_handle]=climada_EDS_DFC(EDS);
            [legend_str,legend_handle]=emdat_barplot(em_data,'','','EM-DAT indexed',legend_str,legend_handle,'southeast');
            xlim([0 250]);
            saveas(fig1,[climada_global_results_dir filesep country_ISO3{country_i} '_FL_test_DFC_combined_emdat'],'png'); % in /cluster/work/climate/dbresch/climada_data/results
        end
        delete(fig1)
    end % ~exist
    
end % country_i

fprintf('job execution took %f sec.\n',etime(clock, t0))
if ~run_on_desktop,delete(pool);exit;end % no need to delete the pool on mac, the cluster appreciates exit