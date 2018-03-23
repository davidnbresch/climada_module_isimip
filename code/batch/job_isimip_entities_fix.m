% batch job for cluster: bsub -W 4:00 -R "rusage[mem=1000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip_entities_fix
% MODULE:
%   isimip
% NAME:
%   job_isimip_entities_fix
% PURPOSE:
%   generate isimip entities for all countries
%
%   some hints to work with the cluster (explicit paths, edit this ;-)
%   copy job to cluster:       scp -r Documents/_GIT/climada_modules/isimip/code/batch/job_isimip_entities_fix.m dbresch@euler.ethz.ch:/cluster/home/dbresch/euler_jobs/.
%   run on cluster:            bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip_entities_fix
%
%   copy results back local:   scp -r dbresch@euler.ethz.ch:/cluster/work/climate/dbresch/climada_data/entities/*.mat Documents/_GIT/climada_data/entities/.
% CALLING SEQUENCE:
%   bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip_entities_fix
% EXAMPLE:
%   bsub -R "rusage[mem=1000]" -n 24 matlab -nodisplay -singleCompThread -r job_isimip_entities_fix
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   to work disk, see PARAMETERS
% MODIFICATION HISTORY:
% David N. Bresch, dbresch@ethz.ch, 20180323, copy from job_isimip_entities
%-


% PARAMETERS
% one should only have to edit this section
%
cd /cluster/home/dbresch/climada % to make sure the cluster finds climada

startup % climada_global exists afterwards
pwd % just to check where the job is running from
% N_pool_workers=24; % for parpool
% climada_global.parfor=1; % for parpool

% pool=parpool(N_pool_workers);

climada_global.entities_dir='/cluster/work/climate/dbresch/climada_data/isimip/entities';
% params.grid_resolution='0150as';
% params.entity_prefix='FL1950_';
% isimip_gdp_entity('all',params,1950,2020);

workdir='/cluster/work/climate/dbresch/climada_data/isimip/entities';
dd=dir([workdir filesep '*0150as*.mat']);
%for file_i=1:length(dd)
for file_i=70:length(dd)
    if ~dd(file_i).isdir
        fprintf('%s:\n',dd(file_i).name);
        entity=climada_entity_load([workdir filesep dd(file_i).name]);
        
        try
            % fix the rounding errors in the centroid coordinates
            if strfind(dd(file_i).name,       '0360as')
                entity = fix_coords_isimip(entity,'0360as');
            elseif strfind(dd(file_i).name,   '0150as')
                entity = fix_coords_isimip(entity,'0150as');
            end
            
            save(entity.assets.filename,'entity',climada_global.save_file_version);
            
        catch ME
            fprintf('ERROR: catch:\n');
            ME
        end

    end % isdir
end % file_i

%delete(pool)

exit % the cluster appreciates this, gives back memory etc.