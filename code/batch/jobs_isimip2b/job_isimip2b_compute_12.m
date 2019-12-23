% batch job for cluster: bsub -W 4:00 -R "rusage[mem=5000]" -n 1 matlab -nodisplay -singleCompThread -r job_isimip2b_compute_all


cd /cluster/home/bguillod/climada % to make sure the cluster finds climada

startup % climada_global exists afterwards
pwd % just to check where the job is running from
climada_global.entities_dir='/cluster/work/climate/bguillod/climada_data/isimip/entities';

years_range_calib = [1992 2010];
RegID_def_folder=[climada_global.data_dir filesep 'isimip'];
NatID_RegID_file=[RegID_def_folder filesep 'NatID_RegID_isimip_flood_filtered_' num2str(years_range_calib(1)) '-' num2str(years_range_calib(2)) '.csv'];
if exist(NatID_RegID_file,'file')
    NatID_RegID_flood = readtable(NatID_RegID_file);
else
    error('unable to load the NatID_RegID_file')
end

all_regions = unique(NatID_RegID_flood.Reg_name);
% remove PIS1 and PIS2
all_regions = setdiff(all_regions,{'PIS1','PIS2',''});
all_scenarios = {'historical','rcp26','rcp60'};
all_damfuns = {'JRC','DFC','YL2','YL2aY'};
entity_years = [0 2005];
params=              struct;


for i=12:12%length(all_regions)
    for j=1:length(all_damfuns)
        for k=1:length(all_scenarios)
            for l=1:length(entity_years)
                if strcmp(all_scenarios{k},'historical') && entity_years(l)==0
                    continue
                end
                fprintf('region %s, damfun %s, scen %s, entity %s\n',all_regions{i},all_damfuns{j},all_scenarios{k},num2str(entity_years(l)));
                params.entity_year=entity_years(l);
                temp = isimip2b_compute_damages_region(all_regions{i},all_scenarios{k},all_damfuns{j},params);
            end
        end
    end
end
print('DONE')
exit
