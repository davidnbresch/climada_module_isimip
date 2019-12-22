function [output_table]=isimip2b_compute_damages_region(RegionID,scenario,damfun_name,params)

if ~isfield(params,'output_folder'),params.output_folder=[climada_global.data_dir filesep 'isimip/results/damages_2b'];end
if ~isdir(params.output_folder),mkdir(params.output_folder);end % create it

% function [status,output_eval_filename,output]=isimip2b_FL_countrydata(RegionID,years_range,params,params_MDR,params_calibration,params_computation)
NatID_RegID_file=[params.RegID_def_folder filesep 'NatID_RegID_isimip_flood_filtered_' num2str(params.years_select_countries(1)) '-' num2str(params.years_select_countries(2)) '.csv'];
if exist(NatID_RegID_file,'file')
    NatID_RegID_flood = readtable(NatID_RegID_file);
else
    [flag,NatID_RegID_flood] = isimip_FL_select_countries(NatID_RegID_file,[],params.years_select_countries);
    if ~flag
        fprintf('** ERROR ** unable to retrieve file %s *****',NatID_RegID_file)
        error('unable to load the NatID_RegID_file')
    end
end
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
if sum(NatID_RegID_flood.Reg_name == RegionID)==0
    error('no country belonging to the given RegionID, perhaps non-existing RegionID?');
end
NatID_RegID_flood = NatID_RegID_flood(NatID_RegID_flood.Reg_name == RegionID,:);

% keep two list of countries: one to be used for calibration, the other one for computation thereafter
% fully excluded countries:
fully_out = ~NatID_RegID_flood.in_extrap;
if sum(fully_out)>0
    fprintf('** WARNING ** these %s countries are entirely excluded: *****\n%s\n\n',num2str(sum(fully_out)),strjoin(NatID_RegID_flood.ISO(fully_out),', '))
end

countries_iso3 = NatID_RegID_flood.ISO(~fully_out);

output_table = table();

for i=1:length(countries_iso3)
    output_table = isimip2b_compute_damages_country(countries_iso3(i),scenario,damfun_name,params);
end

output_file = [RegionID '_' scenario '_' damfun_name 'Entity-Year' num2str(params.entity_year) '-subMATSIRO' num2str(params.subtract_matsiro) '.csv'];
writetable(output_table,output_file);

end