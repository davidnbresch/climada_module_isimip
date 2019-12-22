function [emdat_table] = get_emdat_isimip_FL(RegionID,years_range)

global climada_global
% this function reads all the flood events from EM-DAT for a region and
% writes out a csv file with each individual event. Will allow to count the
% number of events per country/region and range of years

% % RegionID
% fileID=fopen('/cluster/work/climate/bguillod/climada_data/isimip/RegID_names_isimip_flood.txt','r');
% regions_def=strsplit(fscanf(fileID,'%c'),'\n');
% fclose(fileID);
% % skip PIS1 and PIS2 - not enough data there anyway
% regions_def=regions_def(1:15);
% regionIDs = {};
% for i=1:length(regions_def)
% temp = strsplit(regions_def{i},':');
% regionsIDs{i} = temp{1};
% get_emdat_isimip_FL(regionsIDs{i},[1992 2010])
% end

RegID_def_folder=[climada_global.data_dir filesep 'isimip'];

NatID_RegID_file=[RegID_def_folder filesep 'NatID_RegID_isimip_flood_filtered_' num2str(years_range(1)) '-' num2str(years_range(2)) '.csv'];
NatID_RegID_flood = readtable(NatID_RegID_file);
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
if sum(NatID_RegID_flood.Reg_name == RegionID)==0
    error('no country belonging to the given RegionID, perhaps non-existing RegionID?');
end
NatID_RegID_flood = NatID_RegID_flood(NatID_RegID_flood.Reg_name == RegionID,:);
fully_out = ~NatID_RegID_flood.in_extrap;
if sum(fully_out)>0
    fprintf('** WARNING ** these %s countries are entirely excluded: *****\n%s\n\n',num2str(sum(fully_out)),strjoin(NatID_RegID_flood.ISO(fully_out),', '))
end
% countries excluded from calibration but retained for application
calib_out = (NatID_RegID_flood.in_calib < 2) & NatID_RegID_flood.in_extrap;
if sum(calib_out)>0
    fprintf('** WARNING ** these %s countries are excluded from calibration but are included in application (extrapolation): *****\n%s\n\n',num2str(sum(calib_out)),strjoin(NatID_RegID_flood.ISO(calib_out),' , '))
end
% countries retained
countries_keep = ~fully_out & ~calib_out;
if sum(~countries_keep)==0
    fprintf('Region %s: No country has to be skipped\n',RegionID);
end
if sum(countries_keep) == 0
    fprintf('ISSUE: Region %s: No country remains\n',RegionID);
end

% only done for the calibrated countries, for the other ones EM-DAT is not
% relevant!
countries_iso3 = NatID_RegID_flood.ISO(countries_keep)';

% initialize emmpty structure
emdat_table = table('Size',[0,3],'VariableNames',{'year','country','damage'},'VariableTypes',{'double','string','double'});
all_years = years_range(1):years_range(2);
for i=1:length(countries_iso3)
    evalc("[country_emdat,~,cl,em] = emdat_get_country_names(countries_iso3{i},['FL';'F1';'F2'],years_range,0);");%silent
    if cl<0
        % changes_list=-3 or -1 should not happen because these have been filtered out already
        if ismember(cl, [-3 -1]),fprintf('** WARNING ** changes_list=%s for %s - this should NOT HAPPEN',num2str(cl),countries_iso3{i});end
        % changes_list=-2 is ok, but no value to be added
        continue
    elseif cl == 99
        % check for each year whether data is ok to be used (cls is changes_list from each individual year)
        cls = NaN([length(all_years) 1]);
        for yi=1:length(all_years)
            evalc("[country_emdat_yi,~,cl_yi,em_yi] = emdat_get_country_names(countries_iso3{i},['FL';'F1';'F2'],repmat(all_years(yi),[1 2]),0);");%silent
            cls(yi) = cl_yi;
        end
        if sum(cls < 0)
            % there shouldn't be negative changes_list values in any year but just check
            fprintf('** WARNING ** changes_list=%s for %s on year %s - this should NOT HAPPEN',num2str(cls(cls<0)),countries_iso3{i},num2str(all_years(cls<0)))
        end
        if all(cls==99)
            % there is no valid value so nothing to add
            continue
        else
            % some years with changes_list=99, others not, so read data
            emdat_table_temp = emdat_load_events(country_emdat,years_range);
            % use only values where cls is 0,1 or 2
            years_ok = all_years(ismember(cls,[0 1 2]));
            emdat_table_temp = emdat_table_temp(ismember(emdat_table_temp.year,years_ok),:);
            emdat_table = [emdat_table; emdat_table_temp];
            % other values are 99 and hence should not be added in there
        end
    else
        % no changing country issue - simply read in data
        emdat_table_temp = emdat_load_events(country_emdat,years_range);
        emdat_table = [emdat_table; emdat_table_temp];
    end
    % output wanted: country, year, damage, where several country-year
    %   cases can exist if several events in the same country.
end

output_file = [climada_global.data_dir filesep 'isimip/results/emdat' filesep 'emdat_FloodEvents_' RegionID '_' num2str(years_range(1)) '-' num2str(years_range(2)) '.csv'];
writetable(emdat_table, output_file);
end

function emdata = emdat_load_events(country_emdat,years_range)
% function to load EM-DAT data and sum damages per year
peril_ID = ['FL';'F1';'F2'];
exposure_growth = 0;
em_data_i=emdat_read('',country_emdat,peril_ID,exposure_growth,0);
emdata = table('Size',[0,3],'VariableNames',{'year','country','damage'},'VariableTypes',{'double','string','double'});
% keep only years within range
if ~isempty(em_data_i)
    ii=find(em_data_i.year >= years_range(1) & em_data_i.year <= years_range(2));
    if ~isempty(ii)
        emdata = table(em_data_i.year(ii),repmat(string(country_emdat),[length(ii) 1]),em_data_i.damage(ii),...
            'VariableNames',{'year','country','damage'});
        emdata = sortrows(emdata,'year');
    end
end
end
