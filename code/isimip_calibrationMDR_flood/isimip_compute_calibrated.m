function [ status, output_table ] = isimip_compute_calibrated(entity_list, hazard_list, ...
    emdat_list, years_range, damFun, use_YDS, years_i_in)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_compute_calibrated
% PURPOSE:
%   Use after the identification of calibrated parameters. This function
%      (1) calculate damages for each entity and hazard with the calibrated
%          damage function
%      (2) creates a table with yearly damages per country for each
%          hazard set and for EM-DAT
% CALLING SEQUENCE:
%   [ status,output_table ]=isimip_compute_calibrated(entity_list,hazard_list,...
%       emdat_list, years_range, damFun, use_YDS, years_i_in)
% EXAMPLE:
%   see isimip_flood_calibration.m
% INPUTS:
%   entity_list: list of entities, see function isimip_flood_calibration
%   hazard_list: list of hazard, see function isimip_flood_calibration
%   emdat_list: list of observed damages, see function isimip_flood_calibration
%   years_range: range of years for which to create the table (e.g. [1971 2010])
%   damFun: the damagefunction to be used, e.g. as generated by climada_damagefunctions_generate_from_fun
%   years_i_in: logical array indicating which years are kept for which
%      entity (indices of years_range).
% OPTIONAL INPUT PARAMETERS:
%   use_YDS: =1 to use yearly-varying assets. Default =0. Take care to be
%       consistent with how emdat_list was generated (if exposure_growth=0
%       in emdat_read, set use_YDS=1; otherwise use_YDS=0 and the entities
%       in entity_list should correspond to the year used for
%       exposure_growth, e.g. 2005)
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_table: a table containing the damage data.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181107, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181127, use table
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190116, remove params_MDR in call to climada_damagefunctions_generate_from_fun
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190116, change input parameter for easier call from isimip_flood_calibration and reset climada_global.damage_at_centroid
%   
%-

% initialization
global climada_global
status = 0; %init
damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 0;

if ~exist('use_YDS','var'),             use_YDS=0;    end

% retrieve country names
country_list=cellfun(@(x) x.assets.admin0_ISO3,entity_list, 'UniformOutput', 0);

% most contain: country, year, data (obs or model name), damage
output_table_header = {'country', 'year', 'dataset', 'damage', 'used_in_calibration'};
% compute damages
all_years = years_range(1):years_range(2);
% prepare output vectors - numbers of 'models' is +1 because there is also EM-DAT
n_models = length(hazard_list{1})+1;
n_rows = length(entity_list)*n_models*length(all_years);
n_rows_step = length(all_years);
damages_fullmat = NaN([length(entity_list),length(hazard_list{1}),length(all_years)]);
emdat_fullmat = NaN([length(entity_list),length(all_years)]);
hazard_model_name=cell(length(hazard_list{1}),1);
country_col=cell([n_rows 1]);
year_col=NaN([n_rows 1]);
dataset_col=country_col;
damages_col=year_col;
calib_col=year_col;
for i=1:length(entity_list)
    for j = 1:length(hazard_list{i})
        temp = climada_EDS_calc_fast(entity_list{i}, hazard_list{i}{j}, damFun, use_YDS,1,'',0,1);
        % get years IDs
        [~,iis] = ismember(hazard_list{i}{j}.yyyy', all_years);
        damages_fullmat(i,j,iis) = temp.damage*entity_list{i}.assets.currency_unit;
        clear temp;
        if i==1
            % extract model name
            [~,temp,~]=fileparts(hazard_list{i}{j}.filename);
            temp=strsplit(temp,'_');
            hazard_model_name{j}=[temp{2} '_' temp{3}];
        end
        inds = ((i-1)*n_models+(j-1))*n_rows_step+(1:n_rows_step);
        country_col(inds) = country_list(i);
        year_col(inds) = all_years;
        dataset_col(inds) = repmat(hazard_model_name(j),[length(all_years) 1]);
        damages_col(inds) = damages_fullmat(i,j,:);
        calib_col(inds) = years_i_in(:,i);
    end
    [~,iis2] = ismember(emdat_list{i}.year, all_years);
    emdat_fullmat(i,iis2) = emdat_list{i}.values;
    inds = ((i-1)*n_models+(n_models-1))*n_rows_step+(1:n_rows_step);
    country_col(inds) = country_list(i);
    year_col(inds) = all_years';
    dataset_col(inds) = repmat({'EM-DAT'},[length(all_years) 1]);
    damages_col(inds) = emdat_fullmat(i,:);
    calib_col(inds) = years_i_in(:,i);
end

status=1;

% convert to a table containing: country, year, data (obs or model name),
%                                damage, year_used_in_calibration
output_table = table(country_col,year_col,dataset_col,damages_col,calib_col,...
    'VariableNames',output_table_header);

climada_global.damage_at_centroid=damage_at_centroid_temp; % reset

end

% what we want to plot:
% 1) damage function
% 2) For the regional sum (each model + EM-DAT):
%   a) return time plot
% 	b) time series of yearly values
%   c) scatter plot of yearly values (obs vs model)
% 3) At country level, all countries in one plot:
%   a) scatter plot of country ED
%   b) scatter plot of country return values
% 4) At country level, plot by country:
%   a) return time plot
%   b) time series of yearly values
%   c) scatter plot yearly values (obs vs model)