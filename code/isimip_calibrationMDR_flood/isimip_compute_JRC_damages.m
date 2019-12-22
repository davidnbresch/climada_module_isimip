function [status,output_eval_filename,output]=isimip_compute_JRC_damages(RegionID,years_range,params,params_damfun)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_compute_JRC_damages
% PURPOSE:
%   for a given Region (group of countries) and range of years, load the
%   isimip hazard from all models, the EM-DAT damages, the entity, and
%   the JRC damage function, and compute damages per year and country.
%
% CALLING SEQUENCE:
%   [status,output_eval_filename]=isimip_flood_calibration(RegionID,years_range)
% EXAMPLE:
%   RegionID='NAM';
%   years_range=[1990 2000];
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   params_MDR.damFun_xVals=0:0.5:15;
%   [status,output_eval_filename]=isimip_flood_calibration(RegionID,years_range,params,params_MDR)
% INPUTS:
%   RegionID: Region name (full name)
% OPTIONAL INPUT PARAMETERS:
%   years_range: Range of years to be used. Default = [1971 2010].
%   params: a structure with fields:
%     entity_folder: the folder where the entities are located (default:
%        [climada_global.data_dir filesep 'isimip/entities'] ).
%     entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%        entity_prefix='Try1' will result in Try1_DEU_0150as.mat
%     hazard_folder: the folder where the .mat hazard files are located (default:
%        [climada_global.hazards_dir filesep 'isimip'] ).
%     hazard_raw_folder: the folder where the raw .nc hazard files are located (default:
%        [climada_global.data_dir filesep 'isimip'] ).
%     hazard_protection: one of 'flopros' (default), '0', '100'
%     subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%        fraction from the data. Default =1.
%     entity_year: =0 to use yearly-varying assets (default; assumed is 2005
%        USD value of YEARLY assets). Otherwise, specify fixed year of
%        assets to be used (e.g., 2005). EM-DAT data will also be
%        growth-corrected to the year provided, and left uncorrected if
%        entity_year==0.
%     output_folder: the folder where the calibration results are to be
%        saved. Default in isimip/results/damages_JRCdamFun whitin the
%        climada_global.data_dir folder.
%     verbose_excluded_countries: =1 to display the reason for each
%        excluded country (default=0).
%     years_select_countries: range of years (similar to parameter
%        'years_range') over which countries are selected. By default it is
%        [1992 2010] to choose the same countries as in the calibration.
%   params_damfun: a structure used as input to function
%       continent_jrc_damfun, containing:
%     filename_suffix: a suffix to the JRC filename (e.g. for different
%       PAA). Default is 'PAA1'.
%     filepath: the path which contains the damage function file.
%       Default is [climada_global.data_dir filesep 'isimip/entities/damfun'].
% OUTPUTS:
%   status: 2 if successful, 1 if calibration successful but computation and csv file is not, 0 if even calibration failed.
%   output_eval_filename: a file name for the .csv file generated in applying the calibrated damage function.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190405, initial (from isimip_food_calibration)
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% init output
status=0;
output=[];
output_eval_filename=[];

%% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('RegionID','var'),error('Input parameter RegionID is missing');end %
if ~exist('years_range','var'),             years_range=    [1971 2010];end
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_damfun','var'),           params_damfun=       struct;end


%% check for some parameter fields we need
% params
if ~isfield(params,'entity_folder'),    params.entity_folder=   [climada_global.data_dir filesep 'isimip' filesep 'entities'];end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'entity_prefix'),    params.entity_prefix='FL1950';end
if ~isfield(params,'hazard_folder'),    params.hazard_folder =     [climada_global.hazards_dir filesep 'isimip'];end
if ~isfield(params,'hazard_raw_folder'),params.hazard_raw_folder = [climada_global.data_dir    filesep 'isimip'];end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=1;end
if ~isfield(params,'entity_year'),      params.entity_year=0;end
if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end
if ~isfield(params,'output_folder'),params.output_folder=[climada_global.data_dir filesep 'isimip/results/damages_JRCdamFun'];end
if ~isdir(params.output_folder),mkdir(params.output_folder);end % create it
if ~isfield(params,'verbose_excluded_countries'),params.verbose_excluded_countries=0;end
if ~isfield(params,'years_select_countries'),params.years_select_countries=[1992 2010];end
% params_MDR
if ~isfield(params_damfun,'filename_suffix'),params_damfun.filename_suffix='PAA1';end
if ~isfield(params_damfun,'filepath')
    params_damfun.filepath=[climada_global.data_dir filesep 'isimip/entities/damfun'];
end

% pass hazard folders to climada_global
% used in isimip_get_flood_filename and isimip_get_flood_hazard_filename
climada_global.isimip.hazard_folder=params.hazard_folder;
climada_global.isimip.hazard_raw_folder=params.hazard_raw_folder;

%% 1) Get variables and paths, check countries
isimip_simround = '2a';
% get countries that belong to the region, those used for calibration and
% those used only for 'extrapolation'
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

fprintf('Countries selected:\n%s\n\n',strjoin(countries_iso3,', '))
clear fully_out;

%% 2) prepare output: fill in params for calibrate_MDR_steps including file name to be saved
params_step=struct;
%define filename
filename_head = ['damages_JRC_' RegionID '_' num2str(years_range(1)) '-' num2str(years_range(2))];
filename_haz = ['Haz-Prot' params.hazard_protection '-subMATSIRO' num2str(params.subtract_matsiro)];
filename_ent = ['Entity-Year' num2str(params.entity_year)];
output_eval_filename = [params.output_folder filesep filename_head '_' filename_haz '_' filename_ent '_eval.csv'];
fprintf('Results will be saved in: %s\n', output_eval_filename);

%% 3) load entities , hazards and EM-DAT
[entity_list,hazard_list,emdat_list] = load_entity_hazard_emdat(countries_iso3,years_range,params,isimip_simround);


%% 4) get continent
countries_name_for_continent = countries_iso3;
for i=1:length(countries_iso3)
    % identify country name
    country_iso3 = countries_iso3{i};
    countries_name_for_continent{i} = climada_country_name(country_iso3);
    if isempty(countries_name_for_continent{i})
        switch country_iso3
            case 'SCG'
                countries_name_for_continent{i} = 'Serbia';
            case 'ANT' % actually will return an error for ANT earlier on
                countries_name_for_continent{i} = 'Aruba';
            case 'PSE'
                countries_name_for_continent{i} = climada_country_name('PSX');
            otherwise
                fprintf('*** WARNING: country name not recognized at all for %s **\n',country_iso3);
        end
    end
end



%% 7) Compute damages per country and year for each combination using the JRC damage function

damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 0;

% output contains: country, year, data (obs or model name), damage
output_table_header = {'country', 'year', 'dataset', 'damage'};

all_years = years_range(1):years_range(2);
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
for i=1:length(entity_list)
    [~,damfun_file]=continent_jrc_damfun(countries_name_for_continent{i}, params_damfun);
    [damFun,~]=climada_damagefunctions_read(damfun_file);
    if isempty(damFun)
        fprintf('WARNING: damage function file not available for country: %s, setting modelled damages to missing', countries_iso3{i})
        for j = 1:length(hazard_list{i})
            % get years IDs
            damages_fullmat(i,j,:) = NaN;
            clear temp;
            if i==1
                % extract model name
                [~,temp,~]=fileparts(hazard_list{i}{j}.filename);
                temp=strsplit(temp,'_');
                hazard_model_name{j}=[temp{2} '_' temp{3}];
            end
            inds = ((i-1)*n_models+(j-1))*n_rows_step+(1:n_rows_step);
            country_col(inds) = countries_iso3(i);
            year_col(inds) = all_years;
            dataset_col(inds) = repmat(hazard_model_name(j),[length(all_years) 1]);
            damages_col(inds) = damages_fullmat(i,j,:);
        end
        [~,iis2] = ismember(emdat_list{i}.year, all_years);
        emdat_fullmat(i,iis2) = emdat_list{i}.values;
        inds = ((i-1)*n_models+(n_models-1))*n_rows_step+(1:n_rows_step);
        country_col(inds) = countries_iso3(i);
        year_col(inds) = all_years';
        dataset_col(inds) = repmat({'EM-DAT'},[length(all_years) 1]);
        damages_col(inds) = emdat_fullmat(i,:);
        continue
    end

    for j = 1:length(hazard_list{i})
        temp = climada_EDS_calc_fast(entity_list{i}, hazard_list{i}{j}, damFun, ~params.entity_year,1,'',0,1);
        % get years IDs
        damages_fullmat(i,j,:) = temp.damage*entity_list{i}.assets.currency_unit;
        clear temp;
        if i==1
            % extract model name
            [~,temp,~]=fileparts(hazard_list{i}{j}.filename);
            temp=strsplit(temp,'_');
            hazard_model_name{j}=[temp{2} '_' temp{3}];
        end
        inds = ((i-1)*n_models+(j-1))*n_rows_step+(1:n_rows_step);
        country_col(inds) = countries_iso3(i);
        year_col(inds) = all_years;
        dataset_col(inds) = repmat(hazard_model_name(j),[length(all_years) 1]);
        damages_col(inds) = damages_fullmat(i,j,:);
    end
    [~,iis2] = ismember(emdat_list{i}.year, all_years);
    emdat_fullmat(i,iis2) = emdat_list{i}.values;
    inds = ((i-1)*n_models+(n_models-1))*n_rows_step+(1:n_rows_step);
    country_col(inds) = countries_iso3(i);
    year_col(inds) = all_years';
    dataset_col(inds) = repmat({'EM-DAT'},[length(all_years) 1]);
    damages_col(inds) = emdat_fullmat(i,:);
end

% convert to a table containing: country, year, data (obs or model name), damage
output_table = table(country_col,year_col,dataset_col,round(damages_col),...
    'VariableNames',output_table_header);

climada_global.damage_at_centroid=damage_at_centroid_temp; % reset


% check status, print, write out table
fprintf('Evaluation of isimip_compute_calibrated will be saved in: %s\n', output_eval_filename);
writetable(output_table,output_eval_filename);


end


function [entity_list,hazard_list,emdat_list] = load_entity_hazard_emdat(countries_iso3,years_range,params,isimip_simround)
% function to load lists of entity,hazard,emdat data for a list of
% countries over a range of years
all_years = years_range(1):years_range(2);
%% 1) load entities - N entities for N countries
entity_list=cell(length(countries_iso3),1);
for i=1:length(countries_iso3)
    country_iso3 = countries_iso3{i};
    entity_file_isimip_i=[params.entity_folder filesep params.entity_prefix strtrim(country_iso3) '_0150as_entity'];
    entity_isimip_i=climada_entity_load(entity_file_isimip_i,1); % try to load, flag to 1 to avoir overwrite
    if isempty(entity_isimip_i)
        error('*** ERROR: entity file not found %s\n\n',entity_file_isimip_i);
    end
    entity_isimip_i.assets.centroid_index = 1:length(entity_isimip_i.assets.centroid_index);
    % check that all analysed years are included
    if sum(ismember(all_years, entity_isimip_i.assets.Values_yyyy)) ~= diff(years_range)+1
        error('*** ERROR: not all requested years are available in the entities **\n\n');
    end
    % Subset years
    if params.entity_year
        entity_list{i}=climada_subset_years(entity_isimip_i, 'entity', params.entity_year);
        entity_list{i}.assets.Value=entity_list{i}.assets.Values(1,:);
        entity_list{i}.assets.reference_year = params.entity_year;
    else
        entity_list{i}=climada_subset_years(entity_isimip_i, 'entity', all_years);
    end
end


%% 2) load hazards. 46*N hazards, as there are 46 model combinations
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
if years_range(2)>2001
    forcings = {'gswp3', 'princeton', 'wfdei'};
else
    forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};
end
hazard_list=cell(length(countries_iso3),1);
for i=1:length(countries_iso3)
    hazard_list{i} = {};
    ii = 0;
    for j=1:length(ghms)
        for k=1:length(forcings)
            ghm=ghms{j};
            forcing = forcings{k};
            [flddph_filename,~,fld_path] = isimip_get_flood_filename(isimip_simround, ghm, forcing, params.hazard_protection, 'historical');
            flood_filename=[fld_path filesep flddph_filename];
            if ~exist(flood_filename, 'file')
                % no NetCDF file, warn only on first encounter
                if i==1
                    fprintf('     * Warning: FL hazard (netcdf) missing, probably inexistent model combination %s\n',flood_filename);
                end
                % skip that ghm/forcing combination
                continue
            else
                % NetCDF file exists, has the hazard file already been created?
                hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity_list{i},isimip_simround,[0 0],params.subtract_matsiro);
                if ~exist(hazard_FL_file, 'file')
                    % no hazard file (ghm-forcing-entity combination), create it
                    fprintf('     * Warning: hazard file missing, creating now file %s\n',hazard_FL_file);
                    hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity_list{i},0,isimip_simround,[0 0],'nearest',1,params.subtract_matsiro);
                else
                    % hazard file already exists, just load it
                    hazard_FL=climada_hazard_load(hazard_FL_file);
                end
                % post-process hazard
                hazard_FL.yyyy = double(string(hazard_FL.yyyy));
                try
                    hazard_FL=climada_subset_years(hazard_FL, 'hazard', all_years);
                catch ME
                    switch ME.identifier
                        case 'climada_subset_years:missingYear'
                            fprintf('     * Warning: not all requested years are included in the hazard data - skipping %s\n', flood_filename);
                            continue;
                        otherwise
                            error(['     * Error when subsetting the years - skipping ' flood_filename '***']);
%                             continue;
                    end
                end
                ii=ii+1;
                hazard_list{i}{ii}=hazard_FL;
            end
        end
    end
end


%% 3) load EM-DAT
emdat_list=cell(length(countries_iso3),1);
for i=1:length(countries_iso3)
    % pre-fill data with NaN. Only country/years with surely 0 damage in
    %   EM-DAT will be changed
    emdat_list{i}.year=all_years;
    emdat_list{i}.values=NaN([length(all_years) 1]);
    % check changes in countries for the whole time period (using evalc to
    %   avoid printing message from emdat_get_country_names)
    evalc("[country_emdat,~,cl,em] = emdat_get_country_names(countries_iso3{i},['FL';'F1';'F2'],years_range,0);");%silent
    if cl<0
        % changes_list=-3 or -1 should not happen because these have been filtered out already
        if ismember(cl, [-3 -1]),fprintf('** WARNING ** changes_list=%s for %s - this should NOT HAPPEN',num2str(cl),countries_iso3{i});end
        % changes_list=-2 is ok, but in that case leave NaN
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
            % there is no valid value so we can leave the NaN
            continue
        else
            % some years with changes_list=99, others not, so read data
            emdat_list_temp = emdat_load_yearlysum(country_emdat,['FL';'F1';'F2'],params.entity_year,years_range);
            % use only values where cls is 0,1 or 2
            emdat_list{i}.values(ismember(cls,[0 1 2])) = emdat_list_temp.values(ismember(cls,[0 1 2]));
            % other values are 99 and hence NaN should be left in there
        end
    else
        % no changing country issue - simply read in data
        emdat_list{i} = emdat_load_yearlysum(country_emdat,['FL';'F1';'F2'],params.entity_year,years_range);
    end
end
end

function emdata = emdat_load_yearlysum(country_emdat,peril_ID,exposure_growth,years_range)
% function to load EM-DAT data and sum damages per year
all_years=years_range(1):years_range(2);
emdata.year=all_years;
emdata.values=zeros([length(all_years) 1]);
em_data_i=emdat_read('',country_emdat,peril_ID,exposure_growth,0);
% if EM-DAT data available for this country, use, if not, assume zeros
if ~isempty(em_data_i)
    for iy=1:length(all_years)
        ii=find(all_years(iy) == em_data_i.year);
        if ~isempty(ii)
            emdata.values(iy,1) = sum(em_data_i.damage(ii));
        end
    end
end
end
