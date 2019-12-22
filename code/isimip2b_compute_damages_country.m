function [output_table]=isimip2b_compute_damages_country(country_iso3,scenario,damfun_name,params)
% function [status,output_eval_filename,output]=isimip2b_FL_countrydata(RegionID,years_range,params,params_MDR,params_calibration,params_computation)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   isimip_flood_calibration
% PURPOSE:
%   for a given Region (group of countries) and range of years, load the
%   isimip hazard from all models, the EM-DAT damages, the entity, and
%   define a functional shape for the damage function together with input
%   parameters, then call 'calibrate_MDR_steps' to do the calibration.
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
%   years_range: Range of years to be used. Default = [1992 2010].
%   params: a structure with fields:
%     entity_folder: the folder where the entities are located (default:
%        [climada_global.data_dir filesep 'isimip/entities'] ).
%     hazard_folder: the folder where the .mat hazard files are located (default:
%        [climada_global.hazards_dir filesep 'isimip'] ).
%     hazard_raw_folder: the folder where the raw .nc hazard files are located (default:
%        [climada_global.data_dir filesep 'isimip'] ).
%     hazard_protection: one of 'flopros' (default), '0', '100'
%     subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%        fraction from the data. Default =0.
%     entity_year: =0 to use yearly-varying assets (default; assumed is 2005
%        USD value of YEARLY assets). Otherwise, specify fixed year of
%        assets to be used (e.g., 2005). EM-DAT data will also be
%        growth-corrected to the year provided, and left uncorrected if
%        entity_year==0.
%     output_folder: the folder where the calibration results are to be
%        saved. Default in isimip/results/calibration whitin the
%        climada_global.data_dir folder.
%     keep_countries_0emdat: how to deal with countries with 0 EM-DAT
%        damage, as follows (see also emdat_get_country_names, emdat_isdata+1):
%        0 = keep all countries with an EM-DAT entry even if no
%        damage>0 is recorded;
%        1 = keep countries with at least one entry with damage>0 even if
%        these entries are outside of the considered years and perilID;
%        2 (default) = to keep only the countries with at least one entry with
%        damage>0 for the requested peril_ID and years.
%     verbose_excluded_countries: =1 to display the reason for each
%        excluded country (default=0).
% OUTPUTS:
%   status: 2 if successful, 1 if calibration successful but computation and csv file is not, 0 if even calibration failed.
%   output_eval_filename: a file name for the .csv file generated in applying the calibrated damage function.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180711, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181009, many changes
%    including a new input parameter added 'entity_year'
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181015, additional input
%    parameter params.damFun_xVals.
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181107, fixed formatting
%    of parameters in the fprint at the end
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181127, add output,
%    direct definition of MDR_fun and adding call to isimip_compute_calibrated
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190103, filter out countries which do not 'exist' as ISO3
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190114, adding params.keep_countries_0emdat as a country filtering criterion.
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190116, many improvements such as treatment of EM-DAT, application of isimip_compute_calibrated for years not used in calibration, removing old/unused code extracts, removing parallelization of optimization
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190121, new options params_calibration.calib_options and params_computation.do
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190121, use of strcmp to determine params_computation.do
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, adding RTarea as a possible type of cost function, and adding parameter params_calibration.underestimation_factor
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, set default value of params_MDR.damFun_xVals to 0:0.5:15
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190207, bug fix: set default value of params_calibration.calib_options.params.step_tolerance to 0.001
% David N. Bresch, david.bresch@gmail.com, 20190211, climada_global.isimip.*
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% init output
output_table=[];

%% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('RegionID','var'),error('Input parameter RegionID is missing');end %
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_MDR','var'),              params_MDR=          struct;end
if ~exist('params_calibration','var'),      params_calibration=  struct;end
if ~exist('params_computation','var'),      params_computation=  struct;end

years_range_calib=[1992 2010];

%% check for some parameter fields we need
% params
%if ~isfield(params,'entity_folder'),    params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';end
if ~isfield(params,'entity_folder'),    params.entity_folder=   [climada_global.data_dir filesep 'isimip' filesep 'entities'];end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'hazard_folder'),    params.hazard_folder =     [climada_global.hazards_dir filesep 'isimip'];end
if ~isfield(params,'hazard_raw_folder'),params.hazard_raw_folder = [climada_global.data_dir    filesep 'isimip'];end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=1;end
if ~isfield(params,'entity_year'), params.entity_year=0;end

calib_params_folder = [climada_global.data_dir filesep 'isimip/results/calibration'];

% pass hazard folders to climada_global
% used in isimip_get_flood_filename and isimip_get_flood_hazard_filename
climada_global.isimip.hazard_folder=params.hazard_folder;
climada_global.isimip.hazard_raw_folder=params.hazard_raw_folder;



%% 1) Get variables and paths, check countries
isimip_simround = '2b';
% get countries that belong to the region, those used for calibration and
% those used only for 'extrapolation'
NatID_RegID_file=[params.RegID_def_folder filesep 'NatID_RegID_isimip_flood_filtered_' num2str(years_range_calib(1)) '-' num2str(years_range_calib(2)) '.csv'];
if exist(NatID_RegID_file,'file')
    NatID_RegID_flood = readtable(NatID_RegID_file);
else
    error('unable to load the NatID_RegID_file')
%     [flag,NatID_RegID_flood] = isimip_FL_select_countries(NatID_RegID_file,[],years_range);
%     if ~flag
%         fprintf('** ERROR ** unable to retrieve file %s *****',NatID_RegID_file)
%         error('unable to load the NatID_RegID_file')
%     end
end
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
if sum(NatID_RegID_flood.ISO == country_iso3)==0
    error('country not found in NatID_RegID_isimip_flood_filtered_1992-2010.csv');
end
RegionID = NatID_RegID_flood.RegionID(NatID_RegID_flood.ISO == country_iso3);
clear fully_out;

%% 2) load damage function
% create the damage function 'damFun'
if strcmp(damfun_name,'JRC')
    country_name_for_continent = climada_country_name(country_iso3);
    if isempty(country_name_for_continent)
        switch country_iso3
            case 'SCG'
                country_name_for_continent = 'Serbia';
            case 'ANT' % actually will return an error for ANT earlier on
                country_name_for_continent = 'Aruba';
            case 'PSE'
                country_name_for_continent = climada_country_name('PSX');
            otherwise
                fprintf('*** WARNING: country name not recognized at all for %s **\n',country_iso3);
        end
    end
    params_damfun=struct;
    params_damfun.filename_suffix='PAA1';
    params_damfun.filepath=[climada_global.data_dir filesep 'isimip/entities/damfun'];
    [~,damfun_file]=continent_jrc_damfun(country_name_for_continent, params_damfun);
    [damFun,~]=climada_damagefunctions_read(damfun_file);
    clear country_name_for_continent params_damfun damfun_file;
else
    if strcmp(damfun_name,'DFC')
        calib_file_end = '_1992-2010_RTarea-uf2-MMM-patternsearch-reg3-mesh0.25-step0.001_Haz-Protflopros-subMATSIRO1_Entity-Year0_Filters-Country2-emdat0-YDS0_pars0.0001-1-0.0001-5.mat';
    elseif strcmp(damfun_name,'YL2')
        calib_file_end = '_1992-2010_dlog2-uf2-MMM-patternsearch-reg3-mesh0.25-step0.001_Haz-Protflopros-subMATSIRO1_Entity-Year0_Filters-Country2-emdat0-YDS0_Filter-RegYears0_pars0.0001-1-0.0001-5.mat';
    elseif strcmp(damfun_name,'YL2aY')
        calib_file_end = '_1992-2010_dlog2-uf2-MMM-patternsearch-reg3-mesh0.25-step0.001_Haz-Protflopros-subMATSIRO1_Entity-Year0_Filters-Country2-emdat0-YDS0_pars0.0001-1-0.0001-5.mat';
    else
        error('unexpected value in damfun_name');
    end
    calib_file = [calib_params_folder filesep 'calib_' RegionID calib_file_end];
    opt_pars = load(calib_file).optimal_pars;
    MDR_fun=@(x,pars)pars(1)*(1-exp(-pars(2)*x));
    damFun = climada_damagefunctions_generate_from_fun(0:0.5:15, MDR_fun, opt_pars);
    clear calib_file_end calib_file opt_pars MDR_fun;
end



%% 3) load entities and hazards
if strcmp(scenario,'historical')
    years_range = [1950 2005];
else
    years_range = [2006 2099];
end
[entity,hazard_list] = load_entity_hazard(country_iso3,years_range,params,isimip_simround,scenario);


%% 7) Compute damages per country and year for each combination using the identified optimal parameter combination
damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 0;

nyears = length(all_years);
n_rows = length(hazard_list)*nyears;
output_table_header = {'country', 'year', 'dataset', 'damage'};
damages_col = NaN([n_rows 1]);
year_col = NaN([n_rows 1]);
dataset_col = NaN([n_rows 1]);
country_col=repmat(country_iso3,[n_rows 1]);
for j = 1:length(hazard_list)
    temp = climada_EDS_calc_fast(entity, hazard_list{j}, damFun, ~params.entity_year,1,'',0,1);
    inds = (j-1)*nyears+(1:nyears);
    year_col(inds) = all_years;
    damages_col(inds) = temp.damage*entity.assets.currency_unit;
    clear temp;
    % extract model name
    [~,temp,~]=fileparts(hazard_list{j}.filename);
    temp=strsplit(temp,'_');
    hazard_model_name_j = [temp{2} '_' temp{3}];
    dataset_col(inds) = repmat(hazard_model_name_j,[length(all_years) 1]);
    % TODO CLEAN UP THIS MESS
    % columns needed:
    % 1) those given as input (can be created afterwards): country, scenario, damfun_name, entity_year (2005 or transient)
    % 2) others (absolutely needed): year, dataset, damage
    % in other cases we had: country, year, dataset, damage
end

% convert to a table containing: country, year, data (obs or model name), damage
output_table = table(country_col,year_col,dataset_col,round(damages_col),...
    'VariableNames',output_table_header);

climada_global.damage_at_centroid=damage_at_centroid_temp; % reset

end


function [entity,hazard_list] = load_entity_hazard(country_iso3,years_range,params,isimip_simround,scenario)
% function to load lists of entity,hazard,emdat data for a list of
% countries over a range of years
all_years = years_range(1):years_range(2);
%% 1) load entities - N entities for N countries
if strcmp(scenario, 'historical') || params.entity_year==2005
    entity_prefix = 'FL1950';
else
    entity_prefix = 'FL2006-2099';
end
entity_file_isimip=[params.entity_folder filesep entity_prefix strtrim(country_iso3) '_0150as_entity'];
entity_isimip_i=climada_entity_load(entity_file_isimip,1); % try to load, flag to 1 to avoir overwrite
if isempty(entity_isimip_i)
    error('*** ERROR: entity file not found %s\n\n',entity_file_isimip);
end
entity_isimip_i.assets.centroid_index = 1:length(entity_isimip_i.assets.centroid_index);

% Subset years
if params.entity_year
    % check that all analysed years are included
    if ~ismember(params.entity_year, entity_isimip_i.assets.Values_yyyy)
        error('*** ERROR: the requested year is not available in the entity **\n\n');
    end
    entity=climada_subset_years(entity_isimip_i, 'entity', params.entity_year);
    entity.assets.Value=entity.assets.Values(1,:);
    entity.assets.reference_year = params.entity_year;
else
    % check that all analysed years are included
    if sum(ismember(all_years, entity_isimip_i.assets.Values_yyyy)) ~= diff(years_range)+1
        error('*** ERROR: not all requested years are available in the entities **\n\n');
    end
    entity=climada_subset_years(entity_isimip_i, 'entity', all_years);
end


%% 2) load hazards. 46*N hazards, as there are 46 model combinations
ghms = {'H08', 'LPJmL', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'WaterGAP2'};
forcings = {'gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc5'};
hazard_list = {};
ii = 0;
for j=1:length(ghms)
    for k=1:length(forcings)
        ghm=ghms{j};
        forcing = forcings{k};
        [flddph_filename,~,fld_path] = isimip_get_flood_filename(isimip_simround, ghm, forcing, params.hazard_protection, scenario);
        flood_filename=[fld_path filesep flddph_filename];
        if ~exist(flood_filename, 'file')
            % no NetCDF file
            fprintf('     * Warning: FL hazard (netcdf) missing, probably inexistent model combination %s\n',flood_filename);
            % skip that ghm/forcing combination
            continue
        else
            % NetCDF file exists, has the hazard file already been created?
            hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity,isimip_simround,[0 0],params.subtract_matsiro);
            if ~exist(hazard_FL_file, 'file')
                % no hazard file (ghm-forcing-entity combination), create it
                fprintf('     * Warning: hazard file missing, creating now file %s\n',hazard_FL_file);
                hazard_FL=isimip_flood_load(flood_filename,hazard_FL_file,entity,0,isimip_simround,[0 0],'nearest',1,params.subtract_matsiro);
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
            hazard_list{ii}=hazard_FL;
        end
    end
end
end
