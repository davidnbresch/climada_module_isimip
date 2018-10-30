function [status,output_filename]=isimip_flood_calibration(RegionID,years_range,params,params_MDR,params_calibration)
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
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% EXAMPLE:
%   RegionID='NAM';
%   years_range=[1990 2000];
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   params_MDR.damFun_xVals=0:0.5:15;
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% INPUTS:
%   RegionID: Region name (full name)
%   years_range: Range of years to be used (e.g., [1990 2010])
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%     entity_folder: the folder where the entities are located (default:
%        [climada_global.data_dir filesep 'isimip/entities'] ).
%     entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%        entity_prefix='Try1' will result in Try1_DEU_0150as.mat
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
%   params_MDR: parameters to define the MDR function. A structure with fields:
%     remove_years_0emdat: determines how to filter years with zero EM-DAT
%        damages (e.g., no damage).
%        0=no filter (default; all years retained);
%        1=remove years for which the sum of EM-DAT damages over all
%        countries is 0;
%        2=for each country, remove years with zero EM-DAT damages.
%     remove_years_0YDS: list defining whether and how to exclude years
%        with zero simulated damage as follows:
%        do: =1 to exclude years with zero simulated damage. Default=0.
%           If =0, the other fields are not required.
%        threshold: hazard intensity value above which (>) damage is
%           considered to occur to filter out years (required if do==1).
%        what: determines how to choose if several hazards are provided per
%           country (possible choices: 'any', 'all', 'mean', 'median'). For
%           instance, if 'any', then if any hazard set gives 0 damage the
%           year-country event is removed.
%        min_val: damage value (e.g. in USD) below which damage is not
%           considered to occur (default =0).
%     pars_range: range of parameter values (cell array, each element as
%        [min_val max_val]). Default {[0.0001,1],[0.0001,5]}.
%     damFun_xVals: vector of value of hazard intensity to be used when
%        creating the damage function based on MDR_fun (e.g. 0:0.5:10).
%        The second last values will be set to the last value in order to
%        ensure a maximum MDR value.
%   params_calibration: parameters for the calibration. A structure with fields:
%       type (string): cost function, one of:
%           'AED2':  "Annual Expected Damage": result is the squared
%                   difference of mean year damage 
%           'R2':   (DEFAULT) result is R^2 (= the mean of squared differences of
%                   year damages of emdat and climada for each specific historical year).
%           'R4':   result is R^4 (= the sum of ^4 differences of
%                   year damages of emdat and climada for each specific historical year).
%           'R':    result is R (= the sum of the absolute differences of
%                   year damages of emdat and climada for each specific historical year).
%           'dlog2':as R2 but with log: result is R (= the mean of the squared differences of
%                   year log damages of emdat and climada for each specific historical year).
%           'dabslog':as R but with log (= the mean of the absolute differences of
%                   yearly log damages of emdat and climada for each specific historical year).
%           'RP':   "Return Period": as AED but for different return
%                   periods with weights (not implemented yet) - only makes
%                   sense for long time series
%       MM_how (string): how to deal with Multi-Model hazard sets, one of:
%           'MMM':  Multi-Model Mean damage estimate vs observated damages (default).
%           'MMMed':Multi-Model Median damage estimate vs observated damages.
%       step_tolerance: parameter step tolerance for patternsearch
%           algorithm. Default=0.001.
%       parallel: =true to use parellelization for the optimization. Default=false.
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180711, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181009, many changes
%    including a new input parameter added 'entity_year'
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181015, additional input
%    parameter params.damFun_xVals.
%   
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables
status=0;

%% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('RegionID','var'),error('Input parameter RegionID is missing');end %
if ~exist('years_range','var'),             years_range=    [1990 2010];end
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_MDR','var'),              params_MDR=          struct;end
if ~exist('params_calibration','var'),      params_calibration=  struct;end


%% check for some parameter fields we need
% params
if ~isfield(params,'entity_folder'),    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'entity_prefix'),    params.entity_prefix='FL1950';end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=0;end
if ~isfield(params,'entity_year'), params.entity_year=0;end
if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end
if ~isfield(params, 'output_folder'),params.output_folder=[climada_global.data_dir filesep 'isimip/results/calibration'];end
% params_MDR
if ~isfield(params_MDR,'remove_years_0emdat'),params_MDR.remove_years_0emdat=0;end
if ~isfield(params_MDR,'remove_years_0YDS'),params_MDR.remove_years_0YDS.do=0;end
if params_MDR.remove_years_0YDS.do
    if ~isfield(params_MDR.remove_years_0YDS,'threshold')
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.threshold required since params_MDR.remove_years_0YDS.do==1 *****')
    end
    if ~isfield(params_MDR.remove_years_0YDS, 'what') && iscell(hazard_list{1})
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.what required since params_MDR.remove_years_0YDS.do==1 and multiple hazard per country are provided *****')
    end
    if ~isfield(params_MDR.remove_years_0YDS, 'min_val'),params_MDR.remove_years_0YDS.min_val=0;end
end
if ~isfield(params_MDR, 'pars_range')
    params_MDR.pars_range = {};
    params_MDR.pars_range{1} = [0.0001 1];
    params_MDR.pars_range{2} = [0.0001 5];
end
if ~isfield(params_MDR,'damFun_xVals'), warning('** warning ** params_MDR.damFun_xVals not set, intensity steps of the default damage function will be used *****');end
% params_calibration
if ~isfield(params_calibration,'type'),params_calibration.type='R2';end
if ~isfield(params_calibration,'MM_how'),params_calibration.MM_how='MMM';end
if ~isfield(params_calibration,'parallel'),params_calibration.parallel=false;end


%% 0) Get variables and paths
% get countries that belong to the region
NatID_RegID_file = [params.RegID_def_folder filesep 'NatID_RegID_isimip_flood.csv'];
NatID_RegID_flood = readtable(NatID_RegID_file);
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
if sum(NatID_RegID_flood.Reg_name == RegionID)==0
    error('no country belonging to the given RegionID, perhaps non-existing RegionID?');
end
countries=NatID_RegID_flood.ISO(NatID_RegID_flood.Reg_name == RegionID);

% define variables and paths
isimip_simround = '2a';
% isimip_data_dir = [climada_global.data_dir filesep 'isimip'];
countries_iso3=cell(length(countries));
for i=1:length(countries)
    [~,countries_iso3{i}] =  climada_country_name(countries{i});
end

%% 0+) prepare output: fill in params for calibrate_MDR_steps including file name to be saved
params_step=struct;
%define filename
filename_calib = ['calib_' RegionID '_' num2str(years_range(1)) '-' num2str(years_range(2)) '_calib-' params_calibration.type '-' params_calibration.MM_how '-step' num2str(params_calibration.step_tolerance)];
filename_haz = ['Haz-Prot' params.hazard_protection '-subMATSIRO' num2str(params.subtract_matsiro)];
filename_ent = ['_Entity-Year' num2str(params.entity_year)];
filename_filter = ['_Filters-emdat' num2str(params_MDR.remove_years_0emdat) '-YDS' num2str(params_MDR.remove_years_0YDS.do)];
if params_MDR.remove_years_0YDS.do
   filename_filter = [filename_filter '-t' num2str(params_MDR.remove_years_0YDS.threshold) '-w' params_MDR.remove_years_0YDS.what '-m' num2str(params_MDR.remove_years_0YDS.min_val)];
end
filename_pars = ['_pars' num2str(params_MDR.pars_range{1}(1)) '-' num2str(params_MDR.pars_range{1}(2)) '-' num2str(params_MDR.pars_range{2}(1)) '-' num2str(params_MDR.pars_range{2}(2))];
filename = [filename_calib '_' filename_haz '_' filename_ent '_' filename_filter '_' filename_pars '.mat'];
% add to params_step (params in calibrate_MDR_steps)
params_step.savefile=[params.output_folder filesep filename];
output_filename=params_step.savefile;
fprintf('Output file will be: %s\n', output_filename);
[temp1,temp2,temp3]=fileparts(output_filename);
params_calibration.write_outfile=[temp1 filesep temp2 '_steps.dat'];
fprintf('Results from each optimization steps will be saved in: %s\n', params_calibration.write_outfile);

%% 1) load entities - N entities for N countries
entity_list=cell(length(countries),1);
for i=1:length(countries)
    country_iso3 = countries_iso3{i};
    entity_file_isimip_i=[params.entity_folder filesep params.entity_prefix strtrim(country_iso3) '_0150as_entity'];
    entity_isimip_i=climada_entity_load(entity_file_isimip_i,1); % try to load, flag to 1 to avoir overwrite
    if isempty(entity_isimip_i)
        error('*** ERROR: entity file not found %s\n\n',entity_file_isimip_i);
    end
    entity_isimip_i.assets.centroid_index = 1:length(entity_isimip_i.assets.centroid_index);
    % check that all analysed years are included
    if sum(ismember(years_range(1):years_range(2), entity_isimip_i.assets.Values_yyyy)) ~= diff(years_range)+1
        error('*** ERROR: not all requested years are available in the entities **\n\n');
    end
    % Subset years
    if params.entity_year
        entity_list{i}=climada_subset_years(entity_isimip_i, 'entity', params.entity_year);
        entity_list{i}.assets.Value=entity_list{i}.assets.Values(1,:);
        entity_list{i}.assets.reference_year = params.entity_year;
    else
        entity_list{i}=climada_subset_years(entity_isimip_i, 'entity', years_range(1):years_range(2));
    end
end


%% 2) load hazards. 46*N hazards, as there are 46 model combinations
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};
hazard_list=cell(length(countries),1);
for i=1:length(countries)
    hazard_list{i} = {};
    ii = 0;
    for j=1:length(ghms)
        for k=1:length(forcings)
            ghm=ghms{j};
            forcing = forcings{k};
            [flddph_filename,~,fld_path] = isimip_get_flood_filename(isimip_simround, ghm, forcing, params.hazard_protection, 'historical');
            flood_filename=[fld_path filesep flddph_filename];
            hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity_list{i},isimip_simround,[0 0],params.subtract_matsiro);
            if exist(hazard_FL_file, 'file')
                hazard_FL=climada_hazard_load(hazard_FL_file);
                hazard_FL.yyyy = double(string(hazard_FL.yyyy));
                try
                    hazard_FL=climada_subset_years(hazard_FL, 'hazard', years_range(1):years_range(2));
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
            else
                fprintf('     * Warning: FL hazard missing, probably inexistent model combination %s\n',flood_filename);
                continue
            end
        end
    end
end


%% 3) load EM-DAT
emdat_list=cell(length(countries),1);
all_years = years_range(1):years_range(2);
for i=1:length(countries)
    country_iso3 = countries_iso3{i};
    % EM-DAT damages corrected by growth exposure only if fixed entities
    % are used.
    em_data_i=emdat_read('',country_iso3,['FL';'F1';'F2'],params.entity_year,0);
    emdat_damage = zeros([length(all_years) 1]);
    % if EM-DAT data available for this country, use, if not leave zeros
    if ~isempty(em_data_i)
        for iy=1:length(all_years)
            ii=find(all_years(iy) == em_data_i.year);
            if ~isempty(ii)
                emdat_damage(iy,1) = sum(em_data_i.damage(ii));
            end
        end
    end
    emdat_list{i}.values=emdat_damage;
    emdat_list{i}.year=all_years;
end


%% 4) Define MDR function shape and parameters
FunctionHandle = str2func('MDR_functional_shape');
MDR_fun = @(x,pars)FunctionHandle(x,pars);


%% 5) fill in params_MDR including file name to be saved
params_MDR.years_range = years_range;
params_MDR.use_YDS = ~params.entity_year;



%% 6) Call calibrate_MDR_steps (TO DO)
opt_pars = calibrate_MDR_steps(entity_list, hazard_list, emdat_list, ...
    MDR_fun, params_step, params_MDR, params_calibration);
fprintf('best set of parameters identified for region %s: %s %s\n', RegionID, opt_pars(1), opt_pars(2));
% save('input_calibrate_MDR_steps.mat','RegionID', 'entity_list', 'hazard_list', 'emdat_list', 'MDR_fun','params_MDR','params_calibration','-v7.3')

status=1;

end


function mdr = MDR_functional_shape(x,pars) %scale,shape)
% climada isimip create MDR function of x (hazard intensity) as a function
%   of n parameters 'pars' (here, scale and shape)
% MODULE:
%   isimip
% NAME:
%   MDR_functional_shape
% PURPOSE:
%   create MDR function of x (hazard intensity) as a function
%   of scale and shape.
%   next call: isimip...
% CALLING SEQUENCE:
%   mdr_func=MDR_functional_shape([scale,shape])
% EXAMPLE:
%   mdr_func=MDR_functional_shape([0.5,0.5])
% INPUTS:
%   scale: scale parameter of function scale*(1-exp(-shape*x))
%   shape: scale parameter of function scale*(1-exp(-shape*x))
% OPTIONAL INPUT PARAMETERS:
%   none.
% OUTPUTS:
%   mdr: function of x which returns scale*(1-exp(-shape*x))
%
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180911, initial
%-

if length(pars)~=2
    error('** ERROR ** 2 parameters should be provided *****')
end

mdr=pars(1)*(1-exp(-pars(2)*x));

end

