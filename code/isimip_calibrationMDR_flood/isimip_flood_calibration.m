function [status,output_eval_filename,output]=isimip_flood_calibration(RegionID,years_range,params,params_MDR,params_calibration,params_computation)
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
%       calib_options: options on the calibration method to use, a struct with fields:
%           method: one of
%               'patternsearch' (default): use pattern search algorithm
%               'regular_sampling' : regularly sample the parameter space
%           params: parameters for the method.
%               For patternsearch: structure with fields
%                   random: if =1, starting points are chosen randomly,
%                       otherwise they are linearly distributed over the
%                       parameter ranges. Default=0.
%                   nstart: number of starting points (default=1). 
%                       If random==1, the total number of starting points;
%                       if random==0, the number of starting points in each
%                       parameter dimension.
%                   InitialMeshSize: parameter in patternsearch (default=.25)
%                   step_tolerance: parameter step tolerance for patternsearch
%                       algorithm. Default=0.001.
%               For regular_sampling: structure with fields
%                   n_per_dim: number of samples per parameter (total
%                       number of simulations will then be (n_per_dim)^n
%                       where n is the number of parameters). Default=5.
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
%           'RTarea': area of the difference between return time curves in
%                   log10-log10 space, exluding cases with a return value
%                   of 0 in either observed or modelled damages.
%                   Area for underestimated damages is scaled according to
%                   underestimation_factor.
%           'RP':   "Return Period": as AED but for different return
%                   periods with weights (not implemented yet) - only makes
%                   sense for long time series
%       MM_how (string): how to deal with Multi-Model hazard sets, one of:
%           'MMM':  Multi-Model Mean damage estimate vs observated damages (default).
%           'MMMed':Multi-Model Median damage estimate vs observated damages.
%       underestimation_factor: factor by which to scale up cases where climada
%           underestimates damages BEFORE evaluating the cost function
%           (default =1, i.e. no scaling). This allows to account for fact that
%           observated damages are usually rather a lower bound of damages. For
%           instance, a value of 2 means that the contribution of these
%           cases (where simulated damages < observed damages) to the cost
%           function will be worth double of the same difference for other
%           cases. Note that as the factor is applied before evaluating the
%           cost function, for e.g. cost function type 'R2' underestimates
%           would be worth 4 times more with a factor of 2 (use sqrt(2) to
%           actually have an effect of 2).
%   params_computation: A structure with fields:
%       years_range: Range of years to be used for the computation of
%           damages with the calibrated damage function (using
%           isimip_compute_calibrated). Default = [1971 2010]
%       do: if =1, does the computation of damages with calibrated
%           parameters. Otherwise, does not do it. Default=1 if
%           params_calibration.calib_options.method=='patternsearch',
%           otherwise default=0.
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
if ~exist('years_range','var'),             years_range=    [1992 2010];end
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_MDR','var'),              params_MDR=          struct;end
if ~exist('params_calibration','var'),      params_calibration=  struct;end
if ~exist('params_computation','var'),      params_computation=  struct;end


%% check for some parameter fields we need
% params
if ~isfield(params,'entity_folder'),    params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'entity_prefix'),    params.entity_prefix='FL1950';end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=0;end
if ~isfield(params,'entity_year'), params.entity_year=0;end
if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end
if ~isfield(params, 'output_folder'),params.output_folder=[climada_global.data_dir filesep 'isimip/results/calibration'];end
if ~isfield(params,'keep_countries_0emdat'), params.keep_countries_0emdat=2;end
if ~ismember(params.keep_countries_0emdat, 0:2),error('** ERROR ** input parameter params.keep_countries_0emdat should be one of 0,1,2 *****');end
if ~isfield(params,'verbose_excluded_countries'),params.verbose_excluded_countries=0;end
% params_MDR
if ~isfield(params_MDR,'remove_years_0emdat'),params_MDR.remove_years_0emdat=0;end
if ~isfield(params_MDR,'remove_years_0YDS'),params_MDR.remove_years_0YDS.do=0;end
if params_MDR.remove_years_0YDS.do
    if ~isfield(params_MDR.remove_years_0YDS,'threshold')
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.threshold required since params_MDR.remove_years_0YDS.do==1 *****')
    end
    if ~isfield(params_MDR.remove_years_0YDS, 'what')
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.what required since params_MDR.remove_years_0YDS.do==1 and multiple hazard per country are expected *****')
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
if ~isfield(params_calibration,'calib_options'),params_calibration.calib_options=struct;end
if ~isfield(params_calibration.calib_options,'method'),params_calibration.calib_options.method='patternsearch';end
if ~isfield(params_calibration.calib_options,'params'),params_calibration.calib_options.params=struct;end
switch params_calibration.calib_options.method
    case 'patternsearch'
        if ~isfield(params_calibration.calib_options.params,'random'),params_calibration.calib_options.params.random=0;end
        if ~isfield(params_calibration.calib_options.params,'nstart'),params_calibration.calib_options.params.nstart=1;end
        if ~isfield(params_calibration.calib_options.params,'InitialMeshSize'),params_calibration.calib_options.params.InitialMeshSize=0.25;end
        if ~isfield(params_calibration.calib_options.params,'step_tolerance'),params_calibration.calib_options.params.step_tolerance=0.01;end
    case 'regular_sampling'
        if ~isfield(params_calibration.calib_options.params,'n_per_dim'),params_calibration.calib_options.params.n_per_dim=5;end
    otherwise
        error('Input field params_calibration.calib_options.method is not valid')
end
if ~isfield(params_calibration,'underestimation_factor'),params_calibration.underestimation_factor=1;end
% params_computation
if ~isfield(params_computation,'years_range'),params_computation.years_range=[1971 2010];end
if ~isfield(params_computation,'do'),params_computation.do=strcmp(params_calibration.calib_options.method,'patternsearch');end


%% 1) Get variables and paths, check countries
isimip_simround = '2a';
% get countries that belong to the region, those used for calibration and
% those used only for 'extrapolation'
NatID_RegID_file=[params.RegID_def_folder filesep 'NatID_RegID_isimip_flood_filtered_' num2str(years_range(1)) '-' num2str(years_range(2)) '.csv'];
if exist(NatID_RegID_file,'file')
    NatID_RegID_flood = readtable(NatID_RegID_file);
else
    [flag,NatID_RegID_flood] = isimip_FL_select_countries(NatID_RegID_file,[],years_range);
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
% countries excluded from calibration but retained for application
calib_out = (NatID_RegID_flood.in_calib < params.keep_countries_0emdat) & NatID_RegID_flood.in_extrap;
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

countries_iso3_all = NatID_RegID_flood.ISO(~fully_out);
countries_iso3 = NatID_RegID_flood.ISO(countries_keep)';
countries_iso3_extrapolation = NatID_RegID_flood.ISO((~fully_out) & calib_out)';
[~,countries_calib_inds] = intersect(countries_iso3_all,countries_iso3,'stable');

fprintf('Countries selected for calibration:\n%s\n\n',strjoin(countries_iso3,', '))
fprintf('Additional countries selected for application (extrapolation):\n%s\n\n',strjoin(countries_iso3_extrapolation,', '))
clear countries_keep calib_out fully_out;

%% 2) prepare output: fill in params for calibrate_MDR_steps including file name to be saved
params_step=struct;
%define filename
switch params_calibration.calib_options.method
    case 'patternsearch'
        if params_calibration.calib_options.params.random
            filename_calib_method=[params_calibration.calib_options.method '-rand' num2str(params_calibration.calib_options.params.nstart)];
        else
            filename_calib_method=[params_calibration.calib_options.method '-reg' num2str(params_calibration.calib_options.params.nstart)];
        end
        filename_calib_method=[filename_calib_method '-mesh' num2str(params_calibration.calib_options.params.InitialMeshSize) '-step' num2str(params_calibration.calib_options.params.step_tolerance)];
    case 'regular_sampling'
        filename_calib_method=[params_calibration.calib_options.method '-' num2str(params_calibration.calib_options.params.n_per_dim)];
    otherwise
        error('Input field params_calibration.calib_options.method is not valid')
end
filename_calib = ['calib_' RegionID '_' num2str(years_range(1)) '-' num2str(years_range(2)) '_' ...
    params_calibration.type '-uf' num2str(params_calibration.underestimation_factor) '-' params_calibration.MM_how '-' filename_calib_method];
filename_haz = ['Haz-Prot' params.hazard_protection '-subMATSIRO' num2str(params.subtract_matsiro)];
filename_ent = ['Entity-Year' num2str(params.entity_year)];
filename_filter = ['Filters-Country' num2str(params.keep_countries_0emdat) '-emdat' num2str(params_MDR.remove_years_0emdat) '-YDS' num2str(params_MDR.remove_years_0YDS.do)];
if params_MDR.remove_years_0YDS.do
   filename_filter = [filename_filter '-t' num2str(params_MDR.remove_years_0YDS.threshold) '-w' params_MDR.remove_years_0YDS.what '-m' num2str(params_MDR.remove_years_0YDS.min_val)];
end
filename_pars = ['pars' num2str(params_MDR.pars_range{1}(1)) '-' num2str(params_MDR.pars_range{1}(2)) '-' num2str(params_MDR.pars_range{2}(1)) '-' num2str(params_MDR.pars_range{2}(2))];
filename = [filename_calib '_' filename_haz '_' filename_ent '_' filename_filter '_' filename_pars '.mat'];
% add to params_step (params in calibrate_MDR_steps)
params_step.savefile=[params.output_folder filesep filename];
[temp1,temp2,~]=fileparts(params_step.savefile);
params_calibration.write_outfile=[temp1 filesep temp2 '_steps.dat'];
fprintf('Results from each optimization steps will be saved in: %s\n', params_calibration.write_outfile);
output_eval_filename=[temp1 filesep temp2 '_eval.csv'];

%% 3) load entities , hazards and EM-DAT
[entity_list,hazard_list,emdat_list] = load_entity_hazard_emdat(countries_iso3,years_range,params,isimip_simround);


%% 4) Define MDR function shape and parameters
MDR_fun=@(x,pars)pars(1)*(1-exp(-pars(2)*x));


%% 5) fill in params_MDR including file name to be saved
params_MDR.years_range = years_range;
params_MDR.use_YDS = ~params.entity_year;



%% 6) Call calibrate_MDR_steps (TO DO)
[ opt_pars, years_i_in] = calibrate_MDR_steps(entity_list, hazard_list, emdat_list, ...
    MDR_fun, params_step, params_MDR, params_calibration);
fprintf('best set of parameters identified for region %s: scale=%g , shape=%g\n', RegionID, opt_pars(1), opt_pars(2));
status=1;

%% 7) Compute damages per country and year for each combination using the identified optimal parameter combination

% first, load data for all years
if ~params_computation.do
    fprintf('Evaluation of isimip_compute_calibrated not requested\n')
else
    fprintf('Loading data for isimip_compute_calibrated...\n')
    [entity_list_comp,hazard_list_comp,emdat_list_comp] = load_entity_hazard_emdat(countries_iso3_all,params_computation.years_range,params,isimip_simround);
    fprintf('Data loading completed for isimip_compute_calibrated\n')
    
    % prepare input
    % get years_i_in (year/country used for calibration)
    all_years_comp = params_computation.years_range(1):params_computation.years_range(2);
    [~,years_comp_intersect] = intersect(all_years_comp,years_range(1):years_range(2));
    years_i_in_comp = false(length(all_years_comp),length(countries_iso3_all));
    years_i_in_comp(years_comp_intersect,countries_calib_inds)=years_i_in;
    damFun = climada_damagefunctions_generate_from_fun(params_MDR.damFun_xVals, MDR_fun, opt_pars);
    
    % computation
    [ status2, output ] = isimip_compute_calibrated(entity_list_comp, ...
        hazard_list_comp, emdat_list_comp, ...
        params_computation.years_range, damFun, params_MDR.use_YDS, years_i_in_comp);
    status=status+status2;
    
    % check status, print, write out table
    if status2
        fprintf('Evaluation of isimip_compute_calibrated will be saved in: %s\n', output_eval_filename);
        writetable(output,output_eval_filename);
    else
        fprintf('Evaluation of isimip_compute_calibrated failed\n')
    end
end

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
