function result=calibrate_MDR_steps(RegionID, entity_list, hazard_list, emdat_list, MDR_func, params)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   calibrate_MDR_steps
% PURPOSE:
%   for a given Region (group of countries) and the corresponding lists (one
%   element per country) of entities, hazard, EM-DAT damage data, calibrate
%   the parameters of MDR_func to best match EM-DAT damages. Usually called
%   from 'isimip_flood_calibration'.
%
% CALLING SEQUENCE:
%   [status,output_filename]=calibrate_MDR_steps(RegionID, entity_list, hazard_list, emdat_list, MDR_func, params)
% EXAMPLE:
%   RegionID='NAM';
%   years_range=[1990 2000];
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% INPUTS:
%   RegionID: country name (full name or ISO3)
%   entity_list: a list of entities (cell array), one per country
%   hazard_list: a list of hazard (cell array), or several per country. If
%      several, a cell array of cell arrays.
%   emdat_list: a list of EM-DAT damages (cell array) containing the fields
%      'year' and 'values'.
%   MDR_func: a function of x (hazard intensity) and two parameters to be
%      calibrated (e.g., scale and shape).
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%     entity_folder: the folder where the entities are located (default:
%        [climada_global.data_dir filesep 'isimip/entities'] ).
%     RegID_def_folder: the folder where the file listing countries per
%        Region is located.
%     entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%        entity_prefix='Try1' will result in Try1_DEU_0150as.mat
%     hazard_protection: one of 'flopros' (default), '0', '100'
%     subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%        fraction from the data. Default =0.
%     years_range: range of years to be included (e.g. [1990 2010])
%     remove_years_0emdat: determines how to filter years with zero EM-DAT
%        damages (e.g., no damage).
%        0=no filter (all years retained);
%        1=remove years for which the sum of EM-DAT damages over all
%        countries is 0;
%        2=for each year, remove countries with zero EM-DAT damages.
%     remove_years_0YDS: list defining whether to exclude years with zero
%        simulated damage as follows:
%        do: =1 to exclude years with zero simulated damage. Default =0.
%        threshold: value above which (>) damage is considered to occur to
%           filter out years (required if do==1).
%        what: determines how to choose if several hazards are provided per
%           country (possible choices: 'any', 'all', 'mean', 'median'). For
%           instance, if 'any', then if any hazard set gives 0 damage the
%           year-country event is removed.
%        min_val: damage value below which damage is not considered to
%           occur (default =0).
%     use_YDS: =1 to use yearly-varying assets. Default =0.
%     pars_range: range of parameter values (cell array, each element as
%        [min_val max_val].
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180911, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181009, split into
%    sub-functions and progress in implementation
%   
%-

% https://ch.mathworks.com/help/gads/examples/constrained-minimization-using-pattern-search.html

%%


if ~exist('RegionID','var'),error('Input parameter RegionID is missing');end %
if ~exist('entity_list','var'),error('Input parameter entity_list is missing');end
if ~exist('hazard_list','var'),error('Input parameter hazard_list is missing');end
if ~exist('emdat_list','var'),error('Input parameter emdat_list is missing');end
if ~exist('MDR_func','var'),error('Input parameter MDR_func is missing');end
if ~exist('params','var'),                  params=              struct;end

% check for some parameter fields we need
if ~isfield(params,'entity_folder'),    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'entity_prefix'),    params.entity_prefix='FL1950';end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=0;end
if ~isfield(params,'years_range'),params.years_range=[min(emdat_list{1}.year) max(emdat_list{1}.year)];end
if ~isfield(params,'remove_years_0emdat'),params.remove_years_0emdat=0;end
if ~isfield(params,'remove_years_0YDS'),params.remove_years_0YDS.do=0;end
if params.remove_years_0YDS.do
    if ~isfield(params.remove_years_0YDS,'threshold')
        error('** ERROR ** Input parameter params.remove_years_0YDS.threshold required since params.remove_years_0YDS.do==1 *****')
    end
    if ~isfield(params.remove_years_0YDS, 'what') && iscell(hazard_list{1})
        error('** ERROR ** Input parameter params.remove_years_0YDS.what required since params.remove_years_0YDS.do==1 and multiple hazard per country are provided *****')
    end
    if ~isfield(params.remove_years_0YDS, 'min_val'),params.remove_years_0YDS.min_val=0;end
end
if ~isfield(params,'savedir'),          error('Input parameter params.savedir is missing');end
if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end
if ~isfield(params,'use_YDS'),params.use_YDS=0;end
if ~isfield(params, 'pars_range')
    params.pars_range = {};
    params.pars_range{1} = [0 1];
    params.pars_range{2} = [0 5];
end


% determine useful parameters
country_list=cellfun(@(x) x.assets.admin0_ISO3,entity_list, 'UniformOutput', 0);
n_countries = length(emdat_list);

% put hazard_list in double-cell mode in case only 1 hazard was provided.
if ~iscell(hazard_list{1})
    for i = 1:n_countries
        hazard_temp=hazard_list{i};
        hazard_list{i}={};
        hazard_list{i}{1}=hazard_temp;
        clear hazard_temp;
    end
end







% if ~exist('calibrate_countries','var'),calibrate_countries=[];end
% if ~exist('hazard_filename','var'),hazard_filename=[];end
% if ~exist('number_free_parameters','var'),number_free_parameters=[];end
% if ~exist('years_considered','var'),years_considered=[];end



% TCbasins = {'NAT' 'NWP' 'NIN' 'SIN' 'SWP' 'NNN'};
% TCbasinIDs=[ 1     2     3     4     5     0];
% if isempty(TCBasinID),TCBasinID=1;end
% if isempty(resolution),resolution=300;end
% if isempty(calibrate_countries),calibrate_countries=0;end
% if isempty(hazard_filename),hazard_filename=['GLB_0360as_',peril_ID,'_hist'];end
% if isempty(number_free_parameters),number_free_parameters=2;end



delta_shape_parameter = 49; % v_half - v_threshold (m/s)
encode = 0;
optimizerType='R2';
remove_0_YDS_years = 1;
%optimizerType='R';
%optimizerType='logR';

full_parameter_search = 1;
fminconSwitch = 0;
save_output = 1;
force_overwrite_output = 0;

%

savedir = [climada_global.data_dir filesep 'output_sam'];
if ~exist(savedir,'dir') && save_output
    system(['mkdir ' savedir]);
end



% CHECK TIME PERIOD
same_period = check_time_period(params.years_range, emdat_list, entity_list, hazard_list, params.use_YDS);
% 1) if years_range is not the same as all_years_avail, print message
% and abort
if ~same_period
    error('** ERROR ** mismatch in the time period *****')
end
% Now we know that all input have consistent time period:
% params.years_range and the entities, hazard and emdat data lists.


% are years without damage in EM-DAT and/or without simulated damage to be removed?
[entity_list, hazard_list, emdat_list] = filter_years(params.remove_years_0emdat, params.remove_years_0YDS, params.years_range, entity_list, hazard_list, emdat_list, use_YDS)


    


% CODE WRITING HERE


%% set boundaries and starting values
% IDEALLY CHOOSE THE MIDDLE OF THE RANGE OF VALUES, AND DEFINE RANGE
scale_0=0.5;
shape_0=2.5;
x0=[scale_0 shape_0];
bounds.lb=[0.0001 1];
bounds.up=[0.0001 5];
        

% lower bounds are normalized to 1
norm.lb = ones(size(bounds.lb));
% upper bounds are normalized to 2
norm.ub = norm.lb+1;
% normalization of the starting values
norm.x0 = (x0-bounds.lb)./(bounds.ub-bounds.lb) .* (norm.ub-norm.lb) + norm.lb;


%%
% define anonymous function with input factor x (parameters of the damage
% function):
% delta_shape_parameter is specific to TCs, can be removed.
fun = @(x)calibrate_TC_DF_emdat_region(x,delta_shape_parameter, entity, hazard, em_data_region, norm, bounds,optimizerType); % sets all inputvar for the function except for x, use normalized x.


%parpool('local_small')
if full_parameter_search
    
    options = optimoptions('patternsearch','UseParallel',false,...
    'UseCompletePoll', true, 'UseVectorized', false,...
    'MaxFunctionEvaluations',1200,'Display','iter',...
    'Cache','on','InitialMeshSize',.25,...
    'PollMethod','GPSPositiveBasis2N','StepTolerance',0.001);
    tic
    [x_result,fval] = patternsearch(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],options);
    toc
    % convert normalized value of calibrated parameters to their 'real'
    % values
    result.region=(x_result-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb
    fval
    if save_output && ~calibrate_countries
        
        save_file_name=[savedir filesep regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_' peril_ID '_decay_region_calibrate_litpop_gdp_' num2str(number_free_parameters) '-' num2str(value_mode) '-' num2str(resolution) '.mat'];
        
        while (~force_overwrite_output && exist(save_file_name,'file')),save_file_name=strrep(save_file_name,'.mat','_.mat');end % avoid overwriting
        
        save(save_file_name,'result','fval','resolution','years_considered','-v7.3');
        
    end
end
% clear entity hazard

% calibration of single countries with more than or equal to N_min damage years
if calibrate_countries
    N_min = 1;
    
    clear entity
    em_data_country.year = em_data_region.year;
    result.result_c = NaN*ones(length(country_list),length(x0));
    result.fval_c = NaN*ones(length(country_list),1);
    N_damageyears = result.fval_c;
    options = optimoptions('patternsearch','UseParallel',false,...
        'UseCompletePoll', true, 'UseVectorized', false,...
        'MaxFunctionEvaluations',1200,'Display','iter',...
        'Cache','on','InitialMeshSize',.25,...
        'PollMethod','GPSPositiveBasis2N','StepTolerance',0.001);
    
    for i_c = 1:length(country_list)
        %%
        N_damageyears(i_c) = sum(em_data_region.(['damage_' country_list{i_c}])>0);
        disp([country_list{i_c} ', N=' num2str(N_damageyears(i_c))]);
        if N_damageyears(i_c) >= N_min
            entity = climada_entity_load([country_list{i_c} '_GDP_LitPop_BM2016_' num2str(resolution) 'arcsec_ry' num2str(reference_year) '.mat']);
            em_data_country.damage = em_data_region.(['damage_' country_list{i_c}]);
            fun = @(x)calibrate_TC_DF_emdat_region(x,delta_shape_parameter, entity, hazard, em_data_country, norm, bounds,optimizerType);
            if sum(entity.assets.Value)>0
                [x_result,result.fval_c(i_c)] = patternsearch(fun,norm.x0,[],[],[],[],norm.lb,norm.ub,[],options);
            
                result.result_c(i_c,:)=(x_result-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb
            else
                %result.result_c(i_c,:)=[NaN NaN];
                warning('entity contains no values');
            end
        end
    end
    %%
    if save_output
        save_file_name = [savedir filesep regions.mapping.TCBasinName{find(regions.mapping.TCBasinID==TCBasinID,1)}...
            '_' peril_ID '_decay_region_countries_calibrate_litpop_gdp_' num2str(number_free_parameters) '-' num2str(value_mode) '-' num2str(resolution) '.mat'];
        
        while (~force_overwrite_output && exist(save_file_name,'file')),save_file_name=strrep(save_file_name,'.mat','_.mat');end % avoid overwriting
        
        save(save_file_name,'country_list','result','fval','resolution','years_considered','-v7.3');
        save_file_name
    end
end


if on_cluster, cd /cluster/home/eberenzs/ ;exit; end;


end

function same_period = check_time_period(years_range, emdat_list, entity_list, hazard_list, use_YDS)
% CHECK TIME PERIOD
same_period=1;
n_countries = length(emdat_list);
all_years_avail = emdat_list{1}.year';
% 1) if years_range is not the same as all_years_avail, print message
% and abort
if ~isequal(all_years_avail, (years_range(1):years_range(2))')
    if years_range(1) < min(all_years_avail)
        fprintf('years_range(1) not covered by input data - aborted\n');
    elseif years_range(1) > max(all_years_avail)
        fprintf('years_range(2) not covered by input data - aborted\n');
    else
        fprintf('years_range exceeded by input data - aborted\n');
    end
    same_period=0;
    return
end
% 2) Check that entity_list, hazard_list and emdat_list include the same years
for i = 1:n_countries
    if use_YDS
        % entity values and years only checked if YDS is used; otherwise it
        % does not matter (only entity_list{i}.assets.Value will be used)
        if ~isequal(entity_list{i}.assets.Values_yyyy, all_years_avail)
            fprintf('Not all same years in entity_list{%i}\n',i);
            same_period=0;
            return
        end
    end
    if ~isequal(emdat_list{i}.year', all_years_avail)
        fprintf('Not all same years in emdat_list{%i}\n',i);
        same_period=0;
        return
    end
    for j = 1:length(hazard_list{i})
        if ~isequal(hazard_list{i}{j}.yyyy, all_years_avail)
            fprintf('Not all same years in hazard_list{%i}{%j}\n',i,j);
            same_period=0;
            return
        end
    end
end
% Now we know that all input have consistent time period:
% params.years_range and the entities, hazard and emdat data lists.
end



function [entity_list, hazard_list, emdat_list] = filter_years(remove_years_0emdat, remove_years_0YDS, years_range, entity_list, hazard_list, emdat_list, use_YDS)
% function to filter out years, depending on remove_years_0emdat and remove_years_0YDS
country_list=cellfun(@(x) x.assets.admin0_ISO3,entity_list, 'UniformOutput', 0);
n_countries = length(emdat_list);
% array of years to keep per country - initialized as 'keep all'
all_years_in = years_range(1):years_range(2);
years_in = true(length(all_years_in), n_countries);
% remove years with 0 damage in EM-DAT
if remove_years_0emdat
    % find out non-zero EM-DAT damages
    for i = 1:n_countries
        years_in(:,i) = emdat_list{i}.values > 0;
    end
    % if remove years not by country but only years without any EM-DAT
    %     damage in any country
    if remove_years_0emdat == 1
        years_in_temp = sum(years_in, 2)>0;
        for i = 1:n_countries
            years_in(:,i) = years_in_temp;
        end
        fprintf('A total of %i year(s) are considered for calibration (between %i and %i)\n',...
            sum(years_in_temp),min(all_years_in(years_in_temp)),max(all_years_in(years_in_temp)));
    elseif remove_years_0emdat == 2
        list_to_print = ['Total of year(s) considered for calibration for each country (' num2str(min(all_years_in)) '-' num2str(max(all_years_in)) ') :'];
        for i = 1:n_countries
            list_to_print=[list_to_print ' ' country_list{i} ': ' num2str(sum(years_in(:,i)))];
        end
        fprintf([list_to_print, '\n'])
    else
        error('** ERROR ** unexpected value in params.remove_years_0emdat *****')
        return
    end
end
% remove years without simulated damage
if remove_years_0YDS.do
    threshold_0YDS=remove_years_0YDS.threshold;
    % identify peril_ID from the hazard
    peril_ID = hazard_list{1}{1}.peril_ID;
    % empty matrix with computed damages per year, country and hazard set
    EDS_all = zeros(length(all_years_in), n_countries, length(hazard_list{1}));
    for i = 1:n_countries
        DamFunID=unique(entity_list{i}.assets.DamageFunID(~isnan(entity_list{i}.assets.DamageFunID)));
        % check that only 1 DamageFunID is present
        if length(DamFunID)>1,warning('Several DamageFunID are present in entities, only the first one will be changed');end
        % select DamFun to modify
        choose_DF = find(strcmp(entity_list{i}.damagefunctions.peril_ID,peril_ID) .* (entity_list{i}.damagefunctions.DamageFunID==DamFunID(1)));
        % set PPA to 1 for all intensities
        entity_list{i}.damagefunctions.PAA(choose_DF) = 1;
        % select the x value corresponding to the threshold
        df_intensity = entity_list{i}.damagefunctions.Intensity(choose_DF);
        df_intensity_threshold = find(df_intensity >= threshold_0YDS);
        if isempty(df_intensity_threshold) || df_intensity_threshold(1)==length(df_intensity)
            error('** ERROR ** No intensity value below its maximum exceeds the threshold for 0YDS in damagefunctions - needs fixing *****')
        elseif df_intensity_threshold(1)==1
            warning('** WARNING ** intensity threshold for 0YDS in damagefunctions is the minimum of the intensity scale - this will only work correctly if the hazard intensity cannot be below this value (e.g., 0 if negative values do not exist) *****')
        end
        df_intensity_threshold = df_intensity_threshold(1);
        % set all MDD values up to threshold to 0
        entity_list{i}.damagefunctions.MDD(choose_DF(1:df_intensity_threshold)) = 0;
        % set intensity in damage function to threshold at the single point
        % to ensure exact result
        entity_list{i}.damagefunctions.Intensity(choose_DF(df_intensity_threshold)) = threshold_0YDS;
        % set all other MDD values to 1
        entity_list{i}.damagefunctions.MDD(choose_DF((df_intensity_threshold+1):end)) = 1;
        for j=1:length(hazard_list{i})
            % now compute damages for this country/event
            if use_YDS
                yds_params.silent_mode=2;
                EDS_ij = isimip_YDS_calc(entity_list{i}, hazard_list{i}{j},yds_params);
            else
                EDS_ij = climada_EDS_calc(entity_list{i},hazard_list{i}{j},'',0,2); % the damage calculation
            end
            EDS_all(:,i,j) = EDS_ij.damage;
        end
    end
    EDS_all=EDS_all*EDS_ij.currency_unit;
    % now apply criterion
    switch remove_years_0YDS.what
        case 'any' % remove cases with any hazard set giving 0 damage
            EDS_keep = sum(EDS_all<remove_years_0YDS.min_val,3)>0;
        case 'all' % remove cases when all hazard sets give 0 damage
            EDS_keep = sum(EDS_all>=remove_years_0YDS.min_val,3)==0;
        case 'mean' % remove cases when the mean damage from hazard sets is == 0
            EDS_keep = mean(EDS_all,3)>remove_years_0YDS.min_val;
        case 'median' % remove cases when the mean damage from hazard sets is == 0
            EDS_keep = median(EDS_all,3)>remove_years_0YDS.min_val;
        otherwise
            error('** ERROR ** unexpected value in remove_years_0YDS.what *****')
    end
    % now add to years_in
    years_in = years_in & EDS_keep;
    %     % (2) calculate YDS
    %     EDS = climada_EDS_calc(entity,hazard,[],[],2);
    %     % SAM: this is just to sum damages per year - not relevant for me.
    %     [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
    %     [LIA,LOCB] = ismember(em_data_region.year,YDS.yyyy);
    %     % exclude years with 0-damage-entries:
    %     year_i_vector = find(LIA);
    %     year_i_vector(em_data_region.damage(year_i_vector)==0)=[];
    %     % years without any event in climada (hazard set)
    %     year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];
    %     em_data_region.damage = em_data_region.damage(year_i_vector);
    %     em_data_region.year = em_data_region.year(year_i_vector);
    %     for i_c = 1:length(country_list)
    %         em_data_region.(['damage_' country_list{i_c}])=em_data_region.(['damage_' country_list{i_c}])(year_i_vector);
    %     end
    %     % figure(9); hold on;plot(em_data_region.year, em_data_region.damage,'rx');hold off;
    %     N_years = length(find(em_data_region.damage>0));
    %     disp(N_years)
end

% remove years from all data
for i = 1:n_countries
    % entity is subset only if multiple years are available and will be
    % used - otherwise no point subsetting entity
    if use_YDS
        entity_list{i}=climada_subset_years(entity_list{i}, 'entity', all_years_in(years_in(:,i)));
    end
    emdat_list{i}=climada_subset_years(emdat_list{i}, 'obs', all_years_in(years_in(:,i)));
    for j = 1:length(hazard_list{i})
        hazard_list{i}{j}=climada_subset_years(hazard_list{i}{j}, 'hazard', all_years_in(years_in(:,i)));
    end
end

end