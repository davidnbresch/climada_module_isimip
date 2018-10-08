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
%     pars_range: range of parameter values (cell array, each element as
%        [min_val max_val].
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180711, initial
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
if ~isfield(params,'years_range'),      params.years_range=emdat_list{1}.year;end
if ~isfield(params,'remove_years_0emdat'),params.remove_years_0emdat=0;end
if ~isfield(params,'savedir'),          error('Input parameter params.savedir is missing');end
if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end
if ~isfield(params, 'pars_range')
    params.pars_range = {};
    params.pars_range{1} = [0 1];
    params.pars_range{2} = [0 5];
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

% determine useful parameters
country_list=cellfun(@(x) x.assets.admin0_ISO3,entity_list, 'UniformOutput', 0);
n_countries = length(emdat_list);


% CHECK TIME PERIOD
same_period = check_time_period(params.years_range, emdat_list, entity_list, hazard_list);
% 1) if years_range is not the same as all_years_avail, print message
% and abort
if ~same_period
    error('** ERROR ** mismatch in the time period *****')
end
% Now we know that all input have consistent time period:
% params.years_range and the entities, hazard and emdat data lists.


% are years without damage in EM-DAT to be removed?
if params.remove_years_0emdat
    all_years_in = years_range(1):years_range(2);
    years_in = false(length(all_years_in), n_countries);
    % find out non-zero EM-DAT damages
    for i = 1:n_countries
        years_in(:,i) = emdat_list{i}.values > 0;
    end
    % if remove years not by country but only years without any EM-DAT
    %     damage in any country
    if params.remove_years_0emdat == 1
        years_in_temp = sum(years_in, 2)>0;
        for i = 1:n_countries
            years_in(:,i) = years_in_temp;
        end
        fprintf('A total of %i year(s) are considered for calibration (between %i and %i)\n',...
            sum(years_in_temp),min(all_years_in(years_in_temp)),max(all_years_in(years_in_temp)));
    elseif params.remove_years_0emdat == 2
        list_to_print = ['Total of year(s) considered for calibration for each country (' num2str(min(all_years_in)) '-' num2str(max(all_years_in)) ') :'];
        for i = 1:n_countries
            list_to_print=[list_to_print ' ' country_list{i} ': ' num2str(sum(years_in(:,i)))];
        end
        fprintf([list_to_print, '\n'])
    else
        error('** ERROR ** unexpected value in params.remove_years_0emdat *****')
        return
    end
    % remove years from all data
    for i = 1:n_countries
        entity_list{i}=climada_subset_years(entity_list{i}, 'entity', all_years_in(years_in(:,i)));
        emdat_list{i}=climada_subset_years(emdat_list{i}, 'obs', all_years_in(years_in(:,i)));
        for j = 1:length(hazard_list{i})
            hazard_list{i}{j}=climada_subset_years(hazard_list{i}{j}, 'hazard', all_years_in(years_in(:,i)));
        end
    end

end
    


% CODE WRITING HERE


% FILTER OUT YEARS WITHOUT ANY SIMULATED DAMAGE IN EACH REGION (NO EVENT),
% TO DO SO TAKE A SIMPLE DAMAGE FUNCTION AND DETECT YEARS WITH ANY DAMAGE.
% REMOVE OTHER YEARS
if remove_0_YDS_years
    % figure(9); hold on;plot(em_data_region.year, em_data_region.damage,'bo');hold off;
    v_threshold = 15;
    v_half = 20;
    scale = 1;
    choose_DF = strcmp(entity.damagefunctions.peril_ID,'TC') .* (entity.damagefunctions.DamageFunID==1);

    % set PPA to 1 for all intensities
    entity.damagefunctions.PAA(choose_DF==1) = 1;
    % set MDD with x(1)= v_half [m/s] % 
    % based on Elliott et al. 2015, similar also Sealy & Strobl et al. 2017 and Emanuel, 2011
    VV_temp = max((entity.damagefunctions.Intensity(choose_DF==1) - v_threshold),0)/(v_half-v_threshold);
    entity.damagefunctions.MDD(choose_DF==1)= VV_temp.^3 ./ (1+VV_temp.^3);
    % linearly damp MDD by multiplication with scale x(2), 0<x(2)<=1
    entity.damagefunctions.MDD(choose_DF==1)= scale * entity.damagefunctions.MDD(choose_DF==1);
    % (2) calculate YDS
    EDS = climada_EDS_calc(entity,hazard,[],[],2);
    % SAM: this is just to sum damages per year - not relevant for me.
    [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
    [LIA,LOCB] = ismember(em_data_region.year,YDS.yyyy);
    % exclude years with 0-damage-entries:
    year_i_vector = find(LIA);
    year_i_vector(em_data_region.damage(year_i_vector)==0)=[]; 
    % years without any event in climada (hazard set)
    year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];
    em_data_region.damage = em_data_region.damage(year_i_vector);
    em_data_region.year = em_data_region.year(year_i_vector);
    for i_c = 1:length(country_list)
        em_data_region.(['damage_' country_list{i_c}])=em_data_region.(['damage_' country_list{i_c}])(year_i_vector);
    end
    % figure(9); hold on;plot(em_data_region.year, em_data_region.damage,'rx');hold off;
    N_years = length(find(em_data_region.damage>0));
    disp(N_years)
end


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


function same_period = check_time_period(years_range, emdat_list, entity_list, hazard_list)
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
    if ~isequal(entity_list{i}.assets.Values_yyyy, all_years_avail)
        fprintf('Not all same years in entity_list{%i}\n',i);
        same_period=0;
        return
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