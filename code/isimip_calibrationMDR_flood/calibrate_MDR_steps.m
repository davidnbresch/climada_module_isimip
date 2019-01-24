function [ optimal_pars, years_i_in ] = calibrate_MDR_steps(entity_list, hazard_list, ...
    emdat_list, MDR_fun, params, params_MDR, params_calibration)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   calibrate_MDR_steps
% PURPOSE:
%   for a given Region (group of countries) and the corresponding lists (one
%   element per country) of entities, hazard, EM-DAT damage data, calibrate
%   the parameters of MDR_fun to best match EM-DAT damages. Usually called
%   from 'isimip_flood_calibration'.
%
% CALLING SEQUENCE:
%   [status,output_filename]=calibrate_MDR_steps(entity_list, hazard_list, emdat_list, MDR_fun, params_MDR)
% EXAMPLE:
%   [status,output_filename]=calibrate_MDR_steps(entity_list, hazard_list, emdat_list, MDR_fun, params_MDR)
% INPUTS:
%   entity_list: a list of entities (cell array), one per country
%   hazard_list: a list of hazard (cell array), or several per country. If
%      several, a cell array of cell arrays.
%   emdat_list: a list of EM-DAT damages (cell array) containing the fields
%      'year' and 'values'.
%   MDR_fun: a function of x (hazard intensity) and two parameters to be
%      calibrated (e.g., scale and shape).
%   params_MDR: a structure with fields:
%     years_range (optional): range of years to be included (e.g. [1990
%        2010]). If not provided, retrieved from emdat_list{1}.year .
%     use_YDS: =1 to use yearly-varying assets. =0 to use fixed asset values.
%     remove_years_0emdat: determines how to filter years with zero EM-DAT
%        damages (e.g., no damage).
%        0=no filter (all years retained);
%        1=remove years for which the sum of EM-DAT damages over all
%        countries is 0;
%        2=for each country, remove years with zero EM-DAT damages.
%     remove_years_0YDS: list defining whether and how to exclude years
%        with zero simulated damage as follows:
%        do: =1 to exclude years with zero simulated damage. If =0, the
%           other fields are not required.
%        threshold: hazard intensity value above which (>) damage is
%           considered to occur to filter out years (required if do==1).
%        what: determines how to choose if several hazards are provided per
%           country (possible choices: 'any', 'all', 'mean', 'median'). For
%           instance, if 'any', then if any hazard set gives 0 damage the
%           year-country event is removed.
%        min_val: damage value (e.g. in USD) below which damage is not
%           considered to occur (default =0).
%     pars_range: range of parameter values (cell array, each element as
%        [min_val max_val].
%     damFun_xVals: vector of value of hazard intensity to be used when
%        creating the damage function based on MDR_fun (e.g. 0:0.5:10).
%        The damage function is capped above/below the range of the
%        intensity scale given, i.e. using climada_damagefunctions_generate_from_fun
%        with params.IntensityCap=1
%   params_calibration: parameters for the calibration:
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
%           'MMM':  Multi-Model Mean damage estimate vs observated damages.
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
%       exclude_years_0totals: =1 to exclude years where observed
%           damage (total damage for the whole region) is 0. Does not have
%           any effect if params_calibration.type is set to 'RTarea' as
%           those values are not included in the cost function anyway.
%       write_outfile: name of a file where the result from each step
%           should be saved.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%     savefile: file where the output should be saved.
% OUTPUTS:
%   optimal pars: the values of the optimal parameter combination.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180911, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181009, split into
%    sub-functions and progress in implementation
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181127, added output years_i_in
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190116, removing parallel option for patternsearch
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190121, allowing to use patternsearch or regular sampling of parameters, and several initial points for patternsearch (random or uniformly distributed in the parameter space)
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190122, initialize file params_calibration.write_outfile in all cases
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, adding RTarea as a possible type of cost function, and adding parameter params_calibration.underestimation_factor
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, update description of params_MDR.damFun_xVals
%-


% https://ch.mathworks.com/help/gads/examples/constrained-minimization-using-pattern-search.html

% initialization
global climada_global;

%% 0) check input arguments
if ~exist('entity_list','var'),error('Input parameter entity_list is missing');end
if ~exist('hazard_list','var'),error('Input parameter hazard_list is missing');end
if ~exist('emdat_list','var'),error('Input parameter emdat_list is missing');end
if ~exist('MDR_fun','var'),error('Input parameter MDR_fun is missing');end
if ~exist('params','var'),                  params=              struct;end
if ~exist('params_MDR','var'),error('Input parameter params_MDR is missing');end
if ~exist('params_calibration','var'),error('Input parameter params_calibration is missing');end

% check for some parameter fields we need
% params_MDR
if ~isfield(params_MDR,'years_range'),params_MDR.years_range=[min(emdat_list{1}.year) max(emdat_list{1}.year)];end
if ~isfield(params_MDR,'use_YDS'),error('Input parameter params_MDR.use_YDS is missing');end
if ~isfield(params_MDR,'remove_years_0emdat'),error('Input parameter params_MDR.remove_years_0emdat is missing');end
if ~isfield(params_MDR,'remove_years_0YDS'),error('Input parameter params_MDR.remove_years_0YDS is missing');end
if params_MDR.remove_years_0YDS.do
    if ~isfield(params_MDR.remove_years_0YDS,'threshold')
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.threshold required since params_MDR.remove_years_0YDS.do==1 *****')
    end
    if ~isfield(params_MDR.remove_years_0YDS, 'what') && iscell(hazard_list{1})
        error('** ERROR ** Input parameter params_MDR.remove_years_0YDS.what required since params_MDR.remove_years_0YDS.do==1 and multiple hazard per country are provided *****')
    end
    if ~isfield(params_MDR.remove_years_0YDS, 'min_val'),params_MDR.remove_years_0YDS.min_val=0;end
end
if ~isfield(params_MDR, 'pars_range'),error('Input parameter params_MDR.pars_range is missing');end
if ~isfield(params_MDR,'damFun_xVals'), warning('** warning ** params_MDR.damFun_xVals not set, intensity steps of the default damage function will be used *****');end
% params_calibration
if ~isfield(params_calibration,'type'),error('Input parameter params_calibration.type is missing');end
if ~isfield(params_calibration,'MM_how'),error('Input parameter params_calibration.MM_how is missing');end
if ~isfield(params_calibration,'calib_options'),params_calibration.calib_options=struct;end
if ~isfield(params_calibration.calib_options,'method'),params_calibration.calib_options.method='patternsearch';end
if ~isfield(params_calibration.calib_options,'params'),params_calibration.calib_options.params=struct;end
switch params_calibration.calib_options.method
    case 'patternsearch'
        if ~isfield(params_calibration.calib_options.params,'random'),params_calibration.calib_options.params.random=0;end
        if ~isfield(params_calibration.calib_options.params,'nstart'),params_calibration.calib_options.params.nstart=1;end
        if ~isfield(params_calibration.calib_options.params,'InitialMeshSize'),params_calibration.calib_options.params.InitialMeshSize=0.25;end
        if ~isfield(params_calibration.calib_options.params,'step_tolerance'),params_calibration.calib_options.params.step_tolerance=0.001;end
    case 'regular_sampling'
        if ~isfield(params_calibration.calib_options.params,'n_per_dim'),params_calibration.calib_options.params.n_per_dim=5;end
    otherwise
        error('Input field params_calibration.calib_options.method is not valid')
end
if ~isfield(params_calibration,'underestimation_factor'),params_calibration.underestimation_factor=1;end


%% 1) determine useful parameters, formatting
country_list=cellfun(@(x) x.assets.admin0_ISO3,entity_list, 'UniformOutput', 0);
n_countries = length(country_list);

% put hazard_list in double-cell mode in case only 1 hazard was provided.
if ~iscell(hazard_list{1})
    for i = 1:n_countries
        hazard_temp=hazard_list{i};
        hazard_list{i}={};
        hazard_list{i}{1}=hazard_temp;
        clear hazard_temp;
    end
end


%% 2) check time period
same_period = check_time_period(params_MDR.years_range, emdat_list, entity_list, hazard_list, params_MDR.use_YDS);
% 1) if years_range is not the same as all_years_avail, print message
% and abort
if ~same_period
    error('** ERROR ** mismatch in the time period *****')
end
% Now we know that all input have consistent time period:
% params_MDR.years_range and the entities, hazard and emdat data lists.


%% 3) filter years without EM-DAT or simulated damages
% are years without damage in EM-DAT and/or without simulated damage to be removed?
[entity_list, hazard_list, emdat_list, years_i_in] = filter_years(params_MDR.remove_years_0emdat, ...
    params_MDR.remove_years_0YDS, params_MDR.years_range, ...
    entity_list, hazard_list, emdat_list, params_MDR.use_YDS);



%% 4) set boundaries
% retrieve the lower and upper bounds of the parameter
bounds.lb=cellfun(@(x) x(1),params_MDR.pars_range, 'UniformOutput', 1);
bounds.ub=cellfun(@(x) x(2),params_MDR.pars_range, 'UniformOutput', 1);

% normalization of parameter values
% lower bounds are normalized to 1
norm.lb = ones(size(bounds.lb));
% upper bounds are normalized to 2
norm.ub = norm.lb+1;

%% 5) prepare calibration function
% define anonymous function with input factor x (parameters of the damage
% function):
params_MDR2=[];
params_MDR2.use_YDS = params_MDR.use_YDS;
params_MDR2.damFun_xVals = params_MDR.damFun_xVals;
% sets all inputvar for the function except for x, use normalized x.
fun = @(x)calibrate_params_MDR(x,MDR_fun, params_MDR.years_range, ...
    entity_list, hazard_list, emdat_list, norm, bounds,...
    params_MDR2,rmfield(params_calibration,'calib_options'));

%just for testing
% test=calibrate_params_MDR(norm_x0,MDR_fun, params_MDR.years_range, ...
%     entity_list, hazard_list, emdat_list, norm, bounds,...
%     params_MDR2,params_calibration);

%% 6) optimization of parameters

% open file for each steps and overwrite if exists ('w')
if isfield(params_calibration, 'write_outfile')
    fileID=fopen(params_calibration.write_outfile,'w');
    fprintf(fileID,'%s %s %s\n','x(1)','x(2)','cost_value');
    fclose(fileID);
end

n_pars = length(params_MDR.pars_range);
switch params_calibration.calib_options.method
    % ----------------------
    case 'patternsearch'
        
        options = optimoptions('patternsearch','UseParallel',false,...
            'UseCompletePoll', true, 'UseVectorized', false,...
            'MaxFunctionEvaluations',1200,'Display','iter',...
            'Cache','on','InitialMeshSize',params_calibration.calib_options.params.InitialMeshSize,...
            'PollMethod','GPSPositiveBasis2N','StepTolerance',params_calibration.calib_options.params.step_tolerance);
        
        % set starting parameter values: cell array norm_x0
        if params_calibration.calib_options.params.random
            % random start for each starting point
            norm_x0 = cellfun(@(~) random('Uniform',1,2,[1 n_pars]),num2cell(1:params_calibration.calib_options.params.nstart),'UniformOutput',false);
        else
            % equally spaced start points for each dimension (i.e. parameter)
            normx0_vals = linspace(1+1/2/(params_calibration.calib_options.params.nstart),2-1/2/(params_calibration.calib_options.params.nstart),params_calibration.calib_options.params.nstart);
            norm_x0 = cell(repmat(params_calibration.calib_options.params.nstart,[1 length(params_MDR.pars_range)]));
            if params_calibration.calib_options.params.nstart==1
                % only one starting point, take it in the middle (1.5 in normalized space)
                norm_x0{1} = 1.5*ones(size(bounds.lb));
            else
                % several starting points in each dimension (not random)
                for i=1:params_calibration.calib_options.params.nstart
                    for j=1:params_calibration.calib_options.params.nstart
                        norm_x0{i,j}=[normx0_vals(i) normx0_vals(j)];
                    end
                end
            end
        end
        
        % test all starting values, always keep the minimum
        ts=tic;
        fprintf('** lauching pattern search 1 out of %s, norm_x0=(%s,%s) **\n',num2str(numel(norm_x0)),num2str(norm_x0{1}(1)),num2str(norm_x0{1}(2)))
        [x_result,fval] = patternsearch(fun,norm_x0{1},[],[],[],[],norm.lb,norm.ub,[],options);
        for i=2:numel(norm_x0)
            fprintf('** lauching pattern search %s out of %s, norm_x0=(%s,%s) **\n',num2str(i),num2str(numel(norm_x0)),num2str(norm_x0{i}(1)),num2str(norm_x0{i}(2)))
            [x_result_i,fval_i] = patternsearch(fun,norm_x0{i},[],[],[],[],norm.lb,norm.ub,[],options);
            if fval_i < fval
                x_result = x_result_i;
                fval = fval_i;
            end
        end
        fprintf('** all pattern searches completed (%s start values) **\n',num2str(numel(norm_x0)))
        toc(ts)
        
    % ----------------------
    case 'regular_sampling'
        % manually sample many points, regularly spaced, nice for displaying
        % those later on
        % params_calibration.calib_options.n_per_dim values for each parameters
        normx0_vals = linspace(1+1/2/(params_calibration.calib_options.params.n_per_dim),2-1/2/(params_calibration.calib_options.params.n_per_dim),params_calibration.calib_options.params.n_per_dim);
        % indices for looping through
        ndgrid_input = repmat({1:length(normx0_vals)},1,n_pars);
        ngrid_output = cell(1,n_pars);
        [ngrid_output{:}] = ndgrid(ndgrid_input{:});
        for i=1:n_pars,ngrid_output{i}=reshape(ngrid_output{i},1,[]);end
        ind_array = vertcat(ngrid_output{:});
        % output array
        fvals = NaN([1 size(ind_array,2)]);
        % loop through and evaluate cost function
        ts=tic;
        fprintf('looping through %s parameter combinations\n',num2str(size(ind_array,2)))
        for i=1:size(ind_array,2)
            fvals(i) = fun(normx0_vals(ind_array(:,i)));
        end
        fprintf('** all (%s) parameter combinations completed **\n',num2str(size(ind_array,2)))
        toc(ts)
        % pick up the best one
        [fval,i_lowest] = min(fvals);
        x_result = normx0_vals(ind_array(:,i_lowest));
        
    otherwise
        error('Input field params_calibration.calib_options.method is not valid')
end

% convert normalized value of calibrated parameters to their 'real' values
optimal_pars=(x_result-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb;

% computation of damages with the optimal parameter combination
if isfield(params,'savefile')
    save(params.savefile,'optimal_pars','fval','params_MDR','params_calibration','-v7.3');
end

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



function [entity_list, hazard_list, emdat_list, years_in ] = filter_years(remove_years_0emdat, remove_years_0YDS, years_range, entity_list, hazard_list, emdat_list, use_YDS)
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
        error('** ERROR ** unexpected value in params_MDR.remove_years_0emdat *****')
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
