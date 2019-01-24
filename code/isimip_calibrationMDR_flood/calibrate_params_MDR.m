function [ result ] = calibrate_params_MDR(x,MDR_fun,years_range,...
    entity_list, hazard_list, emdat_list,norm,bounds,...
    params_MDR,params_calibration)
% climada isimip flood
% MODULE:
%   isimip
% NAME:
%   calibrate_params_MDR
% PURPOSE:
%   Use within the iteration in the calibration process. This function
%      (1) changes the damagefunction based on x,
%      (2) calculate the YDS based on entity_list, hazard_list, and the new
%          damage function,
%      (3) provides the cost function results (e.g., Rsquared difference)
%          between the (regional total) YDS and emdat_list either based on
%          yearly values or based on specific return times.
% CALLING SEQUENCE:
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% EXAMPLE:
%   RegionID='NAM';
%   years_range=[1990 2000];
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   [result]=calibrate_params_MDR([1.5 1.5],...[to be completed])
% INPUTS:
%   x: double/vector, the two parameters to be calibrated (normalized value
%       if arguments norm and bounds are provided; otherwise assumed
%       to be non-normalized).
%   MDR_fun: function with three arguments: hazard_intensity, and both
%       parameters
%   years_range: range of years included (e.g. [1990 2010])
%   entity_list: list of entities, see function calibrate_MDR_steps
%   hazard_list: list of hazard, see function calibrate_MDR_steps
%   emdat_list: list of observed damages, see function calibrate_MDR_steps
% OPTIONAL INPUT PARAMETERS:
%   norm: structure containing lb and ub, the lower and upper normalization
%       (never set to 0) bounds of both parameters (norm.lb, norm.ub).
%       Usually, norm.lb=[1 1] and norm.ub=[2 2].
%   bounds: structure containing lb and ub, the lower and upper bounds of
%       parameter values (bounds.ub, bounds.lb) expressed as non-normalized
%       values. Goes along with norm: Must be provided if norm is
%       provided, cannot be provided if norm is not provided. These
%       two input parameters are used to compute the raw x value if x is
%       given as a normalized value.
%   params_MDR: parameters for the damage function:
%       damFun_xVals: vector of value of hazard intensity to be used when
%           creating the damage function based on MDR_fun (e.g. 0:0.5:10).
%           The last value will be repeated a second time at
%           max(damFun_xVals)+5 in order to ensure a maximum MDD value.
%       use_YDS: =1 to use yearly-varying assets. Default =0.
%       peril_ID: if not given, is retrieved from hazard_list{1}{1}.peril_ID
%       DamFunID: if not given, is retrieved from
%           entity_list{1}.assets.DamageFunID (must be unique in the entity)
%   params_calibration: parameters for the calibration:
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
%       write_outfile: name of a file where the result from each step
%           should be saved.
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181015, initial
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181022, split into
%    sub-functions, test works
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181106, fast damage computation
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181127, use climada_EDS_calc_fast and climada_damagefunctions_generate_from_fun
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181210, use log10 instead of log
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, adding RTarea as a possible type of cost function, and adding parameter params_calibration.underestimation_factor
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190124, setting cases with 0 to 1 before taking the log (for type dlog2, dabslog)
%-

% initialization
global climada_global
result = 0; %init
damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 0;

% determine peril_ID and DamFunID from entity_list and hazard_list
peril_ID = hazard_list{1}{1}.peril_ID;
DamFunID=unique(entity_list{1}.assets.DamageFunID(~isnan(entity_list{1}.assets.DamageFunID)));

%% 0) check input arguments
if ~exist('x','var'),error('Input parameter x is missing');end %
if ~exist('MDR_fun','var'),error('Input parameter MDR_fun is missing');end %
if ~exist('entity_list','var'),error('Input parameter entity_list is missing');end %
if ~exist('hazard_list','var'),error('Input parameter hazard_list is missing');end %
if ~exist('emdat_list','var'),error('Input parameter emdat_list is missing');end %
if xor(exist('bounds','var'), exist('norm','var')) % check if variables for normalizatiuon are given as input, if no, assume non-normalized input x
    error('Missing input: Function requires either both optional input variables norm and bounds or neither.')
elseif (exist('norm','var') && exist('bounds','var')) % compute x from normalized x using norm and bound:
    norm.x = x; % if x was given normalized, do inverse normalization to get actual value of x
    x=(norm.x-norm.lb).*(bounds.ub-bounds.lb)./(norm.ub-norm.lb)+bounds.lb;
end
if ~exist('params_calibration','var'),params_calibration=[];end
if ~exist('params_MDR','var'),params_MDR=[];end

% check for some parameter fields we need
if ~isfield(params_calibration,'type'),params_calibration.type='R2';end
if ~isfield(params_calibration,'MM_how'),params_calibration.MM_how='MMM';end
if ~ismember(params_calibration.MM_how, {'MMM', 'MMMed'})
    error('Unexpected input value in params_calibration.MM_how')
end
if ~isfield(params_calibration,'underestimation_factor'),params_calibration.underestimation_factor=1;end
if ~isfield(params_MDR,'damFun_xVals')
    % give intensity scale the one of the existing damage function in the entity
    choose_DF = find(strcmp(entity_list{1}.damagefunctions.peril_ID,peril_ID) .* (entity_list{1}.damagefunctions.DamageFunID==DamFunID(1)));
    params_MDR.damFun_xVals=entity_list{1}.damagefunctions.Intensity(choose_DF);
    clear choose_DF;
end
if ~isfield(params_MDR,'peril_ID'),params_MDR.peril_ID=peril_ID;end
if ~isfield(params_MDR,'DamFunID'),params_MDR.DamFunID=DamFunID;end
if ~strcmp(peril_ID, params_MDR.peril_ID),error('** ERROR ** peril_ID does not match *****');end
if DamFunID~=params_MDR.DamFunID,error('** ERROR ** DamFunID does not match *****');end

%% (1) change damage function of entity
damFun = climada_damagefunctions_generate_from_fun(params_MDR.damFun_xVals, MDR_fun, x, params_MDR);

%% (2) calculate damages (YDS or EDS)
all_years = years_range(1):years_range(2);
damages_fullmat = NaN([length(entity_list),length(hazard_list{1}),length(all_years)]);
emdat_fullmat = NaN([length(entity_list),length(all_years)]);
for i=1:length(entity_list)
    for j = 1:length(hazard_list{i})
        temp = climada_EDS_calc_fast(entity_list{i}, hazard_list{i}{j}, damFun, params_MDR.use_YDS,1,'',0,1);
        % get years IDs
        [~,iis] = ismember(hazard_list{i}{j}.yyyy', all_years);
        damages_fullmat(i,j,iis) = temp.damage*entity_list{i}.assets.currency_unit;
        clear temp;
    end
    [~,iis2] = ismember(emdat_list{i}.year, all_years);
    emdat_fullmat(i,iis2) = emdat_list{i}.values;
end

% MM_how (string): how to deal with Multi-Model hazard sets, one of:
%damages_mat = zeros([length(entity_list),length(hazard_list{1})]);
switch params_calibration.MM_how
    case 'MMM'
        damages_mat = squeeze(mean(damages_fullmat,2));
    case 'MMMed'
        damages_mat = squeeze(median(damages_fullmat,2));
    otherwise
        error('** ERROR ** unexpected value in params_calibration.MM_how *****');
end
clear  yds_params damages_fullmat;

% sum damages per year over the whole region. keep only damages for
% selected years (i.e., ignore NaNs)
emdat_yearly = mean(emdat_fullmat, 1, 'omitnan');
damages_yearly = mean(damages_mat, 1, 'omitnan');

% check that NANs are the same; remove those
not_nans = ~isnan(emdat_yearly);
check_nans = (~not_nans) == isnan(damages_yearly);
if ~all(check_nans)
    error('** ERROR ** NANs do not match (EM-DAT vs computed) *****');
end
emdat_yearly=emdat_yearly(not_nans);
damages_yearly=damages_yearly(not_nans);


%% (3) provide difference to em_data
result = cost_function(damages_yearly, emdat_yearly, params_calibration.type, params_calibration.underestimation_factor);

% save result to file
if isfield(params_calibration, 'write_outfile')
    fileID=fopen(params_calibration.write_outfile,'a');
    fprintf(fileID,'%g %g %g\n',x(1),x(2),result);
    fclose(fileID);
end

climada_global.damage_at_centroid = damage_at_centroid_temp;
clear EDS YDS em_data LIA LOCB year_i damage_at_centroid_temp

end

function result = cost_function(damages_yearly, emdat_yearly, type, underestimation_factor)
diff_yearly=damages_yearly-emdat_yearly;
diff_yearly(diff_yearly<0) = diff_yearly(diff_yearly<0)*underestimation_factor;
result=nan;
switch type
    case 'AED2'
        % squared of the difference in yearly mean damage
        result = mean(damages_yearly)-mean(emdat_yearly);
        if result < 0,result=result*underestimation_factor;end
        result = result.^2;
    case 'R2'
        % mean of the yearly squared difference
        result = mean(diff_yearly.^2);
    case 'R4'
        % mean of the yearly ^4 difference
        result = mean(diff_yearly.^4);
    case 'R'
        % mean of the yearly absolute difference
        result = mean(abs(diff_yearly));
    case 'dlog2'
        % mean of the yearly squared difference of the log
        damages_yearly(damages_yearly==0)=1;
        emdat_yearly(emdat_yearly==0)=1;
        result = log10(damages_yearly)-log10(emdat_yearly);
        result(result<0) = result(result<0)*underestimation_factor;
        result = mean(result.^2);
    case 'dabslog'
        % mean of the yearly absolute difference of the log
        % need to set 0 to 1
        damages_yearly(damages_yearly==0)=1;
        emdat_yearly(emdat_yearly==0)=1;
        result = log10(damages_yearly)-log10(emdat_yearly);
        result(result<0) = result(result<0)*underestimation_factor;
        result = mean(abs(result));
    case 'RTarea'
        % indices to take out to remove 0s (rt_n_out lowest values)
        rt_n_out = max([sum(damages_yearly==0) sum(emdat_yearly==0)]);
        i_in = (rt_n_out+1):length(damages_yearly);
        % sort values to easily exclude some
        damages_yearly_sorted=sort(damages_yearly);
        emdat_yearly_sorted=sort(emdat_yearly);
        % get return time and difference in log10 of return values
        rts_in = (length(damages_yearly)+1)./(length(i_in):-1:1);
        rtvals_diff = log10(damages_yearly_sorted(i_in))-log10(emdat_yearly_sorted(i_in));
        % multiply by underestimation factor where needed
        rtvals_diff(rtvals_diff<0) = rtvals_diff(rtvals_diff<0)*underestimation_factor;
        result = trapz(log10(rts_in),abs(rtvals_diff));
    case 'RP'
        % fit GEV
%         pd = fitdist(emdat_yearly', 'GeneralizedExtremeValue');
%         'haha'
%         'hoho'
        error('Type RP not yet implemented, returning squared difference in AED instead');
        damages_yearly(damages_yearly==0)=1;
        emdat_yearly(emdat_yearly==0)=1;
        rts = (length(damages_yearly)+1)./(length(damages_yearly):-1:1);
        plot(log10(rts), sort(damages_yearly));hold on;
        plot(log10(rts), sort(emdat_yearly),'r');
        figure;
        plot(log10(rts), sort(log10(damages_yearly)));hold on;
        plot(log10(rts), sort(log10(emdat_yearly)),'r');
        % area between the curves
        plot(log10(rts),abs(sort(log10(damages_yearly))-sort(log10(emdat_yearly))))
        result = trapz(log10(rts),abs(sort(log10(damages_yearly))-sort(log10(emdat_yearly))));
    otherwise
        error('** ERROR ** unexpected calibration type *****')
end
end
