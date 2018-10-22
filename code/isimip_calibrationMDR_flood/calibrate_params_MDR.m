function [ result ] = calibrate_params_MDR(x,MDR_fun,years_range,entity_list, hazard_list, emdat_list,norm_bnds,bounds,params_MDR,params_calibration)
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
%       if arguments norm_bnds and bounds are provided; otherwise assumed
%       to be non-normalized).
%   MDR_fun: function with three arguments: hazard_intensity, and both
%       parameters
%   years_range: range of years included (e.g. [1990 2010])
%   entity_list: list of entities, see function calibrate_MDR_steps
%   hazard_list: list of hazard, see function calibrate_MDR_steps
%   emdat_list: list of observed damages, see function calibrate_MDR_steps
% OPTIONAL INPUT PARAMETERS:
%   norm_bnds: structure containing lb and ub, the lower and upper normalization
%       (never set to 0) bounds of both parameters (norm.lb, norm.ub)
%   bounds: structure containing lb and ub, the lower and upper bounds of
%       parameter values (bounds.ub, bounds.lb) expressed as normalized
%       values. Goes along with norm_bnds: Must be provided if norm_bnds is
%       provided, cannot be provided if norm_bnds is not provided.
%   params_MDR: parameters for the damage function:
%       damFun_xVals: vector of value of hazard intensity to be used when
%           creating the damage function based on MDR_func (e.g. 0:0.5:10).
%           The second last values will be set to the last value in order to
%           ensure a maximum MDR value.
%       use_YDS: =1 to use yearly-varying assets. Default =0.
%   params_calibration: parameters for the calibration:
%       type (string): cost function, one of:
%           'AED':  "Annual Expected Damage": result is the squared
%                   difference of mean year damage 
%           'R2':   (DEFAULT) result is R^2 (= the sum of squared differences of
%                   year damages of emdat and climada for each specific historical year).
%           'R':    result is R (= the sum of the absolute differences of
%                   year damages of emdat and climada for each specific historical year).
%           'RP':   "Return Period": as AED but for different return
%                   periods with weights (not implemented yet) - only makes
%                   sense for long time series
%       MM_how (string): how to deal with Multi-Model hazard sets, one of:
%           'MMM':  Multi-Model Mean damage estimate vs observated damages.
%           'MMMed':Multi-Model Median damage estimate vs observated damages.
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20181015, initial
%   
%-

% initialization
global climada_global
result = 0; %init
damage_at_centroid_temp = climada_global.damage_at_centroid;
climada_global.damage_at_centroid = 0;

% determine peril_ID and DamFunID from entity_list and hazard_list
peril_ID = hazard_list{1}{1}.peril_ID;
DamFunID=unique(entity_list{1}.assets.DamageFunID(~isnan(entity_list{1}.assets.DamageFunID)));

% check input arguments
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
if ~exist(params_calibration),params_calibration=[];end
if ~isfield(params_calibration,'type'),params_calibration.type='R2';end
if ~isfield(params_calibration,'MM_how'),params_calibration.MM_how='MMM';end
if ~ismember(params_calibration.MM_how, {'MMM', 'MMMed'})
    error('Unexpected input value in params_calibration.MM_how')
end
if ~exist(params_MDR),params_MDR=[];end
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
damFun = create_damagefunction_from_fun(MDR_fun, x, params_MDR);
for i=1:length(entity_list)
    entity_list{i}.damagefunctions=damFun;
end

%% (2) calculate YDS

% commented lines below just for reference, can be deleted later on
%     EDS_FL_2005=climada_EDS_calc(entity_isimip,hazard_FL,'',0,2); % the damage calculation
%     damage_2005(:,i) = EDS_FL_2005.damage*EDS_FL_2005.currency_unit;
%     % for damage, look into isimip_YDS_calc?
%     
%     yds_params.silent_mode=2;
%     YDS_FL=isimip_YDS_calc(entity_isimip,hazard_FL,yds_params);
%     damage(:,i) = YDS_FL.damage*EDS_FL_2005.currency_unit;

clear yds_params;
yds_params.silent_mode=2;
all_years = years_range(1):years_range(2);
EDS_list = {};
damages_fullmat = NaN([length(entity_list),length(hazard_list{1}),length(all_years)]);
emdat_fullmat = NaN([length(entity_list),length(all_years)]);
for i=1:length(entity_list)
    EDS_list{i}={};
    for j = 1:length(hazard_list{i})
        if params_MDR.use_YDS
            EDS_list{i}{j} = isimip_YDS_calc(entity_list{i}, hazard_list{i}{j}, yds_params);
        else
            EDS_list{i}{j} = climada_EDS_calc(entity_list{i}, hazard_list{i}{j},[],[],2);
        end
        % get years IDs
        [~,iis] = ismember(EDS_list{i}{j}.yyyy, all_years);
        damages_fullmat(i,j,iis) = EDS_list{i}{j}.damage*entity_list{i}.assets.currency_unit;
        [~,iis2] = ismember(emdat_list{i}.year, all_years);
        emdat_fullmat(i,iis2) = emdat_list{i}.values;
    end
end

% MM_how (string): how to deal with Multi-Model hazard sets, one of:
damages_mat = zeros([length(entity_list),length(hazard_list{1})]);
switch params_calibration.MM_how
    case 'MMM'
        damages_mat = squeeze(mean(damages_fullmat,2));
    case 'MMMed'
        damages_mat = median(damages_fullmat,2);
    otherwise
        error('** ERROR ** unexpected value in params_calibration.MM_how *****');
end
clear  yds_params;

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
diff_yearly=damages_yearly-emdat_yearly;


%% (3) provide difference to em_data
result = cost_function(damages_yearly, emdat_yearly, params_calibration.type);

climada_global.damage_at_centroid = damage_at_centroid_temp;
clear EDS YDS em_data LIA LOCB year_i damage_at_centroid_temp

end

function result = cost_function(damages_yearly, emdat_yearly, type)
diff_yearly=damages_yearly-emdat_yearly;
result=nan;
switch type
    case 'AED2'
%       squared of the difference in yearly mean damage
        result = (mean(emdat_yearly)-mean(damages_yearly))^2;
    case 'R2'
        % mean of the yearly squared difference
        result = mean(diff_yearly.^2);
    case 'R4'
        result = mean(diff_yearly.^4);
    case 'R'
        result = sum(abs(diff_yearly));
    case 'logR'
                [~,YDS] = evalc('climada_EDS2YDS(EDS,hazard)');
        [LIA,LOCB] = ismember(em_data.year,YDS.yyyy);
        % exclude years with 0-damage-entries:
        year_i_vector = find(LIA);
        year_i_vector(em_data.damage(year_i_vector)==0)=[]; 
     %   year_i_vector(YDS.damage(LOCB(year_i_vector))==0)=[];

        result = sum(log(abs(em_data.damage(year_i_vector) - YDS.damage(LOCB(year_i_vector)))));


    case 'RP'
        warning('Type RP not yet implemented, returning squared difference in AED instead');
        em_data_yyyy_allYears = em_data.year(1):em_data.year(end);
        em_data_damage_allYears = zeros(size(em_data_yyyy_allYears));
        for year_i = em_data_yyyy_allYears
            if max(ismember(em_data.year,year_i))
                em_data_damage_allYears(em_data_yyyy_allYears==year_i)= em_data.damage(em_data.year==year_i);
            end
        end
        result = (mean(em_data_damage_allYears)-EDS.ED)^2;
    otherwise
        error('** ERROR ** unexpected calibration type *****')
end
end

function damFun = create_damagefunction_from_fun(MDR_fun, x, params_MDR)
% function create_damagefunction_from_fun
%   Creates a damage function based on function MDR_fun with parameters
%   x(1) and x(2), with intensity values given by params_MDR.damFun_xVals
%   and peril_ID and DamageFunID given by params_MDR.peril_ID and
%   params_MDR.DamFunID, respectively. PAA is set to 1, and the second last
%   value of MDD is set equal to the last value in order to ensure a
%   maximum value that is not exceeded.

%global climada_global
if ~climada_init_vars,return;end % init/import global variables

% correctly position vectors
if size(params_MDR.damFun_xVals,1)>size(params_MDR.damFun_xVals,2),params_MDR.damFun_xVals=params_MDR.damFun_xVals';end

% create a damage function
damFun=[]; % init output
damFun.filename='create_damagefunction_from_fun';
damFun.Intensity=params_MDR.damFun_xVals;
damFun.DamageFunID=params_MDR.DamFunID*ones(size(damFun.Intensity));
damFun.peril_ID=cellstr(repmat(params_MDR.peril_ID,length(damFun.Intensity),1))';
damFun.Intensity_unit=repmat({'m'}',size(damFun.Intensity));
damFun.peril_ID=cellstr(repmat(params_MDR.peril_ID,length(damFun.Intensity),1))';
damFun.MDD = MDR_func(damFun.Intensity, x(1), x(2));
damFun.MDD(length(damFun.MDD)-1) = damFun.MDD(end);
damFun.PAA = ones(size(damFun.Intensity));
name = 'Calibration damage function';
damFun.name = repmat({name},size(damFun.Intensity));
damFun.datenum = zeros(size(damFun.Intensity))+ now;

end
