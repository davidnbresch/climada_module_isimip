%function res=isimip_tc_calibrate(entity,params)
% climada isimip tropcial cyclone calibrate
% MODULE:
%   isimip
% NAME:
%   isimip_tc_calibrate
% PURPOSE:
%   calibrate the global model for tropical cyclones based in ibtracs data
%
%   previous call: isimip_gdp_entity
%   next call: many...
% CALLING SEQUENCE:
%   res=isimip_tc_calibrate(params)
% EXAMPLE:
%   params.hazard_file='GLB_0360as_TC_hist';
%   entity=isimip_gdp_entity({'PRI','DOM','BRB','CUB'},params);
%   entity=climada_entity_load('BRBCUBDOMPRI') % if repeated
%   isimip_tc_calibrate(entity);
% INPUTS:
%   entity: an isimip entity, output from isimip_gdp_entity
%       in essence a climada entity on isimip NatId grid with additional
%       fields for special isimip use and hazard intensity joined.
%       if ='params', just return default parameters, in res, i.e. the
%           first output, already.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    damage_data_file: the .csv file with damage data to be joined to.
%       Default file is matching_Natcat-damages_ibtracs_1980-2014.csv in
%       the isimip data folder.
%    price_deflator_file: the .csv file with the price deflator.
%       We inflate damages with a factor 1/deflator
%       Default file is GDP_deflator_converted_base2005_1969-2016_source_BEA.csv
%       in the isimip data folder.
%    hazard_file: the filename of a climada hazard set (if not part of
%       entity, as might be the case if entity from isimip_gdp_entity)
%       Does not not need to contain path or .mat extension, since loaded
%       using climada_hazard_load
%    check_plot; =1; plot checkp plots, =0 not (default)
% OUTPUTS:
%   res: the output, empty if not successful
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160218, initial
%-

res=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables


% poor man's version to check arguments
if ~exist('entity','var'),entity=[];end % pass all parameters as structure
if ~exist('params','var'),params=struct;end % pass all parameters as structure

% check for some parameter fields we need
if ~isfield(params,'damage_data_file'),params.damage_data_file='';end
if ~isfield(params,'price_deflator_file'),params.price_deflator_file='';end
if ~isfield(params,'hazard_file'),params.hazard_file='';end
if ~isfield(params,'check_plot'),params.check_plot=[];end

% PARAMETERS
%
verbose=1; % default=1, to suppress output to stdout later
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end
% define default filenames
damage_data_file=[climada_global.data_dir filesep 'isimip' filesep 'matching_Natcat-damages_ibtracs_1980-2014.csv'];
price_deflator_file=[climada_global.data_dir filesep 'isimip' filesep 'GDP_deflator_converted_base2005_1969-2016_source_BEA.csv'];
%
% populate default parameters in params
if isempty(params.damage_data_file),params.damage_data_file=damage_data_file;end
if isempty(params.price_deflator_file),params.price_deflator_file=price_deflator_file;end
if isempty(params.check_plot),params.check_plot=1;end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.damage_data_file,'file'),params.damage_data_file=[isimip_data_dir filesep params.damage_data_file];end
if ~exist(params.price_deflator_file,'file'),params.price_deflator_file=[isimip_data_dir filesep params.price_deflator_file];end

if strcmpi(entity,'params'),res=params;return;end % special case, return the full params structure

if ~isempty(params.damage_data_file)
    if exist(params.damage_data_file,'file')
        damage_data=climada_csvread(params.damage_data_file);
        
        % damage_data:
        % MR_ID: {1x2220 cell}
        % Event: {1x2220 cell}
        % hist_insured_losses: [1x2220 double]
        % hist_tot_losses: [1x2220 double]
        % Long: [1x2220 double]
        % Lat: [1x2220 double]
        % affected: {1x2220 cell}
        % Country: {1x2220 cell}
        % region: {1x2220 cell}
        % CountryID: [1x2220 double]
        % Date: {1x2220 cell}
        % ibtracs_ID: {1x2220 cell}
        % ISO: {1x2220 cell}
        
        % a unique ID, a number for faster comparison, hence convert 1950166N14262 to 1950166.14262
        n_damage_events=length(damage_data.ibtracs_ID);
        damage_data.ID_no=zeros(1,n_damage_events);
        damage_data.year=zeros(1,n_damage_events);
        for event_i=1:n_damage_events
            if length(damage_data.ibtracs_ID{event_i})>9
                damage_data.ID_no(event_i)=str2double(damage_data.ibtracs_ID{event_i}(1:7))+str2double(damage_data.ibtracs_ID{event_i}(9:end))/100000;
            end % length(damage_data.ibtracs_ID{event_i})>9
            year2digit=str2double(damage_data.Date{event_i}(7:end)); % really bad, only 2-digit year stored
            if year2digit<30
                damage_data.year(event_i)=2000+year2digit;
            else
                damage_data.year(event_i)=1900+year2digit;
            end
        end % event_i
        
        inflator_data.year=1859:2100; % init dummy
        inflator_data.inflator_factor=inflator_data.year*0+1; % init dummy
        if ~isempty(params.price_deflator_file)
            if exist(params.price_deflator_file,'file')
                inflator_data=climada_csvread(params.price_deflator_file);
                inflator_data.inflator_factor=1./(inflator_data.GDP_deflator_base2005/100);
            else
                fprintf('WARNING: %s not found, no price inflation\n',params.price_deflator_file)
            end % exist(params.damage_data_file,'file')
        end % ~isempty(params.damage_data_file)
        
        % damage_data:
        % MR_ID: {1x2220 cell}
        % Event: {1x2220 cell}
        % hist_insured_losses: [1x2220 double]
        % hist_tot_losses: [1x2220 double]
        % Long: [1x2220 double]
        % Lat: [1x2220 double]
        % affected: {1x2220 cell}
        % Country: {1x2220 cell}
        % region: {1x2220 cell}
        % CountryID: [1x2220 double]
        % Date: {1x2220 cell}
        % ibtracs_ID: {1x2220 cell}
        % ISO: {1x2220 cell}
        
        unique_years=unique(damage_data.year); % unique list of years
        for year_i=1:length(unique_years)
            inflator_pos=find(inflator_data.year==unique_years(year_i));
            damage_pos=find(damage_data.year==unique_years(year_i));
            
            damage_data.inflator_factor(damage_pos)=inflator_data.inflator_factor(inflator_pos); % to keep track
            
            damage_data.damage_ins(damage_pos)=damage_data.hist_insured_losses(damage_pos).*damage_data.inflator_factor(damage_pos);
            damage_data.damage_tot(damage_pos)=damage_data.hist_tot_losses(damage_pos).*damage_data.inflator_factor(damage_pos);
        end % year_i
        
        if params.check_plot % plot damage data
            figure('Name','reported total damages','Color',[1 1 1]);
            yyaxis left
            plot(damage_data.year,log10(damage_data.hist_tot_losses),'.b');hold on
            plot(damage_data.year,log10(damage_data.damage_tot)     ,'.g');hold on
            xlabel('year');ylabel('log10(damage) [USD 2005]');
            yyaxis right
            plot(damage_data.year,damage_data.inflator_factor,'-r');hold on
            ylabel('inflation factor');
            legend({'historic','inflaed'});
            title('reported total damages');
        end % params.check_plot
        
        % now, we have
        % damage_data.damage_tot(record_i): total damage for record_i
        % damage_data.year(record_i): year
        % damage_data.ID_no(record_i): ibtracs ID for damage record
        % entity.hazard.ID_no(event_i): ibtracs ID for event_i
        
        % match ID_no
        matched=0;
        not_matched=0;
        n_records=length(damage_data.ID_no);
        damage_data.hazard_index  =damage_data.ID_no*0; % init
        damage_data.hazard_matched=damage_data.hazard_index; % init
        for record_i=1:n_records
            pos=find(entity.hazard.ID_no==damage_data.ID_no(record_i));
            if length(pos)==1
                damage_data.hazard_index(record_i)=pos;
                damage_data.hazard_matched(record_i)=1;
                matched=matched+1;
            else
                not_matched=not_matched+1;
            end
        end % record_i
        damage_data.hazard_matched=logical(damage_data.hazard_matched);
        fprintf('%i of %i (%i%%) records matched with hazard events\n',matched,n_records,ceil(matched/n_records*100));
        
        %         if isfield(entity,'hazard')
        %             [YDS,EDS,stats]=isimip_YDS_calc(entity,entity.hazard);
        %         elseif ~isempty(params.hazard_file)
        %             hazard=climada_hazard_load(params.hazard_file);
        %             if isempty(hazard),fprintf('ERROR: hazard %s not found\n',params.hazard_file);return,end
        %             [YDS,EDS,stats]=isimip_YDS_calc(entity,hazard);
        %         end
        
        EDS=climada_EDS_calc(entity,entity.hazard);

        % sum up over countries
        damage_sim=EDS(1).damage;
        for EDS_i=2:length(EDS)
            damage_sim=damage_sim+EDS(EDS_i).damage;
        end % EDS_i
        
        % per event view: match calculated damages with reported ones
        if params.check_plot % plot damage data
            figure('Name','event reported versus simulated damages','Color',[1 1 1]);
            yyaxis left
            plot(damage_data.year(damage_data.hazard_matched),log10(damage_data.damage_tot(damage_data.hazard_matched)),'.g');hold on
            xlabel('year');ylabel('log10(reported damage) [USD 2005]');
            yyaxis right
            plot(damage_data.year(damage_data.hazard_matched),...
                log10(damage_sim(damage_data.hazard_index(damage_data.hazard_matched))),'.r');hold on
            ylabel('log10(simulated damage) [USD]');
            legend({'reported','simulated'});
            title('total damages');
        end % params.check_plot
        
        % per year view: match calculated damages with reported ones
        if params.check_plot % plot damage data
           
            % sum data up over years
            year_mtc=damage_data.year(damage_data.hazard_matched);
            damage_tot_mtc=damage_data.damage_tot(damage_data.hazard_matched);
            damage_sim_mtc=damage_sim(damage_data.hazard_index(damage_data.hazard_matched));
            [year_mtc_uni,uni_pos]=unique(year_mtc);
            damage_tot_mtc_uni = year_mtc_uni*0;
            damage_sim_mtc_uni = year_mtc_uni*0;
            for year_i=1:length(year_mtc_uni)
                pos=find(year_mtc==year_mtc_uni(year_i));
                if ~isempty(pos)
                    damage_tot_mtc_uni(year_i) = sum(damage_tot_mtc(pos));
                    damage_sim_mtc_uni(year_i) = sum(damage_sim_mtc(pos));
                end
            end % year_i
            
            figure('Name','year reported versus simulated damages','Color',[1 1 1]);
            yyaxis left
            plot(year_mtc_uni,log10(damage_tot_mtc_uni),'.g');hold on
            xlabel('year');ylabel('log10(reported damage) [USD 2005]');
            yyaxis right
            plot(year_mtc_uni,log10(damage_sim_mtc_uni),'.r');hold on
            ylabel('log10(simulated damage) [USD]');
            legend({'reported','simulated'});
            title('total damages');
            
            % and the scatter plot
            figure('Name','year reported versus simulated damages','Color',[1 1 1]);
            plot(log10(damage_tot_mtc_uni),log10(damage_sim_mtc_uni),'xr');hold on
            plot(log10(damage_tot_mtc_uni),log10(damage_tot_mtc_uni),'-g');
            log10_max_damage=max(max(log10(damage_tot_mtc_uni)),max(log10(damage_sim_mtc_uni)));
            plot([0 log10_max_damage],[0 log10_max_damage],':g');
            xlabel('log10(reported damage) [USD 2005]');ylabel('log10(simulated damage) [USD]')
            title(sprintf('total damages years %i..%i',min(year_mtc_uni),max(year_mtc_uni)));
        end % params.check_plot
        
    else
        fprintf('ERROR: %s not found\n',params.damage_data_file)
    end % exist(params.damage_data_file,'file')
end % ~isempty(params.damage_data_file)

%end % isimip_tc_calibrate