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
%   entity=climada_entity_load('BRBCUBDOMPRI'); % if repeated
%   isimip_tc_calibrate(entity);
%
%   entity=climada_entity_load('GLB_isimip_entity'); % full globe
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
%    regions_file: the .csv file with the region definition (groups of
%       countries)
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
if ~isfield(params,'regions_file'),params.regions_file='';end
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
regions_file=[climada_global.data_dir filesep 'isimip' filesep 'TC-damage-function-regions.csv'];

%
% populate default parameters in params
if isempty(params.damage_data_file),params.damage_data_file=damage_data_file;end
if isempty(params.price_deflator_file),params.price_deflator_file=price_deflator_file;end
if isempty(params.regions_file),params.regions_file=regions_file;end
if isempty(params.check_plot),params.check_plot=1;end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.damage_data_file,'file'),params.damage_data_file=[isimip_data_dir filesep params.damage_data_file];end
if ~exist(params.price_deflator_file,'file'),params.price_deflator_file=[isimip_data_dir filesep params.price_deflator_file];end
if ~exist(params.regions_file,'file'),params.regions_file=[isimip_data_dir filesep params.regions_file];end

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
        
        if ~isempty(params.regions_file) % assign countries to regions
            if exist(params.regions_file,'file')
                regions_table=climada_csvread(params.regions_file);
                
                % regions_table fields renamed
                if isfield(regions_table,'ISO'),regions_table.ISO3=regions_table.ISO;regions_table=rmfield(regions_table,'ISO');end
                if isfield(regions_table,'ID'),regions_table.NatId=regions_table.ID;regions_table=rmfield(regions_table,'ID');end
                if isfield(regions_table,'Reg_ID'),regions_table.RegId=regions_table.Reg_ID;regions_table=rmfield(regions_table,'Reg_ID');end
                if isfield(regions_table,'Reg_name'),regions_table.RegName=regions_table.Reg_name;regions_table=rmfield(regions_table,'Reg_name');end
                
                % % in case we need region for each asset, use the following
                % if isfield(entity.assets,'NatId')
                %     entity.assets.RegId=entity.assets.NatId*0; % init
                %     unique_assets_NatId=unique(entity.assets.NatId);
                %     for i=1:length(unique_assets_NatId)
                %         NatId_assets_pos=find(entity.assets.NatId==unique_assets_NatId(i));
                %         NatId_region_pos=find(regions_table.NatId==unique_assets_NatId(i));
                %         if ~isempty(NatId_assets_pos) && length(NatId_region_pos)==1
                %             entity.assets.RegId(NatId_assets_pos)=regions_table.RegId(NatId_region_pos);
                %         else
                %             fprintf('WARNING: Region matching error, proceed with caution\n');
                %         end
                %     end % i
                % end
                
            end
        end
        
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
        
        if isfield(entity,'hazard')
            %hazard=entity.hazard;entity=rmfield(entity,'hazard');
            %[~,EDS]=isimip_YDS_calc(entity,hazard);
            [~,EDS]=isimip_YDS_calc(entity,entity.hazard);
        elseif ~isempty(params.hazard_file)
            hazard=climada_hazard_load(params.hazard_file);
            if isempty(hazard),fprintf('ERROR: hazard %s not found\n',params.hazard_file);return,end
            [~,EDS]=isimip_YDS_calc(entity,hazard);
        end
        
        %EDS=climada_EDS_calc(entity,entity.hazard); % FAST for TEST:
        % avoids using isimip_YDS_calc, hence comment section above. Does
        % not support regions, though
        
        n_EDS=length(EDS);
        
        if isfield(EDS(1),'NatId') % EDS has a NatId, can be mapped to region
            RegId_used=zeros(1,n_EDS);
            NatId_used=zeros(1,n_EDS);
            for EDS_i=1:n_EDS
                pos=find(regions_table.NatId==EDS(EDS_i).NatId);
                if length(pos)==1
                    EDS(EDS_i).RegId   = regions_table.RegId(pos);
                    EDS(EDS_i).RegName = regions_table.RegName{pos};
                    RegId_used(EDS_i)  = EDS(EDS_i).RegId;
                    NatId_used(EDS_i)  = EDS(EDS_i).NatId;
                    ISO3_used{EDS_i}   = EDS(EDS_i).comment; % sorry comment for historic reasons
                    RegName_used{EDS_i}= EDS(EDS_i).RegName;
                else
                    fprintf('WARNING: EDS(%i) for %s not mapped to a region\n',EDS_i,EDS(EDS_i).comment);
                end
            end % EDS_i
            [unique_RegId_used,i]=unique(RegId_used);
            unique_RegName_used=RegName_used(i);
            n_regions=length(unique_RegId_used);
        else
            n_regions=0; % just global view
        end
        
        if n_regions>0
            fprintf('processing %i regions: ',n_regions)
            for region_i=1:n_regions;fprintf('%s ',unique_RegName_used{region_i}),end;fprintf('\n')
        end
        
        for region_i=0:n_regions
            
            if region_i==0 % 'global', sum up over all EDSs
                damage_sim=EDS(1).damage;
                region_ISO3_str='all';
                for EDS_i=2:length(EDS)
                    damage_sim=damage_sim+EDS(EDS_i).damage;
                end % EDS_i
            else
                % sum up over countries within region_i
                damage_sim=EDS(1).damage*0; % init empty
                region_ISO3_str=''; % init
                %fprintf('region %s (Id %i): ',unique_RegName_used{region_i},unique_RegId_used(region_i));
                for EDS_i=1:n_EDS
                    if EDS(EDS_i).RegId==unique_RegId_used(region_i)
                        damage_sim=damage_sim+EDS(EDS_i).damage;
                        region_ISO3_str=[region_ISO3_str EDS(EDS_i).comment ' '];
                    end
                end % EDS_i
                region_ISO3_str = [unique_RegName_used{region_i} ': ' deblank(region_ISO3_str)];
                fprintf('%s\n',region_ISO3_str);
            end
            
            %         % global per event view: match calculated damages with reported ones
            %         if params.check_plot % plot damage data
            %             figure('Name',sprintf('%s: event reported versus simulated damages',region_ISO3_str),'Color',[1 1 1]);
            %             yyaxis left
            %             plot(damage_data.year(damage_data.hazard_matched),log10(damage_data.damage_tot(damage_data.hazard_matched)),'.g');hold on
            %             xlabel('year');ylabel('log10(reported damage) [USD 2005]');
            %             yyaxis right
            %             plot(damage_data.year(damage_data.hazard_matched),...
            %                 log10(damage_sim(damage_data.hazard_index(damage_data.hazard_matched))),'.r');hold on
            %             ylabel('log10(simulated damage) [USD]');
            %             legend({'reported','simulated'});
            %             title(region_ISO3_str);
            %         end % params.check_plot
            
            % per year view: match calculated damages with reported ones
            if params.check_plot % plot damage data
                
                % sum data up over years
                year_mtc=damage_data.year(damage_data.hazard_matched);
                damage_tot_mtc=damage_data.damage_tot(damage_data.hazard_matched);
                damage_sim_mtc=damage_sim(damage_data.hazard_index(damage_data.hazard_matched));
                year_mtc_uni=unique(year_mtc);
                damage_tot_mtc_uni = year_mtc_uni*0;
                damage_sim_mtc_uni = year_mtc_uni*0;
                for year_i=1:length(year_mtc_uni)
                    pos=find(year_mtc==year_mtc_uni(year_i));
                    if ~isempty(pos)
                        damage_tot_mtc_uni(year_i) = sum(damage_tot_mtc(pos));
                        damage_sim_mtc_uni(year_i) = sum(damage_sim_mtc(pos));
                    end
                end % year_i
                
                %             figure('Name',sprintf('%s: year reported versus simulated damages',region_ISO3_str),'Color',[1 1 1]);
                %             yyaxis left
                %             plot(year_mtc_uni,log10(damage_tot_mtc_uni),'.g');hold on
                %             xlabel('year');ylabel('log10(reported damage) [USD 2005]');
                %             yyaxis right
                %             plot(year_mtc_uni,log10(damage_sim_mtc_uni),'.r');hold on
                %             ylabel('log10(simulated damage) [USD]');
                %             legend({'reported','simulated'});
                %             title(region_ISO3_str);
                
                % and the scatter plot
                figure('Name',sprintf('%s: year reported versus simulated damages',region_ISO3_str),'Color',[1 1 1]);
                plot(log10(damage_tot_mtc_uni),log10(damage_sim_mtc_uni),'xr');hold on
                plot(log10(damage_tot_mtc_uni),log10(damage_tot_mtc_uni),'-g');
                log10_max_damage=max(max(log10(damage_tot_mtc_uni)),max(log10(damage_sim_mtc_uni)));
                plot([0 log10_max_damage],[0 log10_max_damage],':g');
                xlabel('log10(reported damage) [USD 2005]');ylabel('log10(simulated damage) [USD]')
                title(sprintf('%s years %i..%i',region_ISO3_str,min(year_mtc_uni),max(year_mtc_uni)));
            end % params.check_plot
            
        end % region_i
        
    else
        fprintf('ERROR: %s not found\n',params.damage_data_file)
    end % exist(params.damage_data_file,'file')
end % ~isempty(params.damage_data_file)

%end % isimip_tc_calibrate