%function res=isimip_tc_calibrate(entity,params)
% climada isimip tropcial cyclone calibrate
% MODULE:
%   isimip
% NAME:
%   isimip_tc_calibrate
% PURPOSE:
%   calibrate the global model for tropical cyclones based on ibtracs data
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
%   params.hazard_file='GLB_0360as_TC_hist';
%   entity=isimip_gdp_entity('ALL_IN_ONE',params);
%   entity=climada_entity_load('GLB_0360as_entity'); % full globe
%   isimip_tc_calibrate(entity);
% INPUTS:
%   entity: an isimip entity, output from isimip_gdp_entity
%       in essence a climada entity on isimip NatID grid with additional
%       fields for special isimip use and hazard intensity joined.
%       if ='params', just return default parameters, in res, i.e. the
%           first output, already.
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields (see also ISO3='params' above):
%    damage_data_file: the .csv file with damage data to be joined to.
%       Default file is matching_Natcat-damages_ibtracs_1980-2014.csv in
%       the isimip data folder. See there for required columns
%    price_deflator_file: the .csv file with the price deflator. Needs to
%       contain cloumns labeled year and GDP_deflator_base (or GDP_deflator_base2005)
%       We inflate damages with a factor 1/deflator
%       Default file is GDP_deflator_converted_base2005_1969-2016_source_BEA.csv
%       in the isimip data folder.
%    inflation_reference_year: the reference year for de/inflation, default=2005
%    regions_file: the .csv file with the region definition (groups of
%       countries)
%    hazard_file: the filename of a climada hazard set (if not part of
%       entity, as might be the case if entity from isimip_gdp_entity)
%       Does not not need to contain path or .mat extension, since loaded
%       using climada_hazard_load
%    check_plot; =1; plot checkp plots, =0 not (default)
%       if =-1, only plot inflated damages
%       if =-2, only plot year view of damages (currently disabled)
%       if =-3, only plot scatter plot reported vs simulated
% OUTPUTS:
%   res: the output, empty if not successful
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20160218, initial
% David N. Bresch, david.bresch@gmail.com, 20160224, all done, NatID
% David N. Bresch, david.bresch@gmail.com, 20160304, store_to_entity, DFC added
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
if ~isfield(params,'inflation_reference_year'),params.inflation_reference_year=[];end
if ~isfield(params,'store_to_entity'),params.store_to_entity=[];end

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
if isempty(params.inflation_reference_year),params.inflation_reference_year=2005;end
if isempty(params.store_to_entity),params.store_to_entity=0;end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(params.damage_data_file,'file'),params.damage_data_file=[isimip_data_dir filesep params.damage_data_file];end
if ~exist(params.price_deflator_file,'file'),params.price_deflator_file=[isimip_data_dir filesep params.price_deflator_file];end
if ~exist(params.regions_file,'file'),params.regions_file=[isimip_data_dir filesep params.regions_file];end

if strcmpi(entity,'params'),res=params;return;end % special case, return the full params structure

check_plot=0;if params.check_plot>0 || params.check_plot==-1,check_plot=1;end % plot damage data
damage_data=isimip_damage_read(params.damage_data_file,params.price_deflator_file,entity.hazard.ID_no,check_plot);

if ~isempty(damage_data)
    
           
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
        
        % now, we have
        % damage_data.damage_tot(record_i): total damage for record_i
        % damage_data.year(record_i): year
        % damage_data.ID_no(record_i): ibtracs ID for damage record
        % entity.hazard.ID_no(event_i): ibtracs ID for event_i
        
        inflator_data.year=1859:2100; % init dummy
        inflator_data.inflator_factor=inflator_data.year*0+1; % init dummy
        if ~isempty(params.price_deflator_file)
            if exist(params.price_deflator_file,'file')
                if isempty(strfind(params.price_deflator_file,num2str(params.inflation_reference_year)))
                    fprintf('WARNING: price inflation reference year (%i) might not be matching data in %s\n',...
                        params.inflation_reference_year,params.price_deflator_file)
                end
                inflator_data=climada_csvread(params.price_deflator_file);
                if isfield(inflator_data,'GDP_deflator_base2005')
                    inflator_data.GDP_deflator_base=inflator_data.GDP_deflator_base2005;
                    inflator_data=rmfield(inflator_data,'GDP_deflator_base2005');
                end
                inflator_data.inflator_factor=1./(inflator_data.GDP_deflator_base/100);
            else
                fprintf('WARNING: %s not found, no price inflation\n',params.price_deflator_file)
            end % exist(params.damage_data_file,'file')
        end % ~isempty(params.damage_data_file)
        
        % inlfate/deflate damages to reference year (see params.inflation_reference_year)
        unique_years=unique(damage_data.year); % unique list of years
        for year_i=1:length(unique_years)
            inflator_pos=find(inflator_data.year==unique_years(year_i));
            damage_pos=find(damage_data.year==unique_years(year_i));
            damage_data.inflator_factor(damage_pos)=inflator_data.inflator_factor(inflator_pos); % to keep track
        end % year_i
        damage_data.damage_tot=damage_data.hist_tot_losses.*damage_data.inflator_factor;
        damage_data.damage_ins=damage_data.hist_insured_losses.*damage_data.inflator_factor;
        
        if params.check_plot>0 || params.check_plot==-1 % plot damage data
            figure('Name','reported total damages','Color',[1 1 1]);
            yyaxis left
            % plot(damage_data.year,damage_data.damage_tot     ,'og');hold on % linear plots
            % plot(damage_data.year,damage_data.hist_tot_losses,'.b');
            % xlabel('year');ylabel(['damage [USD ' num2str(params.inflation_reference_year) ']']);
            plot(damage_data.year,log10(damage_data.damage_tot)     ,'og');hold on % log10 plots
            plot(damage_data.year,log10(damage_data.hist_tot_losses),'.b');
            xlabel('year');ylabel(['log10(damage) [USD ' num2str(params.inflation_reference_year) ']']);
            yyaxis right
            plot(damage_data.year,damage_data.inflator_factor,'-r');hold on
            ylabel('inflation factor');
            legend({'inflated','historic'});
            title('reported total damages');
        end % params.check_plot
        
        NatID_RegID=isimip_NatID_RegID; % obtain all country NatIDs and region IDs
        
        % damages, map countries to regions
        n_NatIDs=length(NatID_RegID.NatID);
        damage_data.RegID  =damage_data.NatID*0-1; % init
        damage_data.RegName=cell(1,length(damage_data.RegID)); % init
        for NatID_i=1:n_NatIDs
            pos=find(damage_data.NatID == NatID_RegID.NatID(NatID_i));
            if ~isempty(pos)
                damage_data.RegID(pos)   = NatID_RegID.TCRegID(NatID_i);
                damage_data.RegName(pos) = NatID_RegID.TCRegName(NatID_i);
            end
        end % NatID_i
        pos=find(damage_data.RegID<0);
        if ~isempty(pos),fprintf('WARNING: %i NatIDs not matched to RegIDs\n',length(pos));end
        
        % match ibtracs ID_no
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
        
        fprintf('%i of %i (%i%%) records matched with hazard events (%i%% of total damage value)\n',...
            matched,n_records,ceil(matched/n_records*100),...
            ceil(sum(damage_data.damage_tot(damage_data.hazard_matched))/sum(damage_data.damage_tot)*100));
        
        % only keep matched damage records
        damage_data=climada_subarray(damage_data,damage_data.hazard_matched);
        
        if params.store_to_entity
            fprintf('> matched damage data stored into %s\n',entity.assets.filename);
            entity.damage_data=damage_data;
            save(entity.assets.filename,'entity','-v7.3'); % -v7.3 for size..
        else
            [fP,fN,fE]=fileparts(params.damage_data_file);
            damage_data_save_file=[fP filesep fN '.mat'];
            fprintf('> matched damage data stored as %s\n',damage_data_save_file);
            save(damage_data_save_file,'damage_data');
        end
        
        % % experimental section to replace damage function
        % [damagefunctions,dmf_info_str]=climada_damagefunctions_generate(0:5:120,25,1.5,1,'s-shape','TC',0,100); % fits globally well
        % fprintf('replacing TC damagefunction with: %s\n',dmf_info_str);
        % entity=climada_damagefunctions_replace(entity,damagefunctions);
        
        % calculate simulated damages
        if isfield(entity,'hazard')
            %hazard=entity.hazard;entity=rmfield(entity,'hazard');
            %[~,EDS]=isimip_YDS_calc(entity,hazard);
            [~,EDS]=isimip_YDS_calc(entity,entity.hazard);
        elseif ~isempty(params.hazard_file)
            hazard=climada_hazard_load(params.hazard_file);
            if isempty(hazard),fprintf('ERROR: hazard %s not found\n',params.hazard_file);return,end
            [~,EDS]=isimip_YDS_calc(entity,hazard);
        end
        
        % FAST for TEST, comment section above
        % avoids using isimip_YDS_calc, assets fix to params.inflation_reference_year,
        % does not support regions, though
        %pos=find(entity.assets.Values_yyyy==params.inflation_reference_year);
        %if ~isempty(pos),entity.assets.Value=entity.assets.Values(pos,:);end
        %EDS=climada_EDS_calc(entity,entity.hazard);
        
        n_EDS=length(EDS);
        
        % map EDS countries to regions
        if isfield(EDS(1),'NatID') % EDS has a NatID, can be mapped to region
            RegID_used=zeros(1,n_EDS);
            NatID_used=zeros(1,n_EDS);
            ISO3_used=cell(1,n_EDS);
            RegName_used=cell(1,n_EDS);
            for EDS_i=1:n_EDS
                pos=find(NatID_RegID.NatID==EDS(EDS_i).NatID);
                if length(pos)==1
                    EDS(EDS_i).RegID   = NatID_RegID.TCRegID(pos);
                    EDS(EDS_i).RegName = NatID_RegID.TCRegName{pos};
                    RegID_used(EDS_i)  = EDS(EDS_i).RegID;
                    NatID_used(EDS_i)  = EDS(EDS_i).NatID;
                    ISO3_used{EDS_i}   = EDS(EDS_i).comment; % sorry comment for historic reasons
                    RegName_used{EDS_i}= EDS(EDS_i).RegName;
                else
                    fprintf('WARNING: EDS(%i) for %s not mapped to a region\n',EDS_i,EDS(EDS_i).comment);
                end
            end % EDS_i
            [unique_RegID_used,i]=unique(RegID_used);
            if unique_RegID_used(1)==0 && length(i)>1,unique_RegID_used=unique_RegID_used(2:end);i=i(2:end);end
            fprintf('NOTE: RegID = 0 omitted\n');
            unique_RegName_used=RegName_used(i);
            n_regions=length(unique_RegID_used);
        else
            n_regions=0; % just global view
        end
        
        if n_regions>0
            fprintf('processing %i regions: ',n_regions)
            for region_i=1:n_regions;fprintf('%s ',unique_RegName_used{region_i}),end;fprintf('\n')
        end
        
        % construct the year list for simulated damages
        damage_sim_year=entity.hazard.yyyy(damage_data.hazard_index);
                        
        for region_i=0:n_regions
            
            if region_i==0 % 'global', sum up over all EDSs
                damage_sim=EDS(1).damage;
                for EDS_i=2:length(EDS)
                    damage_sim=damage_sim+EDS(EDS_i).damage;
                end % EDS_i
                damage_tot      = damage_data.damage_tot;
                damage_tot_year = damage_data.year;
                RegName='all';region_ISO3_str='all countries';
            else
                % sum up over countries within region_i
                damage_sim=EDS(1).damage*0; % init empty
                region_ISO3_str=''; % init
                %fprintf('region %s (Id %i): ',unique_RegName_used{region_i},unique_RegID_used(region_i));
                for EDS_i=1:n_EDS
                    if EDS(EDS_i).RegID==unique_RegID_used(region_i)
                        damage_sim=damage_sim+EDS(EDS_i).damage;
                        region_ISO3_str=[region_ISO3_str EDS(EDS_i).comment ' '];
                    end
                end % EDS_i
                region_ISO3_str = deblank(region_ISO3_str);
                RegName=unique_RegName_used{region_i};
                pos=find(damage_data.RegID==unique_RegID_used(region_i));
                damage_tot      = damage_data.damage_tot(pos);
                damage_tot_year = damage_data.year(pos);
            end
            fprintf('%s: %s\n',RegName,region_ISO3_str);
            
            % reduce simulated damage to events we have reported damage data for
            damage_sim=damage_sim(damage_data.hazard_index);
            damage_sim(damage_sim<1)=0; % avoid damages less than 1 USD
            
            % % global per event view: match calculated damages with reported ones
            % if params.check_plot>0 || params.check_plot==-2 % plot damage data
            %
            %     if region_i==0
            %         plot_2_fig=figure('Name',sprintf('ALL: event reported versus simulated damages'),'Color',[1 1 1]);
            %     elseif region_i>0
            %         if region_i==1
            %             plot_2_fig=figure('Name',sprintf('regions: event reported versus simulated damages'),'Color',[1 1 1]);
            %         else
            %             figure(plot_2_fig)
            %         end
            %         subplot(3,3,region_i)
            %     end
            %
            %     % plot per year (does not work, since events not summed up)
            %     % yyaxis left
            %     % plot(damage_tot_year,log10(damage_tot),'.g');hold on
            %     % xlabel('year');ylabel('log10(reported damage) [USD 2005]');
            %     % yyaxis right
            %     % plot(damage_sim_year,log10(damage_sim),'.r');hold on
            %     % ylabel('log10(simulated damage) [USD]');
            %     % legend({'reported','simulated'});
            %     % title(sprintf('%s',RegName));
            %
            %     % scatter plot
            %     plot(log10(damage_tot),log10(damage_sim),'xr');hold on
            %     plot(log10(damage_tot),log10(damage_tot),'-g');
            %     log10_max_damage=max(max(log10(damage_tot)),max(log10(damage_sim)));
            %     plot([0 log10_max_damage],[0 log10_max_damage],':g');
            %     xlabel('log10(reported damage) [USD 2005]');ylabel('log10(simulated damage) [USD]')
            %     title(sprintf('%s years %i..%i',RegName,min(year_unique),max(year_unique)));
            %
            % end % params.check_plot
            
            % per year view: match calculated damages with reported ones
            if params.check_plot>0 || params.check_plot==-3 % plot damage data
                
                if region_i==0
                    plot_3_fig=figure('Name',sprintf('ALL: year reported versus simulated damages'),'Color',[1 1 1]);
                elseif region_i>0
                    if region_i==1
                        plot_3_fig=figure('Name',sprintf('regions: year reported versus simulated damages'),'Color',[1 1 1]);
                    else
                        figure(plot_3_fig)
                    end
                    subplot(3,3,region_i)
                end
                
                % sum data up over years
                year_unique=unique(damage_tot_year); % unique reported damage years
                damage_tot_yearsum = year_unique*0;
                damage_sim_yearsum = year_unique*0;
                for year_i=1:length(year_unique)
                    pos=find(damage_tot_year==year_unique(year_i));
                    if ~isempty(pos)
                        damage_tot_yearsum(year_i) = sum(damage_tot(pos));
                        damage_sim_yearsum(year_i) = sum(damage_sim(pos));
                    end
                end % year_i
                
                %             figure('Name',sprintf('%s: year reported versus simulated damages',region_ISO3_str),'Color',[1 1 1]);
                %             yyaxis left
                %             plot(year_unique,log10(damage_tot_yearsum),'.g');hold on
                %             xlabel('year');ylabel('log10(reported damage) [USD 2005]');
                %             yyaxis right
                %             plot(year_unique,log10(damage_sim_yearsum),'.r');hold on
                %             ylabel('log10(simulated damage) [USD]');
                %             legend({'reported','simulated'});
                %             title(region_ISO3_str);
                
                % and the scatter plot
                plot(log10(damage_tot_yearsum),log10(damage_sim_yearsum),'xr');hold on
                plot(log10(damage_tot_yearsum),log10(damage_tot_yearsum),'-g');
                log10_max_damage=max(max(log10(damage_tot_yearsum)),max(log10(damage_sim_yearsum)));
                plot([0 log10_max_damage],[0 log10_max_damage],':g');
                xlabel(['log10(reported damage) [USD ' num2str(params.inflation_reference_year) ']']);ylabel('log10(simulated damage) [USD]')
                title(sprintf('%s years %i..%i',RegName,min(year_unique),max(year_unique)));
                %title(sprintf('%s years %i..%i',region_ISO3_str,min(year_unique),max(year_unique)));
                
            end % params.check_plot
            
            if params.check_plot>0 || params.check_plot==-4 % plot DFC of damage data
                
                if region_i==0
                    plot_4_fig=figure('Name',sprintf('ALL: DFC'),'Color',[1 1 1]);
                elseif region_i>0
                    if region_i==1
                        plot_4_fig=figure('Name',sprintf('regions: DFC'),'Color',[1 1 1]);
                    else
                        figure(plot_4_fig)
                    end
                    subplot(3,3,region_i)
                end
                
                % and the scatter plot
                year_freq=max(year_unique)-min(year_unique)+1;
                frequency=ones(1,length(damage_tot_yearsum))*1/year_freq;
                
                [sorted_damage_tot,exceedence_freq_tot]=climada_damage_exceedence(damage_tot_yearsum,frequency);
                nonzero_pos         = find(exceedence_freq_tot);
                sorted_damage_tot   = sorted_damage_tot(nonzero_pos);
                exceedence_freq_tot = exceedence_freq_tot(nonzero_pos);
                return_period_tot   = 1./exceedence_freq_tot;
                
                [sorted_damage_sim,exceedence_freq_sim]=climada_damage_exceedence(damage_sim_yearsum,frequency);
                nonzero_pos         = find(exceedence_freq_sim);
                sorted_damage_sim   = sorted_damage_sim(nonzero_pos);
                exceedence_freq_sim = exceedence_freq_sim(nonzero_pos);
                return_period_sim   = 1./exceedence_freq_sim;
                
                %plot(return_period_tot,sorted_damage_tot);hold on
                %plot(return_period_sim,sorted_damage_sim);hold on
                plot(return_period_tot,log10(sorted_damage_tot));hold on
                plot(return_period_sim,log10(sorted_damage_sim));hold on
                legend({'reported','simulated'});
                %ylabel(['damage [USD ' num2str(params.inflation_reference_year) ']']);xlabel('years')
                ylabel(['log10(damage) [USD ' num2str(params.inflation_reference_year) ']']);xlabel('years')
                title(sprintf('%s years %i..%i',RegName,min(year_unique),max(year_unique)));
                
            end % params.check_plot
            
        end % region_i
        
    else
        fprintf('ERROR: %s not found\n',params.damage_data_file)
    end % exist(params.damage_data_file,'file')
end % ~isempty(params.damage_data_file)

%end % isimip_tc_calibrate