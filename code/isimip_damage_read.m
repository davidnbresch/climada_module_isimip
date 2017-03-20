function damage_data=isimip_damage_read(damage_data_file,price_deflator_file,hazard_ID_no,check_plot)
% climada template
% MODULE:
%   isimip
% NAME:
%   isimip_damage_read
% PURPOSE:
%   read isimip reported damage records (from .xlsx)
%
%   Note: use and adapt this code to read other damage data
%
%   previous call: -
%   next call: isimip_tc_calibrate
% CALLING SEQUENCE:
%   damage_data=isimip_damage_read(damage_data_file,price_deflator_file,hazard_ID_no,check_plot)
% EXAMPLE:
%   damage_data=isimip_damage_read; % read default isimip damage table
%   damage_data_file='matching_Natcat-damages_ibtracs_1980-2014.csv';
%   price_deflator_file='GDP_deflator_converted_base2005_1969-2016_source_BEA.csv';
%   hazard=climada_hazard_load('GLB_0360as_TC_hist')
%   damage_data=isimip_damage_read(damage_data_file,price_deflator_file,hazard.ID_no);
% INPUTS:
%   damage_data_file: file with damage data, needs to have some fields, see
%       documentation. Path added, if not provided
%   price_deflator_file: file to deflate/inflate historic damage data
%       Needs to have some fields, see documentation. Path added, if not provided
% OPTIONAL INPUT PARAMETERS:
%   hazard_ID_no(event_i): the ID_no for event i, we match reported damage
%       and events according to this ID_no and only return matched records.
%   check_plot_ if=1, plot raw and inflated damages (default=0).
% OUTPUTS:
%   damage_data: the output, a struct, empty if not successful, key fields are:
%       damage(i): damage for record i
%       year(i): year of record i
%       hazard_index(i): index into hazard.*(hazard_index(i)), if
%           hazard_ID_no has been provided, otherwise not returned
%      optional (in decreasing order of importance):
%       ISO{i}: ISO3 country code of record i
%       ID_no(i): ibtracs ID no, 1950166N14262 converted to 1950166114262
%           (just converted to integer, after N->1, S->0)
%       damage_ins(i): insured damage record i
%       NatID(i): isimip NatID
%       inflator_factor(i): de/inflation factor applied to record i
%       RegID(i): region ID for record i
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20170305, initial
% David N. Bresch, david.bresch@gmail.com, 20170314, new damage file format adopted
% David N. Bresch, david.bresch@gmail.com, 20170315, switched to damage.damage to be fully consistent with climada
% David N. Bresch, david.bresch@gmail.com, 20170320, ID_no as integer
%-

damage_data=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('damage_data_file','var'),damage_data_file='';end
if ~exist('price_deflator_file','var'),price_deflator_file='';end
if ~exist('hazard_ID_no','var'),hazard_ID_no='';end
if ~exist('check_plot','var'),check_plot=0;end

% locate the module's (or this code's) data folder (usually  a folder
% 'parallel' to the code folder, i.e. in the same level as code folder)
%module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% define all parameters here - no parameters to be defined in code below
%
% define the defaut folder for isimip TC track data
isimip_data_dir=[climada_global.data_dir filesep 'isimip'];
if ~isdir(isimip_data_dir)
    mkdir(climada_global.data_dir,'isimip'); % create it
    fprintf('NOTE: store your isimip input data in %s\n',isimip_data_dir);
end
%
if isempty(damage_data_file),damage_data_file=      [climada_global.data_dir filesep 'isimip' filesep ...
        'matching_Natcat-damages_ibtracs_1980-2014.csv'];end
if isempty(price_deflator_file),price_deflator_file=[climada_global.data_dir filesep 'isimip' filesep ...
        'GDP_deflator_converted_base2005_1969-2016_source_BEA.csv'];end

% prepend isimip_data_dir in case only filenames are passed
if ~exist(damage_data_file,'file'),damage_data_file=[isimip_data_dir filesep damage_data_file];end
if ~exist(price_deflator_file,'file'),price_deflator_file=[isimip_data_dir filesep price_deflator_file];end

if ~isempty(damage_data_file)
    if exist(damage_data_file,'file')
        damage_data=climada_csvread(damage_data_file);
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
        if isfield(damage_data,'CountryID'),damage_data.NatID=damage_data.CountryID;damage_data=rmfield(damage_data,'CountryID');end
        
        % a unique ID, a number for faster comparison, hence convert 1950166N14262 to 1950166014262
        n_damage_events=length(damage_data.ibtracs_ID);
        damage_data.ID_no=zeros(1,n_damage_events);
        %damage_data.year=zeros(1,n_damage_events); % OLD, really bad, only 2-digit year stored as dd/mm/yy       
        damage_data.year=fix(damage_data.Date./10000); %yyyymmdd

        for event_i=1:n_damage_events
            if length(damage_data.ibtracs_ID{event_i})>9
                %damage_data.ID_no(event_i)=str2double(damage_data.ibtracs_ID{event_i}(1:7))+str2double(damage_data.ibtracs_ID{event_i}(9:end))/100000;
                damage_data.ID_no(event_i)=str2double(strrep(strrep(damage_data.ibtracs_ID{event_i},'S','1'),'N','0')); % N->0, S->1
            end % length(damage_data.ibtracs_ID{event_i})>9
            % OLD, really bad, only 2-digit year stored as dd/mm/yy
            %year2digit=str2double(damage_data.Date{event_i}(7:end));
            %if year2digit<30
            %    damage_data.year(event_i)=2000+year2digit;
            %else
            %    damage_data.year(event_i)=1900+year2digit;
            %end
        end % event_i
        
        % now, we have
        % damage_data.damage(record_i): total damage for record_i
        % damage_data.year(record_i): year
        % damage_data.ID_no(record_i): ibtracs ID for damage record
        % entity.hazard.ID_no(event_i): ibtracs ID for event_i
        
        inflator_data.year=1859:2100; % init dummy
        inflator_data.inflator_factor=inflator_data.year*0+1; % init dummy
        if ~isempty(price_deflator_file)
            if exist(price_deflator_file,'file')
                inflator_data=climada_csvread(price_deflator_file);
                if isfield(inflator_data,'GDP_deflator_base2005')
                    inflator_data.GDP_deflator_base=inflator_data.GDP_deflator_base2005;
                    inflator_data=rmfield(inflator_data,'GDP_deflator_base2005');
                end
                inflator_data.inflator_factor=1./(inflator_data.GDP_deflator_base/100);
            else
                fprintf('WARNING: %s not found, no price inflation\n',params.price_deflator_file)
            end % exist(price_deflator_file,'file')
        end % ~isempty(price_deflator_file)
        
        % inlfate/deflate damages to reference year (depends on the factors in the file)
        unique_years=unique(damage_data.year); % unique list of years
        for year_i=1:length(unique_years)
            inflator_pos=inflator_data.year==unique_years(year_i);
            damage_pos=damage_data.year==unique_years(year_i);
            damage_data.inflator_factor(damage_pos)=inflator_data.inflator_factor(inflator_pos); % to keep track
        end % year_i
        damage_data.damage=damage_data.hist_tot_losses.*damage_data.inflator_factor;
        damage_data.damage_ins=damage_data.hist_insured_losses.*damage_data.inflator_factor;
        
        if check_plot % plot damage data
            figure('Name','reported total damages','Color',[1 1 1]);
            yyaxis left
            plot(damage_data.year,log10(damage_data.damage)     ,'og');hold on % log10 plots
            plot(damage_data.year,log10(damage_data.hist_tot_losses),'.b');
            xlabel('year');ylabel(['log10(damage) [USD]']);
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
        
        if ~isempty(hazard_ID_no)
            % match ibtracs ID_no
            matched=0;
            not_matched=0;
            n_records=length(damage_data.ID_no);
            damage_data.hazard_index  =damage_data.ID_no*0; % init
            damage_data.hazard_matched=damage_data.hazard_index; % init
            for record_i=1:n_records
                pos=find(hazard_ID_no==damage_data.ID_no(record_i));
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
                ceil(sum(damage_data.damage(damage_data.hazard_matched))/sum(damage_data.damage)*100));
            
            % only keep matched damage records
            damage_data=climada_subarray(damage_data,damage_data.hazard_matched);
                        
        end % ~isempty(hazard_ID_no)
               
        [fP,fN]=fileparts(damage_data_file);
        damage_data_save_file=[fP filesep fN '.mat'];
        fprintf('> matched damage data stored as %s\n',damage_data_save_file);
        save(damage_data_save_file,'damage_data');
    end % isempty(damage_data_file)
end % exist(damage_data_file,'file')

end % isimip_damage_read