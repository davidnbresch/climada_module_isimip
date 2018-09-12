function [status,output_filename]=isimip_flood_calibration(RegionID,years_range,params)
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
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% EXAMPLE:
%   RegionID='NAM';
%   years_range=[1990 2010];
%   params.entity_folder='/cluster/work/climate/dbresch/climada_data/isimip/entities';
%   params.entity_prefix='FL1950';
%   [status,output_filename]=isimip_flood_calibration(RegionID,years_range)
% INPUTS:
%   country: country name (full name or ISO3)
%   ghm: Global hydrological model name
%   forcing: observational forcing name
% OPTIONAL INPUT PARAMETERS:
%   params: a structure with fields:
%     entity_folder: the folder where the entities are located (default:
%        [climada_global.data_dir filesep 'isimip/entities'] ).
%     entity_prefix: if not ='', pre-pend the entity filename with it, e.g.
%        entity_prefix='Try1' will result in Try1_DEU_0150as.mat
%     hazard_protection: one of 'flopros' (default), '0', '100'
%     subtract_matsiro: =1 to subtract the 2-yr return value of MATSIRO flood
%        fraction from the data. Default =0.
% OUTPUTS:
%   status: 1 if successful, 0 if not.
%   output_filename: a file name for the .mat file generated.
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20180711, initial
%   
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('RegionID','var'),error('Input parameter RegionID is missing');end %
if ~exist('years_range','var'),             years_range=    [1990 2010];end
if ~exist('params','var'),                  params=              struct;end

% check for some parameter fields we need
if ~isfield(params,'entity_folder'),    params.entity_folder=[climada_global.data_dir filesep 'isimip/entities'];end
if ~isfield(params,'RegID_def_folder'), params.RegID_def_folder=[climada_global.data_dir filesep 'isimip'];end
if ~isfield(params,'entity_prefix'),    params.entity_prefix='FL1950';end
if ~isfield(params,'hazard_protection'),params.hazard_protection='flopros';end
if ~isfield(params,'subtract_matsiro'), params.subtract_matsiro=0;end

if ~isempty(params.entity_prefix)
    if ~strcmp(params.entity_prefix(end),'_'),params.entity_prefix=[params.entity_prefix '_'];end
end

% get countries that belong to the region
% for now, had oc for testing only
NatID_RegID_file = [params.RegID_def_folder filesep 'NatID_RegID_isimip_flood.csv'];
NatID_RegID_flood = readtable(NatID_RegID_file);
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
if sum(NatID_RegID_flood.Reg_name == RegionID)==0
    error('no country belonging to the given RegionID, perhaps non-existing RegionID?');
end
countries=NatID_RegID_flood.ISO(NatID_RegID_flood.Reg_name == RegionID);

% define variables and paths
isimip_simround = '2a';
isimip_data_dir = [climada_global.data_dir filesep 'isimip'];
countries_iso3={};
for i=1:length(countries)
    [~,countries_iso3{i}] =  climada_country_name(countries{i});
end

%steps:
% 1) load entities - N entities for N countries
entity_list={};
for i=1:length(countries)
    country_iso3 = countries_iso3{i};
    entity_file_isimip_i=[params.entity_folder filesep params.entity_prefix strtrim(country_iso3) '_0150as_entity'];
    entity_isimip_i=climada_entity_load(entity_file_isimip_i,1); % try to load, flag to 1 to avoir overwrite
    if isempty(entity_isimip_i)
        error('*** ERROR: entity file not found %s\n\n',entity_file_isimip_i);
    end
    entity_isimip_i.assets.centroid_index = 1:length(entity_isimip_i.assets.centroid_index);
    % check that all analysed years are included
    if sum(ismember(years_range(1):years_range(2), entity_isimip_i.assets.Values_yyyy)) ~= diff(years_range)+1
        error('*** ERROR: not all requested years are available in the entities **\n\n');
    end
    % hazard
%     hazard_yyyy_pos=find(entity.assets.Values_yyyy>=years_range(1) & entity.assets.Values_yyyy<=years_range(2));
    entity_list{i}=entity_isimip_i;
end


% 2) load hazards. 46*N hazards, as there are 46 model combinations
ghms = {'CLM', 'DBH', 'H08', 'JULES-TUC', 'JULES-UoE', 'LPJmL', 'MATSIRO', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'VEGAS', 'VIC', 'WaterGAP'};
forcings = {'gswp3', 'princeton', 'watch', 'wfdei'};
hazard_list={};
for i=1:length(countries)
    hazard_list{i} = {};
    ii = 0;
    for j=1:length(ghms)
        for k=1:length(forcings)
            ghm=ghms{j};
            forcing = forcings{k};
            [flddph_filename,fldfrc_filename,fld_path] = isimip_get_flood_filename(isimip_simround, ghm, forcing, params.hazard_protection, 'historical');
            flood_filename=[fld_path filesep flddph_filename];
            hazard_FL_file=isimip_get_flood_hazard_filename(flood_filename,entity_list{i},isimip_simround,[0 0],params.subtract_matsiro);
            if exist(hazard_FL_file, 'file')
                hazard_FL=climada_hazard_load(hazard_FL_file);
                hazard_FL.yyyy = double(string(hazard_FL.yyyy));
                hazard_yyyy_pos=find(hazard_FL.yyyy>=years_range(1) & hazard_FL.yyyy<=years_range(2));
                hazard_FL.intensity       =hazard_FL.intensity(hazard_yyyy_pos,:);
                hazard_FL.fraction        =hazard_FL.fraction(hazard_yyyy_pos,:);
                hazard_FL.frequency       =hazard_FL.frequency(hazard_yyyy_pos); % all one (i.e. we sum up)
                hazard_FL.orig_event_count=length(hazard_yyyy_pos);
                hazard_FL.event_count     =length(hazard_yyyy_pos);
                hazard_FL.event_ID        =hazard_FL.event_ID(hazard_yyyy_pos);
                hazard_FL.orig_event_flag =hazard_FL.orig_event_flag(hazard_yyyy_pos);
                hazard_FL.yyyy            =hazard_FL.yyyy(hazard_yyyy_pos);
                hazard_FL.orig_years      =1;
                ii=ii+1;
                hazard_list{i}{ii}=hazard_FL;
            else
                fprintf('     * FL hazard missing, probably inexistent model combination %s\n',flood_filename);
                continue
            end
        end
    end
end
        
% 3) load EM-DAT
emdat_list={};
all_years = hazard_FL.yyyy;
for i=1:length(countries)
    country_iso3 = countries_iso3{i};
    em_data_i=emdat_read('',country_iso3,['FL';'F1';'F2'],2005,0);
    emdat_damage_2005 = zeros([length(all_years) 1]);
    % if EM-DAT data available for this country, use, if not leave zeros
    if ~isempty(em_data_i)
        for iy=1:length(all_years)
            ii=find(all_years(iy) == em_data_i.year);
            if length(ii)>0
                emdat_damage_2005(iy,1) = sum(em_data_i.damage(ii));
            end
        end
    end
    emdat_list{i}.values=emdat_damage_2005;
    emdat_list{i}.year=all_years;
end


% 4) Define MDR function shape and parameters
FunctionHandle = str2func('make_MDR_function_1mExp');
MDR_func = @(x,scale,shape)FunctionHandle(x,scale,shape);
% maybe the range of the parameters should also be provided here?

% 5) Call calibrate_MDR_steps (TO DO)
% calibrate_MDR_steps(entity_list, hazard_list, emdat_list, MDR_func, ...)

end
