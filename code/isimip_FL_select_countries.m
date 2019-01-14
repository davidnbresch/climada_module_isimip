function [flag,output_table]=isimip_FL_select_countries(output_file,output_overwrite,years_range)
% climada isimip select countries to be used for floods in ISIMIP for
%   calibration and extrapolation, respectively
% MODULE:
%   isimip
% NAME:
%   isimip_FL_select_countries
% PURPOSE:
%   Indicates which of the ISIMIP list of countries which ones will be used for
%   calibration and extrapolation, respectively, for flooding. Essentially,
%   it depends on the output of emdat_get_country_names. The output is
%   written in file output_file, which duplicates file NatID_RegID_isimip_flood.csv
%
%   next call: isimip_flood_calibration
% CALLING SEQUENCE:
%   flag=isimip_FL_select_countries(output_file)
% EXAMPLE:
%   flag=isimip_FL_select_countries('',0,[1992 2010]);
% INPUTS:
%   (no requested input)
% OPTIONAL INPUT PARAMETERS:
%   output_file: name of the file to be written out, e.g. 'NatID_RegID_isimip_flood_filtered.csv'
%   output_overwrite: if =0 (default), only write output_file if it does not yet exist. If=1, overwrite anyway.
%   years_range: Range of years to be used. Default = [1992 2010].
% OUTPUTS:
%   flag: =1 if successful; otherwise =0.
%   output_table: table of the full output saved in file output_file, which
%      contains the following columns:
%      ISO,ID,Reg_ID,Reg_name,in_calib,in_extrap
%      in_calib: =-1 if never because country does not exist (or is not
%                    valid)
%                =n if used in calibration for params.keep_countries_0emdat <=n
%                (i.e., 2 means always in;
%                       1 means only in if params.keep_countries_0emdat is 0 or 1;
%                       0 means only if params.keep_countries_0emdat is 0).
%      in_extrap: 0 if the country is not recognized, otherwise =1 (only ANT).
% MODIFICATION HISTORY:
% Benoit P. Guillod, benoit.guillod@env.ethz.ch, 20190114, initial
%-

global climada_global

%% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('years_range','var'),             years_range=    [1992 2010];end
if ~exist('output_file','var'),output_file=['NatID_RegID_isimip_flood_filtered_' num2str(years_range(1)) '-' num2str(years_range(2)) '.csv'];end
if ~exist('output_overwrite','var'),output_overwrite=false;end

% init output
flat=0;
output_table=table;

%% read country/region data
% find file
RegID_def_folder=[climada_global.data_dir filesep 'isimip'];
NatID_RegID_file = [RegID_def_folder filesep 'NatID_RegID_isimip_flood.csv'];
NatID_RegID_flood = readtable(NatID_RegID_file);
NatID_RegID_flood.Reg_name = string(NatID_RegID_flood.Reg_name);
countries_iso3=NatID_RegID_flood.ISO;

%% find filter
% prepare outputs
n_rows=length(countries_iso3);
in_calib=NaN([n_rows 1]);
in_extrap=true([n_rows 1]);
% check each country
for i=1:length(countries_iso3)
    [~,~,cl,em] = emdat_get_country_names(countries_iso3{i},['FL';'F1';'F2'],years_range,0);
    if ismember(cl, [-3 -1 99])
        % country does not exist or is not found
        in_extrap(i) = false;
        in_calib(i) = -1;
    elseif cl == -2
        % country is not in EM-DAT, not used in calibration but can be
        % extrapolated
        in_calib(i) = -1;
    else
        % get em+1, max 2
        in_calib(i) = min([em+1 2]);
    end
end


%% create final output
output_table=NatID_RegID_flood;
output_table.in_calib=in_calib;
output_table.in_extrap=in_extrap;

%% save file (if appropriate)
[fp,fn,ext] = fileparts(output_file);
if isempty(fp),fp=RegID_def_folder;end
output_file=[fp filesep fn ext];
if exist(output_file,'file')
    if ~output_overwrite
        fprintf('** WARNING ** file %s already exists; not overwritten *****',output_file)
    else
        fprintf('** WARNING ** file %s overwritten *****',output_file)
        writetable(output_table,output_file)
    end
else
    writetable(output_table,output_file)
    fprintf('file %s has been written *****',output_file)
end
flag=1;

end
