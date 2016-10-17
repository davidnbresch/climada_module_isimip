function entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot)
% climada isimip entity population
% MODULE:
%   isimip
% NAME:
%   isimip_ssp2_entity
% PURPOSE:
%   create the entity based on isimip ssp (population) data
%
%   Checks the country ID (NatId) on NatId_filename and takes all gridcells
%   within the requested country. If the resolution of the NetId does not
%   match the population data or if there is no NatId_filename provided,
%   the code uses the climada country shape files to select the gridcells
%   within the country (the code notifies).
%
%   The asset values are then scaled by GDP and the income group, see
%   climada_entity_value_GDP_adjust_one.
%
%   Programmer's note: While the present code works for a list or even all
%   countries, it shall later to be upgraded to work more efficiently with
%   a list of countries, in order to avoid loading all files each time...  
%
%   Octave: please install the netCDF package first:
%   pkg install -forge netcdf -local -auto
%   or: pkg install -verbose -forge -auto netcdf
%
%   next call: isimip...
% CALLING SEQUENCE:
%   entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot)
% EXAMPLE:
%   entity=isimip_ssp2_entity('DEU')
% INPUTS:
%   ISO3: the ISO3 country code, ='all': process all countries (be careful,
%       test a few first)
%       > promted for (to select from a list) if not given
% OPTIONAL INPUT PARAMETERS:
%   ssp2_filename: filename of the .nc file with the population values
%   NatId_filename: filename of the .nc file with the national grid IDs (to
%       assign grid cells to countries). Default set in PARAMETERS
%       To avoid using NatId and just use the grid cells within a country
%       shape, set NatId_filename='ignore'.
%   check_plot: whether show a check plot (=1, default), or not (=0)
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
%       PLUS the field entity.assets.Values(n_times,n_centroids) with the
%       variable as on the ssp2 file for each timestep at each centroid.
%       These values are NOT scaled (different from entity.assets.Value,
%       which is scaled by GDP and income group).   
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161017, initial
%-

entity=[]; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('ISO3','var'),          ISO3           = '';end
if ~exist('ssp2_filename','var'), ssp2_filename  = '';end
if ~exist('NatId_filename','var'),NatId_filename = '';end
if ~exist('check_plot','var'),    check_plot     =  0;end

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
%
% define default filenames
if isempty(ssp2_filename),ssp2_filename=[isimip_data_dir filesep 'hyde_ssp2_1860-2100_0.1deg_yearly_2005-2008.nc4'];end
if isempty(NatId_filename),NatId_filename=[isimip_data_dir filesep 'Nat_id_grid_ISIMIP.nc'];end
%
% the entity template to populate a defaiult entity
entity_template=[climada_global.entities_dir filesep 'entity_template'];
%
% admin0 shape file (for fallback selection option):
admin0_shape_file=climada_global.map_border_file;
%
% define the isimip country ISO3 code list
ISO3_list={
    'ABW'    1
    'AFG'    2
    'AGO'    3
    'AIA'    4
    'ALB'    5
    'AND'    6
    'ANT'    7
    'ARE'    8
    'ARG'    9
    'ARM'    10
    'ASM'    11
    'ATG'    12
    'AUS'    13
    'AUT'    14
    'AZE'    15
    'BDI'    16
    'BEL'    17
    'BEN'    18
    'BFA'    19
    'BGD'    20
    'BGR'    21
    'BHR'    22
    'BHS'    23
    'BIH'    24
    'BLR'    25
    'BLZ'    26
    'BMU'    27
    'BOL'    28
    'BRA'    29
    'BRB'    30
    'BRN'    31
    'BTN'    32
    'BWA'    33
    'CAF'    34
    'CAN'    35
    'CHE'    36
    'CHL'    37
    'CHN'    38
    'CIV'    39
    'CMR'    40
    'COD'    41
    'COG'    42
    'COK'    43
    'COL'    44
    'COM'    45
    'CPV'    46
    'CRI'    47
    'CUB'    48
    'CYM'    49
    'CYP'    50
    'CZE'    51
    'DEU'    52
    'DJI'    53
    'DMA'    54
    'DNK'    55
    'DOM'    56
    'DZA'    57
    'ECU'    58
    'EGY'    59
    'ERI'    60
    'ESP'    61
    'EST'    62
    'ETH'    63
    'FIN'    64
    'FJI'    65
    'FLK'    66
    'FRA'    67
    'FRO'    68
    'FSM'    69
    'GAB'    70
    'GBR'    71
    'GEO'    72
    'GGY'    73
    'GHA'    74
    'GIB'    75
    'GIN'    76
    'GLP'    77
    'GMB'    78
    'GNB'    79
    'GNQ'    80
    'GRC'    81
    'GRD'    82
    'GTM'    83
    'GUF'    84
    'GUM'    85
    'GUY'    86
    'HKG'    87
    'HND'    88
    'HRV'    89
    'HTI'    90
    'HUN'    91
    'IDN'    92
    'IMN'    93
    'IND'    94
    'IRL'    95
    'IRN'    96
    'IRQ'    97
    'ISL'    98
    'ISR'    99
    'ITA'    100
    'JAM'    101
    'JEY'    102
    'JOR'    103
    'JPN'    104
    'KAZ'    105
    'KEN'    106
    'KGZ'    107
    'KHM'    108
    'KIR'    109
    'KNA'    110
    'KOR'    111
    'KWT'    112
    'LAO'    113
    'LBN'    114
    'LBR'    115
    'LBY'    116
    'LCA'    117
    'LIE'    118
    'LKA'    119
    'LSO'    120
    'LTU'    121
    'LUX'    122
    'LVA'    123
    'MAC'    124
    'MAR'    125
    'MCO'    126
    'MDA'    127
    'MDG'    128
    'MDV'    129
    'MEX'    130
    'MHL'    131
    'MKD'    132
    'MLI'    133
    'MLT'    134
    'MMR'    135
    'MNG'    136
    'MNP'    137
    'MOZ'    138
    'MRT'    139
    'MSR'    140
    'MTQ'    141
    'MUS'    142
    'MWI'    143
    'MYS'    144
    'MYT'    145
    'NAM'    146
    'NCL'    147
    'NER'    148
    'NFK'    149
    'NGA'    150
    'NIC'    151
    'NIU'    152
    'NLD'    153
    'NOR'    154
    'NPL'    155
    'NRU'    156
    'NZL'    157
    'OMN'    158
    'PAK'    159
    'PAN'    160
    'PCN'    161
    'PER'    162
    'PHL'    163
    'PLW'    164
    'PNG'    165
    'POL'    166
    'PRI'    167
    'PRK'    168
    'PRT'    169
    'PRY'    170
    'PSE'    171
    'PYF'    172
    'QAT'    173
    'REU'    174
    'ROU'    175
    'RUS'    176
    'RWA'    177
    'SAU'    178
    'SCG'    179
    'SDN'    180
    'SEN'    181
    'SGP'    182
    'SHN'    183
    'SJM'    184
    'SLB'    185
    'SLE'    186
    'SLV'    187
    'SMR'    188
    'SOM'    189
    'SPM'    190
    'STP'    191
    'SUR'    192
    'SVK'    193
    'SVN'    194
    'SWE'    195
    'SWZ'    196
    'SYC'    197
    'SYR'    198
    'TCA'    199
    'TCD'    200
    'TGO'    201
    'THA'    202
    'TJK'    203
    'TKL'    204
    'TKM'    205
    'TLS'    206
    'TON'    207
    'TTO'    208
    'TUN'    209
    'TUR'    210
    'TUV'    211
    'TWN'    212
    'TZA'    213
    'UGA'    214
    'UKR'    215
    'URY'    216
    'USA'    217
    'UZB'    218
    'VCT'    219
    'VEN'    220
    'VGB'    221
    'VIR'    222
    'VNM'    223
    'VUT'    224
    'WLF'    225
    'WSM'    226
    'YEM'    227
    'ZAF'    228
    'ZMB'    229
    'ZWE'    230
    };

if strcmpi(ISO3,'all')
    % process all
    if verbose,fprintf('processing %i entities\n',length(ISO3_list));end
    for iso3_i=1:length(ISO3_list)
        ISO3=ISO3_list(iso3_i,1);
        entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot);
    end
    if verbose,fprintf('only last entity returned, see climada_entity_load\n');end
end

if isempty(ISO3)
    % prompt for
    [selection] = listdlg('PromptString','Select one country:',...
        'ListString',ISO3_list(:,1),'SelectionMode','Single');
    pause(0.1)
    if ~isempty(selection)
        ISO3=ISO3_list(selection,1);
    else
        return
    end
end

iso3_pos=strmatch(ISO3,ISO3_list(:,1));

if ~isempty(iso3_pos)
    NatId=ISO3_list{iso3_pos,2}; % get the ID
end % ~isempty(iso3_pos)

% read the population data
if verbose,fprintf('reading lon, lat, time and var1 from %s ...',ssp2_filename);end
% if troubles, use ncinfo(flood_fraction_filename,'var') ...
nc.lon       = ncread(ssp2_filename,'lon');
nc.lat       = ncread(ssp2_filename,'lat');
nc.time      = ncread(ssp2_filename,'time');
nc.var1      = ncread(ssp2_filename,'var1');
if verbose,fprintf(' done\n');end

% create the grid
if verbose,fprintf('creating regular grid (meshgrid) ...');end
[gridlon,gridlat] = meshgrid(nc.lon,nc.lat);
% reform to 1D list
gridlon=gridlon';
gridlat=gridlat';
gridlon=reshape(gridlon,[1 numel(gridlon)]);
gridlat=reshape(gridlat,[1 numel(gridlat)]);
if verbose,fprintf(' done\n');end

NatId_pos=[]; % init
if exist(NatId_filename,'file')
    % read the NatId grid
    nc.NatIdGrid = ncread(NatId_filename,'NatIdGrid');
    
    if sum(abs(size(nc.var1(:,:,1))-size(nc.NatIdGrid))) == 0 % grid sizes the same
        NatIdGrid=reshape(nc.NatIdGrid,[1 numel(nc.NatIdGrid)]);
        NatId_pos=find(NatIdGrid == NatId);
        n_centroids=length(NatId_pos);
        if verbose,fprintf('%i NatId grid cells within country %s\n',n_centroids,ISO3);end
    else
        fprintf('Warning: grid sizes do not match, 2nd approach (shape)\n');
    end
else
    fprintf('Warning: NatIdGrid file not found, 2nd approach (shape)\n');
end

if isempty(NatId_pos) % 2nd try, using shapes
    
    % read admin0 (country) shape file
    admin0_shapes=climada_shaperead(admin0_shape_file);
    
    % figure the country shape
    shape_i=[];
    for shape_ii=1:length(admin0_shapes)
        if strcmpi(admin0_shapes(shape_ii).ADM0_A3,ISO3),shape_i=shape_ii;end
    end % shape_ii
    
    if ~isempty(shape_i)
        NatId_pos=climada_inpolygon(gridlon,gridlat,admin0_shapes(shape_i).X,admin0_shapes(shape_i).Y,0);
        n_centroids=sum(NatId_pos);
        if n_centroids>0
            if verbose,fprintf('%i grid cells within country %s shape\n',n_centroids,ISO3);end
        else
            fprintf('ERROR: no gridpoints within shape, aborted\n');
            return
        end
    else
        fprintf('ERROR: %s not in shape file, aborted\n',ISO3);
    end
end

if ~isempty(NatId_pos)
    
    % get template entity
    entity=climada_entity_load(entity_template);
    entity_filename=[climada_global.entities_dir filesep ISO3 '_entity'];
    
    entity.assets.admin0_ISO3=ISO3;
    
    entity.assets.lon=gridlon(NatId_pos);
    entity.assets.lat=gridlat(NatId_pos);
    
    n_times=length(nc.time);
    if verbose,fprintf('extracting %i centroids at %i times ...',n_centroids,n_times);end
    entity.assets.Values=zeros(n_times,n_centroids);
    for time_i=1:n_times
        temp_var1=reshape(nc.var1(:,:,time_i),[1 numel(nc.var1(:,:,time_i))]); % get time slice
        entity.assets.Values(time_i,:)=temp_var1(NatId_pos);
    end  % time_i
    entity.assets.Value=entity.assets.Values(1,:); % copy first one
    if verbose,fprintf(' done\n');end
    
    % complete entity
    entity.assets.Cover      =entity.assets.Value;
    entity.assets.Deductible =entity.assets.Value*0;
    entity.assets.DamageFunID=entity.assets.Deductible+1;
    entity.assets.Category_ID=entity.assets.DamageFunID;
    entity.assets.Region_ID  =entity.assets.DamageFunID;
    
    entity.assets = rmfield(entity.assets,'centroid_index');
    entity.assets = rmfield(entity.assets,'Value_unit');
    entity.assets = rmfield(entity.assets,'hazard');
    entity.assets = climada_assets_complete(entity.assets);
    
    % scale up asset values based on a country's estimated total asset value
    entity=climada_entity_value_GDP_adjust_one(entity,verbose);
    entity.assets.isimip_comment='only entity.assets.Value scaled, entity.assets.Values not modified';
    
    if verbose,fprintf('saving entity as %s\n',entity_filename);end
    %save(entity_filename,'entity');
    
    if check_plot,climada_entity_plot(entity);end
    
end % ~isempty(NatId_pos)

end % isimip_ssp2_entity