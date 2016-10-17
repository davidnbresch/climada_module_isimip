function entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot)
% climada isimip entity population
% MODULE:
%   isimip
% NAME:
%   isimip_ssp2_entity
% PURPOSE:
%   create the entity based on isimip ssp (population) data
%
%   Programmer's note: later to be upgraded to work with a list of
%   countries, in order to avoid loading all files each time...
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
%   check_plot: whether show a check plot (=1, default), or not (=0)
%       Note that plotting might often take longer than the full
%       conversion...
%       if =2, also show the centroids as red dots
%       if =3, also show the original data grid as blue dots (might take time...)
% OUTPUTS:
%   entity: a climada entity structure, see climada_entity_read for a full
%       description of all fields
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
if ~exist('check_plot','var'),    check_plot     =  1;end

% PARAMETERS
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
    fprintf('processing %i entities\n',length(ISO3_list));
    for iso3_i=1:length(ISO3_list)
        ISO3=ISO3_list(iso3_i,1);
        entity=isimip_ssp2_entity(ISO3,ssp2_filename,NatId_filename,check_plot);
    end
    fprintf('only last entity returned, see climada_entity_load\n');
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

% read the NatId grid
nc.NatIdGrid = ncread(NatId_filename,'NatIdGrid');

NatId_pos=find(nc.NatIdGrid == NatId);

if ~isempty(NatId_pos) % there are gridcells within this country
    
    % rewad the population data
    fprintf('reading lon, lat, time and var1 from %s ...',ssp2_filename);
    % if troubles, use ncinfo(flood_fraction_filename,'var') ...
    nc.lon      = ncread(ssp2_filename,'lon');
    nc.lat      = ncread(ssp2_filename,'lat');
    nc.time     = ncread(ssp2_filename,'time');
    nc.var1     = ncread(ssp2_filename,'var1');
    fprintf(' done\n');
    
    if sum(size(nc.var1(:,:,1))-size(nc.NatIdGrid)) == 0
        % grid sizes the same
        
        % get template entity
        entity=climada_entity_load(entity_template);
        
        entity_filename=[climada_global.entities_dir filesep ISO3 '_entity'];
        
        % create the grid
        [gridlon,gridlat] = meshgrid(nc.lon,nc.lat);
        
        entity.assets.lon=gridlon(NatId_pos);
        entity.assets.lat=gridlat(NatId_pos);
        for time_i=1:length(nc.time)
            temp_var1=var1(NatId_pos,time_i);
            entity.assets.Value=temp_var1(NatId_pos);
        end  % time_i
        
        % complete entity
        entity.assets.Cover=entity.assets.Value;
        entity.assets.Deductible=entity.assets.Value*0;
        entity.assets.DamageFunID=entity.assets.Deductible+1;
        
        %fprintf('saving entity as %s\n',entity_filename)
        %save(entity_filename,'entity');
        
        if check_plot,climada_entity_plot(entity);end
        
    else
        fprintf('ERROR: grid sizes do not match, aborted\n');
    end
    
end % ~isempty(NatId_pos)

end % isimip_ssp2_entity