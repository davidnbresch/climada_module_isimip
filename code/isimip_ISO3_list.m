function ISO3_list=isimip_ISO3_list(ISO3_code,isimip_number)
% climada template
% MODULE:
%   isimip
% NAME:
%   isimip_ISO3_list
% PURPOSE:
%   return the isimip standard list of ISO3 country codes and the isimip
%   internal country number.
%
% CALLING SEQUENCE:
%   ISO3_list=isimip_ISO3_list(ISO3_code,isimip_number)
% EXAMPLE:
%   ISO3_list=isimip_ISO3_list; % return whole list
%   isimip_number=isimip_ISO3_list('DEU'); % return isimip number for DEU
%   ISO3_code=isimip_ISO3_list('',52); % return ISO3 code for 52
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   ISO3_code: ISO3 country code(s)
%   isimip_number: isimip internal country number(s)
% OUTPUTS:
%   ISO3_list: the list of ISO3 country codes and the isimip
%       internal country number
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20161117, initial
%-

% define the isimip country ISO3 code list (kept simple, i.e entries on
% same line for easy check/update)
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

if ~exist('ISO3_code','var'),ISO3_code = '';end
if ~exist('isimip_number','var'),isimip_number = [];end

if ~isempty(ISO3_code)
    iso3_pos=strmatch(ISO3_code,ISO3_list(:,1));
    if ~isempty(iso3_pos)
        ISO3_list=ISO3_list{iso3_pos,2}; % get the ID
    else
        ISO3_list=[];
    end % ~isempty(iso3_pos)
    return
end

if ~isempty(isimip_number)
    pos=0;
    for i=1:length(ISO3_list(:,2))
        if ISO3_list{i,2}==isimip_number
            pos=i;
        end
    end
    if pos>0
        ISO3_list=ISO3_list(pos,1);
    else
        ISO3_list=[];
    end
    return
end


end % isimip_ISO3_list