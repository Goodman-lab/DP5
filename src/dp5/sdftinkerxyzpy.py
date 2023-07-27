#################################################################
#
# source activate my-rdkit-env
#
#!/user/bin/python

import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import ChemicalForceFields

###############################################################################
# atom      1     1   CR       "ALKYL CARBON, SP3"         6      12.000    4
# atom      2     2  C=C       "VINYLIC CARBON, SP2"       6      12.000    3
# atom      3     2 CSP2       "GENERIC SP2 CARBON"        6      12.000    3
# atom      4     3  C=O       "GENERAL CARBONYL CARBON"   6      12.000    3
# atom      5     3  C=N       "SP2 CARBON IN C=N"         6      12.000    3
# atom      6     3  CGD       "GUANIDINE C=N"             6      12.000    3
# atom      7     3 C=OR       "KETONE/ALD CARBONYL C"     6      12.000    3
# atom      8     3 C=ON       "AMIDE CARBONYL C"          6      12.000    3
# atom     10     3 CONN       "UREA CARBONYL C"           6      12.000    3
# atom     11     3  COO       "CARBOX AC/ESTER CARB C"    6      12.000    3
# atom     12     3 COON       "CARBAMATE CARBONYL C"      6      12.000    3
# atom     13     3 COOO       "CARBONIC AC/EST CARB C"    6      12.000    3
# atom     14     3 C=OS       "THIOESTER CARB C=O"        6      12.000    3
# atom     15     3  C=S       "THIOESTER C=S"             6      12.000    3
# atom     16     3 C=SN       "THIOAMIDE, C=S"            6      12.000    3
# atom     17     3 CSO2       "CARBON IN >C=SO2"          6      12.000    3
# atom     18     3 CS=O       "C IN >C=S=O"               6      12.000    3
# atom     19     3  CSS       "THIOCARB.AC/EST.CARB.C"    6      12.000    3
# atom     20     3  C=P       "C=P"                       6      12.000    3
# atom     21     4  CSP       "ACETYLENIC CARBON"         6      12.000    2
# atom     22     4  =C=       "ALLENIC CARBON"            6      12.000    2
# atom     23     5   HC       "H ATTACHED TO C"           1       1.008    1
# atom     24     5  HSI       "H ATTACHED TO SI"          1       1.008    1
# atom     25     6   OR       "O ALCOHOL/ETHER"           8      15.995    2
# atom     26     6 OC=O       "ESTER/CARBOX ACID -O-"     8      15.995    2
# atom     27     6 OC=C       "ENOLIC OR PHENOLIC O"      8      15.995    2
# atom     28     6 OC=N       "DIVALENT OXYGEN"           8      15.995    2
# atom     29     6 OC=S       "THIOESTER/THIOACID -O-"    8      15.995    2
# atom     30     6 ONO2       "DIVALENT NITRATE O"        8      15.995    2
# atom     31     6 ON=O       "DIVALENT NITRITE O"        8      15.995    2
# atom     32     6 OSO3       "DIVALENT O-SO3"            8      15.995    2
# atom     33     6 OSO2       "DIVALENT O-SO2"            8      15.995    2
# atom     34     6  OSO       "DIVALENT O-S"              8      15.995    2
# atom     35     6 OS=O       "DIVALENT O-SULFOXIDE S"    8      15.995    2
# atom     36     6  -OS       "GENERAL DIVALENT OX-S"     8      15.995    2
# atom     37     6 OPO3       "DIVALENT O-PO3"            8      15.995    2
# atom     38     6 OPO2       "DIVALENT O-PO2"            8      15.995    2
# atom     39     6  OPO       "DIVALENT O-P"              8      15.995    2
# atom     40     6  -OP       "DIVALENT O-P"              8      15.995    2
# atom     41     6  -O-       "GENERAL DIVALENT O"        8      15.995    2
# atom     42     7  O=C       "GENERAL C=O"               8      15.995    1
# atom     43     7 O=CN       "CARBONYL O, AMIDES"        8      15.995    1
# atom     44     7 O=CR       "CARBONYL O,ALD/KETONES"    8      15.995    1
# atom     45     7 O=CO       "CARBONYL O,CARB.AC/EST"    8      15.995    1
# atom     46     7  O=N       "NITROSO OXYGEN"            8      15.995    1
# atom     47     7  O=S       "O=S IN SULFOXIDES"         8      15.995    1
# atom     48     7 O=S=       "O=S ON S= TO,E.G.,C"       8      15.995    1
# atom     49     8   NR       "N ALIPHATIC AMINES"        7      14.003    3
# atom     50     9  N=C       "N=C IMINES"                7      14.003    2
# atom     51     9  N=N       "N=N AZO COMPOUNDS"         7      14.003    2
# atom     52    10 NC=O       "N AMIDE N-C=O"             7      14.003    3
# atom     53    10 NC=S       "N THIOAMIDE N-C=S"         7      14.003    3
# atom     54    10 NN=C       "NITROGEN IN N-N=C"         7      14.003    3
# atom     55    10 NN=N       "NITROGEN IN N-N=N"         7      14.003    3
# atom     56    11    F       "FLUORINE"                  9      18.998    1
# atom     57    12   CL       "CHLORINE"                 17      35.453    1
# atom     58    13   BR       "BROMINE"                  35      79.904    1
# atom     59    14    I       "IODINE"                   53     126.904    1
# atom     60    15    S       "S THIOETHER/MERCAPTAN"    16      31.972    2
# atom     61    16  S=C       "TERMINAL S=C"             16      31.972    1
# atom     62    17  S=O       "SULFUR IN SULFOXIDES"     16      31.972    3
# atom     63    17 >S=N       "S,TRICOORD,S=N"           16      31.972    3
# atom     64    18  SO2       "S SULFONES -SO2-"         16      31.972    4
# atom     65    18 SO2N       "S IN SULFONAMIDES"        16      31.972    4
# atom     66    18  SO3       "SULFONATE SULFUR"         16      31.972    4
# atom     68    18  SO4       "SULFATE SULFUR"           16      31.972    4
# atom     69    18 =SO2       "SULFONE S=C"              16      31.972    3
# atom     70    18  SNO       "S SULFONE ANALOGS"        16      31.972    4
# atom     71    19   SI       "SILICON"                  14      28.086    4
# atom     72    20 CR4R       "C IN 4-MEMBERED RINGS"     6      12.000    4
# atom     73    21  HOR       "H-O ALCOHOLS"              1       1.008    1
# atom     74    21   HO       "H-O GENERAL OXYGEN"        1       1.008    1
# atom     75    21  HOM       "H HYDROXIDE ANION"         1       1.008    1
# atom     76    22 CR3R       "C 3-MEMBERED RING"         6      12.000    4
# atom     77    23  HNR       "H-N(SP3)"                  1       1.008    1
# atom     78    23  H3N       "H-N(SP3), AMMONIA"         1       1.008    1
# atom     79    23 HPYL       "H-N PYRROLE"               1       1.008    1
# atom     80    23 HNOX       "H-N N-OXIDE"               1       1.008    1
# atom     81    23  HNM       "H ON DICOORD, N(-)"        1       1.008    1
# atom     82    23   HN       "H-N GENERAL NITROGEN"      1       1.008    1
# atom     83    24 HOCO       "H-O CARBOXYLIC ACIDS"      1       1.008    1
# atom     84    24  HOP       "H-O ATTACHED TO P"         1       1.008    1
# atom     85    25  PO4       "P PO4/PHOSPHODIESTER"     15      30.974    4
# atom     86    25  PO3       "TETRACOORDINATE PO3"      15      30.974    4
# atom     87    25  PO2       "TETRACOORDINATE PO2"      15      30.974    4
# atom     88    25   PO       "TETRACOORDINATE PO"       15      30.974    4
# atom     89    25 PTET       "GENERALTETRACOORD P"      15      30.974    4
# atom     90    26    P       "TRICOORD P"               15      30.974    3
# atom     91    27 HN=N       "AZO HYDROGEN"              1       1.008    1
# atom     92    27 HN=C       "IMINE HYDROGEN"            1       1.008    1
# atom     93    28 HNCO       "AMIDE HYDROGEN"            1       1.008    1
# atom     94    28 HNCS       "THIOAMIDE HYDROGEN"        1       1.008    1
# atom     95    28 HNCC       "H-N ENAMINES"              1       1.008    1
# atom     96    28 HNCN       "H-N H-N-C=N"               1       1.008    1
# atom     97    28 HNNC       "H-N H-N-N=C"               1       1.008    1
# atom     98    28 HNNN       "H-N H-N-N=N"               1       1.008    1
# atom     99    28 HNSO       "H-N SULFONAMIDE"           1       1.008    1
# atom    100    28 HNPO       "H-N PHOSPHONAMIDE"         1       1.008    1
# atom    101    28 HNC%       "H-N TRIPLY BONDED C"       1       1.008    1
# atom    102    28 HSP2       "GENERAL H ON SP2 N"        1       1.008    1
# atom    103    29 HOCC       "H-O ENOLS/PHENOLS"         1       1.008    1
# atom    104    29 HOCN       "H-O HO-C=N"                1       1.008    1
# atom    105    30 CE4R       "C=C 4-RING OLEFIN"         6      12.000    3
# atom    106    31  HOH       "H-OH WATER"                1       1.008    1
# atom    107    32 O2CM       "O CARBOX ANION"            8      15.995    1
# atom    108    32  OXN       "N-OXIDE OXYGEN"            8      14.003    1
# atom    109    32  O2N       "NITRO OXYGEN"              8      15.995    1
# atom    110    32 O2NO       "NITRO O IN NITRATE"        8      15.995    1
# atom    111    32  O3N       "NITRATE ANION O"           8      15.995    1
# atom    112    32  O-S       "TERMIN OonTETRAC S"        8      15.995    1
# atom    113    32  O2S       "TERM O-S SULFONE"          8      15.995    1
# atom    114    32  O3S       "TERM O IN SULFONATES"      8      15.995    1
# atom    115    32  O4S       "TERM O IN SO4(-3)"         8      15.995    1
# atom    116    32 OSMS       "TERM O(-)THIOSULFINATE"    8      15.995    1
# atom    117    32   OP       "TERM O IN PHOSPHOXIDES"    8      15.995    1
# atom    118    32  O2P       "TERM O IN PHOSPHINATES"    8      15.995    1
# atom    119    32  O3P       "TERM O IN PHOSPHONATES"    8      15.995    1
# atom    120    32  O4P       "TERM O IN PO4/DIESTER"     8      15.995    1
# atom    121    32 O4CL       "O IN CLO4(-)"              8      15.995    1
# atom    122    33  HOS       "H ON O ATTACHED TO S"      1       1.008    1
# atom    123    34  NR+       "QUATERNARY NSP3(+)"        7      14.003    4
# atom    124    35   OM       "ALKOXIDE O(-)"             8      15.995    1
# atom    125    35  OM2       "OXIDE O ON CSP2(-)"        8      15.995    1
# atom    126    36 HNR+       "H ON QUATERNARY N"         1       1.008    1
# atom    127    36 HIM+       "H ON IMIDAZOLIUM N"        1       1.008    1
# atom    128    36 HPD+       "H ON PYRIDINE N(+)"        1       1.008    1
# atom    129    36 HNN+       "H ON AMIDINIUM N"          1       1.008    1
# atom    130    36 HNC+       "H ON IMINE N(+)"           1       1.008    1
# atom    131    36 HGD+       "H-NH-R GUANIDINIUM"        1       1.008    1
# atom    132    36 HN5+       "H ON N5+,N5A+ORN5B+"       1       1.008    1
# atom    133    37   CB       "C BENZENE/AROMATIC"        6      12.000    3
# atom    134    38 NPYD       "N PYRIDINE"                7      14.003    2
# atom    135    39 NPYL       "N PYRROLE"                 7      14.003    3
# atom    136    40 NC=C       "NITROGEN ON N-C=C"         7      14.003    3
# atom    137    40 NC=N       "NITROGEN IN N-C=N"         7      14.003    3
# atom    138    40 NC=P       "NITROGEN IN N-C=P"         7      14.003    3
# atom    139    40 NC%C       "N-C-C TRIPLE BOND"         7      14.003    3
# atom    140    41 CO2M       "CARBOXYLATE ANION C"       6      12.000    3
# atom    141    41 CS2M       "C THIOCARBOXYLATE"         6      12.000    3
# atom    142    42  NSP       "N  TRIPLE BONDED"          7      14.003    1
# atom    143    43 NSO2       "N IN SULFONAMIDES"         7      14.003    2
# atom    144    43 NSO3       "N SULFONAMIDES SO3"        7      14.003    3
# atom    145    43 NPO2       "N IN PHOSPHONAMIDES"       7      14.003    3
# atom    146    43 NPO3       "N PHOSPHONAMIDES,PO3"      7      14.003    3
# atom    147    43 NC%N       "N-(CYANO GROUP)"           7      14.003    3
# atom    148    44 STHI       "S THIOPHENE"              16      31.972    2
# atom    149    45  NO2       "NITRO GROUP N"             7      14.003    3
# atom    150    45  NO3       "NITRATE GROUP N"           7      14.003    3
# atom    151    46  N=O       "NITROSO NITROGEN"          7      14.003    2
# atom    152    47 NAZT       "TERM N IN AZIDO/DIAZO"     7      14.003    1
# atom    153    48  NSO       "DIVALNforMONOVALOinSO2"    7      14.003    2
# atom    154    49   O+       "OXONIUM(TRICOORD)O"        8      15.995    3
# atom    155    50  HO+       "H ON O+ O"                 1       1.008    1
# atom    156    51  O=+       "OXENIUM (DICOORD) O"       8      15.995    2
# atom    157    52 HO=+       "H ON OXENIUM OXYGEN"       1       1.008    1
# atom    158    53  =N=       "N IN C=N=N OR -N=N=N"      7      14.003    2
# atom    159    54 N+=C       "IMINIUM NITROGEN"          7      14.003    3
# atom    160    54 N+=N       "N(+)=N"                    7      14.003    3
# atom    161    55 NCN+       "N IN +N=C-N"               7      14.003    3
# atom    162    56 NGD+       "N GUANIDINIUM"             7      14.003    3
# atom    163    57 CGD+       "C GUANIDINIUM"             6      12.000    3
# atom    164    57 CNN+       "C IN +N=C-N "              6      12.000    3
# atom    165    58 NPD+       "PYRIDINIUM-TYPE N"         7      14.003    3
# atom    166    59 OFUR       "O FURAN"                   8      15.995    2
# atom    167    60   C%       "ISONITRILE CARBON"         6      12.000    1
# atom    168    61  NR%       "ISONITRILE N/DIAZO N"      7      14.003    2
# atom    169    62   NM       "DEPROT.SULFONAMIDE N-"     7      14.003    2
# atom    170    63  C5A       "A-Cin5MEMB HETEROAROM"     6      12.000    3
# atom    171    64  C5B       "B-Cin5MEMB HETEROAROM"     6      12.000    3
# atom    172    65  N5A       "A-AROMHETEROCYC5RING N"    7      14.003    2
# atom    173    66  N5B       "B-AROMHETEROCYC5RING N"    7      14.003    2
# atom    174    67 N2OX       "NSP2-OXIDE N"              7      14.003    3
# atom    175    68 N3OX       "NSP3-OXIDE N"              7      14.003    4
# atom    176    69 NPOX       "PYRIDINE N-OXIDE N"        7      14.003    3
# atom    177    70  OH2       "OXYGEN ON WATER"           8      15.995    2
# atom    178    71   HS       "H ATACHtoDIVAL,DICOR S"    1       1.008    1
# atom    179    71 HS=N       "H ATACHto4VAL,3COR S=N"    1       1.008    1
# atom    180    71   HP       "H ATACHto3/4CORD P"        1       1.008    1
# atom    181    72  S-P       "TERM S BONDED TO P"       16      31.972    1
# atom    182    72 S2CM       "TER SinTHIOCARBOXYLATE"   16      31.972    1
# atom    183    72   SM       "TERMINAL SULFUR"          16      31.972    1
# atom    184    72 SSMO       "TER S IN THIOSULFINATE"   16      31.972    1
# atom    185    73 SO2M       "S IN SO2(-)"              16      31.972    3
# atom    186    73 SSOM       "3COR SinTHIOSULFINATE"    16      31.972    3
# atom    187    74 =S=O       "SULFINYL S,EG.inC=S=O"    16      31.972    2
# atom    188    75 -P=C       "P=C"                      15      30.974    2
# atom    189    76  N5M       "N(-)in,E.G,3-/4AZOLE"      7      14.003    2
# atom    190    77 CLO4       "CL IN CLO4(-)"            17      35.453    4
# atom    191    78   C5       "GENERAL C 5MEMBHETAR"      6      12.000    3
# atom    192    79   N5       "GENERAL N 5MEMHETERCYC"    7      14.003    2
# atom    193    80 CIM+       "C IMIDAZOLIUM N-C-N"       6      12.000    3
# atom    194    81 NIM+       "IMIDAZOLIUM N"             7      14.003    3
# atom    195    81 N5A+       "POSITIVE N5A N"            7      14.003    3
# atom    196    81 N5B+       "POSITIVE N5B N"            7      14.003    3
# atom    197    81  N5+       "POSITIVE N5 N"             7      14.003    3
# atom    198    82 N5AX       "N-OXIDEN 5-RING A-POS"     7      14.003    3
# atom    199    82 N5BX       "N-OXIDEN 5-RING B-POS"     7      14.003    3
# atom    200    82 N5OX       "N-OXIDN GENER5RINGPOS"     7      14.003    3
# atom    201    87 FE+2       "IRON +2 CATION"           26      55.845    0
# atom    202    88 FE+3       "IRON +3 CATION"           26      55.845    0
# atom    203    89   F-       "FLUORIDE ANION"            9      18.998    0
# atom    204    90  CL-       "CHLORIDE ANION"           17      35.453    0
# atom    205    91  BR-       "BROMIDE ANION"            35      79.904    0
# atom    206    92  LI+       "LITHIUM CATION"            3       6.941    0
# atom    207    93  NA+       "SODIUM CATION"            11      22.990    0
# atom    208    94   K+       "POTASSIUM CATION"         19      39.098    0
# atom    209    95 ZINC       "DIPOSITIVE ZINC"          30      65.390    0
# atom    210    95 ZN+2       "DIPOSITIVE ZINC"          30      65.390    0
# atom    211    96 CA+2       "CALCIUM(+2) CATION"       20      40.078    0
# atom    212    97 CU+1       "COPPER(+1) CATION"        29      63.546    0
# atom    213    98 CU+2       "COPPER(+2) CATION"        29      63.546    0
# atom    214    99 MG+2       "MAGNESIUM(+2) CATION"     12      24.305    0
###############################################################################

def getMMFF_large_atom_type(mmff_props , atom, m):
  # Small to large numbers
  small_to_large_list = [[[]],
      [[1]],
      [[3, "C"], [2,"C=C"]],
      [[4,"C=O"], [5,"C=N"], [6,"NC(N)=N"], [7,"CC=O"], [8,"NC=O"], [10,"NC(=O)N"], [11,"OC=O"], [12,"NC(=O)O"], [13,"NC(=O)O"], [14,"OC(=O)O"], [15,"SC=O"], [16,"NC=S"], [17,"C=S(O)O"], [18,"C=S=O"], [19,"SC=S"], [20,"C=P"]],
      [[21,"C#[C,N]"], [22,"[C,N,O]=C=[C,N,O]"]],
      [[23,"C[H]"], [24,"[Si][H]"]],
      [[41,"O"], [25,"OC"], [26,"OC=O"], [27,"OC=C"], [27,"Occ"], [28,"OC=N"], [29,"OC=S"], [31,"ON=O"], [30,"O[N+]([O-])=O"], [36,"OS"], [34,"OSO"], [35,"OS=O"], [33,"OS(O)=O"], [32,"OS(O)(=O)=O"], [40,"OP"], [39,"OPO"], [38,"OP(=O)O"], [37,"OP(=O)(=O)O"]],
      [[42,"C=O"], [44,"CC=O"], [43,"NC=O"], [45,"OC=O"], [46,"O=N"], [47,"S=O"], [48,"[C,N]=S=O"]],
      [[49]],
      [[50,"C=N"], [51,"N=N"]],
      [[52,"NC=O"], [53,"NC=S"], [54,"NN=C"], [55,"NN=N"]],
      [[56]],
      [[57]],
      [[58]],
      [[59]],
      [[60]],
      [[61]],
      [[62,"S=O"], [63,"S=N"]],
      [[64,"O=S=O"], [70,"OSN"], [65,"N-S(=O)=O"], [66,"OS(O)O"], [67,"C"], [68,"OS(O)(O)O"], [69,"CS(O)(O)C"]],
      [[71]],
      [[72]],
      [[74,"[H]O"], [73,"[H]OC"], [75,"[H][O-]"]],
      [[76]],
      [[82,"[H]N"], [77,"[H]N(C)C"], [78,"[H]N([H])[H]"], [79,"[H]n1cccc1"], [80,"[H]NO"], [81,"[H][N-]"]],
      [[83,"[H]OC=O"], [84,"[H]OP"]],
      [[89,"P"], [88,"PO"], [87,"OPO"], [86,"OP(O)O"], [89,"OP(O)(O)O"]],
      [[90]],
      [[91,"[H]N=N"], [92,"[H]N=C"]],
      [[102,"[H]N"], [93,"[H]NC=O"], [94,"[H]NC=S"], [95,"[H]NC=C"], [96,"[H]NC=N"], [97,"[H]NN=C"], [98,"[H]NN=N"], [99,"[H]NS=O"], [100,"[H]NP=O"], [101,"HN#[C,N]"]],
      [[103,"[H]OC=C"], [104,"[C]OC=N"]],
      [[105]],
      [[106]],
      [[107,"[O-]C=O"], [108,"NO"], [109,"ON=O"], [110,"O[N+]([O-])=O"], [111,"[O-][N+]([O-])=O"], [112,"OS"], [113,"OS=O"], [114,"OS(=O)=O"], [115,"OS(=O)(=O)O"], [116,"O=[S-]S"], [117,"OP"], [118,"OPO"], [119,"OP(=O)O"], [120,"OP(=O)(=O))"], [121,"OCl(=O)(=O)[O-]"]],
      [[122]],
      [[123]],
      [[124,"[O-]"], [125,"[O-]C=[C,N]"]],
      [[126,"[H][N+][H,C][H,C][H,C]"], [127,"C1=[NH+]C=CN1"], [128,"C1=C[NH+]=CC=C1"], [129,"CC(N)=[NH2+]"], [130,"C=[NH2+]"], [131,"NC(N)=[NH2+]"], [132,"[H]N([H])([H])([H])[H]"]],
      [[133]],
      [[134]],
      [[135]],
      [[136,"NC=C"], [137,"NC=N"], [138,"NC=P"], [139,"NC#C"]],
      [[140,"[O-]C=O"], [141,"[S-]C=S"]],
      [[142]],
      [[143,"NS(=O)O"], [144,"NS(=O)(=O)O"], [145,"NP(=O)O"], [146,"NP(=O)(=O)O"], [147,"NC#N"]],
      [[148]],
      [[149,"ON=O"], [150,"O[N+][O-]=O"]],
      [[151]],
      [[152]],
      [[153]],
      [[154]],
      [[155]],
      [[156]],
      [[157]],
      [[158]],
      [[159,"[N+]=C"], [160,"[N+=N]"]],
      [[161]],
      [[162]],
      [[163,"NC(N)=[NH2+]"], [164,"[N+]=CN"]],
      [[165]],
      [[166]],
      [[167]],
      [[168]],
      [[169]],
      [[170]],
      [[171]],
      [[172]],
      [[173]],
      [[174]],
      [[175]],
      [[176]],
      [[177]],
      [[178,"[H]S"], [179,"[H]S=N"], [180,"[H]P"]],
      [[181,"SP"], [183,"[S-]"], [182,"[S-]C=S"], [184,"[S-]S(=O)"]],
      [[185,"[O-]S=O"], [186,"[O-]S=S"]],
      [[187]],
      [[188]],
      [[189]],
      [[190]],
      [[191]],
      [[192]],
      [[193]],
      [[194,"N[N+]1=CNC=C1"], [195,"[H][N+]([H])([H])([H])[H]"], [196,"[H][N+]([H])([H])([H])[H]"], [197,"[H][N+]([H])([H])([H])[H]"]],
      [[198,"[H][N+]([H])([H])([H])[H]"],[199,"[H][N+]([H])([H])([H])[H]"],[200,"[H][N+]([H])([H])([H])[H]"]],
      [[-1]],
      [[-1]],
      [[-1]],
      [[-1]],
      [[201]],
      [[202]],
      [[203]],
      [[204]],
      [[205]],
      [[206]],
      [[207]],
      [[208]],
      [[209,"[Zn]"], [210,"Zn++"]],
      [[211]],
      [[212]],
      [[213]],
      [[214]]]

  MMFF_small_atom_type = mmff_props.GetMMFFAtomType(atom.GetIdx())
  MMFF_large_atom_type = small_to_large_list[MMFF_small_atom_type][0][0]
  if len(small_to_large_list[MMFF_small_atom_type]) > 1:
    for atom_info in small_to_large_list[MMFF_small_atom_type]:
      substructure = Chem.MolFromSmarts(atom_info[1])
      for substructure_match in m.GetSubstructMatches(substructure):
        if substructure_match.count(atom.GetIdx()) > 0:
          MMFF_large_atom_type = atom_info[0]
  return MMFF_large_atom_type






###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#

def main(sdf_file):

    # Tinker xyz format
    # number of atoms, title
    # number, element symbol, x, y, z, MMFF atom type, connectivity
    m = Chem.MolFromMolFile(sdf_file + ".sdf", removeHs=False)
    #print(sys.argv)
    # smi_string=Chem.MolToSmiles(m)
    # print(smi_string)

    conf = m.GetConformer()
    mmff_props = AllChem.MMFFGetMoleculeProperties(m)

    xyz = open(sdf_file + ".xyz", "w+")

    xyz.write("{:>6}  {}\n".format(m.GetNumAtoms(),sdf_file))

    for atom in m.GetAtoms():

        bond_list = []
        attached_atoms = ""
        bond_type = []
        connection_type = []

        for connection in atom.GetNeighbors():
            bond_list.append(connection.GetIdx()+1)
        bond_list.sort()

        for connection in bond_list:
            attached_atoms += "{:>5}".format(str(connection))+" "



        xyz.write("{:>6} {:>2} {:13.6f} {:11.6f} {:11.6f} {:>5} {}\n".format(
        atom.GetIdx()+1,
        atom.GetSymbol(),
        list(conf.GetAtomPosition(atom.GetIdx()))[0],
        list(conf.GetAtomPosition(atom.GetIdx()))[1],
        list(conf.GetAtomPosition(atom.GetIdx()))[2],
        getMMFF_large_atom_type(mmff_props , atom, m ),
        attached_atoms) )

    xyz.close()


