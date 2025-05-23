const NUMBER_OF_SHELLS = [
  1,  1,  2,  2,  3,    3,  4,  4,  3,  4, ##   1 - 10
  5,  5,  6,  6,  6,    6,  6,  7,  8,  8, ##  11 - 20
  9,  9,  9,  9,  9,    9,  9, 10, 10, 10, ##  21 - 30
 11, 11, 11, 11, 11,   12, 13, 13, 14, 14, ##  31 - 40
 14, 14, 14, 14, 14,   15, 15, 15, 16, 16, ##  41 - 50
 16, 16, 16, 17, 18,   18, 19, 19, 19, 19, ##  51 - 60
 19, 19, 19, 20, 19,   19, 19, 19, 19, 20, ##  61 - 70
 21, 21, 21, 21, 21,   21, 21, 21, 22, 22, ##  71 - 80
 23, 23, 23, 23, 24,   24, 25, 25, 26, 26, ##  81 - 90
 27, 27, 27, 26, 26,   27, 27, 26, 26, 26, ##  91 - 100
 27, 27, 28, 28                            ## 101 - 104
]

const INDEX_OF_SHELLS = [1; cumsum(NUMBER_OF_SHELLS) .+ 1]

const ATOMIC_SHELLS = [
  ##  H  ---------------------------------------------------------1
  13.60,

  ##  He ---------------------------------------------------------2
  24.59,

  ##  Li  --------------------------------------------------------3
  58.0,  5.39,

  ##  Be  --------------------------------------------------------5
  115.0, 9.32,

  ##  B   --------------------------------------------------------7
  192.0, 12.93, 8.3,

  ##  C  ---------------------------------------------------------10
  288.0, 16.59, 11.26,

  ##  N  ---------------------------------------------------------
  403.0, 37.3, 20.33, 14.53,

  ##  O  ---------------------------------------------------------
  543.1, 41.6, 28.48, 13.62,

  ##  F  ---------------------------------------------------------
  696.7, 37.85, 17.42,

  ##  Ne ---------------------------------------------------------
  870.1, 48.47, 21.66, 21.56,

  ##  Na ---------------------------------------------------------
  1075.0, 66.0, 34.0, 34.0, 5.14,

  ##  Mg ---------------------------------------------------------
  1308.0, 92.0, 54.0, 54.0, 7.65,

  ##  Al ---------------------------------------------------------
  1564.0, 121., 77.0, 77.0, 10.62, 5.99,

  ##  Si ---------------------------------------------------------
  1844.0, 154.0, 104.0, 104.0, 13.46, 8.15,

  ##  P  ---------------------------------------------------------
  2148.0, 191.0, 135.0, 134.0, 16.15, 10.49,

  ##  S  ---------------------------------------------------------
  2476.0, 232.0, 170.0, 168.0, 20.20, 10.36,

  ##  Cl ---------------------------------------------------------
  2829.0, 277.0, 208.0, 206.0, 24.54, 12.97,

  ##  Ar ---------------------------------------------------------
  3206.3, 326.5, 250.6, 248.5, 29.24, 15.94, 15.76,

  ##  K  ---------------------------------------------------------
  3610.0, 381.0, 299.0, 296.0, 37.0, 19.0, 18.7, 4.34,

  ##  Ca ---------------------------------------------------------
  4041.0, 441.0, 353.0, 349.0, 46.0, 28.0, 28.0, 6.11,

  ##  Sc ---------------------------------------------------------
  4494.0, 503.0, 408.0, 403.0, 55.0, 33.0, 33.0, 8.0, 6.54,

  ##  Ti ---------------------------------------------------------
  4966.0, 567.0, 465.0, 459.0, 64.0, 39.0, 38.0, 8.0, 6.82,

  ##  V  ---------------------------------------------------------
  5465.0, 633.0, 525.0, 518.0, 72.0, 44.0, 43.0, 8.0, 6.74,

  ##  Cr ---------------------------------------------------------
  5989.0, 702.0, 589.0, 580.0, 80.0, 49.0, 48.0, 8.25, 6.77,

  ##  Mn ---------------------------------------------------------
  6539.0, 755.0, 656.0, 645.0, 89.0, 55.0, 53.0, 9.0, 7.43,

  ##  Fe ---------------------------------------------------------
  7112.0, 851.0, 726.0, 713.0, 98.0, 61.0, 59.0, 9.0, 7.87,

  ##  Co ---------------------------------------------------------
  7709.0, 931.0, 800.0, 785.0, 107.0, 68.0, 66.0, 9.0, 7.86,

  ##  Ni ---------------------------------------------------------
  8333.0, 1015.0, 877.0, 860.0, 117.0, 75.0, 73.0, 10.0, 10.0,
     7.64,

  ##  Cu ---------------------------------------------------------
  8979.0, 1103.0, 958.0, 938.0, 127.0, 82.0, 80.0, 11.0, 10.4,
     7.73,

  ##  Zn ---------------------------------------------------------
  9659.0, 1198.0, 1047.0, 1024.0, 141.0, 94.0, 91.0, 12.0, 11.2,
     9.39,

  ##  Ga ---------------------------------------------------------
  10367.0, 1302.0, 1146.0, 1119.0, 162.0, 111.0, 107.0, 21.0,
     20.0,   11.0,    6.0,

  ##  Ge ---------------------------------------------------------
  11103.0, 1413.0, 1251.0, 1220.0, 184.0, 130.0, 125.0, 33.0,
     32.0,   14.3,    7.9,

  ##  As ---------------------------------------------------------
  11867.0, 1531.0, 1362.0, 1327.0, 208.0, 151.0, 145.0, 46.0,
     45.0,   17.0,    9.81,

  ##  Se ---------------------------------------------------------
  12658.0, 1656.0, 1479.0, 1439.0, 234.0, 173.0, 166.0, 61.0,
     60.0,   20.15,   9.75,

  ##  Br ---------------------------------------------------------
  13474.0, 1787.0, 1602.0, 1556.0, 262.0, 197.0, 189.0, 77.0,
     76.0,   23.8,   11.85,

  ##  Kr ---------------------------------------------------------
  14326.0, 1924.6, 1730.9, 1678.4, 292.8, 222.2, 214.4, 95.0,
     93.8,   27.51,  14.65,  14.0,

  ##  Rb ---------------------------------------------------------
  15200.0, 2068.0, 1867.0, 1807.0, 325.0, 251.0, 242.0, 116.0,
    114.0,   32.0,   16.0,   15.3,   4.18,

  ##  Sr ---------------------------------------------------------
  16105.0, 2219.0, 2010.0, 1943.0, 361.0, 283.0, 273.0, 139.0,
    137.0,   40.0,   23.0,   22.0,   5.69,

  ##  Y  ---------------------------------------------------------
  17038.0, 2375.0, 2158.0, 2083.0, 397.0, 315.0, 304.0, 163.0,
    161.0,   48.0,   30.0,   29.0,   6.48,   6.38,  

  ##  Zr ---------------------------------------------------------
  17998.0, 2536.0, 2311.0, 2227.0, 434.0, 348.0, 335.0, 187.0,
    185.0,   56.0,   35.0,   33.0,   8.61,  6.84,

  ##  Nb ----------------------------------------------------------
  18986.0, 2702.0, 2469.0, 2375.0, 472.0, 382.0, 367.0, 212.0,
    209.0,   62.0,   40.0,   38.0,   7.17,  6.88,

  ##  Mo -----------------------------------------------------------
  20000.0, 2872.0, 2632.0, 2527.0, 511.0, 416.0, 399.0, 237.0,
    234.0,   68.0,   45.0,   42.0,   8.56,  7.10,

  ##  Tc -----------------------------------------------------------
  21044.0, 3048.0, 2800.0, 2683.0, 551.0, 451.0, 432.0, 263.0,
    259.0,   74.0,   49.0,   45.0,   8.6,   7.28,

  ##  Ru -----------------------------------------------------------
  22117.0, 3230.0, 2973.0, 2844.0, 592.0, 488.0, 466.0, 290.0,
    286.0,   81.0,   53.0,   49.0,   8.5,   7.37,

  ##  Rh -----------------------------------------------------------
  23220.0, 3418.0, 3152.0, 3010.0, 634.0, 526.0, 501.0, 318.0,
    313.0,   87.0,   58.0,   53.0,   9.56,  7.46,

  ##  Pd -----------------------------------------------------------
  24350.0, 3611.0, 3337.0, 3180.0, 677.0, 565.0, 537.0, 347.0,
    342.0,   93.0,   63.0,   57.0,   8.78,  8.34,  7.52,

  ##  Ag -----------------------------------------------------------
  25514.0, 3812.0, 3530.0, 3357.0, 724.0, 608.0, 577.0, 379.0,
    373.0,  101.0,   69.0,   63.0,  11.0,  10.0,   7.58,

  ##  Cd ------------------------------------------------------------
  26711.0, 4022.0, 3732.0, 3542.0, 775.0, 655.0, 621.0, 415.0,
    408.0,  112.0,   78.0,   71.0,  14.0,  13.0,   8.99,

  ##  In ------------------------------------------------------------
  27940.0, 4242.0, 3943.0, 3735.0, 830.0, 707.0, 669.0, 455.0,
    447.0,  126.0,   90.0,   82.0,  21.0,  20.0,  10.0,   5.79,

  ##  Sn ------------------------------------------------------------
  29200.0, 4469.0, 4160.0, 3933.0, 888.0, 761.0, 719.0, 497.0,
    489.0,  141.0,  102.0,   93.0,  29.0,  28.0,  12.0,   7.34,

  ##  Sb ------------------------------------------------------------
  30419.0, 4698.0, 4385.0, 4137.0, 949.0, 817.0, 771.0, 542.0,
    533.0,  157.0,  114.0,  104.0,  38.0,  37.0,  15.0,   8.64,

  ##  Te ------------------------------------------------------------
  31814.0, 4939.0, 4612.0, 4347.0, 1012.0, 876.0, 825.0, 589.0,
    578.0,  174.0,  127.0,  117.0,   48.0,  46.0,  17.84,  9.01,

  ##  I  ------------------------------------------------------------
  33169.0, 5188.0, 4852.0, 4557.0, 1078.0, 937.0, 881.0, 638.0,
    626.0,  193.0,  141.0,  131.0,   58.0,  56.0,  20.61, 10.45,

  ##  Xe ------------------------------------------------------------
  34570.0, 5460.0, 5110.0, 4790.0, 1148.7, 1002.1, 940.6, 689.0,
    676.4,  213.2,  157.0,  145.5,   69.5,   67.5,  23.39, 13.43,
     12.13,

  ##  Cs ------------------------------------------------------------
  35985.0, 5714.0, 5359.0, 5012.0, 1220.0, 1068.0, 1000.0, 742.0,
    728.0,  233.0,  174.0,  164.0,   81.0,   79.0,   25.0,  14.0,
     12.3,    3.89,

  ##  Ba -------------------------------------------------------------
  37441.0, 5989.0, 5624.0, 5247.0, 1293.0, 1138.0, 1063.0, 797.0,
    782.0,  254.0,  193.0,  181.0,   94.0,   92.0,   31.0,  18.0,
     16.0,    5.21,

  ##  La -------------------------------------------------------------
  38925.0, 6266.0, 5891.0, 5483.0, 1365.0, 1207.0, 1124.0, 851.0,
    834.0,  273.0,  210.0,  196.0,  105.0,  103.0,   36.0,  22.0,
     19.0,    5.75,   5.58,

  ##  Ce -------------------------------------------------------------
  40443.0, 6548.0, 6164.0, 5723.0, 1437.0, 1275.0, 1184.0, 903.0,
    885.0,  291.0,  225.0,  209.0,  114.0,  111.0,   39.0,  25.0, 
     22.0,    6.0,    5.65,

  ##  Pr -------------------------------------------------------------
  41991.0, 6835.0, 6440.0, 5964.0, 1509.0, 1342.0, 1244.0, 954.0,
    934.0,  307.0,  238.0,  220.0,  121.0,  117.0,   41.0,
     27.0,   24.0,    6.0,  5.42 ,

  ##  Nd -------------------------------------------------------------
  43569.0, 7126.0, 6722.0, 6208.0, 1580.0, 1408.0, 1303.0, 1005.0,
    983.0,  321.0,  250.0,  230.0,  126.0,  122.0,   42.0,
     28.0,   25.0,    6.0,  5.49,

  ##  Pm -------------------------------------------------------------
  45184.0, 7428.0, 7013.0, 6459.0, 1653.0, 1476.0, 1362.0, 1057.0,
   1032.0,  325.0,  261.0,  240.0,  131.0,  127.0,   43.0,
     28.0,   25.0,    6.0,   5.55,

  ##  Sm -------------------------------------------------------------
  46834.0, 7737.0, 7312.0, 6716.0, 1728.0, 1546.0, 1422.0, 1110.0,
   1083.0,  349.0,  273.0,  251.0,  137.0,  132.0,   44.0,
     29.0,   25.0,    6.0,    5.63,

  ##  Eu -------------------------------------------------------------
  48519.0, 8052.0, 7617.0, 6977.0, 1805.0, 1618.0, 1484.0, 1164.0,
   1135.0,  364.0,  286.0,  262.0,  143.0,  137.0,   45.0,
     30.0,   26.0,    6.0,    5.68,

  ##  Gd -------------------------------------------------------------
  50239.0, 8376.0, 7930.0, 7243.0, 1884.0, 1692.0, 1547.0, 1220.0,
   1189.0,  380.0,  300.0,  273.0,  150.0,  143.0,   46.0,
     31.0,   27.0,    6.16,   6.0,    6.0,

  ##  Tb -------------------------------------------------------------
  51996.0, 8708.0, 8252.0, 7514.0, 1965.0, 1768.0, 1612.0, 1277.0,
   1243.0,  398.0,  315.0,  285.0,  157.0,  150.0,   48.0,
     32.0,   28.0,    6.0,    5.85,

  ##  Dy -------------------------------------------------------------
  53789.0, 9046.0, 8581.0, 7790.0, 2048.0, 1846.0, 1678.0, 1335.0,
   1298.0,  416.0,  331.0,  297.0,  164.0,  157.0,   50.0,
     33.0,   28.0,    6.0,    5.93,

  ##  Ho -------------------------------------------------------------
  55618.0, 9394.0, 8918.0, 8071.0, 2133.0, 1926.0, 1746.0, 1395.0,
   1354.0,  434.0,  348.0,  310.0,  172.0,  164.0,   52.0,
     34.0,   29.0,    6.02,   6.0,

  ##  Er -------------------------------------------------------------
  57486.0, 9751.0, 9264.0, 8358.0, 2220.0, 2008.0, 1815.0, 1456.0,
   1412.0,  452.0,  365.0,  323.0,  181.0,  172.0,   54.0,
     35.0,   30.0,    6.10,   6.0, 

  ##  Tu -------------------------------------------------------------
  59390.0, 10116.0, 9617.0, 8648.0, 2309.0, 2092.0, 1885.0, 1518.0,
   1471.0,   471.0,  382.0,  336.0,  190.0,  181.0,   56.0,
     36.0,    30.0,    7.0,    6.18,

  ##  Yb -------------------------------------------------------------
  61332.0, 10486.0, 9978.0, 8944.0, 2401.0, 2178.0, 1956.0, 1580.0,
   1531.0,   490.0,  399.0,  349.0,  200.0,  190.0,  
     58.0,    37.0,   31.0,    8.0,    7.0,    6.25,

  ##  Lu -------------------------------------------------------------
  63314.0, 10870.0, 10349.0, 9244.0, 2499.0, 2270.0, 2032.0, 1647.0,
   1596.0,   514.0,   420.0,  366.0,  213.0,  202.0,
     62.0,    39.0,    32.0,   13.0,   12.0,    7.0,    6.6,

  ##  Hf -------------------------------------------------------------
  65351.0, 11271.0, 10739.0, 9561.0, 2604.0, 2369.0, 2113.0, 1720.0,
   1665.0,   542.0,   444.0,  386.0,  229.0,  217.0,
     68.0,    43.0,    35.0,   21.0,   20.0,    7.5,    7.0,

  ##  Ta -------------------------------------------------------------
  67416.0, 11682.0, 11136.0, 9881.0, 2712.0, 2472.0, 2197.0, 1796.0,
   1737.0,   570.0,   469.0,  407.0,  245.0,  232.0,
     74.0,    47.0,    38.0,   30.0,   28.0,    8.3,    7.9,

  ##  W  -------------------------------------------------------------
  69525.0, 12100.0, 11544.0, 10207.0, 2823.0, 2577.0, 2283.0, 1874.0,
   1811.0,   599.0,   495.0,   428.0,  261.0,  248.0,
     80.0,    51.0,    41.0,    38.0,   36.0,    9.0,    8.0,

  ##  Re -------------------------------------------------------------
  71676.0, 12527.0, 11959.0, 10535.0, 2937.0, 2686.0, 2371.0, 1953.0,
   1887.0,   629.0,   522.0,   450.0,  278.0,  264.0,
     86.0,    56.0,    47.0,   45.0,    45.0,     9.6,    7.9,

  ##  Os -------------------------------------------------------------
  73871.0, 12968.0, 12385.0, 10871.0, 3054.0, 2797.0, 2461.0, 2035.0,
   1964.0,   660.0,   551.0,   473.0,  295.0,  280.0,
     92.0,    61.0,    56.0,    54.0,   49.0,    9.6,    8.5,

  ##  Ir -------------------------------------------------------------
  76111.0, 13419.0, 12824.0, 11215.0, 3175.0, 2912.0, 2554.0, 2119.0,
   2044.0,   693.0,   581.0,   497.0,  314.0,  298.0,
     99.0,    67.0,    66.0,    64.0,   53.0,    9.6,    9.1,

  ##  Pt -------------------------------------------------------------
  78395.0, 13880.0, 13273.0, 11564.0, 3300.0, 3030.0, 2649.0, 2206.0,
   2126.0,   727.0,   612.0,   522.0,  335.0,  318.0,
    106.0,    78.0,    75.0,    71.0,   57.0,     9.6,    9.0,

  ##  Au -------------------------------------------------------------
  80725.0, 14353.0, 13734.0, 11919.0, 3430.0, 3153.0, 2748.0, 2295.0,
   2210.0,   764.0,   645.0,   548.0,  357.0,  339.0,
    114.0,    91.0,    87.0,    76.0,   61.0,    12.5,  11.1,    9.23,

  ##  Hg -------------------------------------------------------------
  83102.0, 14839.0, 14209.0, 12284.0, 3567.0, 3283.0, 2852.0, 2390.0,
   2300.0,   806.0,   683.0,   579.0,  382.0,  363.0,
    125.0,   107.0,   103.0,    85.0,   68.0,   14.0,   12.0,   10.4,

  ##  Tl -------------------------------------------------------------
  85530.0, 15347.0, 14698.0, 12658.0, 3710.0, 3420.0, 2961.0, 2490.0,
   2394.0,   852.0,   726.0,   615.0,  411.0,  391.0,
    139.0,   127.0,   123.0,    98.0,   79.0,   21.0,   19.0,    8.0,
      6.11,

  ##  Pb -------------------------------------------------------------
  88005.0, 15861.0, 15200.0, 13055.0, 3857.0, 3560.0, 3072.0, 2592.0,
   2490.0,   899.0,   769.0,   651.0,  441.0,  419.0,
    153.0,   148.0,   144.0,   111.0,   90.0,   27.0,   25.0,   10.0,
      7.42,

  ##  Bi -------------------------------------------------------------
  90526.0, 16388.0, 15711.0, 13419.0, 4007.0, 3704.0, 3185.0, 2696.0,
   2588.0,   946.0,   813.0,   687.0,  472.0,  448.0,  170.0,
    167.0,   165.0,   125.0,   101.0,   34.0,   32.0,   12.0,    7.29,

  ##  Po -------------------------------------------------------------
  93105.0, 16939.0, 16244.0, 13814.0, 4161.0, 3852.0, 3301.0, 2802.0,
   2687.0,   994.0,   858.0,   724.0,  503.0,  478.0,  193.0,  187.0,
    181.0,   139.0,   112.0,    41.0,   38.0,   15.0,    8.43,

  ##  At -------------------------------------------------------------
  95730.0, 17493.0, 16785.0, 14214.0, 4320.0, 4005.0, 3420.0, 2910.0,
   2788.0,  1044.0,   904.0,   761.0,  535.0,  508.0,  217.0,  211.0,
    196.0,   153.0,   123.0,    48.0,   44.0,   19.0,   11.0,    9.3,

  ##  Rn -------------------------------------------------------------
  98404.0, 18049.0, 17337.0, 14619.0, 4483.0, 4162.0, 3452.0, 3109.0,
   2890.0,  1096.0,   951.0,   798.0,  567.0,  538.0,  242.0,  235.0,
    212.0,   167.0,   134.0,    55.0,   51.0,   24.0,   14.0,   10.7,

  ##  Fr -------------------------------------------------------------
 101137.0, 18639.0, 17907.0, 15031.0, 4652.0, 4324.0, 3666.0, 3134.0,
   2998.0,  1153.0,  1003.0,   839.0,  603.0,  572.0,  268.0,  260.0,
    231.0,   183.0,   147.0,    65.0,   61.0,   33.0,   19.0,   14.0,
      4.0,

  ##  Ra -------------------------------------------------------------
 103922.0, 19237.0, 18484.0, 15444.0, 4822.0, 4491.0, 3793.0, 3254.0,
   3111.0,  1214.0,  1060.0,   884.0,  642.0,  609.0,  296.0,  287.0,
    253.0,   201.0,   161.0,    77.0,   73.0,   40.0,   25.0,   19.0,
      5.28,

  ##  Ac -------------------------------------------------------------
 106755.0, 19840.0, 19083.0, 15871.0, 5002.0, 4656.0, 3921.0, 3374.0,
   3223.0,  1274.0,  1116.0,   928.0,  680.0,  645.0,  322.0,  313.0,
    274.0,   218.0,   174.0,    88.0,   83.0,   45.0,   29.0,   22.0,
      6.3,     5.7,

  ##  Th -------------------------------------------------------------
 109651.0, 20472.0, 19693.0, 16300.0, 5182.0, 4830.0, 4049.0, 3494.0,
   3335.0,  1333.0,  1171.0,   970.0,  717.0,  679.0,  347.0,  338.0,
    293.0,   233.0,   185.0,    97.0,   91.0,   50.0,   33.0,   25.0,
      6.0,     6.0,

  ##  Pa -------------------------------------------------------------
 112601.0, 21105.0, 20314.0, 16733.0, 5367.0, 5001.0, 4178.0, 3613.0,
   3446.0,  1390.0,  1225.0,  1011.0,  752.0,  712.0,  372.0,  362.0,
    312.0,   248.0,   195.0,   104.0,   97.0,   50.0,   32.0,
     24.0,     6.0,     6.0,     6.0,

  ##  U  -------------------------------------------------------------
 115606.0, 21757.0, 20948.0, 17166.0, 5548.0, 5182.0, 4308.0, 3733.0,
   3557.0,  1446.0,  1278.0,  1050.0,  785.0,  743.0,  396.0,  386.0,
    329.0,   261.0,   203.0,   110.0,  101.0,   52.0,   34.0,
     24.0,     6.1,     6.0,     6.0,

  ##  Np -------------------------------------------------------------
 118678.0, 22426.0, 21600.0, 17610.0, 5723.0, 5366.0, 4440.0, 3854.0,
   3669.0,  1504.0,  1331.0,  1089.0,  819.0,  774.0,  421.0,  410.0,
    346.0,   274.0,   211.0,   116.0,  106.0,   54.0,   35.0,
     25.0,     6.0,     6.0,     6.0,

  ##  Pu -------------------------------------------------------------
 121818.0, 23097.0, 22266.0, 18056.0, 5933.0, 5541.0, 4557.0, 3977.0,
   3783.0,  1563.0,  1384.0,  1128.0,  853.0,  805.0,  446.0,  434.0,
    356.0,   287.0,   219.0,   122.0,  111.0,   53.0,   34.0,
     23.0,     6.0,     6.0,

  ##  Am -------------------------------------------------------------
 125027.0, 23773.0, 22944.0, 18504.0, 6121.0, 5710.0, 4667.0, 4102.0,
   3898.0,  1623.0,  1439.0,  1167.0,  887.0,  836.0,  467.0,  452.0,
    355.0,   301.0,   220.0,   123.0,  112.0,   54.0,   44.0,
     36.0,     6.0,     6.0,

  ##  Cm -------------------------------------------------------------
 128220.0, 24460.0, 23779.0, 18930.0, 6288.0, 5895.0, 4797.0, 4236.0,
   4014.0,  1664.0,  1493.0,  1194.0,  919.0,  864.0,  494.0,  479.0,
    384.0,   314.0,   239.0,   126.0,  119.0,   60.0,   39.0,
     27.0,    11.0,     5.0,     6.0,

  ##  Bk -------------------------------------------------------------
 131590.0, 25275.0, 24385.0, 19452.0, 6556.0, 6147.0, 4977.0, 4366.0,
   4133.0,  1729.0,  1554.0,  1236.0,  955.0,  898.0,  520.0,  504.0,
    401.0,   329.0,   248.0,   142.0,  124.0,   63.0,   41.0,
     27.0,    12.0,     6.0,     4.0,

  ##  Cf -------------------------------------------------------------
 135960.0, 26110.0, 25250.0, 19930.0, 6754.0, 6359.0, 5109.0, 4492.0,
   4247.0,  1789.0,  1610.0,  1273.0,  987.0,  925.0,  546.0,  529.0,
    412.0,   338.0,   251.0,   142.0,  129.0,   61.0,   39.0,
     25.0,     9.0,     6.0,

  ##  Es -------------------------------------------------------------1488
 139490.0, 26900.0, 26020.0, 20410.0, 6977.0, 6754.0, 5252.0, 4630.0,
   4369.0,  1857.0,  1674.0,  1316.0, 1024.0,  959.0,  573.0,  554.0,
    429.0,   353.0,   260.0,   148.0,  135.0,   63.0,   40.0,
     25.0,     9.0,     6.0,

  ##  Fm -------------------------------------------------------------1514
 143090.0, 27700.0, 26810.0, 20900.0, 7205.0, 6793.0, 5397.0, 4766.0,
   4498.0,  1933.0,  1746.0,  1366.0, 1068.0, 1000.0,  606.0,  587.0,
    453.0,   375.0,   275.0,   160.0,  145.0,   69.0,   45.0,
     29.0,    15.0,     7.0,

  ## Md --------------------------------------------------------------1540
 146526.0, 28387.0, 27438.0, 21356.0, 7440.0, 7001.0, 5552.0, 4889.0,
   4615.0,  2024.0,  1816.0,  1424.0, 1105.0, 1034.0,  618.0,  597.0,
    471.0,   389.0,   272.0,   154.0,  137.0,   12.9,   10.5,   61.0,
     37.0,    17.0,     5.9,

  ## No --------------------------------------------------------------1567
 149208.0, 29221.0, 28255.0, 21851.0, 7678.0, 7231.0, 5702.0, 5028.0,
   4741.0,  2097.0,  1885.0,  1469.0, 1145.0, 1070.0,  645.0,  624.0,
    490.0,   406.0,   280.0,   161.0,  142.0,   13.6,   11.1,   63.0,
     38.0,    18.0,     6.0,

  ## Lr --------------------------------------------------------------1594
 152970.0, 30083.0, 29103.0, 22359.0, 7930.0, 7474.0, 5860.0, 5176.0,
   4876.0,  2180.0,  1963.0,  1523.0, 1192.0, 1112.0,  680.0,  658.0,
    516.0,   429.0,   296.0,   174.0,  154.0,   19.9,   17.0,   71.0,
     44.0,    21.0,     3.9,     6.9,

  ## Rf --------------------------------------------------------------1622
 156288.0, 30881.0, 29986.0, 22907.0, 8161.0, 7738.0, 6009.0, 5336.0,
   5014.0,  2237.0,  2035.0,  1554.0, 1233.0, 1149.0,  725.0,  701.0,
    535.0,   448.0,   319.0,   190.0,  171.0,   26.0,   22.8,   82.0,
     55.0,    33.0,     5.0,     7.5
]



function binding_energies(Z)
    i = INDEX_OF_SHELLS[Z]
    n = NUMBER_OF_SHELLS[Z]

    return ATOMIC_SHELLS[i:i + n - 1] .* co.eV
end
