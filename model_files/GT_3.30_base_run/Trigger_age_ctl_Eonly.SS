#V3.30.13-safe;_2019_03_09;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.0
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#C estimate either M or h
#_data_and_control_files: Trigger_age_datEonly.SS // Trigger_age_ctl_var_adjEonly.SS
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
5 #_Nblock_Patterns
 2 1 1 1 2 #_blocks_per_pattern
# begin and end years of blocks
 1945 1998 1999 2008
 1945 2012
 2009 2017 #new_block_start
 2009 2017
 1945 2008 2009 2012
#
# controls for all timevary parameters
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  21-24 keep last dev for rest of years
#
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
 #_Age_natmort_by sex x growthpattern
 0.387 0.50841 0.40535 0.35185 0.32007 0.29966 0.28592 0.27635 0.26955 0.26463 0.26104 0.25839 0.25643 0.25496 0.25387 0.25305
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.75 #_Age(post-settlement)_for_L1;linear growth below this
9 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
3 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
#_Age_Maturity by growth pattern
 0 0 0.79 0.91 0.98 0.99 1 1 1 1 1 1 1 1 1 1
2 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
# Sex: 1  BioPattern: 1  Growth
 10 40 17.6 19.9094 1 0 -4 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 20 70 46.3317 44.8 6.5 0 -4 0 0 0 0 0.5 0 0 # L_at_Amax_Fem_GP_1
 0.01 0.6 0.2811 0.14 0.18771 0 -4 0 0 0 0 0.5 0 0 # VonBert_K_Fem_GP_1
 0.05 0.4 0.1967 0.2 99 0 -5 0 0 0 0 0.5 0 0 # CV_young_Fem_GP_1
 0.05 0.2 0.1667 0.2 99 0 -5 0 0 0 0 0.5 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 0 1 2.16e-05 2.16e-05 0.8 0 -2 0 0 0 0 0.5 0 0 # Wtlen_1_Fem_GP_1
 0 4 3.007 3.007 0.8 0 -2 0 0 0 0 0.5 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 25 100 31 31 0.8 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -1 0 -0.065 -0.065 0.8 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 40 60 51.357 51.357 0 0 -3 0 0 0 0 0 0 0 # Eggs_scalar_Fem_GP_1
 1 4 2.8538 2.8538 0 0 -3 0 0 0 0 0 0 0 # Eggs_exp_len_Fem_GP_1
# Hermaphroditism
#  Recruitment Distribution
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -1 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             6            12       8.11776       9.80476            99             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.1             1          0.99       0.45679            99             0         -1          0          0          0          0          0          0          0 # SR_BH_steep
           0.1             1.5      0.679168           0.6            99             0          4          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             0         -4          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0        -99          0          0          0          0          0          0          0 # SR_autocorr
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1982 # first year of main recr_devs; early devs can preceed this era
2016 # last year of main recr_devs; forecast devs start in following year
6 #_recdev phase
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 4 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
1975.1   #_last_early_yr_nobias_adj_in_MPD
1997.4   #_first_yr_fullbias_adj_in_MPD
2017.0   #_last_yr_fullbias_adj_in_MPD
2018.0   #_first_recent_yr_nobias_adj_in_MPD
0.8414   #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 36 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017F 2018F 2019F 2020F 2021F 2022F 2023F 2024F 2025F 2026F 2027F 2028F 2029F 2030F 2031F 2032F 2033F 2034F 2035F 2036F 2037F
#-0.96344 -0.227156 0.414004 0.0617242 0.148486 0.334103 1.48058 1.4009 -1.80909 1.50191 -2.23955 0.697059 0.443105 -0.160053 0.235944 0.0651691 -2.06274 0.00988645 0.667072 0.916548 0.374461 0.399839 0.0465852 -0.0504453 -0.297027 -0.115674 0.184981 -0.392917 -0.185627 -0.409684 -0.492147 -0.127945 0.37616 0.432567 -0.657564 -0.55477 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1982	-0.96344
1983	-0.227156
1984	0.414004
1985	0.0617242
1986	0.148486
1987	0.334103
1988	1.48058
1989	1.4009
1990	-1.80909
1991	1.50191
1992	-2.23955
1993	0.697059
1994	0.443105
1995	-0.160053
1996	0.235944
1997	0.0651691
1998	-2.06274
1999	0.00988645
2000	0.667072
2001	0.916548
2002	0.374461
2003	0.399839
2004	0.0465852
2005	-0.0504453
2006	-0.297027
2007	-0.115674
2008	0.184981
2009	-0.392917
2010	-0.185627
2011	-0.409684
2012	-0.492147
2013	-0.127945
2014	0.37616
2015	0.432567
2016	-0.657564
2017	-0.55477


# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info
0.2 # F ballpark
-2001 # F ballpark year (neg value to disable)
2 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
1.5 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
 0.0005 3 0 # overall start F value; overall phase; N detailed inputs to read
#Fleet Yr Seas F_value se phase (for detailed setup of F_Method=2; -Yr to fill remaining years)

#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#2037 2038
# F rates by fleet
# Yr:  1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036 2037
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# Rec_E_1 0.000281814 0.00100093 0.00542386 0.00899757 0.0126523 0.0179337 0.0222039 0.0266943 0.0314103 0.036336 0.0415135 0.0451284 0.04892 0.0529928 0.0574763 0.0623298 0.0644804 0.0658207 0.0670947 0.0690227 0.0713532 0.07391 0.0753017 0.0766054 0.0785151 0.0809381 0.0830228 0.0875677 0.0932235 0.0965552 0.101189 0.108507 0.119692 0.128118 0.136787 0.143732 0.258412 0.81058 0.227275 0.0958943 0.13946 0.332104 0.59498 0.572308 0.656308 0.714405 0.589476 0.657766 0.495487 0.493039 0.64873 0.277673 0.281574 0.451026 0.530461 0.514654 0.561979 0.957556 1.04314 0.989838 0.825422 0.574426 0.566562 0.389168 0.396661 0.68475 1.14863 0.59902 1.00108 0.423909 0.136189 0.572832 0.349521 0.327897 0.62616 0.62007 0.607396 0.594078 0.582037 0.571332 0.561841 0.553448 0.546037 0.539494 0.533718 0.528616 0.524107 0.520121 0.516595 0.513476 0.510714 0.508269 0.506104
# Com_E_2 1.70511e-06 1.36482e-05 2.7369e-05 3.95784e-05 5.38282e-05 6.85968e-05 8.40886e-05 0.000100451 0.000117795 0.000134263 0.000153741 0.000174359 0.000196099 0.000219292 0.000242134 0.000269247 0.000297618 0.000325436 0.000352065 0.00182728 0.00208236 0.00105385 0.00151557 0.00107703 0.00184487 0.00204666 0.00394412 0.0061753 0.00526223 0.00534255 0.00843636 0.00968908 0.00712442 0.00707787 0.00986345 0.00987789 0.0100624 0.0119292 0.0114951 0.00829041 0.0109761 0.0140597 0.0171968 0.0388401 0.0540926 0.0688815 0.0681487 0.0452682 0.0802601 0.0585465 0.0560726 0.0397722 0.0280331 0.0302022 0.0417257 0.0272913 0.0449986 0.0511698 0.054275 0.0499414 0.0431003 0.0234178 0.0191083 0.0152212 0.0399567 0.0286758 0.0719794 0.0570235 0.0573375 0.0331105 0.0322243 0.0263289 0.0257487 0.0228001 0.0498613 0.0493748 0.0483766 0.0473182 0.0463573 0.0455017 0.0447426 0.0440711 0.0434781 0.0429546 0.0424923 0.042084 0.0417232 0.0414042 0.041122 0.0408723 0.0406513 0.0404556 0.0402822
# Com_LL_2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00523495 0.00518128 0.0119238 0.0331938 0.0440028 0.0377447 0.0519311 0.0643057 0.0932944 0.0802177 0.271207 0.110767 0.000313231 1.04186 0.0900173 0.121817 0.0419479 0.0439565 0.0611645 0.0291655 0.0576895 0.0376484 0.0483956 0.0330835 0.105332 0.23547 0.102654 0.106183 0.0941962 0.15973 0.0779542 0.0200701 0.00989288 0.00801534 0.0125313 0.00505201 0.00609665 0.00822188 0.00420999 0.00323203 0.0109631 0.0109244 0.0107703 0.0104877 0.010158 0.00984044 0.00955218 0.00929564 0.00906937 0.00887075 0.00869672 0.00854425 0.00841055 0.00829311 0.00818982 0.00809891 0.0080188 0.00794816 0.0078858
# Shrimp_Eff_4 0.00170965 0.00636804 0.032522 0.0852559 0.13744 0.243955 0.308555 0.363537 0.375014 0.486048 0.481158 0.614521 0.714684 0.918694 0.985791 0.985253 0.616457 1.04695 1.17959 1.38121 0.908527 0.770423 0.919644 1.07263 1.1704 0.831683 0.677321 1.08979 1.06874 1.05366 1.04186 0.950063 1.22712 1.0384 0.996949 0.628019 1.07429 1.00268 1.18579 1.5 1.30065 1.20484 1.04572 1.09651 1.2722 1.12419 1.03352 1.10836 0.940192 1.08837 1.27291 1.39964 1.38375 1.5 1.13652 0.969917 1.03275 1.29226 0.956431 0.882353 0.717931 0.528283 0.456225 0.38233 0.511422 0.269358 0.336648 0.352362 0.365902 0.295813 0.333301 0.348207 0.363565 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
# HB_E_6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         2         1         0         0         0         0  #  Com_E_2
         4         1         0         0         0         0  #  Shrimp_Eff_4
         5         1         0         0         0         0  #  MRFSS_5
         6         1         0         0         0         0  #  HB_E_6
         7         1         0         0         0         0  #  LarvalGW_7
         8         1         0         0         0         0  #  SEAMAP_Trawl_E_9
         9         1         0         0         0         0  #  Video_GW_10
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            0      -7.04108             0             1             0          1          0          0          0          0          0          0          0  #  LnQ_base_Com_E_2(2)
           10000        100000     34363.8             0             1             0          1          0          0          0          0          0          0          0  #  Q_base_Shrimp_Eff_4(4)
           -15            0      -7.52947             0             1             0          1          0          0          0          0          0          0          0  #  LnQ_base_MRFSS_5(5)
           -15            0      -7.19321             0             1             0          1          0          0          0          0          0          3          1  #  LnQ_base_HB_E_6(6)
           -35            0      -21.3884             0             1             0          3          0          0          0          0          0          0          0  #  LnQ_base_LarvalGW_7(7)
           -15            0       -4.4746             0             1             0          1          0          0          0          0          0          0          0  #  LnQ_base_SEAMAP_Trawl_E_8(8)
           -15            0     -7.27326             0             1             0          1          0          0          0          0          0          0          0  #  LnQ_base_Video_GW_9(9)
# timevary Q parameters
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type     PHASE  #  parm_name
           -15            0      -1.21546             0             1             0      3  # LnQ_base_HB_E_6(6)_BLK3add_2008
# info on dev vectors created for Q parms are reported with other devs after tag parameter section
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 1 2 0 0 # 1 Rec_E_1
 1 2 0 0 # 2 Com_E_2
 1 0 0 0 # 3 Com_LL_2
 24 3 0 0 # 4 Shrimp_Eff_4
 15 0 0 1 # 5 MRFSS_5
 15 0 0 1 # 6 HB_E_6
 0 0 0 0 # 7 LarvalGW_7
 15 0 0 4 # 8 SEAMAP_Trawl_E_8
 1 0 0 0 # 9 Video_GW_9
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 0 0 0 0 # 1 Rec_E_1
 0 0 0 0 # 2 Com_E_2
 0 0 0 0 # 3 Com_LL_2
 0 0 0 0 # 4 Shrimp_Eff_4
 0 0 0 1 # 5 MRFSS_5
 0 0 0 1 # 6 HB_E_6
 0 0 0 0 # 7 LarvalGW_7
 0 0 0 0 # 8 SEAMAP_Trawl_E_8
 0 0 0 0 # 9 Video_GW_9
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   Rec_E_1 LenSelex
            10            50       27.0438            30            99             0          3          0          0          0          0          0          4          1  #  Size_inflection_Rec_E_1(1)
           -10            10       6.43329             5            99             0          3          0          0          0          0          0          0          0  #  Size_95%width_Rec_E_1(1)
             10           73       35.9583            15            99             0          4          0          0          0          0        0.5          1          1  #  Retain_L_infl_Rec_E_1(1)
             0            10       1.05045             1            99             0          3          0          0          0          0        0.5          1          1  #  Retain_L_width_Rec_E_1(1)
           -10            10       1.99012            10            99             0          3          0          0          0          0        0.5          1          1  #  Retain_L_asymptote_logit_Rec_E_1(1)
            -1             2             0             0            99             0         -4          0          0          0          0        0.5          0          0  #  Retain_L_maleoffset_Rec_E_1(1)
           -10            10            -1           -10            99             0         -2          0          0          0          0        0.5          0          0  #  DiscMort_L_infl_Rec_E_1(1)
            -1             2             1             1            99             0         -4          0          0          0          0        0.5          0          0  #  DiscMort_L_width_Rec_E_1(1)
            -1             2          0.25          0.05            99             0         -2          0          0          0          0        0.5          0          0  #  DiscMort_L_level_old_Rec_E_1(1)
            -1             2             0             0            99             0         -4          0          0          0          0        0.5          0          0  #  DiscMort_L_male_offset_Rec_E_1(1)
# 2   Com_E_2 LenSelex
             10           50        28.092            30            99             0          3          0          0          0          0          0          4          1  #  Size_inflection_Com_E_2(2)
            -5            10       5.41644             5            99             0          4          0          0          0          0          0          0          0  #  Size_95%width_Com_E_2(2)
             10           60            35            15            99             0         -3          0          0          0          0        0.5          5          1  #  Retain_L_infl_Com_E_2(2)
           -10            10       2.12681       2.99532            99             0          3          0          0          0          0        0.5          0          0  #  Retain_L_width_Com_E_2(2)
           -10            10       1.84361            10            99             0          3          0          0          0          0        0.5          5          1  #  Retain_L_asymptote_logit_Com_E_2(2)
            -1             2             0             0            99             0         -6          0          0          0          0        0.5          0          0  #  Retain_L_maleoffset_Com_E_2(2)
           -10            10            -1           -10            99             0         -2          0          0          0          0        0.5          0          0  #  DiscMort_L_infl_Com_E_2(2)
            -1             2             1             1            99             0         -4          0          0          0          0        0.5          0          0  #  DiscMort_L_width_Com_E_2(2)
            -1             2          0.25          0.05            99             0         -2          0          0          0          0        0.5          0          0  #  DiscMort_L_level_old_Com_E_2(2)
            -1             2             0             0            99             0         -4          0          0          0          0        0.5          0          0  #  DiscMort_L_male_offset_Com_E_2(2)
# 3   Com_LL_2 LenSelex
            20            73       56.6336           104            99             0          3          0          0          0          0          0          0          0  #  Size_inflection_Com_LL_2(3)
          0.01            30       10.7383            30            99             0          3          0          0          0          0          0          0          0  #  Size_95%width_Com_LL_2(3)
# 4   Shrimp_Eff_4 LenSelex
             0            20         0.008            30            99             0         -3          0          0          0          0          0          0          0  #  Size_DblN_peak_Shrimp_Eff_4(4)
           -15            10      -13.5601            -3         0.001             0         -3          0          0          0          0          0          0          0  #  Size_DblN_top_logit_Shrimp_Eff_4(4)
           -25            10             2             5            99             0         -3          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_Shrimp_Eff_4(4)
           -10            10       4.00318             6            99             0          3          0          0          0          0          0          0          0  #  Size_DblN_descend_se_Shrimp_Eff_4(4)
           -20            15           999             0            99             0         -3          0          0          0          0          0          0          0  #  Size_DblN_start_logit_Shrimp_Eff_4(4)
           -10            20      -5.56368             0            99             0          6          0          0          0          0          0          0          0  #  Size_DblN_end_logit_Shrimp_Eff_4(4)
# 5   MRFSS_5 LenSelex
# 6   HB_E_6 LenSelex
# 7   LarvalGW_7 LenSelex
# 8   SEAMAP_Trawl_E_9 LenSelex
# 9   Video_GW_10 LenSelex
             10           40        28.052       28.8061             7             6          3          0          0          0          0          0          0          0  #  Size_inflection_Video_GW_10(10)
          0.01            20       13.0584       12.8744             1             6          3          0          0          0          0          0          0          0  #  Size_95%width_Video_GW_10(10)
# 1   Rec_E_1 AgeSelex
# 2   Com_E_2 AgeSelex
# 3   Com_LL_2 AgeSelex
# 4   Shrimp_Eff_4 AgeSelex
# 5   MRFSS_5 AgeSelex
# 6   HB_E_6 AgeSelex
# 7   LarvalGW_7 AgeSelex
# 8   SEAMAP_Trawl_E_8 AgeSelex
# 9   Video_GW_9 AgeSelex
# timevary selex parameters
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
           -20            10     -0.894743       20.0922            99             0      5  # Size_inflection_Rec_E_1(1)_BLK4add_2009
           -20           10      -10.8574         30.48            99             0      5  # Retain_L_infl_Rec_E_1(1)_BLK1add_1945
           -20           10      -7.97135         35.56            99             0      5  # Retain_L_infl_Rec_E_1(1)_BLK1add_1999
             0            10      0.996545             2             1             6      5  # Retain_L_width_Rec_E_1(1)_BLK1add_1945
             0            10             0             0             0             0      -5  # Retain_L_width_Rec_E_1(1)_BLK1add_1999
           -10            15      0.618624       20.0922            99             0      5  # Retain_L_asymptote_logit_Rec_E_1(1)_BLK1add_1945
           -10            10     -0.397759            10            99             0      5  # Retain_L_asymptote_logit_Rec_E_1(1)_BLK1add_1999
           -20            10        5.6638            10            99             0      5  # Size_inflection_Com_E_2(2)_BLK4add_2009
           -30            15      -18.1686       20.0922            99             0      5  # Retain_L_infl_Com_E_2(2)_BLK5add_1945
           -20            10      -1.02998            10            99             0      5  # Retain_L_infl_Com_E_2(2)_BLK5add_2009
           -10            15       3.12177       20.0922            99             0      5  # Retain_L_asymptote_logit_Com_E_2(2)_BLK5add_1945
           -10            10             0            10            99             0      -5  # Retain_L_asymptote_logit_Com_E_2(2)_BLK5add_2009
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      3     4     1     3     1     0     0     0     0     0     0     0
#      5     1     2     4     1     0     0     0     0     0     0     0
#      5     3     3     1     1     0     0     0     0     0     0     0
#      5     4     5     1     1     0     0     0     0     0     0     0
#      5     5     7     1     1     0     0     0     0     0     0     0
#      5    11     9     4     1     0     0     0     0     0     0     0
#      5    13    10     5     1     0     0     0     0     0     0     0
#      5    15    12     5     1     0     0     0     0     0     0     0
     #
# Input variance adjustments factors:
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      4      1     0.976
      4      2    0.6074
      4      3    0.8607
      4      8     0.611
      4      9     0.762
 -9999   1    0  # terminator
#
4 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 2 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch;
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 7 1 1 0 1
 7 9 1 0 1
#
# lambdas (for info only; columns are phases)
#first set are for indices component=1
#  1 2  3  0 1   #_COM_E_2 
#  1 4  3  0 1    #_SHRIMP_EFF_4
#  1 5  3  0 1   #_MRFSS_5
#  1 6  3  0 1   #_HB_E_6
#  1 7  3  0 1   #_LARVAL_7
#  1 8  3  0 1   #_SEAMAP_E_8
#  1 9 3  0 1   #_VIDEO_9
#  1 1 1 1 #_discard:_1
#  1 1 1 1 #_discard:_2
#  0 0 0 0 #_discard:_3
#  1 1 1 1 #_discard:_4
#  0 0 0 0 #_discard:_5
#  0 0 0 0 #_discard:_6
#  0 0 0 0 #_discard:_7
#  0 0 0 0 #_discard:_8
#  0 0 0 0 #_discard:_9
#  0 0 0 0 #_discard:_10
#  1 1 1 1 #_lencomp:_1
#  1 1 1 1 #_lencomp:_2
#  1 1 1 1 #_lencomp:_3
#  0 0 0 0 #_lencomp:_4
#  0 0 0 0 #_lencomp:_5
#  0 0 0 0 #_lencomp:_6
#  0 0 0 0 #_lencomp:_7
#  1 1 1 1 #_lencomp:_8
#  1 1 1 1 #_lencomp:_9
#  1 1 1 1 #_lencomp:_10
#  0 0 0 0 #_size-age:_1
#  0 0 0 0 #_size-age:_2
#  0 0 0 0 #_size-age:_3
#  0 0 0 0 #_size-age:_4
#  0 0 0 0 #_size-age:_5
#  0 0 0 0 #_size-age:_6
#  0 0 0 0 #_size-age:_7
#  0 0 0 0 #_size-age:_8
#  0 0 0 0 #_size-age:_9
#  0 0 0 0 #_size-age:_10
#  1 1 1 1 #_init_equ_catch
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
#  0 0 0 0 # F_ballpark_lambda
-9999  1  1  1  1  #  terminator
0 # (0/1) read specs for more stddev reporting
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999


