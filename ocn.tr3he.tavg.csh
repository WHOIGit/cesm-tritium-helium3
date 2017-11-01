#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/pop2conf/tr3he_tavg_contents << EOF
$s1  SHUM
$s1  TR3HE_IFRAC
$s1  TR3HE_ATM_PRESS
$s1  HELIUM3_KS
$s1  HELIUM3_KB
$s1  HELIUM3_SURF_SAT
$s1  STF_HELIUM3
$s1  STF_HELIUM3_DGE
$s1  STF_HELIUM3_CTB
$s1  STF_HELIUM3_PTB
$s1  INERT_HELIUM3_KS
$s1  INERT_HELIUM3_KB
$s1  INERT_HELIUM3_SURF_SAT
$s1  STF_INERT_HELIUM3
$s1  STF_INERT_HELIUM3_DGE
$s1  STF_INERT_HELIUM3_CTB
$s1  STF_INERT_HELIUM3_PTB
$s1  HELIUM4_KS
$s1  HELIUM4_KB
$s1  HELIUM4_SURF_SAT
$s1  STF_HELIUM4
$s1  STF_HELIUM4_DGE
$s1  STF_HELIUM4_CTB
$s1  STF_HELIUM4_PTB
$s1  NEON_KS
$s1  NEON_KB
$s1  NEON_SURF_SAT
$s1  STF_NEON
$s1  STF_NEON_DGE
$s1  STF_NEON_CTB
$s1  STF_NEON_PTB
$s1  ARGON_KS
$s1  ARGON_KB
$s1  ARGON_SURF_SAT
$s1  STF_ARGON
$s1  STF_ARGON_DGE
$s1  STF_ARGON_CTB
$s1  STF_ARGON_PTB
$s1  STF_TRITIUM_EVAP
$s1  STF_TRITIUM_EVAP_UP
$s1  STF_TRITIUM_EVAP_DOWN
$s1  STF_TRITIUM_PREC
$s1  STF_TRITIUM
$s1  HQ
$s1  HELIUM3
$s1  INERT_HELIUM3
$s1  HELIUM4
$s1  NEON
$s1  ARGON
$s1  TRITIUM
$s1  TRITIUM_DECAY
EOF

if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/pop2conf/tr3he_tavg_contents << EOF
$s1  KPP_SRC_TRITIUM
$s1  KPP_SRC_HELIUM3
$s1  DIA_IMPVF_TRITIUM
$s1  DIA_IMPVF_HELIUM3
$s1  HDIFE_TRITIUM
$s1  HDIFE_HELIUM3
$s1  HDIFN_TRITIUM
$s1  HDIFN_HELIUM3
$s1  HDIFB_TRITIUM
$s1  HDIFB_HELIUM3
$s1  UE_TRITIUM
$s1  UE_HELIUM3
$s1  VN_TRITIUM
$s1  VN_HELIUM3
$s1  WT_TRITIUM
$s1  WT_HELIUM3
EOF
endif

#===============================================================================
# The following are fields computed by the CFC modules that are not placed in
# the tavg file by default.
#
#1  HELIUM3_SCHMIDT
#1  HELIUM3_PV
#1  HELIUM3_surf_sat
#===============================================================================
