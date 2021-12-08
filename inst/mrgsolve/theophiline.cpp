[ PROB ]
HMCTMD::OneCptIV.Amikacin

[ SET ] 
req = "CP", end = 30, rtol=1e-5, delta=0.5

[ GLOBAL ]
#define CP (CENT/VD)

[ CMT ] @annotated
CENT : Central compartment (mg)
AUC  : AUC


[ PARAM ] @annotated
WT       : 70      : Body weight (kg)
LBW      : 70      : LBW
CLCR     : 30      : clcr
CL_NR    : 0.0417  : Non-Renal Clearance
CL_SLOPE : 0.815   : Renal Clearnace Slope
VD_NR    : 0.27    : Non-Renal Volume
VD_SLOPE : 0       : Non-Renal Volume Slope
UNIT_CL  : 0.06    : Unit Clearance (60 / 1000)

[ MAIN ]
double CL    = (CL_NR * LBW + CL_SLOPE * CLCR) * UNIT_CL;
double VD    = VD_NR * LBW + (VD_SLOPE * CLCR);

[ ODE ]
dxdt_CENT   =  -CL * CP;

dxdt_AUC = CENT/VD;
if(SS_ADVANCE) dxdt_AUC=0;


[TABLE]
double IPRED = CP;

[ capture ] @annotated
CL      : clearance
VD      : Vd
IPRED   : Individual-predicted concenrtation (mg/L)