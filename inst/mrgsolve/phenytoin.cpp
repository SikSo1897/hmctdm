[ PROB ]
OneCpt.Oral.Phenytoin(Non_linear)

[ SET ] 
req = "VD IPRED", end = 30, rtol=1e-5, delta=0.01

[ GLOBAL ]
#define CP (CENT/VD)

[CMT] @annotated
CENT : Central compartment (mg)

[PARAM] @annotated
CrCL          : 1    : CrCL
WT            : 60   : WT
LBW           : 60   : lbw
Vmax          : 500  : Maximum velocity (mg/day)
Km            : 5    : Michaelis constant (mg/L)
VD_NR         : 0.8  : V1 (L)
CL_SLOPE      : 0.01 : cl_slope

[MAIN]
double Norm_Vmax = Vmax * pow(WT/70, 0.6);
double VD = VD_NR * LBW;
double CL = (CL_SLOPE * CrCL) * 60 / 1000;

[ODE]
dxdt_CENT    = -(CL + (Norm_Vmax/24) / (Km + CP)) * CP;

[TABLE]
double IPRED = CP;

[ capture ] @annotated
Vmax       : Vmax
Km         : Km
VD         : VD
CL         : CL
IPRED      : Individual-predicted concenrtation (mg/L)