// This file mainly contains the definitions of types and variables.

typedef struct {
  int **l, c, siz;
} ListInt2;
typedef struct {
  int *l, c, siz;
} ListInt;
typedef struct {
  double *l;
  int c;
} ListDbl;
typedef struct {
  int *cntU, *cntB, *cntW;
} AbpDyn;
// Duration
typedef struct {
  long long dur;
  double durR, facK;
} Duration;
// Rheological measurements
typedef struct {
  double dsp, acc, accDsp, accDspUS, lim[2], accDspFl;
} Strain;
typedef struct {
  double *acc, accShUS, accShErr, curr, goal;
} Stress;
typedef struct {
  int dur, tgl;
  double rate, magDsp, mag, facK;
} PreStraStre;
typedef struct {
  int prd;
  double prdR, amp, magDsp;
} SinuStraStre;
// Dynamic behaviors on boundaries
typedef struct {
  int gTgl, gTglPa, gTglSpr, *cntMe;
  double *p, *facK, dep, depR;
} BndReb;
typedef struct {
  int gTgl, *cntMe, maxF;
  double *p, *facK0, *facX;
} BndUnb;
typedef struct {
  int gTgl, stfUnit;
  double stf[3*2], thk[6];
  ListInt rankF;
} BndMv;
// Dynamic behaviors of actins
typedef struct { 
  int tgl, gTgl, cntMe, gTglCap[2], gTglBst[2];
  double k, *p, *pWA, facKWA, facX, critF;
} Sever;
typedef struct { 
  int tgl, gTgl, cntMe;
  double k, ang, p;
} Anneal;
typedef struct {
  int tgl, gTgl, gTglFN, cntMe, cntFNme;
  double k, facP;
} Nucle;
typedef struct {
  int tgl, gTgl, cntMe;
  double k[2], p[2], pWA[2], facP[2], facKWA;
} AssDisCapUncAge;
typedef struct {
  int tgl, gTgl, cntMe, *l, c;
  double k, p, facP, dist, *l2, x1, x2;
} Degrade;
typedef struct {
  int gTgl, *cntMe2, maxF, maxNFA, gTglPa;
  double *p, facK0, facX;
} Mature;
// Dynamic behaviors of ABPs
typedef struct {
  int tgl, gTgl, cntMe;
  double facP, k, facK;
} InaUnbMbind;
typedef struct {
  int tgl, gTgl, gTglPa, gTglCrsAng, gTglOppDir, cntMe;
  double p, facK, por;
} Reb;
typedef struct {
  int tgl, gTgl, cntMe, maxF;
  double *p, facK0, facX, facK0c, facXc;
} Unb;
typedef struct {
  int tgl, gTgl, cntMe, maxF[2];
  double *p, facK0, facX;
} MotUnbWalk;
typedef struct {
  int *cnt, nL;
  double *len, *allF, *sprF, *bendF;
} RecLenForce;
typedef struct {
  int prd, tgl, gTgl, mode;
  double prdR;
} FuncCont;
typedef struct {
  int prd, tgl, gTgl, tgl2, gTglBnd, gTglInfo, mode;
  double prdR, minF, maxF;
} RecConfVmd;
typedef struct {
  int prd, tgl, tgl2, gTgl, gTglCho, nActL, nAbpL;
  double prdR, porActL, porAbpL;
  ListInt actMe, abpMe, act, abp;
} RecTraj;
typedef struct {
  double *r, *v, *f, *fBr, *maxFbr, *rPrev, *bdRFix, *bndRFix, *distMb;
  int *ch, *id, *fix, *iF, *mbReb, *nFA, *bdFix, *cap, *age, *len;
  ListInt cyl;
} MoleAct;
typedef struct {
  double *r, *v, *f, *fBr, *maxFbr, *distMb;
  int *ch, *id, *mId, *age;
} MoleAbp;
typedef struct {
  double n, inv;
} Double;
typedef struct {
  int gTgl;
  double stf, stf2, eq, eq2, lo, hi, facStf, drag, *val;
} Force;
typedef struct {
  int gTglActMv, gTglRnd, gotF;
  double stfRep, radRnd, *drag, *r, *f, *stfSpr;
} Boundary;
typedef struct {
  double wid[3], base[3];
  int n[3], *l;
} Cell;
typedef struct {
  double *dia, *repDia, maxSiz;
  Double *len, *dtRepAct, *dtSprAct, *drag;
  Force *a90, *cr, *bend, *spr, rep;
} AbpF;
typedef struct {
  int gTgl, maxF, mode;
  double *p, pf, k, k0, x;
} MotTurnover;
typedef struct {
  Force bend, spr;
  MotTurnover to;
  int tgl, gTgl, gTglConSiz, *id;
  int cntNucMe, cntAssMe, cntTurnMe, nMotPerTF, nNucMe, nMotPerSide;
  double cenDist, kAss, pAss;
} MotSelfAss;
typedef struct {
  int nHead;
  double k01, k10, k12, k21, k20;
} MotMechChem;
typedef struct {
  double dia, dragR;
  Force bend, spr, rep;
} ActF;
// Actin bursting depolymerization
typedef struct {
  int tgl, gTgl, cntMe;
  double k[2], p[2], facKWA;
  ListInt fil;
} Burst;
// Membrane
typedef struct {
  int tgl, gTgl, gTglPa, cntMe;
  double p, facK, por, facK0, facX;
} MembNucReb;
typedef struct {
  int tgl, gTgl, cntMe, dur, gTglCT;
  double p, facK, f, durR, rCT[3];
} MembPro;
typedef struct {
  double stf, diff, *conc;
} MembMyoCyto;
typedef struct {
  int *ch, *id, *reb, *fix, *pro, *nFA, *idx, *bias;
  double *r, *v, *f, *fBr, *nDir, maxFbr, *rPrev, *rFix, *areaL;
  double thk, dragR, len, *eqLen, *unitEqAr, *unitNDir;
  MembMyoCyto myo, cyto;
  ListInt unit, unitCp;
  Force bend, spr, spr2, spr3, area, vol;
  double nucThk, nucDragR, nucLen, nucMaxFbr, tensF;
  Force nucBend, nucSpr, nucSpr2, nucArea, nucVol;
} MembNuc;
typedef struct {
  int *l, *reb, *side, nReb, nRebMe, cntMe, gTgl;
  double *info, por;
} MembSldList;
typedef struct {
  int gTgl, prd;
  double critDist, eqDist, stf, drag;
  MembSldList act, abp;
} MembSld;

typedef struct {
  int n, rheoDir, rheoPrd, rheoWay, rheoSig, rheoUpdPrd;
  int rheoPreTgl;
  double thkRep, stfRep, facStfRep, rad2, *maxFbr, *rad, *r, *f, *drag;
  double rheoAmp, rheoPrdR, rheoAmp2, *rheoFAcc;
  double rheoPreRate, rheoPreMag, rheoPreMagDsp;
  double rheoFGoal, *rheoFAccErr, *rheoFCurr;
} Bead;
typedef struct {
  int *cntMe, gTgl;
} BeadBind;

#include <mpi.h>
#include "header.h"
#include "bead.c"
#include "boundary.c"
#include "calForce.c"
#include "dataMgr.c"
#include "error.c"
#include "gatPrint.c"
#include "init.c"
#include "membrane.c"
#include "paraProc.c"
#include "process.c"
#include "record.c"
#include "rheology.c"
#include "rng.c"
#include "tools.c"
#include "update.c"

// Variables related to time steps
long long currTimeStep;
double dt, dtReal;
time_t initTime;
Duration netForm, rheo, motActiv;
// Given concentration
int nAbpGoal, nAbpDet[NK_ABP], nAbpMall[NK_ABP], nAbpGoalDet[NK_ABP];
double cAct, RAbp[NK_ABP];
// Number of particles
int nAct, nAbp, nMot; 
int nActMe, nActCp, nActMeCp, nActC;
int nAbpMe, nAbpCp, nAbpMeCp, nAbpC;
int nAcpInaMe, nMotInaMe, nAcpMme, nMotMme;
int nActGoal, nActMall, nActFilaMe;
int nActMin, nAbpMin, nMbMin; 
// Chain, position, forces, and lists of particles
ListInt actM, acpM, motM, iFilaP;
MoleAct act;
MoleAbp abp;
ActF actF;
AbpF abpF;
MotSelfAss motSA;
int *chAct, *chAbp, cntAbpDC, tglNeiAbpSC, tglNeiAbpDC;
int nChAc, nChAcX, nChAcY, nChAb, nActPerSeg, confNChAcX, confNChAcY;
double *rAct, *rAbp;
// Update neighboring list
ListInt neigh;
Cell cell;
double dispHeuSq, maxNeiLenSq;
// Arrays related to lists
ListInt sendAbpDyn, noAbpDyn;
ListInt sendActDyn, noActDyn;
ListInt sendMbDyn, noMbDyn;
int *fixAct, *abpMotId;
// Domain (width, periodic boundary conditions, or repulsive force)
int pbc[NDIM], neiPbc[NDIM], neiPbcMb[NDIM], confPbc[NDIM], dir2D;
double dimDom[NDIM], dimDomH[NDIM], minDimDomC;
Boundary bnd;
// For boundaries
FuncCont recBndLoc, updSubdSize, recBndUnbReb, recBndActMat;
FuncCont recBndTracF;
double *rGridInit, volMe;
Mature bndMat;
BndReb bndReb;
BndUnb bndUnb;
Force bndVol;
BndMv bndMv;

/*-------------------------- For parallel processing -------------------------*/
int mpiMethod, nCpu;
int *iAct, *iAbp;
// Ranks of CPUs
int rank, *adjRank, *iRank, *cntAdjRank;
// For boundaries and indicies of subdomains
int nCell[NDIM], iCell[NDIM], nGrid[NDIM];
double edge[NDIM*2], **rGrid, neiEdge;
// Used in CopyParticles and MoveParticles
int *cntCpPar, *cntMvPar, modeActCh;
ListInt *cpPar, *mvPar, insNeiPar; //, mvParAll;
// Varilables related to messages
int sizeBufMsg, *mpiTestSendFlag, *mpiTestRecvFlag;
char **bufSendMsg, **bufRecvMsg;
MPI_Status status;
MPI_Request *sReq, *rReq;
// Variables related to longCh
ListInt longCh, longChExtMsg, longChIntMsg;
double *longChDist, maxDisp, maxActCh;
/*----------------------------------------------------------------------------*/

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/
FuncCont updMono;
// Dynamics of actin
AssDisCapUncAge actAss, actDis, actCap, actUnc, actAge;
Burst actBst;
Sever actSev;
Anneal actAnn;
Nucle actNuc;
Degrade actDgd;
int tglActMoDyn, tglActDyn, tglActFormDyn;
int tglActDisBstSev, tglActDisBstSevAbp;
int gTglActDynPres, gTglActDynNF;
int gTglActSevCap, gTglActSevBst, gTglActCapAll;
int gTglActTherm, gTglAcpTherm, gTglMotTherm;
int durNoActDyn, cntDisFila;
// Dynamics of ACPs and motors
AbpDyn abpDyn;
Reb acpReb, acp, motReb;
Unb acpUnb;
MotUnbWalk motUnb, motWalk;
InaUnbMbind acpInaUnb, motInaUnb, acpMoBind, motMoBind;
FuncCont updMotN, abpAge;
int tglAbpAcInaDyn, tglAbpInaMoDyn, tglAcpAcInaDyn, tglMotAcInaDyn;
int gTglAcpDynPres, gTglMotUnbRebPres, gTglMotWalkPres;
int gTglAcpDynNF, gTglMotUnbRebNF, gTglMotWalkNF, gTglMotWalkSld;
int gTglImpAcpM, gTglImpMotM;
int durNoAcpUnbReb, durNoMotUnbReb, durNoMotWalk;
ListDbl unbLog, bindLog, toLog, sevLog;
MotMechChem motMC;
// Dynamics of membrane
int tglMbDyn, tglMbNucDyn, durNoMbDyn, cntNucAssAnn;
int tglMbNucLocAre; 
/*----------------------------------------------------------------------------*/
 
/*----------------------- For rheological measurements -----------------------*/
// For both ways
FuncCont recStre, recProg, recConf;
int rheoWay;
// Segment rheology
RecTraj recTraj;
// Bulk rheology
ListInt meaStrePar, appStraPar, meaStreParMe, appStraParMe, rankMeaStre;
Strain stra;
Stress stre;
PreStraStre pres;
SinuStraStre sinuStr;
int dirStr, dirNoPBC, dirOther, prdUpdSinuStre;
int bulkRheoWay, bulkRheoType, signStr, gotMeaStre;
/*----------------------------------------------------------------------------*/

/*----------------------------- For microrheology ----------------------------*/
int gTglBead;
Bead bead;
BeadBind beadBind;
/*----------------------------------------------------------------------------*/

/*---------------------------- For data recording ----------------------------*/
// For general records
char fnOut[80], dataFold[80];
// Record the accumulated chain lengths and forces
RecLenForce recAct, recAbp, recMb;
// Record configuration for VMD
RecConfVmd recConfVmd;
FuncCont recConfPV;
// Record configuration for Matlab
FuncCont recConfMlb;
// Record internal elastic/viscous stresses
FuncCont recSecStre;
int nSecStreDiv[NDIM];
double *viscStre, *elasStre;
// Record the percolation of networks
FuncCont recPerc;
int *iPerc, dirRecPerc;
// Find supportive framework
FuncCont findSupp;
int kindFindSupp, *iSupp;
double porFindSupp;
// Record the turnover of ACPs or motors 
int tglRecAbpTurn;
double *abpTurn;
// Record (longitudinal) forces of ACPs or motors
FuncCont recLongF;
double *recLongSprFabp, *recInstSprFabp;
// Record etc
FuncCont confVmdInfo, recFilaL, recE, recInfo;
FuncCont recMotSize, recMotPos, recCrsDist, recAbpDyn, recConn;
FuncCont recPoreSize, recAbpUnb, recAbpBind, recAbpTurn, recActSev;
FuncCont recMbCen, recMbDim, recBeadLoc;
/*----------------------------------------------------------------------------*/

/*--------------------------- Membrane & nucleus------------------------------*/
int gTglLoadMbNucData, dirNormMbNuc, dimMbNuc;
int nMb, nMbMe, nMbCp, nMbC, nUnitMb, nChMb, nObjMbNuc;
int *chMb, *iMb, *unitMb, *rebMb, *actRebMb, dirMbNuc[NDIM], *idxMb;
double *rMb, *nDirMb, *eqLenMb, *cenMb, *nDirUnitMb, *eqArUnitMb;
double stfRepMbNuc;
Cell cellMb;
ListInt neighMb;
/*---------------------------------- Membrane --------------------------------*/
int gTglMb, gTglMbTherm, gTglMbActNuc, gTglMbCont;
int nMbAct, nMbPerObj;
int levMb, sideMb;
double radMb, radMbInit, thkMbActNuc;
MembNuc memb;
Unb mbUnb, mbDef;
Reb mbReb, mbFix;
Mature mbMat;
MembPro mbPro;
FuncCont mbVol, mbAre;
MembSld mbSld;
/*----------------------------------- Nucleus --------------------------------*/
int gTglNuc, gTglNucTherm, gTglNucActNuc;
int nNucAct, nNucPerObj;
int levNuc;
double radNuc, radNucInit;
Unb nucUnb;
MembNucReb nucReb;
FuncCont nucVol, nucAre;
/*----------------------------------------------------------------------------*/

// Check errors
int stopSig;
double magUnstF;
// Misc.
int gTglLoadNetData, gTglLoadNetDataFix;
int seed, *allIntL;
double *arrAcos, *allDblL;

int typeBias, tglBias;
double thkBias, degBias;
