// This file contains the definitions of parameters, function prototypes
// variables, and functions.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
/******************************************************************************/
/************************** Definition of Parameters **************************/
/******************************************************************************/

// General parameters and constants
#define NDIM 					3		
#define NK_ABP					3
#define KT						1.0			// Boltzmann Energy
#define KT_IN_J					4.142e-21	// Boltzmann Energy in J
#define VISCOSITY				0.8599e1	// Viscosity of water
#define PI						3.141593	
#define N_AVO					6.022e23
#define POS_LARGE_VALUE			1e9
#define NEG_LARGE_VALUE			-1e9
#define POS_SMALL_VALUE			1e-10		
// For neighboring list
#define DT_NL_UPDATE			0.4
#define DT_NL_UPDATE_BUF		0.1
// Forces
#define MAG_UNSTABLE_FORCE		10000.0e-12
// For stress-controlled bulk rheology
#define KPI_STRESS_P			1.0e-5*0.01
//#define KPI_STRESS_P			1.0e-5*0.1
#define KPI_STRESS_I			1.0e-7*0.1
#define PRD_UPD_SINU_STRESS		1.0e-1
// For active microrheology using force
#define KPI_FORCE_P				1.0e-5
#define KPI_FORCE_I				1.0e-7
#define PRD_UPD_SINU_FORCE		1.0e-4
// Record general data
#define OUTPUT_FILE				"Output"
#define DELETE_FILE				1		
// Extent of pre-assigned arrays 
#define DEG_ARR_ACOS			10000
#define DEG_ARR_ACTSEV			10000
// For parallelization (MoveParticles and CopyParticles)
#define OFFSET_LIST															\
	{{0,1,2,3}, {0,1,2,3,4}, {0,3,4}, {0,1,2,3}, {0,1,2,3,4}, {0,3,4}, 		\
	{0,1}, {0,1}, {0}, {0,1,2,3,5,6,7,8}, {0,1,2,3,4,5,6,7,8,9,10}, 		\
	{0,3,4,5,8,9,10}, {0,1,2,3,5,6,7,8,12,13},					 			\
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13}, {0,3,4,5,8,9,10,11,12}, 		 		\
	{0,1,5,6,12,13}, {0,1,5,6,10,11,12,13}, {0,5,10,11,12}}
#define OFFSET_LEN	{4,5,3,4,5,3,2,2,1,8,11,7,10,14,9,6,8,5}
#define OFFSET_CPU_LIST														\
	{{0,1,3,4,9,10,12}, {1,4,10}, {1,2,4,5,10,11,14}, {3,4,12}, {4}, 		\
	{4,5,14}, {3,4,6,7,12,15,16}, {4,7,16}, {4,5,7,8,14,16,17}, {9,10,12},	\
	{10}, {10,11,14}, {12}, {-1}, {14}, {12,15,16}, {16}, {14,16,17},		\
	{9,10,12,18,19,21,22}, {10,19,22}, {10,11,14,19,20,22,23}, {12,21,22},	\
	{22}, {14,22,23}, {12,15,16,21,22,24,25}, {16,22,25}, 					\
	{14,16,17,22,23,25,26}}
#define OFFSET_CPU_LEN {7,3,7,3,1,3,7,3,7,3,1,3,1,1,1,3,1,3,7,3,7,3,1,3,7,3,7}

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/

// For severing of actins
#define MAX_ACT_SEV_FORCE	400e-12
#define X_ACT_SEV			2e-9
// For annealing of actins
#define DTL_ACT_ANN			0.95
#define DTH_ACT_ANN			1.05
// For maturation of actins
#define MAX_ACT_MAT_FORCE	500e-12
#define K0_ACT_MAT			0.01
#define X_ACT_MAT			400.0e-12
// For unbinding of filaments on boundaries
#define MAX_BND_UNB_FORCE	500e-12
#define K0_BND_UNB			1
#define X_BND_UNB			100.0e-12
#define K_BND_REB			1
// Unbinding/binding/binding of ACPs
#define MAX_ACP_UNB_FORCE   500.0e-12
#define K0_ACP_UNB			0.115
#define X_ACP_UNB			104.0e-12
#define K_ACP_BIND			1.0e2
// Binding of motors
#define K_MOT_BIND			1.0e8
// Turnover of motor filaments
#define MAX_MOT_TURN_FORCE   500.0e-12
// Lower and upper limits of distances allowing binding of ACPs or motors
#define DTL_ACP_REB			0.9	
#define DTH_ACP_REB			1.1	
#define DTL_MOT_REB			0.7037	
#define DTH_MOT_REB			1.2963	
// Time for which a dynamic behavior is prohibited after a previous one occurs
#define DUR_ACP_NO_UNB		2e-9    
#define DUR_MOT_NO_UNB		2e-9   
#define DUR_MOT_NO_WALK		2e-9    
#define DUR_ACT_NO_DYN		2e-9    
#define DUR_MEMB_NO_DYN		2e-9    
/*----------------------------------------------------------------------------*/

/*----------------------------------- Membranes ------------------------------*/
// Volume conservation
#define STF_MEMB_VOL		0.		// in Pa
//#define STF_MEMB_VOL2		2250.	// in Pa
#define STF_MEMB_VOL2		1000.	// in Pa
#define VISC_MEMB_VOL		100.	// in Ns/m^2
// Area conservation
//#define STF_MEMB_AREA		8.284e1
#define STF_MEMB_AREA		1e-3	 // in N/m
// Bending and extensional stiffness and relevant parameters
#define STF_MEMB_BEND		2.4e-19		// in Nm
#define STF_MEMB_SPR		1e-4
#define STF_MEMB_ACT_SPR	1e-3
#define STF_MEMB_FIX_SPR	1e-3
#define DTL_MEMB_SPR		1
#define DTH_MEMB_SPR		1
// Repulsive force
#define STF_REP_MEMB		(1.691e-3*3)
// Parameters involved with actomyosin contractility
#define STF_MEMB_SPR_MYO	720e-6		// in N/m
#define DIFF_MEMB_SPR_MYO	0.0987		// in 1/s
#define STF_MEMB_SPR_CYTO	240e-6		// in 1/s
#define DIFF_MEMB_SPR_CYTO	0.025		// in 1/s
// Unbinding
#define MAX_MEMB_UNB_FORCE	1000.0e-12	// in N
#define K0_MEMB_UNB			1.			// in 1/s
#define X_MEMB_UNB			500e-12		// in m
// Binding
#define K_MEMB_REB			1000.		// in 1/s
#define DTL_MEMB_REB		0.5			// relative to eq length
#define DTH_MEMB_REB		1.5			// relative to eq length
// Fixation
#define K_MEMB_FIX			100.		// in 1/s
// Release
#define K0_MEMB_DEF			0.115		// in 1/s
#define X_MEMB_DEF			104.0e-12	// in m
#define MAX_MEMB_DEF_FORCE	400.0e-12	// in N
// Protrusion
#define K_MEMB_PROT			4e-4		// in 1/s
#define MAX_V_MEMB_PROT		0.2e-6		// in m/s
#define DUR_MEMB_PROT		7.5			// in s
// Maturation
#define K0_MEMB_MAT			0.01		// in 1/s
#define X_MEMB_MAT			400.0e-12	// in m
#define MAX_MEMB_MAT_FORCE	500.0e-12	// in N
/*----------------------------------------------------------------------------*/

/*------------------------------------ Nucleus -------------------------------*/
// Volume conservation
#define STF_NUC_VOL			0.
#define STF_NUC_VOL2		2250.		// in Pa
// Area conservation
#define STF_NUC_AREA		8.284e1
// Bending and extensional stiffness
#define STF_NUC_BEND		0.215e-18
#define STF_NUC_ACT_SPR		1e-3
#define STF_NUC_SPR			5e-3
#define DTL_NUC_SPR			1.0
#define DTH_NUC_SPR			1.0
// Repulsive forces
#define STF_REP_BEAD		1.691e-4
// Unbinding
#define MAX_NUC_UNB_FORCE	1000.0e-12
#define K0_NUC_UNB			1.
#define X_NUC_UNB			0.5e-9
// Binding
#define K_NUC_REB			1000.
#define DTL_NUC_REB			0.5
#define DTH_NUC_REB			1.5
/*----------------------------------------------------------------------------*/

/*-- Mechanical stiffness and geometric parameters of actin, ACP, and motor --*/

// Mechanical stiffness for actin
#define STF_ACT_SPR			1.691e-1
#define STF_ACT_BEND		3.697e-26
#define STF_ACT_REP			1.691e-3
#define STF_ABP_REP			1.691e-3
// Geometric parameters for actin
#define ANGL_ACT_BEND		10.
#define DIA_CYL_ACT			7.0e-9
#define DIA_CYL_ACT_UM		(DIA_CYL_ACT * 1e6)
#define DIA_CYL_ACT_NM		(DIA_CYL_ACT * 1e9)
// Mechanical stiffness for ACP crosslinking filaments at right angle
#define STF_ACPC_SPR		2e-3	
#define STF_ACPC_BEND		1.036e-18	 
#define STF_ACPC_90			1.036e-19
#define STF_ACPC_TOR		4.142e-18
// Geometric parameters for ACP crosslinking filaments at right angle
#define L_ACPC_ARM			20.0e-9
#define ANG_ACPC_BEND		0
#define ANGL_ACPC_BEND		15.
#define ANGL_ACPC_90		10.
#define ANGL_ACPC_TOR		10.
#define DIA_CYL_ACPC		10.0e-9
// Mechanical stiffness for ACP bundling filaments in parallel
#define STF_ACPB_SPR		2.0e-3
#define STF_ACPB_BEND		1.036e-19
#define STF_ACPB_90			1.036e-19
#define STF_ACPB_TOR		4.142e-18
// Geometric parameters for ACP bundling filaments in parallel
#define L_ACPB_ARM			20.0e-9
#define ANG_ACPB_BEND		0.
#define ANGL_ACPB_BEND		15.
#define ANGL_ACPB_90		10.
#define ANGL_ACPB_TOR		10.
#define DIA_CYL_ACPB		10.0e-9
// Mechanical stiffness for motors
#define STF_MOT_SPR			1.0e-3
#define STF_MOT_SPR2		1.0e-3
#define STF_MOT_BEND		1.036e-18	 
#define STF_MOT_90			0
#define STF_MOT_TOR			4.142e-18
#define MOT_STALLF			5.04e-12
// Geometric parameters for motors
#define L_MOT_ARM			10.0e-9	
#define ANG_MOTBACK_TOR		180.
#define ANG_MOTBACK_BEND	0.
#define ANG_MOT_TOR			180.
#define ANG_MOT_BEND		0
#define ANGL_MOT_BEND		15.
#define ANGL_MOT_90			10.
#define ANGL_MOT_TOR		10.
#define DIA_CYL_MOT			10.0e-9
// Mechanical stiffness for motors
#define STF_MOTBACK_SPR		1.691e-2
#define STF_MOTBACK_BEND	5.072e-18
// Geometric parameters for motors
#define L_MOTBACK_DIST		42.0e-9	
#define L_MOTBACK_CEN_DIST	42.0e-9	
#define ANG_MOT_BEND		0
// Boundaries
#define STF_REP_BND			1e-3
#define STF_SPR_BND			1e-3
#define STF_BND_VOL         1e19
// Degree of coarse-graining using cylindrical segments and 
// the corresponding length of one actin segments
#define L_SCALE_IN_M			(nActPerSeg * DIA_CYL_ACT)
#define L_SCALE_IN_UM			(nActPerSeg * DIA_CYL_ACT_UM)
#define L_SCALE_IN_NM			(nActPerSeg * DIA_CYL_ACT_NM)
/*----------------------------------------------------------------------------*/

/******************************************************************************/
/***************************** Definition of Tools ****************************/
/******************************************************************************/
// Alternative forms of arrays
#define P2(a,b,c)			a[(b)*NDIM+(c)]
#define P2A(a,b,c,d)		a[(b)*(d)+(c)]
// Frequently used expressions 
#define SELECT(a, b, c, d, e)											\
	if ((a) > (b)) { (c)=(d); } else { (c)=(e); }						
#define SELECT_LOHI(a, b, c, d)											\
	if ((a) > (b)) { (c)=(b); (d)=(a); } else { (c)=(a); (d)=(b); }
#define BREAK(a)			if (a) { break; }
#define CONT(a)             if (a) { continue; }
#define ISACTM(a)		    (P2A(act.ch,(a),0,nChAc) < 0  				\
		&& P2A(act.ch,(a),1,nChAc) < 0)
#define ISABPM(a)		    (P2A(abp.ch,(a),0,nChAb) < 0 				\
		&& P2A(abp.ch,(a),1,nChAb) < 0 && (K_ABP(a) != 2 				\
		|| motSA.gTgl == 0 || (P2A(abp.ch,(a),3,nChAb) < 0 				\
		&& P2A(abp.ch,(a),4,nChAb) < 0 && K_ABP(a) == 2 				\
		&& motSA.gTgl != 0)))
#define ISABPIM(a)		    (P2A(abp.ch,(a),0,nChAb) < 0 				\
		&& P2A(abp.ch,(a),1,nChAb) < 0 && ((K_ABP(a) != 2 				\
		&& gTglImpAcpM != 0) || (K_ABP(a) == 2 && gTglImpMotM != 0 		\
		&& motSA.gTgl == 0) || (P2A(abp.ch,(a),3,nChAb) < 0 			\
		&& P2A(abp.ch,(a),4,nChAb) < 0 && K_ABP(a) == 2 				\
		&& motSA.gTgl != 0)))
#define ISNUC(a)			(gTglNuc != 0 && (a) >= nObjMbNuc / 2)
#define ISMTF(a)			(motSA.gTgl != 0 && (a) == 2)
#define K_ABP(a)			abp.ch[(a)*nChAb+2]
#define SPRING(a,b,c)		REVSIGN(a) * ((b) - (c))
#define CYL_DRAG(a,b)		(3.*PI*VISCOSITY*(a)*(3.+2.*(b)/(a))/5.)
// Frequently used for-loops
#define FOR_ACT(a)			for((a) = 0; (a) < nAct; (a)++)
#define FOR_ABP(a)			for((a) = 0; (a) < nAbp; (a)++)
#define FOR_NDIM(a)			for((a) = 0; (a) < NDIM; (a)++)
#define FOR_ACTME(a)		for((a) = 0; (a) < nActMe; (a)++)
#define FOR_ABPME(a)		for((a) = 0; (a) < nAbpMe; (a)++)
#define FOR_ACTCP(a)		for((a) = 0; (a) < nActCp; (a)++)
#define FOR_ABPCP(a)		for((a) = 0; (a) < nAbpCp; (a)++)
#define FOR_ACTMECP(a)		for((a) = 0; (a) < nActMe + nActCp; (a)++)
#define FOR_ABPMECP(a)		for((a) = 0; (a) < nAbpMe + nAbpCp; (a)++)
#define FOR_ACTC(a)			for((a) = 0; (a) < nActC; (a)++)
#define FOR_ABPC(a)			for((a) = 0; (a) < nAbpC; (a)++)
#define FOR_MB(a)			for((a) = 0; (a) < nMb; (a)++)
#define FOR_MBME(a)			for((a) = 0; (a) < nMbMe; (a)++)
#define FOR_MBCP(a)			for((a) = 0; (a) < nMbCp; (a)++)
#define FOR_MBMECP(a)		for((a) = 0; (a) < nMbMe + nMbCp; (a)++)
#define FOR_MBC(a)			for((a) = 0; (a) < nMbC; (a)++)
// Algebra
#define AVG2(a,b)			(0.5*((a)+(b)))
#define REVSIGN(a)			(-1.*(a))
#define INV(a)				(1./(a))
#define SQR(a)				((a)*(a))
#define CUBE(a)				((a)*(a)*(a))
// Alternative forms of memory allocation
#define MALLOC(a,b,c)		(a)=(b *)malloc(sizeof(b)*(c))
#define MALLOC2(a,b,c)		(a)=(b **)malloc(sizeof(b *)*(c))
#define MALLOC3(a,b,c)		(a)=(b ***)malloc(sizeof(b **)*(c))
#define MALLOC4(a,b,c)		(a)=(b ****)malloc(sizeof(b ***)*(c))
// Packing and unpacking data in meassages
#define MPI_PACK_DBL(a, b, c)		MPI_Pack((a), (b), MPI_DOUBLE,		\
	bufSendMsg[(c)], sizeBufMsg, &posi, MPI_COMM_WORLD);
#define MPI_PACK_INT(a, b, c)			MPI_Pack((a), (b), MPI_INT,		\
	bufSendMsg[(c)], sizeBufMsg, &posi, MPI_COMM_WORLD);
#define MPI_UNPACK_DBL(a, b, c)		MPI_Unpack(bufRecvMsg[(c)],			\
	sizeBufMsg,&posi,(a),(b), MPI_DOUBLE, MPI_COMM_WORLD);
#define MPI_UNPACK_INT(a, b, c)			MPI_Unpack(bufRecvMsg[(c)],		\
	sizeBufMsg,&posi,(a),(b), MPI_INT, MPI_COMM_WORLD);
// Convert between indices and numbers
#define V2IND_ASSIGN_INT(a,b,c,d)										\
	(b)=(a)/(d), (c)=(a)%(d)
#define V3IND_ASSIGN_DBL(a,b,c,d)										\
	(d)[0]=(double)(((a)/(b)[2])/(b)[1])*(c),							\
	(d)[1]=(double)(((a)/(b)[2])%(b)[1])*(c),							\
	(d)[2]=(double)((a)%(b)[2])*(c)							
#define V3IND_ASSIGN_INT(a,b,c,d)										\
	(d)[0]=(((a)/(b)[2])/(b)[1])*(c),									\
	(d)[1]=(((a)/(b)[2])%(b)[1])*(c),									\
	(d)[2]=((a)%(b)[2])*(c)							
#define V3IND_BACK_INT(a,b,c)											\
	(a)=((b)[0]*(c)[1]+(b)[1])*(c)[2]+(b)[2]
#define V3IND_ASSIGN_CONST_INT(a,b,c)									\
	(c)[0]=(((a)/b)/b),													\
	(c)[1]=(((a)/b)%b),													\
	(c)[2]=((a)%b)							
#define V3IND_BACK_CONST_INT(a,b,c)		(a)=((b)[0]*c+(b)[1])*c+(b)[2]
// Unit conversion between conventional units and simulation units
#define T_SEC2TS(a)			((int)((a) / dtReal))
#define T_S2SEC(a)			((a)*(CYL_DRAG(DIA_CYL_ACT, DIA_CYL_ACT		\
	* nActPerSeg)*SQR(L_SCALE_IN_M)/KT_IN_J))
#define T_SEC2S(a)			((a)/(CYL_DRAG(DIA_CYL_ACT, DIA_CYL_ACT 	\
	* nActPerSeg)*SQR(L_SCALE_IN_M)/KT_IN_J))
#define F_PN2S(a)			((a)/KT_IN_J*L_SCALE_IN_M/1.0e12)
#define F_S2PN(a)			((a)*KT_IN_J/L_SCALE_IN_M*1.0e12)
#define F_N2S(a)			((a)/KT_IN_J*L_SCALE_IN_M)
#define F_S2N(a)			((a)*KT_IN_J/L_SCALE_IN_M)
#define F_N2S(a)			((a)/KT_IN_J*L_SCALE_IN_M)
#define E_S2J(a)			((a)*KT_IN_J)
#define L_S2M(a)			((a)*L_SCALE_IN_M)
#define L_S2UM(a)			((a)*L_SCALE_IN_UM)
#define L_S2NM(a)			((a)*L_SCALE_IN_NM)
#define L_M2S(a)			((a)/L_SCALE_IN_M)
#define L_UM2S(a)			((a)/L_SCALE_IN_UM)
#define L_NM2S(a)			((a)/L_SCALE_IN_NM)
#define KS_S2NPM(a)			((a)*KT_IN_J/SQR(L_SCALE_IN_M))
#define KS_NPM2S(a)			((a)/KT_IN_J*SQR(L_SCALE_IN_M))
#define KB_S2NM(a)			((a)*KT_IN_J)
#define KB_NM2S(a)			((a)/KT_IN_J)
#define PA2S(a)				((a)/KT_IN_J*CUBE(L_SCALE_IN_M))
#define S2PA(a)				((a)*KT_IN_J/CUBE(L_SCALE_IN_M))
#define NPV2S(a)			((a)/KT_IN_J*pow(L_SCALE_IN_M, 4.))
#define S2NPV(a)			((a)*KT_IN_J/pow(L_SCALE_IN_M, 4.))
// Unit conversion of angles
#define DEG2RAD(a)			((a)*PI/180.)
#define RAD2DEG(a)			((a)/PI*180.)
// Probabilities and rates
#define K2P(a) 				(1.-exp(-1.*(a)*dtReal))
#define P2K(a) 				(log(1.-(a))/(-1.*dtReal))
#define FAC_P(a) 			(-1.*(a)/N_AVO/CUBE(L_SCALE_IN_M)*1.0e3*dtReal)
/*---------------------------- Vector calculation ----------------------------*/

// for 3-D vector
#define V3REVSIGN(a)		(a)[0]=-1.*(a)[0], (a)[1]=-1.*(a)[1],		\
	(a)[2]=-1.*(a)[2]
#define V3DOT(a,b)			((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define V3CROSS(a,b,c)		(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1], 		\
	(a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2], (a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0]
#define V3ADD(a,b,c)													\
	(a)[0]=(b)[0]+(c)[0],(a)[1]=(b)[1]+(c)[1],(a)[2]=(b)[2]+(c)[2]
#define V3SUB(a,b,c)													\
	(a)[0]=(b)[0]-(c)[0],(a)[1]=(b)[1]-(c)[1],(a)[2]=(b)[2]-(c)[2]
#define V3LEN_SQ(a)			V3DOT(a,a)
#define V3LEN(a)			sqrt(V3DOT(a,a))
#define V3SET(a,b,c,d)		(a)[0]=(b), (a)[1]=(c), (a)[2]=(d)
#define V3SET_ALL(a,b)		(a)[0]=(b), (a)[1]=(b), (a)[2]=(b)
#define V3ZERO(a)			V3SET_ALL(a,0)
#define V3MUL(a,b,c)		(a)[0]=(b)[0]*(c)[0], (a)[1]=(b)[1]*(c)[1], \
	(a)[2]=(b)[2]*(c)[2]
#define V3DIV(a,b,c)		(a)[0]=(b)[0]/(c)[0], (a)[1]=(b)[1]/(c)[1], \
	(a)[2]=(b)[2]/(c)[2]
#define V3DIV_INT(a,b,c)	(a)[0]=(int)((b)[0]/(c)[0]), 				\
							(a)[1]=(int)((b)[1]/(c)[1]),				\
							(a)[2]=(int)((b)[2]/(c)[2])
#define V3SCALE(a,b)		(a)[0]*=(b), (a)[1]*=(b), (a)[2]*=(b);
#define VV3ADD(a,b)			V3ADD(a,a,b)
#define VV3SUB(a,b)			V3SUB(a,a,b)
#define VV3DIV(a,b)			V3DIV(a,a,b)
#define V3PROD(a)			((a)[0]*(a)[1]*(a)[2])
#define V3SUM(a)			((a)[0]+(a)[1]+(a)[2])
#define VS3ADD(a,b,c,d)													\
	(a)[0]=(b)[0]+(c)[0]*(d),(a)[1]=(b)[1]+(c)[1]*(d),					\
	(a)[2]=(b)[2]+(c)[2]*(d)
#define VSV3ADD(a,b,c,d)												\
	(a)[0]=(b)[0]+(c)[0]*(d)[0],(a)[1]=(b)[1]+(c)[1]*(d)[1],			\
	(a)[2]=(b)[2]+(c)[2]*(d)[2]
#define V3AVG(a,b,c)													\
	V3ADD(a,b,c), V3SCALE(a,0.5)
#define VVS3ADD(a,b,c)		VS3ADD(a,a,b,c)
#define VSS3ADD(a,b,c,d,e)												\
	(a)[0]=(b)[0]*(d)+(c)[0]*(e),(a)[1]=(b)[1]*(d)+(c)[1]*(e),			\
	(a)[2]=(b)[2]*(d)+(c)[2]*(e)
#define VS3SUB(a,b,c,d)													\
	(a)[0]=(b)[0]-(c)[0]*(d),(a)[1]=(b)[1]-(c)[1]*(d),					\
	(a)[2]=(b)[2]-(c)[2]*(d)
#define VVS3SUB(a,b,c)		VS3SUB(a,a,b,c)
#define VVSS3SUB(a,b)		(a)[0]-=(b), (a)[1]-=(b), (a)[2]-=(b)
#define VVSS3ADD(a,b)		(a)[0]+=(b), (a)[1]+=(b), (a)[2]+=(b)
#define V3COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2]
#define VS3COPY(a,b,c)													\
	(a)[0]=(b)[0]*(c),(a)[1]=(b)[1]*(c),(a)[2]=(b)[2]*(c)
#define V3COS(a,b)			(V3DOT((a),(b))/(V3LEN(a)*V3LEN(b)))
#define V3ANG(a,b)			Acos(V3COS(a,b))
// for 2-D vector
#define V2SET(a,b,c)		(a)[0]=(b), (a)[1]=(c)
#define V2SET_ALL(a,b)		(a)[0]=(b), (a)[1]=(b)
#define V2ZERO(a)			V2SET_ALL(a,0)
#define V2COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1]
// for 4-D vector
#define V4SET(a,b,c,d,e)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d), (a)[3]=(e)
#define V4SET_ALL(a,b)		V4SET(a,b,b,b,b)
#define V4COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3]
// for 5-D vector
#define V5SET(a,b,c,d,e,f)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d),			\
	(a)[3]=(e), (a)[4]=(f)
#define V5SET_ALL(a,b)		V4SET(a,b,b,b,b,b)
#define V5COPY(a,b)			(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4]
#define V6SET(a,b,c,d,e,f,g)	(a)[0]=(b), (a)[1]=(c), (a)[2]=(d),			\
	(a)[3]=(e), (a)[4]=(f), (a)[5]=(g)
#define V6COPY(a,b)				(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4], (a)[5]=(b)[5] 
#define V7COPY(a,b)				(a)[0]=(b)[0],(a)[1]=(b)[1],(a)[2]=(b)[2],  \
	(a)[3]=(b)[3],(a)[4]=(b)[4], (a)[5]=(b)[5], (a)[6]=(b)[6]

/******************************************************************************/
/***************************** Function Prototypes ****************************/
/******************************************************************************/
/*-------------------------------- process.c ---------------------------------*/
// Main process 
void MainProcess(void);
// Initialize and excute single steps & post-process after single steps
void SingleProcess(void), TaskAfterSingleProcess(void);
void InitSingleProcess(void);
inline void InitSingleProcessSubroutine(double *arr, int cnt) ;
// Switch toggles 
void InitToggleParameters(int mode);
// Procedures necessary without pre-assembled network data
void PrepareStateWoNetworkData(void);
/*---------------------------------- init.c ----------------------------------*/
// Initialize the number of cells based on the given number of cores
void InitCellNumber(void);
// Initialize process
void InitRun(void), InitOutput(void);
// Define arrays and assign values for them
void AllocPreArrays(void), AllocArrays(void), AssignArrayValues(void);
// Assign initial values of variables
void AssignInitValues(void), CheckPeriodToggleParameter(FuncCont *);
// Initialize record files and seeds for generating random numbers
void InitRecFiles(void), InitRandSeed(void);
// Check whether 'condition' and 'parallel' exist
void InitFileCheck(void);
// Initialize the unbinding and walking rates based on the number of heads
void InitMotorWalkUnbindRates(int);
/*-------------------------------- calForce.c --------------------------------*/
// Brownian forces (thermal fluctuation)
void CalcBrownForces(void);
inline void CalcBrownForcesSubroutine(double *, double, int, int, double *);
// Repulsive forces between cylindrical segments
void CalcRepulsiveForces(void);
void CalcRepulsiveForcesSubroutine(int *, int *, double [][NDIM],
        double *[4], double *, int);
double CalcRepulsiveForcesSubSubroutine(double [][NDIM], double *, 
		double *, double, int);
void CalcRepulsiveForcesSubSubSubroutine(double [][NDIM], int, int);
// Spring forces
void CalcSpringForces(void), CalcSpringForcesSubroutine(int, int, int, int *);
// Bending forces and necessary subroutines
void CalcFilaBendForces(void), CalcFilaBendForcesSubroutine(int, int, int);
void CalcAbpBendForces(void), CalcAbpBendForcesSubroutine1(int, int);
void CalcAbpBendForcesSubroutine2(int, int *);
// For motor thick filaments
void CalcMotorBackboneForces(void);
// Tools for force calculation
void CalcCosine(double*,double *,double *,double *,double *,double *,double *);
double CalcBendForceSubroutine(double *, double *, double, double, 
		double *, double *);
double CalcBendForce(double *, double *, double *, double *, 
		double *, double, double, double *, double *);
void AddSpringForce(double, double, double *, double *, double *);
/*---------------------------------- update.c --------------------------------*/
// Update locations of particles reflecting forces
void UpdateNewLocation(void);
// Update neighboring list 
void UpdatePrevLocationForNL(void), MeasureDisplacementForNL(void);
void UpdateNeighborListSubroutine(double [][NDIM], double *, int);
int UpdateNeighborListSubroutine2(int, int *, int *);
int UpdateNeighborListSubroutine3(int, int, int *, int *);
void UpdateNeighborList(void);
void DeleteActinSegmentInNeighborList(int *);
void DeleteElementInNeighborList(int, int);
void InsertElementInNeighborList(int, int, int);
// Dynamic behaviors of actins
void UpdateActinNucleation(void), UpdateActinAssembly(void);
void UpdateActinDisassembly(void), UpdateActinDisassemblySubroutine(int, int);
void UpdateActinDisassemblySubroutine2(int, int);
void UpdateActinBurstDisassembly(void);
void UpdateActinSevering(void), UpdateActinSeveringSubroutine(int *);
void UpdateActinSeverAnnealSubroutine(int, int, ListInt *, int);
void UpdateActinAnnealing(void), UpdateActinAnnealingSubroutine(int *);
void UpdateActinCapUncap(void), UpdateActinDegradation(int *, double, int);
double CalcActinDegradationRate(int, int, int);
void UpdateNoActinDynamicsList(void);
int CheckActinAvailability(int, int);
// Dynamic behaviors (unbinding/binding/walking) of ACPs and motors
void UpdateAbpBinding(void), UpdateAbpBindingSubroutine(int *, int);
void UpdateActiveAbpUnbinding(void);
void UpdateActiveAbpUnbindingSubroutine(int, int);
void UpdateInactAbpUnbinding(void), UpdateInactAbpUnbindingSubroutine(int ,int);
void UpdateAbpMonomerBinding(void), UpdateMotorWalking(void);
void UpdateMotorAssembly(void), UpdateMotorTurnover(void);
void UpdateNoAbpUnbindList(void);
double UpdateMotorWalkingSubroutine(int, int);
double CalcUnbindingRate(int, int, int, int);
// Updating chains and information
void UpdateChainList(void), SortChainList(void);
/*------------------------------- dataManage.c -------------------------------*/
// Load and extract network configuration
void LoadConfig(int), ExtractConfig(void);
void ExtractConfigSubroutine(double *, double *, double, int, int *, int *);
void CheckAnswerYesNo(char *, const char *, int *);
// Load initial parameters
void LoadInitParameter(void);
void LoadInitParameterSubroutine(FILE *, const char *, int *);
void LoadInitParameterSubroutine2(char *, int, int *, int *);
// Eliminate inactive or free components
void AddFreeAbp(void), AddFreeActin(void);
void AddFreeActinAbpSubroutine(int, int *, int *);
void DeleteAbp(void), DeleteInactiveAbp(void), DeleteFreeActin(void);
// Sever filaments for each purpose
void SeverFilaOnPlane(int, int), SeverLongFilament(void);
void SeverFilaOnPlaneSubroutine(int, int);
// Break chains between actins and ABPs
void DetachAbpOnActin(int), DetachActinOnAbp(int);
// Inspect the integrity of chain data
void InspectChainListAndLength(void);
int InspectChainList(void), InspectChainLength(void);
// Pack information
int BinaryPackActinChainArray(int), BinaryPackAbpChainArray(int);
// Find the supportive framework of a network
void FindSupportiveFramework(void);
// Find the percolated structures of a network
void FindPercolatedFilaments(int, int);
void FindPercolatedFilamentsSubroutine(ListInt *, ListInt *, ListInt *, 
		ListInt2 *, int *, int *,int *, int, int, int *);
/*--------------------------------- record.c ---------------------------------*/
// Progress 
void RecordProgressSubroutine(int *, FILE *);
void RecordProgressSubroutine2(int *, int *, int, FILE *, int);
void RecordProgressSubroutine3(int *, int *, FILE *);
void RecordProgress(void);
// Record initial parameters
void RecordInitParameterSubroutine(const char *, int, int, int, FILE *); 
void RecordInitParameterSubroutine2(const char *, FuncCont *, FILE *);
void RecordInitParameter(void);
// Record counters
void RecordCounters(void), RecordCountersSubroutine(int, int, FILE *);
// Network configuration
void RecordConfig(char *, int), RecordConfigVmd(int), RecordConfigMatlab(void);
int RecordConfigVmdSubroutine(double *, double *, int *, int *);
int RecordConfigInformation(double *, double *, double *, int *, int *, int *,
		int);
void RecordConfigParaview(void);
// Record information related to the unbinding event of ABPs
void RecordAbpDynamics(int), RecordActinSeverEvent(void);
void RecordAbpUnbindEvent(int, int, int), RecordAbpBindEvent(int, int, int);
void RecordAbpTurnover(int, int, int, int);
// Motor thick filament
void RecordMotorFilamentSize(void), RecordMotorPosition();
// Morphological properties
void RecordFilamentLength(int), RecordConnectivity(int);
void RecordConnectivitySubroutine1(int *, int *, int *);
void RecordConnectivitySubroutine2(int **, int *, int *, int, int);
void RecordCrosslinkDistance(void), RecordCrosslinkAngle(void);
double RecordPoreSizeSubroutine1(double *, double *, double *, double *);
void RecordPoreSizeSubroutine2(ListDbl *, int *);
void RecordPoreSizeSubroutine3(double *, double *, int *, double *, int *);
void RecordPoreSize(void), FilamentEndPositionForceInBundle(void);
// Energy 
void RecordMechEnergy(int);
// Internal stress measurement
void RecordElasViscStressSubroutine1(double *, double *, int *, 
        double, int *, int *, int);
void RecordElasViscStressSubroutine2(double, double, int *, double *, double *, 
		double *, int *, double, int *, int);
void RecordElasViscStress(int, int *, int);
// Instantaneous information
void RecordInstantForces(void), RecordInstantChainLength(void);
void RecordActinAbpInstantOrient(int);
// Accumulated information
void RecordAccuLengthForces(void), ResetAccuLengthForces(void);
void RecordAccuLongSpringForces(int);
void RecordChainList(void);
// Information in unit of individual elements, filament segments, filaments
void RecordIndvSegFilaInformation(int);
/*--------------------------------- gatPrint.c -------------------------------*/
// Tools for recording
void Printf(const char *, ...), Printf0(const char *, ...);
void Fprintf1dFillerInt(FILE *, int, int, int);
void Fprintf1dArrayInt(FILE *, int *, int, int);
void Fprintf1dArrayIntWFil(FILE *, int *, int, int, int);
void Fprintf1dArrayIntWoRet(FILE *, int *, int);
void Fprintf1dFillerDouble(FILE *, double, int, int);
void Fprintf1dArrayDouble(FILE *, double *, int, int);
void Fprintf1dArrayDoubleWFil(FILE *, double *, int, int, int);
void Fprintf1dArrayDoubleWoRet(FILE *, double *, int);
void Fprintf2dArrayIntDouble(FILE *, double *, int, int, int *);
void Fprintf2dArrayInt(FILE *, int *, int, int, int, int);
void Fprintf2dArrayIntWFil(FILE *, int *, int, int, int, int, int);
void Fprintf2dArrayDouble(FILE *, double *, int, int, int, int);
void Fprintf2dArrayDoubleWFil(FILE *, double *, int, int, int, int, int);
void RecordChainArrayInt(FILE *, const char *, int *, int, int, int);
void Gather1dArrayIntWoInd(int, int *, int *, int *);
void RecordGather1dArrayIntWoInd(int, int *, int *, FILE *);
void RecordGather1dArrayIntWoIndWoCnt(int, int *, FILE *, const char *);
void Gather2dArrayInt(int, int *, int, int *, int *, int *);
void Gather2dArrayIntWoIndWoSort(int, int *, int *, int *, int);
void RecordGather2dArrayInt(int, int *, int, int *, int *, FILE *, int);
void RecordGather2dArrayIntWoCnt(int, int, int *, int *, FILE *, int);
void RecordChainArrayDouble(FILE *, const char *, double *, int, int, int);
void Gather1dArrayDoubleWoInd(int, int *, double *, double *);
void RecordGather1dArrayDoubleWoInd(int, int *, double *, FILE *);
void RecordGather1dArrayDoubleWoIndWoCnt(int, double *, FILE *, 
		const char *);
void Gather2dArrayDouble(int, int *, int, int *, double *, double *);
void Gather2dArrayDoubleWoIndWoSort(int, int *, double *, double *, int);
void RecordGather2dArrayDouble(int, int *, int, int *, double *, FILE *, int);
void RecordGather2dArrayDoubleWoCnt(int, int, int *, double *, FILE *, int);
void GatherActChainPosition(int *, int *, double *);
void GatherAbpChainPosition(int *, int *, double *);
void GatherMembChainPosition(int *, int *, double *);
/*---------------------------------- error.c ---------------------------------*/
// Print errors if any 
void RecordError(int);
void RecordErrorArrayInt(FILE *, const char *, int, int *, int);
void RecordErrorArrayDouble(FILE *, const char *, int, double *, int);
void RecordErrorSubroutine1(FILE *, int), RecordErrorSubroutine2(FILE *, int);
void RecordErrorSubroutine3(FILE *, int);
void RecordErrorSpringForce(int, int, double, double, int);
void RecordErrorRepForce(int *, int *, double (*)[NDIM], double, 
		double, double *, int);
void RecordErrorBendingForce(int, int, int, double, int);
void RecordErrorTotalForceSubroutine1(FILE *, int);
void RecordErrorTotalForceSubroutine2(FILE *, int);
void RecordErrorTotalForceSubroutine3(FILE *, int);
void RecordErrorTotalForce(int, int), RecordErrorElement(int, int);
int CheckLargeForce(double, int);
void CheckLargeTotalForceAll(int), CheckNanForce(int);
/*----------------------------------- rng.c ----------------------------------*/
// Initialize the generation of rundom numbers
void init_genrand(unsigned long);
// Generate a Gaussian random number
void genrand_gauss(double *, double *);
// Generate a double (uniform) random number
double genrand_real3(void);
/*--------------------------------- boundary.c -------------------------------*/
// Calculate forces
void CalcBoundRepulsiveForces(void);
// Update information
void UpdateBoundaryActinUnbindMature(void), UpdateBoundaryActinBinding(void);
void UpdateBoundaryLocation(void);
// Handle boundary conditions
void CheckCrossBound(int *, double *, double *);
void ApplyBoundCondVector(double *, int, int);
void ApplyBoundCondVecDiff(double *), ApplyBoundCondAll(void);
void ConvertRectDomainVector(double *, int);
// Interactions between boundaries and others
int CheckParticleInDomain(double *), CheckActinAbpOverlapBoundary(double *);
double HowManyInOutBound(void);
// Record information
void RecordBoundaryLocation(void), RecordBoundaryActinUnbindBind(void);
void RecordBoundaryActinMature(void);
/*--------------------------------- paraProc.c -------------------------------*/
// Handle the list of long chains between actins and ABPs
void UpdateLongChainNormal(void), UpdateLongChainPlympton(void);
void UpdateLongChainNormalSubroutine(int, int, double);
void UpdateLongChainNormalSubSubroutine1(double *, ListInt *, double);
void DeleteLongChain(int, int);
int InsertLongChain(int, int, double);
// Adjust subdomain size to maintain loads in CPUs at a similar level
void UpdateSubdomSectionLocation(int), RecordSubdomSectionLocation(void);
// Move particles between subdomains
void MoveParticles(void), MoveParticlesSubroutine2(int, int);
void MoveParticlesSubroutine1(ListInt *, double *, int);
void MoveParticlesSubroutine3(ListInt *, double *, int, int);
// Copy particles between subdomains
void CopyParticles(void), CopyParticlesSubSubroutine(int, int *, int);
void CopyParticlesNormalSubroutine(int, int *, int *, int);
void CopyParticlesPlymptonSubroutine(int, int, int *, int);
int CopyParticlesPlymptonSubroutine2(int);
// Process the messages regarding the dynamic events of actins and ABPs
// which will be transferred between subdomains
void UpdateActinAbpDynamicsEvents(int); 
void UpdateActinAbpDynamicsEventsSubroutine(ListInt *, int);
void UpdateActinAbpDynamicsEventsSubSubroutine(ListInt *, int, int, int, 
		ListInt *, int *);
void UpdateActinSeverEvents(void); 
void UpdateAbpUnbRebLists(int, int, int, int);
// Manage actin and ABP monomers
void UpdateActinAbpMonomerList(void);
void UpdateActinAbpMonomerListSubroutine(int, int, int);
void UpdateActinAbpMonomerListSubroutine2(int, int, int);
void UpdateMotorNucleCounter(void);
// Collect arrays from adjacent subdomains
void CollectArrayIntFromAdjacentSubdomain(ListInt *, int);
void CollectArrayDblFromSubdomainList(double *, double *, int, ListInt *, int);
/*--------------------------------- membrane.c -------------------------------*/
// Calculate force 
void CalcMembraneRepulsiveForces(void), CalcMembraneSpringForces(void);
void Calc2dMembraneBendForces(void), Calc3dMembraneBendForces(void);
void CalcMembraneVolumeForces(void), CalcMembraneAreaForces(void);
void CalcMembraneProtrusiveForces(void), CalcMembraneSlideForces(void);
void CalcMembraneSubstrateProtrusiveForces(void);
// Update information
void UpdateMembraneFixation(void), UpdateMembraneRelease(void);
void UpdateMembraneBinding(void), UpdateMembraneUnbindMature(void);
void UpdateMembraneBindingSubroutine(int, int, int, int);
int UpdateMembraneUnbindMatureSubroutine(int, int, int, int);
void UpdateMembraneCenterPosition(void);
void UpdateMembraneProtrusion(void), UpdateMembraneProteinRecover(void);
void UpdateMembraneVolume(void), UpdateMembraneTotalArea(void);
void UpdateMembraneUnitArea(void), UpdateMembranePressure(void);
void UpdateMembraneSlideList(void), UpdateMembraneSlidePoint(void);
// Record information
void RecordMembraneDimension(int), RecordMembraneCenter(int);
void RecordMembraneAdhesion(void);
// Initialization
void InitMembraneInformation(void), Init3dMembraneInformation(void);
// Process the messages regarding the dynamic events of membranes
// which will be transferred between subdomains
void UpdateMembraneDynamicsEvents(); 
void UpdateMembraneDynamicsEventsSubroutine(ListInt *);
void UpdateNoMembraneDynamicsList(void);
// Interactions between networks and membrane
void TrimNetworkForMembrane(void);
int TrimNetworkForMembraneSubroutine(double);
int CheckActinAbpOverlapMembrane(double *, double *, int);
double HowManyInOutMembrane(double *);
// Check error
void CheckMembraneError(int);
/*----------------------------------- bead.c ---------------------------------*/
// Calculate force
void CalcBeadRepulsiveForces(void), CalcBeadActinBindingForces(void);
double CalcBeadRepulsiveForcesSubroutine(double, double);
// Update information
void UpdateBeadNewLocation(void);
void ApplyBeadSinuSignal(void), UpdateBeadActinBinding(void);
//Record information
void RecordBeadForceLocation(void), RecordBeadActinBinding(void);
// Initialization
void InitBeadInformation(void);
int CheckActinAbpOverlapBead(double *, double *, int);
/*--------------------------------- rheology.c -------------------------------*/
// Procedures before bulk rheology measurements
void PrepareBulkRheology(void), PrepareNotBulkRheology(void);
// Apply stress and strain
void ApplyStressStrain(void); 
// Choose the particles which will be traced for segment-rheology
int ChooseTrajectoryListSubroutine(int, int);
void ChooseTrajectoryList(void);
void ReplaceElementInTrajList(int, int);
// Rheological data
void RecordStress(int), RecordTrajectory(int);
void RecordTrajectorySubroutine(double *, int);
/*---------------------------------- tools.c ---------------------------------*/
// Find maximum and minimum values of values in an array
int FindMaxMin1dArrayInt(int *, int, int);
double FindMaxMin1dArrayDbl(double *, int, int);
// Set values for an array
void SetAllValue1dArrayInt(int *, int, int);
void SetAllValue1dArrayDouble(double *, int, double);
// Copy values in an array to other array
void Copy1dArrayInt(int *, int *, int);
void Copy1dArrayDouble(double *, double *, int);
// Find, insert, or delete elements in an array
int FindElementArray(int *, int, int, int, int);
int Find2ElementArray(int *, int, int, int, int, int);
int Find3ElementArray(int *, int, int, int, int, int, int);
void InsertElement1dArrayWoChk(int *, int *, int);
int InsertElement1dArrayWChk(int *, int *, int);
void InsertElementArrayByIndex(int *, int *, int *, int, int);
void DeleteElement1dArray(int *, int *, int);
void DeleteElementArrayByIndex(int *, int *, int, int);
int DeleteTwoElementArray(int *, int *, int, int, int, int);
// Sum or average values in array
double AvgArrDbl(double *, int), AvgArrInt(int *, int);
double SumArrDbl(double *, int);
int SumArrInt(int *, int);
// Check the size of array
void CheckArraySize(ListInt *, int *, int, int);
void Check2dArraySize(ListInt2 *, int *, int, int);
// Vector and geometry calculation (general)
int CheckTriaPntDist(double [][NDIM], double *, double);
int CheckTriaSegDist(double [][NDIM], double [][NDIM], double);
int CheckTriaTriaDist(double [][NDIM], double [][NDIM], double);
double CalcTriaPntDist(double [][NDIM], double *, double *, double *);
double CalcTriaSegDist(double [][NDIM], double [][NDIM], double *, double *);
double CalcTriaTriaDist(double [][NDIM], double [][NDIM], double *, double *);
double CalcSegPntDist(double [][NDIM], double *, double *, double *);
double CalcSegSegDist(double [][NDIM], double [][NDIM], double *, 
		double *, int);
double CalcTriaCentroidRadius(double [][NDIM], double *);
void CalcTriaCircumCenter(double [][NDIM], double *);
double CalcTriaArea(double *, double *, double *);
double CalcTetraVol(double *, double *, double *, double *);
void DistributeForceOnTria(double (*)[NDIM], double **, double *, double *);
void NormVec(double *), CalcVec(double *, double *, double *);
double CalcUnitVec(double *, double *, double *);
double CalcDist(double *, double *, int);
double CalcVecDist(double *, double *, double *, int);
// Vector calculation for Specific for elements
double CalcVecActinAbp(double *, int, int, int);
void CalcUnitVecActinAbp(double *, int, int, int);
double CalcDistActinAbp(int, int, int);
double CalcVecDistActinAbp(double *, int, int, int);
void CalcPosOnActSeg(double *, double *, double *, double);
void CalcPosOnActSegSide(double *, double *, double *, int);
void OffsetSegEndPnt(double *, double *, double *);
void CalcInactAbpPosition(double *, double *, double *, double *, int, int);
void CalcAbpArmEndPos(double *, int, int);
void CalcMembUnitNormalDirec(int *, double *);
// Find information related to actin and ABP
int FindAbpActinChain(int, int, int), HowManyAbpActinChain(int, int);
// Generate random information
void GenRandDirecVec(double *), GenRandPosSubdom(double *);
int GenRandIntIndex(int);
// Calculate the rank or iCell of particles depending on location
void CalcIndMolecule(double *, int *);
int CalcRankMolecule(double *), CalcRankIndMolecule(double *, int *);
// Adjustment for rates of dynamic behaviors of actins and ABPs
double AdjustDynamicsRate(double);
// Matrix calculation
void MultiplyMatrixIntSerial(int *, int *, int *, int, int, int);
void MultiplySqMatrixIntParallel(int *, int *, int *, int);
// Misc.
char *GenFileName(const char *);
void SwapInt(int *, int *), SwapDbl(double *, double *);
int TrimIntVal(int, int, int), SetKind(int, int, int), SignDbl(double);
double TrimDblVal(double, double, double), Acos(double), Sqrt(double);
int CompInt(const void *, const void *);
int CompDbl(const void *, const void *);
char* IntToStr(int, int);

void UpdateCenterSever(void);

/******************************************************************************/
/****************************** Global Variables ******************************//******************************************************************************/
 
// Variables related to time steps
extern long long currTimeStep;
extern double dt, dtReal;
extern time_t initTime;
extern Duration netForm, rheo, motActiv;
// Given concentration
extern int nAbpGoal, nAbpDet[NK_ABP], nAbpMall[NK_ABP], nAbpGoalDet[NK_ABP];
extern double cAct, RAbp[NK_ABP];
// Number of particles
extern int nAct, nAbp, nMot; 
extern int nActMe, nActCp, nActMeCp, nActC;
extern int nAbpMe, nAbpCp, nAbpMeCp, nAbpC;
extern int nAcpInaMe, nMotInaMe, nAcpMme, nMotMme;
extern int nActGoal, nActMall, nActFilaMe;
extern int nActMin, nAbpMin, nMbMin;
// Chain, position, forces, and lists of particles
extern ListInt actM, acpM, motM, iFilaP;
extern MoleAct act;
extern MoleAbp abp;
extern ActF actF;
extern AbpF abpF;
extern MotSelfAss motSA;
extern int *chAct, *chAbp, cntAbpDC, tglNeiAbpSC, tglNeiAbpDC;
extern int nChAc, nChAcX, nChAcY, nChAb, nActPerSeg, confNChAcX, confNChAcY;
extern double *rAct, *rAbp;
// Update neighboring list
extern ListInt neigh;
extern Cell cell;
extern double dispHeuSq, maxNeiLenSq;
// Arrays related to lists
extern ListInt sendAbpDyn, noAbpDyn;
extern ListInt sendActDyn, noActDyn;
extern ListInt sendMbDyn, noMbDyn;
extern int *fixAct, *abpMotId;
// Domain (width, periodic boundary conditions, or repulsive force)
extern int pbc[NDIM], neiPbc[NDIM], neiPbcMb[NDIM], confPbc[NDIM], dir2D;
extern double dimDom[NDIM], dimDomH[NDIM], minDimDomC;
extern Boundary bnd;
// For boundaries
extern FuncCont recBndLoc, updSubdSize, recBndUnbReb, recBndActMat;
extern FuncCont recBndTracF;
extern double *rGridInit, volMe;
extern Mature bndMat;
extern BndReb bndReb;
extern BndUnb bndUnb;
extern Force bndVol;
extern BndMv bndMv;

/*-------------------------- For parallel processing -------------------------*/
extern FuncCont updSubdSize;
extern int mpiMethod, nCpu;
extern int *iAct, *iAbp;
// Ranks of CPUs
extern int rank, *adjRank, *iRank, *cntAdjRank;
// For boundaries and indicies of subdomains
extern int nCell[NDIM], iCell[NDIM], nGrid[NDIM];
extern double edge[NDIM*2], **rGrid, neiEdge;
// Used in CopyParticles and MoveParticles
extern int *cntCpPar, *cntMvPar, modeActCh;
extern ListInt *cpPar, *mvPar, insNeiPar; //, mvParAll;
// Varilables related to messages
extern int sizeBufMsg, *mpiTestSendFlag, *mpiTestRecvFlag;
extern char **bufSendMsg, **bufRecvMsg;
extern MPI_Status status;
extern MPI_Request *sReq, *rReq;
// Variables related to longCh
extern ListInt longCh, longChExtMsg, longChIntMsg;
extern double *longChDist, maxDisp, maxActCh;
/*----------------------------------------------------------------------------*/

/*---------------- Dynamic behaviors of actin, ACP, and motor ----------------*/
extern FuncCont updMono;
// Dynamics of actin
extern AssDisCapUncAge actAss, actDis, actCap, actUnc, actAge;
extern Burst actBst;
extern Sever actSev;
extern Anneal actAnn;
extern Nucle actNuc;
extern Degrade actDgd;
extern int tglActMoDyn, tglActDyn, tglActFormDyn;
extern int tglActDisBstSev, tglActDisBstSevAbp;
extern int gTglActDynPres, gTglActDynNF;
extern int gTglActSevCap, gTglActSevBst, gTglActCapAll;
extern int gTglActTherm, gTglAcpTherm, gTglMotTherm;
extern int prdUpdActM, durNoActDyn, cntDisFila;
// Dynamics of ACPs and motors
extern AbpDyn abpDyn;
extern Reb acpReb, motReb;
extern Unb acpUnb; 
extern MotUnbWalk motUnb, motWalk; 
extern InaUnbMbind acpInaUnb, motInaUnb, acpMoBind, motMoBind;
extern FuncCont updMotN, abpAge;
extern int tglAbpAcInaDyn, tglAbpInaMoDyn, tglAcpAcInaDyn, tglMotAcInaDyn;
extern int gTglAcpDynPres, gTglMotUnbRebPres, gTglMotWalkPres;
extern int gTglAcpDynNF, gTglMotUnbRebNF, gTglMotWalkNF, gTglMotWalkSld;
extern int durNoAcpUnbReb, durNoMotUnbReb, durNoMotWalk;
extern int gTglImpAcpM, gTglImpMotM;
extern ListDbl unbLog, bindLog, toLog, sevLog;
extern MotMechChem motMC;
// Dynamics of membrane
extern int tglMbDyn, tglMbNucDyn, durNoMbDyn, cntNucAssAnn;
extern int tglMbNucLocAre; 
/*----------------------------------------------------------------------------*/
 
/*----------------------- For rheological measurements -----------------------*/
// For both ways
extern FuncCont recStre, recProg, recConf;
extern int rheoWay;
// Segment rheology
extern RecTraj recTraj;
// Bulk rheology
extern ListInt meaStrePar, appStraPar, meaStreParMe, appStraParMe, rankMeaStre;
extern Strain stra;
extern Stress stre;
extern PreStraStre pres;
extern SinuStraStre sinuStr;
extern int dirStr, dirNoPBC, dirOther, prdUpdSinuStre;
extern int bulkRheoWay, bulkRheoType, signStr, gotMeaStre;
/*----------------------------------------------------------------------------*/

/*----------------------------- For microrheology ----------------------------*/
extern int gTglBead;
extern Bead bead;
extern BeadBind beadBind;
/*----------------------------------------------------------------------------*/

/*---------------------------- For data recording ----------------------------*/
// For general records
extern char fnOut[80], dataFold[80];
// Record the accumulated chain lengths and forces
extern RecLenForce recAct, recAbp, recMb;
// Record configuration for VMD
extern RecConfVmd recConfVmd;
extern FuncCont recConfPV;
// Record configuration for Matlab
extern FuncCont recConfMlb;
// Record internal elastic/viscous stresses
extern FuncCont recSecStre;
extern int nSecStreDiv[NDIM];
extern double *viscStre, *elasStre;
// Record the percolation of networks
extern FuncCont recPerc;
extern int *iPerc, dirRecPerc;
// Find supportive framework
extern FuncCont findSupp;
extern int kindFindSupp, *iSupp;
extern double porFindSupp;
// Record the turnover of ACPs or motors 
extern int tglRecAbpTurn;
extern double *abpTurn;
// Record (longitudinal) forces of ACPs or motors
extern FuncCont recLongF;
extern double *recLongSprFabp, *recInstSprFabp;
// Record etc
extern FuncCont confVmdInfo, recFilaL, recE, recInfo;
extern FuncCont recMotSize, recMotPos, recCrsDist, recAbpDyn, recConn;
extern FuncCont recPoreSize, recAbpUnb, recAbpBind, recAbpTurn, recActSev;
extern FuncCont recMbCen, recMbDim, recBeadLoc;
/*----------------------------------------------------------------------------*/

/*--------------------------- Membrane & nucleus------------------------------*/
extern int gTglLoadMbNucData, dirNormMbNuc, dimMbNuc;
extern int nMb, nMbMe, nMbCp, nMbC, nUnitMb, nChMb, nObjMbNuc;
extern int *chMb, *iMb, *unitMb, *rebMb, *actRebMb, dirMbNuc[NDIM], *idxMb;
extern double *rMb, *nDirMb, *eqLenMb, *cenMb, *nDirUnitMb, *eqArUnitMb;
extern double stfRepMbNuc;
extern Cell cellMb;
extern ListInt neighMb;
/*---------------------------------- Membrane --------------------------------*/
extern int gTglMb, gTglMbTherm, gTglMbActNuc, gTglMbCont;
extern int nMbAct, nMbPerObj;
extern int levMb, sideMb, maxNFAMb;
extern double radMb, radMbInit, thkMbActNuc;
extern MembNuc memb;
extern Unb mbUnb, mbDef;
extern Mature mbMat;
extern Reb mbReb, mbFix;
extern MembPro mbPro;
extern FuncCont mbVol, mbAre;
extern MembSld mbSld;
/*----------------------------------- Nucleus --------------------------------*/
extern int gTglNuc, gTglNucTherm, gTglNucActNuc;
extern int nNucAct, nNucPerObj;
extern int levNuc;
extern double radNuc, radNucInit;
extern Unb nucUnb;
extern MembNucReb nucReb;
extern FuncCont nucVol, nucAre;
/*----------------------------------------------------------------------------*/

// Check errors
extern int stopSig;
extern double magUnstF;
// Misc.
extern int gTglLoadNetData, gTglLoadNetDataFix;
extern int seed, *allIntL;
extern double *arrAcos, *allDblL;

extern int typeBias, tglBias;
extern double thkBias, degBias;
