// This file contains functions related to initializations.

#include <float.h>

void InitCellNumber(void) {
  int ind;
  FILE *fIn;

  ind = -1;
  fIn = fopen("parallel", "r");
  while((fscanf(fIn, "%d\t%d\t%d\t%d\n", &ind, &nCell[0], &nCell[1], &nCell[2]))
		== 4) { 
	BREAK(ind == nCpu);
  }
  if (ind == -1) {
	Printf0("Error: the given number of cores doesn't exist "
			"in 'parallel'!\n");
	exit(-1);
  }
  fclose(fIn);
  if (nCpu != V3PROD(nCell)) {
	Printf0("Error: the given number of cores doesn't match the product of "
			"the cell numbers!\n");
	exit(-1);
  }
  if (gTglMb != 0 && dimMbNuc == 2) {
	if (nCell[dirNormMbNuc] > 1) {
		Printf0("Error: in a current version, the number of cell in the "
				"direction of axis of 2D membrane should be one!\n");
		exit(-1);
	}
 }
}

void InitRun(void) {
  // Put a random seed into the engine of random numbers
  InitRandSeed(); 
  AllocPreArrays();
  // Assign initial values for parameters
  AssignInitValues();
  // Allocate arrays
  AllocArrays();  
  // Assign initial values in the arrays
  AssignArrayValues();
}

void InitOutput(void) { 
  Printf0("\n===================================== Start "
		"====================================\n");
  Printf0("currTimeStep ");
  if (nAct > 0) {
	if (actAss.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ Printf0("nActM/nActA/"); }
	Printf0(" nAct ");
	if (actNuc.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ Printf0("nFila "); }
	if (actNuc.gTgl != 0) { Printf0("cAcNu "); }
	if (actAss.gTgl != 0) { Printf0("cAcAs "); }
	if (actDis.gTgl != 0 || actBst.gTgl != 0) { Printf0("cAcDs "); }
	if (actBst.gTgl != 0) { Printf0("cAcBD "); }
	if (actSev.gTgl != 0) { Printf0("cAcSv "); }
	if (actAnn.gTgl != 0) { Printf0("cAcAn "); }
	if (actCap.gTgl != 0) { Printf0("cAcCp "); }
	if (actUnc.gTgl != 0) { Printf0("cAcUc "); }
  }
  if (nAbp - nMot > 0) { 
	if (acpInaUnb.gTgl != 0 || acpMoBind.gTgl != 0) { Printf0("nAcpM/"); }
    if (acpUnb.gTgl != 0 || acpReb.gTgl != 0) { Printf0("nAcpI/"); }
	Printf0("nAcpA/ nAcp "); 
	if (acpInaUnb.gTgl != 0) { Printf0("cApIU "); }
	if (acpMoBind.gTgl != 0) { Printf0("cApMB "); }
	if (acpUnb.gTgl != 0) { Printf0("cApUb "); }
	if (acpReb.gTgl != 0) { Printf0("cApRb "); }
  }
  if (nMot > 0) { 
	if (motSA.gTgl == 0) {
		if (motInaUnb.gTgl != 0 || motMoBind.gTgl != 0) { Printf0("nMotM/"); }
	    if (motUnb.gTgl != 0 || motReb.gTgl != 0) { Printf0("nMotI/"); }
		Printf0("nMotA/ nMot "); 
	}
	else {
		if (motInaUnb.gTgl != 0 || motMoBind.gTgl != 0) { Printf0("nMoHM/"); }
	    if (motUnb.gTgl != 0 || motReb.gTgl != 0) { Printf0("nMoHU/"); }
		Printf0("nMoHB/ nMoH "); 
	}
	if (motSA.gTgl == 0) {
		if (motInaUnb.gTgl != 0) { Printf0("cMoIU "); }
		if (motMoBind.gTgl != 0) { Printf0("cMoMB "); }
	}
	if (motUnb.gTgl != 0) { Printf0("cMoUb "); }
	if (motReb.gTgl != 0) { Printf0("cMoRb "); }
	if (motWalk.gTgl != 0) { Printf0("cMoWk "); }
    if (motSA.gTgl != 0) {
		Printf0("cMoNu cMoAs ");
			if (motSA.to.gTgl != 0) {
				Printf0("cMoTn ");
			}
		}
	}

  if (gTglMb != 0) {
	Printf(" nMb  ");
	if (mbSld.act.gTgl != 0) { Printf0("nMbSA "); }
	if (mbSld.abp.gTgl != 0) { Printf0("nMbSB "); }
	if (mbUnb.gTgl != 0) { Printf0("cMbUb "); }
	if (mbReb.gTgl != 0) { Printf0("cMbRb "); }
	if (mbDef.gTgl != 0) { Printf0("cMbDf "); }
	if (mbPro.gTgl != 0) { Printf0("cMbPr "); }
	if (mbFix.gTgl != 0) { Printf0("cMbFx "); }
	if (mbMat.gTgl != 0) { Printf0("cMbMa cMbMu "); 
	}
  }
  if (gTglNuc != 0) {
	if (nucUnb.gTgl != 0) { Printf0("cNuUb "); }
	if (nucReb.gTgl != 0) { Printf0("cNuRb "); }
  }
  Printf0("\n");
}

// Define dynamic arrays
void AllocPreArrays(void) {
  MALLOC(bnd.r,double,NDIM*2);
  MALLOC(abpF.dtRepAct,Double,NK_ABP);
  MALLOC(abpF.dtSprAct,Double,NK_ABP);
  MALLOC(abpF.len,Double,NK_ABP);
  MALLOC(abpF.drag,Double,NK_ABP);
  MALLOC(abpF.dia,double,NK_ABP);
  MALLOC(abpF.repDia,double,NK_ABP);
}

// Allocate arrays
void AllocArrays (void) {
  int n, sizeMpiArr, sizeArr;

  sizeArr = (nActGoal + nAbpGoal) * 2;
  sizeArr += nMb * ((nChMb - nMbAct) > nMbAct ? nChMb - nMbAct : nMbAct);
  MALLOC(allIntL,int,sizeArr);
  MALLOC(allDblL,double,nActGoal + nAbpGoal 
		+ ((nMb > nUnitMb) ? nMb : nUnitMb));

  sizeArr = nActC + nAbpC;
  if (gTglMb != 0) { sizeArr += nMbC; }

  neigh.siz = 2 * (nActC + nAbpC) * 10;
  if (gTglMb != 0) { neigh.siz += 2 * nMbC * 10 * 10; }
  // Neighboring list
  MALLOC(neigh.l,int,neigh.siz);

  if (gTglMb != 0) {
	// Arrays for membrane
	MALLOC(rMb,double,NDIM*nMb);
	MALLOC(nDirMb,double,NDIM*nMb);
	MALLOC(chMb,int,nChMb*nMb);
	MALLOC(iMb,int,nMb);
	MALLOC(idxMb,int,nMb);
	MALLOC(unitMb,int,nUnitMb*dimMbNuc);
	MALLOC(nDirUnitMb,double,nUnitMb*NDIM);
	MALLOC(memb.id,int,nMbC);
	MALLOC(memb.idx,int,nMbC);
	MALLOC(memb.r,double,NDIM*nMbC);
	MALLOC(memb.rPrev,double,NDIM*nMbC);
	MALLOC(memb.nDir,double,NDIM*nMbC);
	MALLOC(memb.f,double,NDIM*nMbC);
	MALLOC(memb.fBr,double,NDIM*nMbC);
	MALLOC(memb.ch,int,nChMb*nMbC);
	MALLOC(memb.areaL,double,nMbC);
	neighMb.siz = 2 * sizeArr * 10;
	MALLOC(neighMb.l,int,neighMb.siz);
	MALLOC(memb.unit.l,int,nUnitMb*dimMbNuc);
	MALLOC(memb.unitNDir,double,nUnitMb*NDIM);
	MALLOC(memb.unitCp.l,int,nUnitMb*dimMbNuc);
	MALLOC(cenMb,double,nObjMbNuc*NDIM);
	MALLOC(memb.vol.val,double,nObjMbNuc);
	MALLOC(memb.area.val,double,nObjMbNuc);
if (tglBias != 0 && typeBias == 2) {
	MALLOC(memb.bias,int,nMbC);
}
	if (dimMbNuc == 3) { 
		MALLOC(eqLenMb,double,(nChMb-nMbAct)*nMb);
		MALLOC(memb.eqLen,double,(nChMb-nMbAct)*nMbC);
		if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
			MALLOC(eqArUnitMb,double,nUnitMb);
			MALLOC(memb.unitEqAr,double,nUnitMb);
		}
	}
	if (mbReb.gTgl != 0 || mbUnb.gTgl != 0 
			|| nucReb.gTgl != 0 || nucUnb.gTgl != 0) {
		MALLOC(rebMb,int,nMb);
		MALLOC(memb.reb,int,nMbC);
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		MALLOC(memb.fix,int,nMbC);
		MALLOC(memb.rFix,double,NDIM*nMbC);
	}
	if (mbPro.gTgl != 0) {
		MALLOC(memb.pro,int,nMbC*2);
		if (gTglMbCont != 0) {
			MALLOC(memb.myo.conc,double,nMbC);
			MALLOC(memb.cyto.conc,double,nMbC);
		}
	}
	if (mbMat.gTgl != 0) {
		MALLOC(memb.nFA,int,nMbC * nMbAct);
	}
	if (gTglLoadNetData != 0) {
		MALLOC(actRebMb,int,nAct);
	}
	MALLOC(act.mbReb,int,nActC);
  }
  // Cylindrical segments of actin
  MALLOC(act.cyl.l,int,nActC*2);
  // Forces and positions of actin
  MALLOC(act.r,double,NDIM*nActC);
  MALLOC(act.rPrev,double,NDIM*nActC);
  MALLOC(act.f,double,NDIM*nActC);
  MALLOC(act.fBr,double,NDIM*nActC);
  MALLOC(act.maxFbr,double,1);
  if (gTglMb != 0) {
	MALLOC(act.distMb,double,nActC);
	MALLOC(abp.distMb,double,nAbpC);
  }
  // Forces and positions of ACPs and motors
  MALLOC(abp.r,double,NDIM*nAbpC);
  MALLOC(abp.f,double,NDIM*nAbpC);
  MALLOC(abp.fBr,double,NDIM*nAbpC);
  MALLOC(abp.maxFbr,double,NK_ABP);
  // Chains of actin, ACP, and motor
  MALLOC(act.ch,int,nActC*nChAc);
  MALLOC(abp.ch,int,nAbpC*nChAb);
  MALLOC(actM.l,int,nActC);
  if (gTglImpAcpM != 0) { MALLOC(acpM.l,int,nAbpC); }
  if (gTglImpMotM != 0) { MALLOC(motM.l,int,nAbpC); }
  // Global index of actin, ACP, and motor
  MALLOC(act.id,int,nActC);
  MALLOC(abp.id,int,nAbpC);
  // Local index of actin, ACP, and motor
  MALLOC(iAct, int, nActGoal);
  MALLOC(iAbp, int, nAbpGoal);
  // List of actin segments clamped on boundaries
  if (gTglLoadNetData != 0) {
	MALLOC(fixAct,int,nAct);
  }
  MALLOC(act.fix,int,nActC);
  // Local index of actin filaments
  MALLOC(act.iF,int,nActC);
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) {
	MALLOC(iFilaP.l,int,nActC);
  }
  if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
    MALLOC(act.nFA,int,nActC);
  }
  if (bndReb.gTgl != 0 && bndReb.gTglSpr == 1) {
	MALLOC(act.bndRFix,double,nActC*NDIM);
  }
  if (gTglBead != 0 && beadBind.gTgl != 0) {
	MALLOC(act.bdFix,int,nActC);	
	MALLOC(act.bdRFix,double,nActC*NDIM);
  }

  // Chains and positions of actin, ACP, and motor, which are used only when
  // a network is initially loaded and when information about chains and 
  // positions are collected from all subdomains.
  if (gTglLoadNetData != 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
  }

  // Sectioning lines
  MALLOC2(rGrid,double,NDIM);
  FOR_NDIM(n) {
	MALLOC(rGrid[n],double,nGrid[n]);
  }

  /*----------------------- Dynamic behaviors of actins ----------------------*/
  if (actSev.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0 
		|| actAnn.gTgl != 0) {
	sendActDyn.siz = 3 * nActC;
	MALLOC(sendActDyn.l,int,sendActDyn.siz);
  }
  if (actSev.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0 
		|| actNuc.gTgl != 0 || actAss.gTgl != 0 || actAnn.gTgl != 0) {
	noActDyn.siz = 2 * nActC;
	MALLOC(noActDyn.l,int,noActDyn.siz);
  }
  if (actSev.gTgl != 0) {
	MALLOC(actSev.p,double,DEG_ARR_ACTSEV + 1);
	MALLOC(actSev.pWA,double,DEG_ARR_ACTSEV + 1);
  }
  if (gTglActCapAll != 0) {
	MALLOC(act.cap,int,nActC);
  }
  if (actAge.gTgl != 0) {
	MALLOC(act.age,int,nActC);
  }
  MALLOC(act.len,int,nActC);
  if (actBst.gTgl != 0) {
	MALLOC(actBst.fil.l,int,nActC*3);
  }
  if (actDgd.gTgl != 0) {
	MALLOC(actDgd.l,int,nActC);
	MALLOC(actDgd.l2,double,nActC*3);
  }
  if (bndMat.gTgl != 0 || bndUnb.gTgl != 0) {
	MALLOC(bndMat.p,double,bndMat.maxF);
	MALLOC(bndMat.cntMe2,int,12);
  }
  /*--------------------------------------------------------------------------*/

  /*------------ Dynamic behaviors of ACP, motor, and boundaries -------------*/
  abpAge.gTgl = (motSA.to.gTgl != 0) ? 1 : 0;
  if (abpAge.gTgl != 0) {
	MALLOC(abp.age,int,nAbpC);
  }
  if (recAbpDyn.tgl != 0) {
	MALLOC(abpDyn.cntU,int,nAbpC);
	MALLOC(abpDyn.cntB,int,nAbpC);
	MALLOC(abpDyn.cntW,int,nAbpC);
  }
  // Unbinding of ACP
  if (acpUnb.gTgl != 0 || acpInaUnb.gTgl != 0) {
	MALLOC(acpUnb.p,double,acpUnb.maxF);
  }
  if (motSA.to.gTgl != 0 && motSA.to.mode == 1) {
	MALLOC(motSA.to.p,double,motSA.to.maxF);
  }
  // Arrays related to unbinding of ACP and motor
  if (acpUnb.gTgl != 0 || acpReb.gTgl != 0 || motUnb.gTgl != 0 
		|| motReb.gTgl != 0 || motWalk.gTgl != 0 
		|| acpInaUnb.gTgl != 0 || motInaUnb.gTgl != 0 
		|| acpMoBind.gTgl != 0 || motMoBind.gTgl != 0
        || (actDis.facKWA > 0. && actDis.gTgl != 0)
        || (actSev.facKWA > 0. && actSev.gTgl != 0) || actBst.gTgl != 0) {
	sendAbpDyn.siz = 3 * nAbpC;
	MALLOC(sendAbpDyn.l,int,sendAbpDyn.siz);
	noAbpDyn.siz = 2 * 3 * nAbpC;
	MALLOC(noAbpDyn.l,int,noAbpDyn.siz);
  }
  if (recAbpUnb.tgl != 0) {
	MALLOC(unbLog.l,double,nAbpC * 29 * 10);
  }
  if (recAbpBind.tgl != 0) {
	MALLOC(bindLog.l,double,nAbpC * 26 * 10);
  }
  if (recActSev.tgl != 0) {
	MALLOC(sevLog.l,double,nActC * 12 * 10);
  }
  if (recAbpTurn.tgl != 0) {
	MALLOC(toLog.l,double,nAbpC * 14 * 10);
  }
  // Unbinding of filaments from boundaries
  if (bndUnb.gTgl != 0) {
	MALLOC(bndUnb.p,double,bndUnb.maxF*NDIM*2);
  }
  if (bndReb.gTgl != 0) {
	MALLOC(bndReb.p,double,NDIM*2);
  }
  if (bndUnb.gTgl != 0 || bndReb.gTgl != 0) {
	MALLOC(bndUnb.cntMe,int,NDIM*2);
	MALLOC(bndReb.cntMe,int,NDIM*2);
  }
  // Moving boundaries
  if (bndMv.gTgl != 0) {
	MALLOC(bndMv.rankF.l,int,nCpu);
  }
  if (bndMv.gTgl != 0 || (rheoWay > 0 && bulkRheoType == 1)) {
	MALLOC(rGridInit,double,NDIM*2);
  }
  MALLOC(bnd.f,double,NDIM*2*NDIM);
  // Self-assembly of motors
  if (motSA.gTgl != 0) {
	MALLOC(abp.mId,int,nAbpC);
	if (gTglLoadNetData != 0) {
		MALLOC(abpMotId,int,nAbp);
		memset(abpMotId, -1, sizeof(int) * nAbp);
	}
  }
  // Membranes
  if (gTglMb != 0) {
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0 || mbReb.gTgl != 0 
			|| mbUnb.gTgl != 0) {
		sendMbDyn.siz = 3 * nMbC;
		MALLOC(sendMbDyn.l,int,sendMbDyn.siz);
		noMbDyn.siz = 3 * nMbC;
		MALLOC(noMbDyn.l,int,noMbDyn.siz);
	}
	if (mbUnb.gTgl != 0) {
		MALLOC(mbUnb.p,double,mbUnb.maxF);
	}
	if (mbDef.gTgl != 0) {
		MALLOC(mbDef.p,double,mbDef.maxF);
	}
	if (mbMat.gTgl != 0) {
		MALLOC(mbMat.p,double,mbMat.maxF);
		MALLOC(mbMat.cntMe2,int,2);
	}
	if (mbSld.act.gTgl != 0) {
	  	MALLOC(act.v,double,NDIM*nActC);
		MALLOC(mbSld.act.l,int,nActC*dimMbNuc);
		MALLOC(mbSld.act.reb,int,nActC);
		MALLOC(mbSld.act.side,int,nActC);
		MALLOC(mbSld.act.info,double,nActC*7);
	}
	if (mbSld.abp.gTgl != 0) {
	  	MALLOC(abp.v,double,NDIM*nAbpC);
		MALLOC(mbSld.abp.l,int,nAbpC*dimMbNuc);
		MALLOC(mbSld.abp.reb,int,nAbpC);
		MALLOC(mbSld.abp.info,double,nAbpC*7);
	}
	MALLOC(memb.v,double,NDIM*nMbC);
  }
  // Nucleus
  if (gTglNuc != 0) {
	if (nucUnb.gTgl != 0) {
		MALLOC(nucUnb.p,double,nucUnb.maxF);
	}
  }
  /*--------------------------------------------------------------------------*/

  /*---------------------  Information of long chains ------------------------*/
  // Information of long chains of ACP and motor
  // Normal 
  if (mpiMethod == 0) {
	MALLOC(longChDist,double,sizeArr * 2);
	MALLOC(longCh.l,int,sizeArr * 2);
	longChExtMsg.siz = sizeArr * 3;
	MALLOC(longChExtMsg.l,int,longChExtMsg.siz);
	longChIntMsg.siz = sizeArr * 2;
	MALLOC(longChIntMsg.l,int,longChIntMsg.siz);
  }
  // Plympton
  else {
	MALLOC(longChIntMsg.l,int,sizeArr);
  }
  maxDisp = 0.1;
  maxActCh = 1.5;
  /*--------------------------------------------------------------------------*/

  /*------------------------- Rheological measurements -----------------------*/
  if (recTraj.tgl != 0) {
	if (recTraj.gTglCho == 0 && gTglLoadNetData != 0) {
		MALLOC(recTraj.act.l,int,nAct);
	}
	MALLOC(recTraj.actMe.l,int,nActC);
  }
  if (recTraj.tgl2 != 0) {
	if (recTraj.gTglCho == 0 && gTglLoadNetData != 0) {
		MALLOC(recTraj.abp.l,int,nAbp);
	}
	MALLOC(recTraj.abpMe.l,int,nAbpC);
  }
  // Bulk rheology
  if (rheoWay > 0) { 
    MALLOC(stre.acc,double,NDIM);
	if (gTglLoadNetData != 0) {
		MALLOC(meaStrePar.l,int,nAct);
		MALLOC(appStraPar.l,int,nAct);
	}
	MALLOC(meaStreParMe.l,int,nActC);
	MALLOC(appStraParMe.l,int,nActC);
	MALLOC(rankMeaStre.l,int,nCpu);
  }
  /*--------------------------------------------------------------------------*/

  /*------------------------------ Record data ------------------------------*/
  // Record the instantaneous forces of ACP or motor
  MALLOC(recInstSprFabp,double,2*2*nAbpC);
  // Record the longitudinal forces of ACP and motor
  if (recLongF.tgl != 0) {
	MALLOC(recLongSprFabp,double,nAbpC*4);
  }
  // Coloring methods in VMD
  if (confVmdInfo.tgl != 0) {
	MALLOC(recAct.len,double,nActC);
	MALLOC(recAct.sprF,double,nActC);
	MALLOC(recAct.bendF,double,nActC);
	MALLOC(recAct.allF,double,nActC * (NDIM + 1));
	MALLOC(recAct.cnt,int,nActC);

	MALLOC(recAbp.len,double,nAbpC * recAbp.nL);
	MALLOC(recAbp.sprF,double,nAbpC);
	MALLOC(recAbp.bendF,double,nAbpC);
	MALLOC(recAbp.allF,double,nAbpC * (NDIM + 1));
	MALLOC(recAbp.cnt,int,nAbpC);

	MALLOC(recMb.len,double,nMbC);
	MALLOC(recMb.sprF,double,nMbC);
	MALLOC(recMb.bendF,double,nMbC);
	MALLOC(recMb.allF,double,nMbC * (NDIM + 1));
	MALLOC(recMb.cnt,int,nMbC);
  }	
  // Record the turnover of ACP and motor
  if (recAbpTurn.tgl != 0 && (acpUnb.gTgl != 0 || motUnb.gTgl != 0)) { 
	MALLOC(abpTurn,double,nAbpC*7);
  }
  MALLOC(iPerc, int, nActGoal + nAbpGoal);
  MALLOC(iSupp, int, nActGoal + nAbpGoal);
  if (recSecStre.tgl != 0) {
	MALLOC(viscStre,double,V3PROD(nSecStreDiv) * NDIM);
	MALLOC(elasStre,double,V3PROD(nSecStreDiv) * NDIM);
  }
  /*--------------------------------------------------------------------------*/

  /*-------------------------- For parallelization ---------------------------*/
  sizeMpiArr = (mpiMethod == 0) ? cntAdjRank[0] : 2;
  MALLOC(sReq, MPI_Request, sizeMpiArr);
  MALLOC(rReq, MPI_Request, sizeMpiArr);
  MALLOC(cntCpPar,int,sizeMpiArr*3);
  MALLOC(cntMvPar,int,sizeMpiArr*3);
  MALLOC(cpPar,ListInt,sizeMpiArr);
  MALLOC(mvPar,ListInt,sizeMpiArr);
  MALLOC2(bufSendMsg,char,sizeMpiArr);
  MALLOC2(bufRecvMsg,char,sizeMpiArr);
  MALLOC(mpiTestSendFlag,int,sizeMpiArr);
  MALLOC(mpiTestRecvFlag,int,sizeMpiArr);
  for(n = 0; n < sizeMpiArr; n++) {
	MALLOC(cpPar[n].l,int,sizeArr);
	MALLOC(mvPar[n].l,int,sizeArr);
	MALLOC(bufSendMsg[n],char,sizeBufMsg);
	MALLOC(bufRecvMsg[n],char,sizeBufMsg);
  }
  MALLOC(insNeiPar.l,int,nActGoal + ((gTglMb != 0) ? nUnitMb : 0));
  /*--------------------------------------------------------------------------*/

  // Misc.
  MALLOC(arrAcos,double,2*DEG_ARR_ACOS+1); 
}

// Assign values to dynamic arrays. 
void AssignArrayValues(void) {
  int n, k, nMax, indCell[NDIM], direc, sizeArr;
  double k0, ang, angMax, cosValue, lenR, dimDomC[NDIM], *pD;

  sizeArr = nActC + nAbpC;
  if (gTglMb != 0) { sizeArr += nMbC; }

  // Maximum Brownian force (this satisfies fluctuation-dissipation theorem)
  act.maxFbr[0] = sqrt(2.0 / dt);
  for(n = 0; n < NK_ABP; n++) {
	abp.maxFbr[n] = act.maxFbr[0] * sqrt(abpF.drag[n].n);
  }
  memb.maxFbr = sqrt(2.0 / dt) * 2.;
  memb.nucMaxFbr = sqrt(2.0 / dt) * 2.;
  // List of actin segments clamped on boundaries
  if (gTglLoadNetData != 0) {
	memset(fixAct, -1, sizeof(int) * nAct);
  }
  memset(act.fix, -1, sizeof(int) * nActC);

  if (gTglBead != 0 && beadBind.gTgl != 0) {
	memset(act.bdFix, -1, sizeof(int) * nActC);
  }

  // Local index of actin filaments
  memset(act.iF, -1, sizeof(int) * nActC);

  if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
    memset(act.nFA, 0, sizeof(int) * nActC);
  }

  // Sectioning lines
  V3DIV(dimDomC, dimDom, nCell);
  FOR_NDIM(k) {
	rGrid[k][0] = 0.;
	rGrid[k][nGrid[k] - 1] = dimDom[k];
	for(n = 1; n < nGrid[k] - 1; n++) {
		rGrid[k][n] = dimDomC[k] * (double)n;
	}
  }
  /*--------------------- Dynamic behaviors of actin and ABP -----------------*/
  CheckPeriodToggleParameter(&updMono);
  /*--------------------------------------------------------------------------*/

  /*----------------------- Dynamic behaviors of actins ----------------------*/
  if (actNuc.gTgl != 0) {
	actNuc.facP = FAC_P(actNuc.k * 2. * nActPerSeg);
  }
  if (actSev.gTgl != 0) {
	for(k = 0; k < 2; k++) {
		pD = (k == 0) ? actSev.p : actSev.pWA;
		k0 = actSev.k * (k == 0 ? 1. : actSev.facKWA);
		angMax = log(log(1e-10) / (-1. * k0 * dtReal)) * actSev.facX;
	    for (n = 0; n < DEG_ARR_ACTSEV + 1; n++) {
	        ang = (double)n / DEG_ARR_ACTSEV * 180.;
    	    pD[n] = (ang < angMax) ? K2P(k0 * exp(ang / actSev.facX)) : 1.;
		}
	}
  }
  if (actAnn.gTgl != 0) {
	actAnn.p = K2P(actAnn.k);
  }
  if (actAss.gTgl != 0) {
	for(n = 0; n < 2; n++) {
		actAss.facP[n] = FAC_P(actAss.k[n]);
	}
  }
  if (actDis.gTgl != 0) {
	for(n = 0; n < 2; n++) {
		actDis.p[n] = K2P(actDis.k[n] / (2. * nActPerSeg));
        actDis.pWA[n] = K2P(actDis.k[n] * actDis.facKWA 
				/ (2. * nActPerSeg));
	}
  }
  if (actBst.gTgl != 0) {
	for(n = 0; n < 2; n++) {
		actBst.p[n] = K2P(actBst.k[n]);
	}
  }
  // Capping 
  if (gTglActCapAll != 0) {
	memset(act.cap, -1, sizeof(int) * nActC);
  }
  if (actCap.gTgl != 0) {
	for(n = 0; n < 2; n++) {
		actCap.p[n] = K2P(actCap.k[n]);
	}
  }
  // Uncapping
  if (actUnc.gTgl != 0) {
	for(n = 0; n < 2; n++) {
		actUnc.p[n] = K2P(actUnc.k[n]);
	}
  }
  // Count the aging of actin filaments
  if (actAge.gTgl != 0) {
	memset(act.age, 0, sizeof(int) * nActC);
  }
  memset(act.len, 0, sizeof(int) * nActC);
  // Maturation of actins  on boundaries
  if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
	for(n = 0; n < bndMat.maxF; n++) {
		bndMat.p[n] = K2P(K0_ACT_MAT * bndMat.facK0 * 
				exp(bndMat.facX * X_ACT_MAT / KT_IN_J * F_S2N((double)n)));
	}
	memset(bndMat.cntMe2, 0, sizeof(int) * 12);
  }
  /*--------------------------------------------------------------------------*/

  /*------- Dynamic behaviors of ACP, motor, boundaries, and membranes -------*/
  // Count the aging of ABPs
  if (abpAge.gTgl != 0) {
	memset(abp.age, 0, sizeof(int) * nAbpC);
  }
  if (recAbpDyn.tgl != 0) {
	memset(abpDyn.cntU, 0, sizeof(int) * nAbpC);
	memset(abpDyn.cntB, 0, sizeof(int) * nAbpC);
	memset(abpDyn.cntW, 0, sizeof(int) * nAbpC);
  }
  CheckPeriodToggleParameter(&updMotN);
  updMotN.tgl = (motSA.gTgl != 0) ? 1 : 0;
  // Unbinding of ACP
  if (acpUnb.gTgl != 0 || acpInaUnb.gTgl != 0) {
	for(n = 0; n < acpUnb.maxF; n++) {
		acpUnb.p[n] = K2P(acpUnb.facK0 * exp(acpUnb.facX / KT_IN_J 
				* F_S2N((double)n))	+ acpUnb.facK0c * 
				exp(acpUnb.facXc / KT_IN_J * F_S2N((double)n)));
	}
  }
  // Walking of motor
  if (motWalk.gTgl != 0 && nAbpGoalDet[2] > 0) {
	InitMotorWalkUnbindRates(0);
  }
  // Unbinding of motor
  if ((motUnb.gTgl != 0 || motInaUnb.gTgl != 0) && nAbpGoalDet[2] > 0) {
	InitMotorWalkUnbindRates(1);
  }
  if (motSA.gTgl != 0) {
	motSA.pAss = K2P(motSA.kAss);
	// Turnover of motor filaments
	if (motSA.to.gTgl != 0) {
		motSA.to.pf = K2P(motSA.to.k);
		if (motSA.to.mode == 1) {
			for(n = 0; n < motSA.to.maxF; n++) {
				motSA.to.p[n] = K2P(motSA.to.k0 *
						exp(motSA.to.x / KT_IN_J * F_S2N((double)n)));
			}
		}
	}
  }
  // Arrays related to unbinding of ACP and motor
  if (acpUnb.gTgl != 0 || acpReb.gTgl != 0 || motUnb.gTgl != 0 
		|| motReb.gTgl != 0 || motWalk.gTgl != 0 || acpInaUnb.gTgl != 0 
		|| motInaUnb.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0 
		|| actSev.gTgl != 0) {
    memset(noAbpDyn.l, -1, sizeof(int) * noAbpDyn.siz);
  }
  
  if (actSev.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0 
		|| actNuc.gTgl != 0 || actAss.gTgl != 0 || actAnn.gTgl != 0) {
	memset(noActDyn.l, -1, sizeof(int) * noActDyn.siz);
  }
  // Unbinding of filaments from boundaries
  if (bndUnb.gTgl != 0) {
	for(n = 0; n < NDIM * 2; n++) {
	    for(k = 0; k < bndUnb.maxF; k++) {
			P2A(bndUnb.p,n,k,bndUnb.maxF) = K2P(K0_BND_UNB * bndUnb.facK0[n] 
					* exp(bndUnb.facX[n] * X_BND_UNB / KT_IN_J 
					* F_S2N((double)k)));
		}
	}
  }
  if (bndReb.gTgl != 0) {
	for(n = 0; n < NDIM * 2; n++) {
		bndReb.p[n] = K2P(K_BND_REB * bndReb.facK[n]);
	}
  }
  // Unbinding of filaments from boundaries
  if (bndUnb.gTgl != 0 || bndReb.gTgl != 0) {
	memset(bndUnb.cntMe, 0, sizeof(int) * NDIM * 2);
	memset(bndReb.cntMe, 0, sizeof(int) * NDIM * 2);
  }
  // Moving boundaries
  if (bndMv.gTgl != 0 || (rheoWay > 0 && bulkRheoType == 1)) {
	// Initial location of sectioning lines should be memorized to calculate
	// stress-strain relationship.
	FOR_NDIM(k) {
		P2A(rGridInit,0,k,NDIM) = rGrid[k][0];
		P2A(rGridInit,1,k,NDIM) = rGrid[k][nGrid[k] - 1];
	}
  }
  if (bndMv.gTgl != 0 || bndUnb.gTgl != 0 || bndReb.gTgl != 0) {
	// "bndMv.rankF" is the list of CPUs having at least one boundary
	bnd.gotF = 0;
	if (bndMv.gTgl != 0) { bndMv.rankF.c = 0; }
	for(n = 0; n < nCpu; n++) {
		V3IND_ASSIGN_INT(n, nCell, 1, indCell);
		FOR_NDIM(k) {
			CONT(pbc[k] != 0);
			CONT(!(indCell[k] == 0 || indCell[k] == nCell[k] - 1));
			if (bndMv.gTgl != 0) {
				InsertElement1dArrayWoChk(bndMv.rankF.l, &bndMv.rankF.c, n);
			}
			if (n == rank) { 
				bnd.gotF = 1; 
			}
			break;
		}
	}
  }
  // If bulk rheology
  if (rheoWay > 0) {
	// Shear or normal
	direc = (bulkRheoType == 0) ? dirNoPBC : dirStr;
	gotMeaStre = 0;
	rankMeaStre.c = 0;
	for(n = 0; n < nCpu; n++) {
		V3IND_ASSIGN_INT(n, nCell, 1, indCell);
		CONT(!(indCell[direc] == nCell[direc] - 1));
		InsertElement1dArrayWoChk(rankMeaStre.l, &rankMeaStre.c, n);
		if (n == rank) { 
			gotMeaStre = 1; 
		}
	}	
  }
  // Membrane
  if (gTglMb != 0) {
	if (mbUnb.gTgl != 0) {
		for(n = 0; n < mbUnb.maxF; n++) {
			mbUnb.p[n] = K2P(K0_MEMB_UNB * mbUnb.facK0 * 
					exp(mbUnb.facX * X_MEMB_UNB / KT_IN_J * F_S2N((double)n)));
		}
	}
	if (mbDef.gTgl != 0) {
		for(n = 0; n < mbDef.maxF; n++) {
			mbDef.p[n] = K2P(K0_MEMB_DEF * mbDef.facK0 * 
					exp(mbDef.facX * X_MEMB_DEF / KT_IN_J * F_S2N((double)n)));
		}
	}
	if (mbMat.gTgl != 0) {
		for(n = 0; n < mbMat.maxF; n++) {
			mbMat.p[n] = K2P(K0_MEMB_MAT * mbMat.facK0 * 
					exp(mbMat.facX * X_MEMB_MAT / KT_IN_J * F_S2N((double)n)));
		}
	}
	if (mbReb.gTgl != 0 || mbUnb.gTgl != 0 
			|| nucReb.gTgl != 0 || nucUnb.gTgl != 0) {
		memset(memb.reb, -1, sizeof(int) * nMbC);
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		memset(memb.fix, -1, sizeof(int) * nMbC);
	} 
	if (mbPro.gTgl != 0) {
		memset(memb.pro, -1, sizeof(int) * nMbC * 2);
		if (gTglMbCont != 0) {
			for(n = 0; n < nMbC; n++) {
				memb.myo.conc[n] = 1.;
				memb.cyto.conc[n] = 1.;
			}
		}
	}
	if (mbMat.gTgl != 0) {
		memset(memb.nFA, 0, sizeof(int) * nMbC * nMbAct);
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0 || mbReb.gTgl != 0 
			|| mbUnb.gTgl != 0) {
		memset(noMbDyn.l, -1, sizeof(int) * noMbDyn.siz);
	}
	if (gTglLoadNetData != 0) {
		memset(actRebMb, -1, sizeof(int) * nAct);
	}
	memset(act.mbReb, -1, sizeof(int) * nActC);
	for(n = 0; n < nObjMbNuc; n++) {
		memb.vol.val[n] = 0.;
		memb.area.val[n] = 0.;
	}
	if (mbSld.act.gTgl != 0) {
		memset(mbSld.act.l, -1, sizeof(int) * nActC * dimMbNuc);
		memset(mbSld.act.reb, -1, sizeof(int) * nActC);
	}
	if (mbSld.abp.gTgl != 0) {
		memset(mbSld.abp.l, -1, sizeof(int) * nAbpC * dimMbNuc);
		memset(mbSld.abp.reb, -1, sizeof(int) * nAbpC);
	}
	MALLOC(memb.v,double,NDIM*nMbC);
  }
  // Nucleus
  if (gTglNuc != 0) {
	if (nucUnb.gTgl != 0) {
		for(n = 0; n < nucUnb.maxF; n++) {
			nucUnb.p[n] = K2P(K0_NUC_UNB * nucUnb.facK0 * 
					exp(nucUnb.facX * X_NUC_UNB / KT_IN_J * F_S2N((double)n)));
		}
	}
  }
  // Self-assembly of motors
  if (motSA.gTgl != 0) { 
	memset(abp.mId, -1, sizeof(int) * nAbpC);
	if (gTglLoadNetData != 0) {
		memset(abpMotId, -1, sizeof(int) * nAbp);
	}
  }
  /*--------------------------------------------------------------------------*/

  /*---------------------  Information of long chains ------------------------*/
  if (mpiMethod == 0) {
	longCh.c = 0;
	longChExtMsg.c = 0;
	memset(longCh.l, -1, sizeof(int) * sizeArr * 2);
	memset(longChExtMsg.l, -1, sizeof(int) * sizeArr * 3);
	memset(longChIntMsg.l, -1, sizeof(int) * sizeArr * 2);
	for(n = 0; n < sizeArr * 2; n++) {
		longChDist[n] = 0.;
	}
  }
  else {
	memset(longChIntMsg.l, -1, sizeof(int) * sizeArr);
  }
  longChIntMsg.c = 0;
  /*--------------------------------------------------------------------------*/

  /*------------------------------- Record data ------------------------------*/
  // Record the longitudinal forces of ACP and motor
  if (recLongF.tgl != 0) {
	for(n = 0; n < nAbpC * 4; n++) {
		recLongSprFabp[n] = 0; 
	}
  }
  // Coloring methods in VMD
  if (confVmdInfo.tgl != 0) {
	for(n = 0; n < nActC; n++) {
		recAct.len[n] = 0.;		
		recAct.sprF[n] = 0.;	
		recAct.bendF[n] = 0.;
		V4SET_ALL(&P2A(recAct.allF,n,0,NDIM + 1),  0.);
		recAct.cnt[n] = 0;
	}	
	for(n = 0; n < nAbpC; n++) {
		for(k = 0; k < recAbp.nL; k++) {	
			P2A(recAbp.len,n,k,recAbp.nL) = 0.;		
		}
		recAbp.sprF[n] = 0.;	
		recAbp.bendF[n] = 0.;
		V4SET_ALL(&P2A(recAbp.allF,n,0,NDIM + 1),  0.);
		recAbp.cnt[n] = 0;
	}
	for(n = 0; n < nMbC; n++) {
		recMb.len[n] = 0.;		
		recMb.sprF[n] = 0.;	
		recMb.bendF[n] = 0.;
		V4SET_ALL(&P2A(recMb.allF,n,0,NDIM + 1),  0.);
		recMb.cnt[n] = 0;
	}	
  }
  /*--------------------------------------------------------------------------*/

  /*------------------------- Rheological measurements -----------------------*/
  stra.accDsp = 0.;
  if (rheoWay > 0) {
    V3ZERO(stre.acc);
    if (bulkRheoWay == 0) { 
		stre.curr = 0.;
		stre.goal = 0.;
		stre.accShUS = 0.;
		stra.accDspUS = 0.; 
		stre.accShErr = 0.;
	}
	if (gTglLoadNetData != 0) {
		memset(meaStrePar.l, -1, sizeof(int) * nAct);
		memset(appStraPar.l, -1, sizeof(int) * nAct);
	}
	memset(meaStreParMe.l, -1, sizeof(int) * nActC);
	memset(appStraParMe.l, -1, sizeof(int) * nActC);
  }
  if (recTraj.tgl != 0) { recTraj.actMe.c = 0; }
  if (recTraj.tgl2 != 0) { recTraj.abpMe.c = 0; }
  memset(iPerc, 0, sizeof(int) * (nActGoal + nAbpGoal));
  memset(iSupp, -1, sizeof(int) * (nActGoal + nAbpGoal));
  if (recSecStre.tgl != 0) {
	for(n = 0; n < V3PROD(nSecStreDiv) * NDIM; n++) {
		viscStre[n] = 0.;
		elasStre[n] = 0.;
  	}
  }
  /*--------------------------------------------------------------------------*/
  // Pre-assigned array for acos(x)
  for (n = 0; n < 2 * DEG_ARR_ACOS + 1; n++) {
    cosValue = (1. / (double)DEG_ARR_ACOS) * (double)n - 1.;
	cosValue = TrimDblVal(cosValue, -1., 1.);
    arrAcos[n] = acos(cosValue);
  }
}

void CheckPeriodToggleParameter(FuncCont *cont) {
  if (cont->prdR < 0) {
	cont->tgl = 0;
	memset(&cont->prd, 1, sizeof(int));
  }
  else {
	cont->tgl = 1;
	cont->prd = T_SEC2TS(cont->prdR);
  }
}

void AssignInitValues(void) {
  int n, k, sizeRankArr, ind, adjCpu, fac, nEle, nUnit, pbcCond[NDIM];
  int oft[NDIM], iAdjCell[NDIM], *chkCpu, *lev, *nPerObj;
  double sumCos, ratio, minDim, dimDomC[NDIM], mbMaxSiz, ang;
  double volDom, volMb, volNuc;
  double *len, rad;
 
  nChAc = 2 + nChAcX * nChAcY;
  /*-------------------------- Actin concentration ---------------------------*/
  if (cAct > 0) {
	if (gTglMb != 0) {
		if (dimMbNuc == 2) {
			volDom = dimDom[dirMbNuc[0]] * dimDom[dirMbNuc[1]];
			volMb = PI * SQR(L_UM2S(radMb));
			volNuc = (gTglNuc != 0) ? PI * SQR(L_UM2S(radNuc)) : 0.;
		}
		else {
			volDom = V3PROD(dimDom);
			volMb = 4 * PI / 3. * CUBE(L_UM2S(radMb));
			volNuc = (gTglNuc != 0) ? 4 * PI / 3. * CUBE(L_UM2S(radNuc)) : 0.;
		}
		if (sideMb == 0) {
			cAct *= (volMb - volNuc) / volDom;
		}
		else {
			cAct *= (volDom - volMb) / volDom;
		}
	}
    nActGoal = (int)(cAct / (2. * nActPerSeg / N_AVO
            / (V3PROD(dimDom) * CUBE(L_SCALE_IN_M)) * 1.0e3));
    if (nActGoal < nAct) {
		Printf0("Error: the number of actins specified in 'condition' "
				"(%d) is smaller than that in 'Config' (%d)!!\n", 
				nActGoal, nAct);
        exit(-1);
    }
	if (gTglMb != 0 && gTglLoadNetData == 0) {
		nActMall = nActGoal - nAct;
	}
	else {
		nActMall = (int)((double)(nActGoal - nAct) / (double)nCpu) * nCpu;
	}
	nActGoal = nAct + nActMall;
  } 
  else { 
  	if (gTglLoadNetData == 0 && gTglMb == 0) {
		Printf0("Error: actin concentration should be specified to begin "
				"without a pre-assembled network !!\n");
        exit(-1);
	}
	nActMall = 0;
	nActGoal = nAct; 
  }
  /*--------------------------------------------------------------------------*/

  /*--------------------------- ABP concentration ----------------------------*/
  for(n = 0 ; n < NK_ABP; n++) {
	if (RAbp[n] >= 0) {
		nAbpGoalDet[n] = (int)(RAbp[n] * nActPerSeg * (double)nActGoal);
		nAbpMall[n] = nAbpGoalDet[n] - nAbpDet[n];
	}
	else {
		nAbpGoalDet[n] = nAbpDet[n];
		nAbpMall[n] = 0;
	}
  }
  nAbpGoal = nAbp + V3SUM(nAbpMall);
  /*--------------------------------------------------------------------------*/

  /*------------- Geometric and mechanical properties of membrane ------------*/
  if (gTglMb != 0) {
	if (dimMbNuc == 2) {
		for(n = 0; n < 2; n++) {
			dirMbNuc[n] = (dirNormMbNuc + n + 1) % NDIM;
		}
	}
	else {
		V3SET(dirMbNuc, 0, 1, 2);
	}
	radMb = L_UM2S(radMb);
	memb.len = L_UM2S(memb.len);
	memb.thk = L_UM2S(memb.thk);
	nChMb = (dimMbNuc == 2) ? 2 : 6;
	nChMb += nMbAct;

	radNuc = L_UM2S(radNuc);
	memb.nucLen = L_UM2S(memb.nucLen);
	memb.nucThk = L_UM2S(memb.nucThk);

	nMb = 0;
	nUnitMb = 0;
	fac = (gTglNuc != 0) ? 2 : 1;
	for(n = 0; n < fac; n++) {
		if (n == 0) {
			len = &memb.len;
			rad = radMb;
			lev = &levMb;
			nPerObj = &nMbPerObj;
		}
		else {
			len = &memb.nucLen;
			rad = radNuc;
			lev = &levNuc;
			nPerObj = &nNucPerObj;
		}
		// 2D membrane
		if (dimMbNuc == 2) {
			nEle = (int)(PI / asin(0.5 * (*len) / rad)) + 1;
			ang = 2. * PI / (double)nEle;
			*len = 2. * rad * sin(0.5 * ang);
			nUnit = nEle;
		}
		// 3D membrane
		else {
			*lev = (int)ceil(log(1.05146 * rad / (*len)) / log(2));
			nEle = 12;
			for(k = 0; k < (*lev); k++) {
				nEle += 30 * (int)pow(4., k);
			}
			nUnit = 20 * (int)pow(4., (*lev));
		}
		*nPerObj = nEle;
		nEle *= nObjMbNuc / fac;
		nUnit *= nObjMbNuc / fac;
		nMb += nEle;
		nUnitMb += nUnit;
	}
	nMbC = nMb * 2;

	memb.vol.stf = PA2S(STF_MEMB_VOL) * memb.vol.facStf;
	memb.vol.stf2 = PA2S(STF_MEMB_VOL2) * memb.vol.facStf;
	memb.area.stf = KS_NPM2S(STF_MEMB_AREA) * memb.area.facStf;
	memb.bend.stf = KB_NM2S(STF_MEMB_BEND) * memb.bend.facStf;
	memb.bend.eq = 0.;
	memb.spr2.stf = KS_NPM2S(STF_MEMB_ACT_SPR) * memb.spr2.facStf;
	memb.spr2.eq = AVG2(L_M2S(DIA_CYL_ACT), memb.thk) + 0.2;
	memb.spr2.lo = DTL_MEMB_REB * memb.spr2.eq;
	memb.spr2.hi = DTH_MEMB_REB * memb.spr2.eq;
	memb.spr3.stf = KS_NPM2S(STF_MEMB_FIX_SPR);
	memb.spr3.eq = 0.;
	memb.spr.stf = KS_NPM2S(STF_MEMB_SPR) * memb.spr.facStf;
	memb.spr.lo = DTL_MEMB_SPR;
	memb.spr.hi = DTH_MEMB_SPR;
	memb.myo.stf = KS_NPM2S(STF_MEMB_SPR_MYO);
	memb.myo.diff = T_S2SEC(DIFF_MEMB_SPR_MYO / (mbPro.durR / DUR_MEMB_PROT));
	memb.cyto.stf = KS_NPM2S(STF_MEMB_SPR_CYTO);
	memb.cyto.diff = T_S2SEC(DIFF_MEMB_SPR_CYTO / (mbPro.durR / DUR_MEMB_PROT));

	memb.nucVol.stf = PA2S(STF_NUC_VOL) * memb.nucVol.facStf;
	memb.nucVol.stf2 = NPV2S(STF_NUC_VOL2) * memb.nucVol.facStf;
	memb.nucArea.stf = PA2S(STF_NUC_AREA) * memb.nucArea.facStf;
	memb.nucBend.stf = KB_NM2S(STF_NUC_BEND) * memb.nucBend.facStf;
	memb.nucBend.eq = 0.;
	memb.nucSpr2.stf = KS_NPM2S(STF_NUC_ACT_SPR) * memb.nucSpr2.facStf;
	memb.nucSpr2.eq = AVG2(L_M2S(DIA_CYL_ACT), memb.nucThk) + 0.2;
	memb.nucSpr2.lo = DTL_NUC_REB * memb.nucSpr2.eq;
	memb.nucSpr2.hi = DTH_NUC_REB * memb.nucSpr2.eq;
	memb.nucSpr.stf = KS_NPM2S(STF_NUC_SPR);
	memb.nucSpr.lo = DTL_NUC_SPR;
	memb.nucSpr.hi = DTH_NUC_SPR;
  }
  else {
	nMb = 0;
	nMbMe = 0;
	nMbCp = 0;
	nMbC = 0;
  }
  nMbMin = (nMbC < nMb) ? nMbC : nMb;
  /*--------------------------------------------------------------------------*/

  // Stiffness for repulsive forces between elements and boundaries
  bnd.stfRep = KS_NPM2S(STF_REP_BND);
  stfRepMbNuc = KS_NPM2S(STF_REP_MEMB);
  // Stiffness of surrounding medium
  for(k = 0 ; k < NDIM * 2; k++) {
	if (bndMv.stfUnit == 0) { bndMv.stf[k] = PA2S(bndMv.stf[k]); }
	else { bndMv.stf[k] = KS_NPM2S(bndMv.stf[k]); }
	bndMv.thk[k] = L_UM2S(bndMv.thk[k]);
	bnd.stfSpr[k] = bnd.stfSpr[k] * KS_NPM2S(STF_SPR_BND);
  }
  if (bndVol.gTgl != 0) {
	bndVol.stf = STF_BND_VOL * bndVol.facStf / (KT_IN_J
			/ pow(L_SCALE_IN_M, 6.));
	bndVol.eq = V3PROD(dimDom);
  }
  bnd.radRnd = L_UM2S(bnd.radRnd);
  bnd.gTglRnd = (bnd.radRnd < 0.) ? 0 : 1;
  if (bnd.radRnd > 0.5 * dimDom[0]) {
    Printf0("Error: the radius of the circular boundary is too large"
            " compared to the domain size!\n");
    exit(-1);
  }

  /*------ Geometric and mechanical properties of actin, ACP, and motor ------*/

  // Diameter of cylindrical segments
  actF.dia = L_M2S(DIA_CYL_ACT);
  abpF.dia[0] = L_M2S(DIA_CYL_ACPC);
  abpF.dia[1] = L_M2S(DIA_CYL_ACPB);
  abpF.dia[2] = L_M2S(DIA_CYL_MOT);
  // Diameter for calculation of repulsive forces
  for (k = 0; k < NK_ABP; k++) {
	abpF.repDia[k] = abpF.dia[k];
  }
  // Size of ACP and motor
  abpF.len[0].n = L_M2S(L_ACPC_ARM);
  abpF.len[1].n = L_M2S(L_ACPB_ARM);
  abpF.len[2].n = L_M2S(L_MOT_ARM);
  // From the size and diameter, various distances calculated.
  for (k = 0; k < NK_ABP; k++) {
	abpF.len[k].inv = INV(abpF.len[k].n);
	abpF.dtRepAct[k].n = AVG2(actF.dia, abpF.dia[k]);
	abpF.dtSprAct[k].n = 0.5 * actF.dia + abpF.len[k].n; 
	abpF.dtRepAct[k].inv = INV(abpF.dtRepAct[k].n);
	abpF.dtSprAct[k].inv = INV(abpF.dtSprAct[k].n);
  }
  // actin
  actF.rep.stf = actF.rep.facStf * KS_NPM2S(STF_ACT_REP);
  actF.bend.stf = actF.bend.facStf * STF_ACT_BEND / L_SCALE_IN_M / KT_IN_J;
  actF.bend.eq = 0.;
  actF.bend.lo = cos(0. + DEG2RAD(ANGL_ACT_BEND));
  actF.spr.stf = actF.spr.facStf * KS_NPM2S(STF_ACT_SPR / nActPerSeg 
		/ (nActPerSeg / 10.)) * 4.;
  actF.spr.eq = 1.;
  actF.spr.lo = DTL_ACT_ANN * actF.spr.eq;
  actF.spr.hi = DTH_ACT_ANN * actF.spr.eq;
  // ABP
  abpF.rep.stf = abpF.rep.facStf * KS_NPM2S(STF_ABP_REP);
  // ACP^C
  abpF.bend[0].stf = abpF.bend[0].facStf * KB_NM2S(STF_ACPC_BEND);
  abpF.bend[0].eq = DEG2RAD(ANG_ACPC_BEND);
  abpF.bend[0].lo = cos(DEG2RAD(TrimDblVal(ANG_ACPC_BEND 
		+ ANGL_ACPC_BEND, 0., 180.)));
  abpF.bend[0].hi = cos(DEG2RAD(TrimDblVal(ANG_ACPC_BEND 
		- ANGL_ACPC_BEND, 0., 180.)));
  abpF.a90[0].stf = abpF.a90[0].facStf * KB_NM2S(STF_ACPC_90);
  abpF.a90[0].eq = PI * 0.5;
  abpF.a90[0].lo = cos(PI / 2. - DEG2RAD(ANGL_ACPC_90));
  abpF.spr[0].stf = abpF.spr[0].facStf * KS_NPM2S(STF_ACPC_SPR); 
  abpF.spr[0].eq = abpF.dtSprAct[0].n;
  abpF.spr[0].lo = DTL_ACP_REB * abpF.dtSprAct[0].n;
  abpF.spr[0].hi = DTH_ACP_REB * abpF.dtSprAct[0].n;
  abpF.cr[0].stf = KB_NM2S(STF_ACPC_TOR);
  abpF.cr[0].eq = PI * 0.5;
  abpF.cr[0].lo = cos(TrimDblVal(abpF.cr[0].eq 
		+ DEG2RAD(ANGL_ACPC_TOR), 0., PI));
  abpF.cr[0].hi = cos(TrimDblVal(abpF.cr[0].eq 
		- DEG2RAD(ANGL_ACPC_TOR), 0., PI));
  // ACP^B
  abpF.bend[1].stf = abpF.bend[1].facStf * KB_NM2S(STF_ACPB_BEND);
  abpF.bend[1].eq = DEG2RAD(ANG_ACPB_BEND);
  abpF.bend[1].lo =cos(DEG2RAD(TrimDblVal(ANG_ACPB_BEND 
		+ ANGL_ACPB_BEND, 0., 180.)));
  abpF.bend[1].hi =cos(DEG2RAD(TrimDblVal(ANG_ACPB_BEND 
		- ANGL_ACPB_BEND, 0., 180.)));
  abpF.a90[1].stf = abpF.a90[1].facStf * KB_NM2S(STF_ACPB_90);
  abpF.a90[1].eq = PI * 0.5;
  abpF.a90[1].lo = cos(PI / 2. - DEG2RAD(ANGL_ACPB_90));
  abpF.spr[1].stf = abpF.spr[1].facStf * KS_NPM2S(STF_ACPB_SPR); 
  abpF.spr[1].eq = abpF.dtSprAct[1].n;
  abpF.spr[1].lo = DTL_ACP_REB * abpF.dtSprAct[1].n;
  abpF.spr[1].hi = DTH_ACP_REB * abpF.dtSprAct[1].n;
  abpF.cr[1].stf = KB_NM2S(STF_ACPB_TOR);
  abpF.cr[1].eq = 0.;
  abpF.cr[1].lo = cos(TrimDblVal(abpF.cr[1].eq 
		+ DEG2RAD(ANGL_ACPB_TOR), 0., PI));
  abpF.cr[1].hi = cos(TrimDblVal(abpF.cr[1].eq 
		- DEG2RAD(ANGL_ACPB_TOR), 0., PI));
  // motor
  // Equilibrium angle of "Bend" depends on the size of motor.
  abpF.bend[2].stf = abpF.bend[2].facStf * KB_NM2S(STF_MOT_BEND);
  abpF.bend[2].eq = DEG2RAD(ANG_MOT_BEND);
  abpF.bend[2].lo = cos(DEG2RAD(TrimDblVal(ANG_MOT_BEND 
		+ ANGL_MOT_BEND, 0., 180.)));
  abpF.bend[2].hi = cos(DEG2RAD(TrimDblVal(ANG_MOT_BEND
		- ANGL_MOT_BEND, 0., 180.)));
  abpF.a90[2].stf = abpF.a90[2].facStf * KB_NM2S(STF_MOT_90);
  abpF.a90[2].eq = PI * 0.5;
  abpF.a90[2].lo = cos(PI / 2. - DEG2RAD(ANGL_MOT_90));
  abpF.spr[2].stf = abpF.spr[2].facStf * KS_NPM2S(STF_MOT_SPR2); 
  abpF.spr[2].stf2 = abpF.spr[2].facStf * KS_NPM2S(STF_MOT_SPR);
  abpF.spr[2].eq = abpF.dtSprAct[2].n;
  abpF.spr[2].lo = DTL_MOT_REB * abpF.dtSprAct[2].n;
  abpF.spr[2].hi = DTH_MOT_REB * abpF.dtSprAct[2].n;
  abpF.cr[2].stf = KB_NM2S(STF_MOT_TOR);
  abpF.cr[2].eq = DEG2RAD(ANG_MOT_TOR);
  abpF.cr[2].lo = cos(TrimDblVal(abpF.cr[2].eq 
		+ DEG2RAD(ANGL_MOT_TOR), 0., PI));
  abpF.cr[2].hi = cos(TrimDblVal(abpF.cr[2].eq 
		- DEG2RAD(ANGL_MOT_TOR), 0., PI));

  if (motSA.gTgl != 0) {
	motSA.bend.stf = KB_NM2S(STF_MOTBACK_BEND);
	motSA.bend.eq = DEG2RAD(ANG_MOTBACK_BEND);
	motSA.spr.stf = KS_NPM2S(STF_MOTBACK_SPR) * abpF.spr[2].facStf;
	motSA.spr.eq = L_M2S(L_MOTBACK_DIST);
	motSA.cenDist = L_M2S(L_MOTBACK_CEN_DIST);
	motSA.nMotPerSide = motSA.nMotPerTF / 2;
  }
  /*--------------------------------------------------------------------------*/

  // Drag coeffcients of actin, ACP, and motors
  actF.dragR = CYL_DRAG(DIA_CYL_ACT, DIA_CYL_ACT * nActPerSeg);
  abpF.drag[0].n = CYL_DRAG(DIA_CYL_ACPC, L_ACPC_ARM) / actF.dragR;
  abpF.drag[1].n = CYL_DRAG(DIA_CYL_ACPB, L_ACPB_ARM) / actF.dragR;
  abpF.drag[2].n = (motSA.gTgl != 0) ? 1. : 
		CYL_DRAG(DIA_CYL_MOT, L_MOT_ARM) / actF.dragR;
  for (k = 0; k < NK_ABP; k++) {
	abpF.drag[k].inv = INV(abpF.drag[k].n);
  }

  // Time variables
  dt = 0.5e-5; 
  dtReal = T_S2SEC(dt);
  time(&initTime);    
  currTimeStep = 0;    
  rheo.dur = T_SEC2TS(rheo.durR);
  netForm.dur = T_SEC2TS(netForm.durR);
  motActiv.dur = T_SEC2TS(motActiv.durR);

  /*-------------------------- Parallel processing ---------------------------*/
  V3IND_ASSIGN_INT(rank, nCell, 1, iCell);
  VS3COPY(dimDomH, dimDom, 0.5);
  V3DIV(dimDomC, dimDom, nCell);
  minDimDomC = POS_LARGE_VALUE;
  FOR_NDIM(k) {
	CONT(k == dir2D);
	if (dimDomC[k] < minDimDomC) { minDimDomC = dimDomC[k]; }
  }
  FOR_NDIM(k) {
	nGrid[k] = nCell[k] + 1; 
	P2A(bnd.r,0,k,NDIM) = dimDomC[k] * iCell[k];
	P2A(bnd.r,1,k,NDIM) = dimDomC[k] * (iCell[k] + 1.);
	if (dimDomC[k] < 2. * neiEdge) {
		Printf0("Error: one of the subdomain widths is too small. It seems "
				"that too many CPUs are used, or the way to divide a "
				"domain is inappropriate. Check the 'parallel' file.\n");
		exit(-1);
	}
  }
  // modeActCh decides which way is more efficient between sending a whole
  // actin chain and sending only filled parts..
  modeActCh = ((double)nAbp / (double)((nChAc - 2) * nAct) < 0.5) ? 1 : 0; 
  /*--------------------------------------------------------------------------*/

  /*---- Dynamic behaviors of actin, ACP, motor, boundaries, and membranes ---*/

  // Put zeros to counters.
  if (actNuc.gTglFN != 0) { actNuc.cntFNme = 0; }
  actNuc.cntMe = 0;
  actAss.cntMe = 0;
  actDis.cntMe = 0;
  actBst.cntMe = 0;
  actSev.cntMe = 0;
  actAnn.cntMe = 0;
  actCap.cntMe = 0;
  actUnc.cntMe = 0;

  acpUnb.cntMe = 0;
  acpInaUnb.cntMe = 0;
  acpReb.cntMe = 0;
  acpMoBind.cntMe = 0;

  motUnb.cntMe = 0;
  motInaUnb.cntMe = 0;
  motReb.cntMe = 0;
  motMoBind.cntMe = 0;
  motWalk.cntMe = 0;
  motSA.cntNucMe = 0;
  motSA.cntAssMe = 0;
  motSA.cntTurnMe = 0;

  mbUnb.cntMe = 0;
  mbReb.cntMe = 0;
  mbFix.cntMe = 0;
  mbDef.cntMe = 0;
  mbPro.cntMe = 0;
  if (mbSld.act.gTgl != 0) { mbSld.act.cntMe = 0; }
  if (mbSld.abp.gTgl != 0) { mbSld.abp.cntMe = 0; }

  nucUnb.cntMe = 0;
  nucReb.cntMe = 0;

  noAbpDyn.c = 0;
  noActDyn.c = 0;
  noMbDyn.c = 0;

  // How long an event is prohibited after a previous event.
  // This reflects time required for a protein to change a state.
  durNoAcpUnbReb = T_SEC2TS(DUR_ACP_NO_UNB);
  durNoMotUnbReb = T_SEC2TS(DUR_MOT_NO_UNB);
  durNoMotWalk = T_SEC2TS(DUR_MOT_NO_WALK);
  durNoActDyn = T_SEC2TS(DUR_ACT_NO_DYN);
  durNoMbDyn = T_SEC2TS(DUR_MEMB_NO_DYN);

  actDgd.dist = L_NM2S(actDgd.dist);
  if (actBst.gTgl != 0) {
	actBst.fil.c = 0;
  }

  actAnn.ang = DEG2RAD(actAnn.ang);

  if (bndMat.gTgl != 0 || bndUnb.gTgl != 0) {
	bndMat.maxF = (int)(F_N2S(MAX_ACT_MAT_FORCE));
	bndMat.maxF++;
  }

  // Binding and binding of ACP and motor
  acpReb.p = K2P(K_ACP_BIND * acpReb.facK);
  motReb.p = K2P(40. * motMC.nHead * motReb.facK);
  
  // Binding of monomers of ACPs and motors
  acpMoBind.facP = FAC_P(K_ACP_BIND * acpMoBind.facK);
  motMoBind.facP = FAC_P(40 * motMC.nHead * motMoBind.facK);
  // Unbinding of ACP
  if (acpUnb.gTgl != 0 || acpInaUnb.gTgl != 0) {
	acpUnb.maxF = (int)(F_N2S(MAX_ACP_UNB_FORCE));
	acpUnb.maxF++;
  }
  // Self-assembly of motor
  if (motSA.gTgl != 0) {
	motSA.nNucMe = nAbpGoalDet[2] / (motSA.nMotPerSide * 2);
  }
  if (motSA.to.gTgl != 0 && motSA.to.mode == 1) {
	motSA.to.maxF = (int)(F_N2S(MAX_MOT_TURN_FORCE));
	motSA.to.maxF++;
  }
  CheckPeriodToggleParameter(&recAbpUnb);
  if (recAbpUnb.tgl != 0) {
	unbLog.c = 0;  
  }
  CheckPeriodToggleParameter(&recAbpBind);
  if (recAbpBind.tgl != 0) {
	bindLog.c = 0;  
  }
  CheckPeriodToggleParameter(&recActSev);
  if (recActSev.tgl != 0) {
	sevLog.c = 0;  
  }
  CheckPeriodToggleParameter(&recAbpTurn);
  if (recAbpTurn.tgl != 0) {
	toLog.c = 0;  
  }
  CheckPeriodToggleParameter(&updSubdSize);
  // Unbinding of filaments from boundaries
  CheckPeriodToggleParameter(&recBndUnbReb);
  if (bndUnb.gTgl != 0) {
	bndUnb.maxF = (int)(F_N2S(MAX_BND_UNB_FORCE));
	bndUnb.maxF++;
  }
  bndReb.dep = L_NM2S(bndReb.depR);

  // Moving boundaries
  CheckPeriodToggleParameter(&recBndLoc);
  // Maturation on boundaries
  CheckPeriodToggleParameter(&recBndActMat);
  // Traction forces
  CheckPeriodToggleParameter(&recBndTracF);

  // Membrane
  if (gTglMb != 0) {
	if (mbReb.gTgl != 0) {
		mbReb.p = K2P(K_MEMB_REB * mbReb.facK);
	}
	if (mbFix.gTgl != 0) {
		mbFix.p = K2P(K_MEMB_FIX * mbFix.facK);
	}
	if (mbUnb.gTgl != 0) {
		mbUnb.maxF = (int)(F_N2S(MAX_MEMB_UNB_FORCE));
		mbUnb.maxF++;
	}
	if (mbDef.gTgl != 0) {
		mbDef.maxF = (int)(F_N2S(MAX_MEMB_DEF_FORCE));
		mbDef.maxF++;
	}
	if (mbMat.gTgl != 0) {
		mbMat.maxF = (int)(F_N2S(MAX_MEMB_MAT_FORCE));
		mbMat.maxF++;
	}
	if (mbPro.gTgl != 0) {
		mbPro.p = K2P(K_MEMB_PROT * mbPro.facK);
		mbPro.dur = T_SEC2TS(mbPro.durR);
	}
	mbPro.f = F_PN2S(mbPro.f);
	if (mbSld.gTgl != 0) {
		mbSld.prd = (updSubdSize.prd > 0) ? updSubdSize.prd : recProg.prd;
		mbSld.critDist = L_NM2S(mbSld.critDist);
		mbSld.eqDist = L_NM2S(mbSld.eqDist);
		mbSld.stf = KS_NPM2S(mbSld.stf);
		mbSld.act.nReb = (mbSld.act.gTgl != 0) 
				? (int)(mbSld.act.por * nActGoal) : 0;

		mbSld.abp.nReb = (mbSld.abp.gTgl != 0) 
				? (int)(mbSld.abp.por * nAbpGoalDet[1]) : 0;
	}
  }
  if (gTglNuc != 0) {
	// Nuclues
	if (nucReb.gTgl != 0) {
		nucReb.p = K2P(K_NUC_REB * nucReb.facK);
	}
	if (nucUnb.gTgl != 0) {
		nucUnb.maxF = (int)(F_N2S(MAX_NUC_UNB_FORCE));
		nucUnb.maxF++;
	}
  }
  /*--------------------------------------------------------------------------*/

  /*------------------------- Rheological measurement ------------------------*/

  // If bulk rheology
  if (rheoWay > 0) {
	CheckPeriodToggleParameter(&recStre);
	sinuStr.prd = T_SEC2TS(sinuStr.prdR);
	if (bulkRheoWay == 0) { 
		prdUpdSinuStre = T_SEC2TS(PRD_UPD_SINU_STRESS);
		if (prdUpdSinuStre == 0) { prdUpdSinuStre = 1; }
	}
	V3COPY(pbcCond, pbc);
	if (bulkRheoType == 0) {
		dirOther = 3 - dirStr - dirNoPBC;
		pbc[dirStr] = 1;
		pbc[dirOther] = 1;
		pbc[dirNoPBC] = 0;
	}
	else {
		pbc[dirStr] = 0;
		pbc[(dirStr + 1) % NDIM] = 1;
		pbc[(dirStr + 2) % NDIM] = 1;
		if (dir2D > -1) {
			pbc[dir2D] = 0;
		}
	}
	FOR_NDIM(k) {
		CONT(pbcCond[k] == pbc[k]);
		Printf0("Warning: periodic boundary condition in the condition file "
				"is ignored for bulk rheology!\n\n");
		break;
	}
	if (pres.mag == 0.) { 
		pres.magDsp = 0.; 
		pres.dur = 0; 
		if (bulkRheoWay == 0) { pres.tgl = 0; }
	}
	else { 
		if (bulkRheoWay == 1) {
			pres.magDsp = pres.rate * dimDom[(bulkRheoType == 0) 
					? dirNoPBC : dirStr] * dtReal; 
			pres.dur = T_SEC2TS(pres.mag / pres.rate);
		}
		else {
			pres.magDsp = pres.rate * dtReal;
			pres.dur = (int)POS_LARGE_VALUE; 
			pres.tgl = 1;
		}
	}
	sumCos = 0.;
	for(n = 0; n < sinuStr.prd / 4; n++) {
		sumCos += cos(PI * 0.5 / ((double)sinuStr.prd * 0.25) * (double)n);
	}
	sinuStr.magDsp = (bulkRheoWay == 1 ? dimDom[(bulkRheoType == 0) ? 
			dirNoPBC : dirStr] : 1.) * sinuStr.amp / sumCos;
  }
  else {
	pres.dur = 0;
  }

  if (rheoWay == 0 || rheoWay == 2) {
	// Tracking trajectory
	if(recTraj.gTglCho != 0) {
		if (gTglMb == 0) {
			recTraj.nActL = (int)((int)(nActGoal * recTraj.porActL) / nCpu) 
					* nCpu;
			recTraj.nAbpL = (int)((int)(nAbpGoal * recTraj.porAbpL) / nCpu) 
					* nCpu;
		}
		else {
			recTraj.nActL = (int)(nActGoal * recTraj.porActL);
			recTraj.nAbpL = (int)(nAbpGoal * recTraj.porAbpL);
		}
		if (recTraj.nActL > 0) {
			recTraj.tgl = 1;
		}
		else {
			recTraj.tgl = 0;
			Printf0("Warning: there are not a sufficient number of actins "
					"to track! No actin will be traced\n\n");
		}
		if (recTraj.nAbpL > 0) {
			recTraj.tgl2 = 1;
		}
		else {
			recTraj.tgl2 = 0;
			Printf0("Warning: there are not a sufficient number of ABPs "
					"to track! No ABP will be traced\n\n");
		}
	}
	recTraj.prd = T_SEC2TS(recTraj.prdR);
  }
  else {
	recTraj.tgl = 0;
	recTraj.tgl2 = 0;
  }
  /*--------------------------------------------------------------------------*/

  /*--------------------------- Neighboring list  ----------------------------*/
  // Find the maximum size of particles
  abpF.maxSiz = NEG_LARGE_VALUE;
  for(k = 0; k < NK_ABP; k++) {
	if (abpF.maxSiz < abpF.len[k].n * 2.) { abpF.maxSiz = abpF.len[k].n * 2.; }
  }
  abpF.maxSiz = (abpF.maxSiz > 1.) ? abpF.maxSiz : 1.;
  dispHeuSq = pow(0.5 * (DT_NL_UPDATE - DT_NL_UPDATE_BUF), 2);

  neiEdge = 2.0;
  FOR_NDIM(k) {
	if (nCell[k] > 1) {
		if (iCell[k] == 0 && pbc[k] == 0) {
			P2A(edge,0,k,NDIM) = 0.;		
			P2A(edge,1,k,NDIM) = neiEdge;
		}
		else if (iCell[k] == nCell[k] - 1 && pbc[k] == 0) {
			P2A(edge,0,k,NDIM) = neiEdge;		
			P2A(edge,1,k,NDIM) = 0.;
		}
		else { 
			P2A(edge,0,k,NDIM) = neiEdge;		
			P2A(edge,1,k,NDIM) = neiEdge;
		}
	}
	else { 
		P2A(edge,0,k,NDIM) = 0.;	
		P2A(edge,1,k,NDIM) = 0.;
	}
	cell.n[k] = (int)((dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM)) 
			/ ((abpF.maxSiz + DT_NL_UPDATE) * 1.2));
	if (cell.n[k] < 1) { cell.n[k] = 1; }
	cell.wid[k] = (dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM)) 
			/ (double)cell.n[k]; 
	cell.base[k] = P2A(bnd.r,0,k,NDIM) - P2A(edge,0,k,NDIM);
	neiPbc[k] = (nCell[k] == 1 && pbc[k] != 0) ? 1 : 0;
	if (cell.n[k] == 1) { neiPbc[k] = 0; }
  }
  /*--------------------------------------------------------------------------*/

  /*------------------ Neighboring list related to membrane  -----------------*/
  if (gTglMb != 0) {
	mbMaxSiz = (abpF.maxSiz > memb.len) ? abpF.maxSiz : memb.len;
	mbMaxSiz = (memb.spr.hi > mbMaxSiz) ? memb.spr.hi : mbMaxSiz;
	mbMaxSiz = (memb.nucSpr.hi > mbMaxSiz) ? memb.nucSpr.hi : mbMaxSiz;
	FOR_NDIM(k) {
		cellMb.n[k] = (k == dirNormMbNuc && dimMbNuc == 2) 
				? 1 : (int)((dimDomC[k] + P2A(edge,0,k,NDIM) 
				+ P2A(edge,1,k,NDIM)) / ((mbMaxSiz + DT_NL_UPDATE) * 1.2));
		if (cellMb.n[k] < 1) { cellMb.n[k] = 1; }
		cellMb.wid[k] = (dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM)) 
				/ (double)cellMb.n[k]; 
		cellMb.base[k] = cell.base[k];
		neiPbcMb[k] = (nCell[k] == 1 && pbc[k] != 0) ? 1 : 0;
		if (cellMb.n[k] == 1) { neiPbcMb[k] = 0; }
	}
  }
  /*--------------------------------------------------------------------------*/

  /*-------------------------- Parallel processing ---------------------------*/
  // It is very important to set these two numbers at proper levels. If they 
  // are too large, many arrays are defined with unnecessarily large size.
  // On the contrary, if they are too small, the accumulation of particles in
  // one CPU can cause errors due to the access to out of range.
  // Without any motor, networks relatively maintain homogenous morphology, so 
  // they don't need to be as large as nAct or nAbp. 
  if (nCpu == 1) {
	nActC = nActGoal;
	nAbpC = nAbpGoal;
  }
  else {
  	if (motWalk.gTgl == 0 || nAbpGoalDet[2] == 0) {
		ratio = 1.;
		FOR_NDIM(k) {
			if (nCell[k] > 1) {
				ratio *= (dimDomC[k] + 2. * neiEdge) / dimDomC[k];
			}
		}
		nActC = (int)((double)nActGoal / nCpu * ratio * 2.0);
		nAbpC = (int)((double)nAbpGoal / nCpu * ratio * 2.0);
	}
	else {
		nActC = nActGoal;
		nAbpC = nAbpGoal;
	}
  }
  nActC = nActGoal * 2;
  nAbpC = nAbpGoal * 2;
  nActMin = (nActC < nActGoal) ? nActC : nActGoal;
  nAbpMin = (nAbpC < nAbpGoal) ? nAbpC : nAbpGoal;

  // adjRank: the first column of this array has the absolute CPU index of 
  //          adjacent subdomains. (Some rows of this array might have the 
  //          value of -1 the subdomain is located at a certain boundary.) 
  //          The second column has the relative CPU index which starts from 0 
  //          and increases by 1.
  // iRank  : this array has information about the absolute CPU index 
  //          corresponding to relative CPU index.
  sizeRankArr = (mpiMethod == 0) ? 27 : NDIM * 2;
  MALLOC(adjRank,int,sizeRankArr * 2);
  MALLOC(iRank,int,sizeRankArr);
  MALLOC(cntAdjRank,int,(mpiMethod == 0) ? 1 : NDIM);
  memset(adjRank, -1, sizeof(int) * sizeRankArr * 2);
  memset(iRank, -1, sizeof(int) * sizeRankArr);
  memset(cntAdjRank, 0, sizeof(int) * ((mpiMethod == 0) ? 1 : NDIM));

  // Normal method 
  // In this method, the subdomain communicates with all adjacent subdomains
  // for moving and copying the information of particles. In the worst case,
  // the subdomain located at the mid of a large domain needs 26 communications.
  // Communications between subdomains mean those between CPUs, which are 
  // time-lagging processes compared to internal computation in a CPU.
  if (mpiMethod == 0) {
	MALLOC(chkCpu, int, nCpu);
	memset(chkCpu, -1, sizeof(int) * nCpu);
	FOR_NDIM(oft[0]) {
	    CONT(pbc[0] == 0 && ((iCell[0] == 0 && oft[0] == 0) ||
	            (iCell[0] == nCell[0] - 1 && oft[0] == 2))); 
	    FOR_NDIM(oft[1]) { 
	        CONT(pbc[1] == 0 && ((iCell[1] == 0 && oft[1] == 0) ||
	                (iCell[1] == nCell[1] - 1 && oft[1] == 2)));
	        FOR_NDIM(oft[2]) { 
	            CONT(pbc[2] == 0 && ((iCell[2] == 0 && oft[2] == 0) ||
	                    (iCell[2] == nCell[2] - 1 && oft[2] == 2)));
	            CONT(oft[0] == 1 && oft[1] == 1 && oft[2] == 1);
	            V3ADD(iAdjCell, iCell, oft);
	            VVSS3SUB(iAdjCell, 1);
	            FOR_NDIM(k) {
	                if (iAdjCell[k] == nCell[k]) {
	                    iAdjCell[k] = 0;
	                }
	                else if (iAdjCell[k] == -1 ) {
	                    iAdjCell[k] = nCell[k] - 1;
	                }
	            }
	            V3IND_BACK_INT(adjCpu, iAdjCell, nCell);
	            ind = ((oft[0] * NDIM) + oft[1]) * NDIM + oft[2];
	            if (adjCpu != rank) {
	                P2A(adjRank,ind,0,2) = adjCpu;
	                if (chkCpu[adjCpu] < 0) {
	                    chkCpu[adjCpu] = cntAdjRank[0];
	                    P2A(adjRank,ind,1,2) = cntAdjRank[0];
	                    iRank[cntAdjRank[0]] = adjCpu;
	                    cntAdjRank[0]++;
	                }
	                else {
	                    P2A(adjRank,ind,1,2) = chkCpu[adjCpu];
	                }
				}
	        }
	    }
	}
	free(chkCpu);
  }
  // Plympton method
  // This method is very effective in a large domain. In this method, the
  // subdomain communicates only with 6 adjacent subdomains. Let say that 
  // all subdomains communicate with adjacent subdomains located in +x and -x
  // directions. After that, the communication occurs in +y and -y direction, 
  // and then it does in +z and -z direction. As a result, all subdomains 
  // have information from all adjacent subdomains. While the amount of 
  // communication can be reduced, the internal computation increases, but
  // its cost is much less than communication between CPUs.
  else {	
	FOR_NDIM(n) {
		for(k = 0; k < 2; k++) {
			V3SET_ALL(oft,0);
			oft[n] += (k == 0) ? -1 : 1;
			V3ADD(iAdjCell, iCell, oft);
			if (iAdjCell[n] == nCell[n]) { 
				iAdjCell[n] = (pbc[n] == 0) ? nCell[n] - 1 : 0; 
			}
			else if (iAdjCell[n] == -1) {
				iAdjCell[n] = (pbc[n] == 0) ? 0 : nCell[n] - 1; 
			}
	
			V3IND_BACK_INT(adjCpu, iAdjCell, nCell);
			ind = n * 2 + k;
			if (adjCpu != rank) {
				P2A(adjRank,ind,0,2) = adjCpu;
				if (k == 1) {
					if (P2A(adjRank,ind,0,2) == P2A(adjRank,ind - 1,0,2)) { 
						cntAdjRank[n]--; 
					}
				}
				P2A(adjRank,ind,1,2) = cntAdjRank[n];
				iRank[n * 2 + cntAdjRank[n]] = adjCpu;
				cntAdjRank[n]++;
			}
		}
	  }
  }

  // Determine the size of the buffer for messages, depending on the number
  // of particles and measurement method
  minDim = POS_LARGE_VALUE; 
  FOR_NDIM(k) {
	if (dimDomC[k] < minDim) { minDim = dimDomC[k]; }
  }
  sizeBufMsg = nActC * (sizeof(double) * 3 + sizeof(int) 
		* ((nChAc > 10 ? 10 : nChAc) + 3)); 
  if (gTglActCapAll != 0) { sizeBufMsg += nActC * sizeof(int); }
  if (actAge.gTgl != 0) { sizeBufMsg += nActC * sizeof(int); }
  if (gTglMb != 0) { sizeBufMsg += nActC * sizeof(int); }
  if (confVmdInfo.tgl != 0) {
	  sizeBufMsg += nActC * (sizeof(double) * (4 + NDIM) + sizeof(int)); 
  }
  sizeBufMsg += nAbpC * (sizeof(double) * (NDIM * 2 + 4) 
		+ sizeof(int) * (1 + ((motSA.gTgl != 0) ? nChAb : nChAb - 2))); 
  if (motSA.gTgl != 0) { sizeBufMsg += nAbpC * sizeof(int); }
  if (abpAge.gTgl != 0) { sizeBufMsg += nAbpC * sizeof(int); }
  if (recAbpDyn.tgl != 0) { sizeBufMsg += nAbpC * sizeof(int) * 3; }
  if (recTraj.tgl2 != 0) { sizeBufMsg += nAbpC * sizeof(int); }
  if (confVmdInfo.tgl != 0) {
	sizeBufMsg += nAbpC * (sizeof(double) * (3 + NDIM + recAbp.nL) 
			+ sizeof(int)); 
  }
  if (recLongF.tgl != 0) { sizeBufMsg += nAbpC * sizeof(double) * 4; }
  if (gTglMb != 0) {
	sizeBufMsg += nMbC * (sizeof(int) * (nChMb + 2) + sizeof(double) * 12);
	if (dimMbNuc == 3) {
		sizeBufMsg += nMbC * sizeof(double) * (nChMb - nMbAct);
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		sizeBufMsg += nMbC * (sizeof(int) + sizeof(double) * NDIM);
	}
	if (mbReb.gTgl != 0 || mbUnb.gTgl != 0 
			|| nucReb.gTgl != 0 || nucUnb.gTgl != 0) {
		sizeBufMsg += nMbC * sizeof(int);
	}
	if (mbPro.gTgl != 0) {
		sizeBufMsg += nMbC * (sizeof(int) * 2 + sizeof(double) * 2);
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0 || mbReb.gTgl != 0 
				|| mbUnb.gTgl != 0) {
		sizeBufMsg += nMbC * sizeof(int) * nMbAct;
	}
  }
  if (rheoWay > 0) {
	sizeBufMsg += (nActC + nAbpC) * sizeof(double) * 4;
  }
  ratio = 2.;
  sizeBufMsg *= (int)ratio;
  /*--------------------------------------------------------------------------*/

  /*------------------------------ Record data -------------------------------*/
  // general data
  CheckPeriodToggleParameter(&recProg);
  CheckPeriodToggleParameter(&recConf);
  CheckPeriodToggleParameter(&recE);
  CheckPeriodToggleParameter(&recCrsDist);
  CheckPeriodToggleParameter(&recMotSize);
  CheckPeriodToggleParameter(&recMotPos);
  if (motSA.gTgl == 0) { 
	recMotSize.tgl = 0; 
	recMotPos.tgl = 0; 
  }
  CheckPeriodToggleParameter(&recConn);
  CheckPeriodToggleParameter(&recPoreSize);
  CheckPeriodToggleParameter(&recAbpDyn);
  // membrane-related
  CheckPeriodToggleParameter(&recMbCen);
  CheckPeriodToggleParameter(&recMbDim);
  // networks via VMD
  CheckPeriodToggleParameter(&findSupp);
  CheckPeriodToggleParameter(&recConfMlb);
  if (recConfVmd.prdR < 0) {
	recConfVmd.tgl = 0;
	memset(&recConfVmd.prd, 1, sizeof(int));
  }
  else {
	recConfVmd.tgl = 1;
	recConfVmd.prd = T_SEC2TS(recConfVmd.prdR);
  }
  if (recConfPV.prdR < 0) {
	recConfPV.tgl = 0;
	memset(&recConfPV.prd, 1, sizeof(int));
  }
  else {
	recConfPV.tgl = 1;
	recConfPV.prd = T_SEC2TS(recConfPV.prdR);
  }

  if (pbc[dirRecPerc] != 0 && recPerc.tgl != 0) {
	Printf0("Warning: periodic boundary condition exists in the direction of "
			"percolation measurement, so percolation won't be measured.\n\n");
	recPerc.tgl = 0;
	memset(&recPerc.prd, 1, sizeof(int));
  }
  else {
	CheckPeriodToggleParameter(&recPerc);
  }
  tglRecAbpTurn = (recAbpTurn.tgl != 0 && (acpUnb.gTgl != 0 
			|| motUnb.gTgl != 0)) ? 1 : 0;

  CheckPeriodToggleParameter(&recInfo);
  if (recInfo.prd != recConfVmd.prd && recInfo.tgl != 0 
		&& recConfVmd.tgl != 0 && recConfVmd.gTglInfo != 0) {
	Printf0("Error: period of recording filament-based information should be "
			"the same as that of recording a network via VMD if both are "
			"needed at the same time.");
	exit(-1);
  }
  confVmdInfo.tgl = ((recConfVmd.tgl != 0 && recConfVmd.gTglInfo != 0) 
			|| recInfo.tgl != 0) ? 1: 0;
  confVmdInfo.prd = (recInfo.tgl != 0) ? recInfo.prd : recConfVmd.prd;
  CheckPeriodToggleParameter(&recSecStre);
  CheckPeriodToggleParameter(&recFilaL);
  CheckPeriodToggleParameter(&recLongF);
  // Number of lengths to record for ABPs
  recAbp.nL = (motSA.gTgl != 0) ? 2 : 1;
  /*--------------------------------------------------------------------------*/
  // Variables related to error process
  stopSig = 0;  
  magUnstF = F_N2S(MAG_UNSTABLE_FORCE);
  // Variables related to chain list
  actM.c = 0;
  acpM.c = 0;
  motM.c = 0;
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) 
  { iFilaP.c = 0i; }
  if (actSev.gTgl != 0) {
    if (actSev.critF > 0.) {
        actSev.critF = F_PN2S(actSev.critF);
    }
  }
  sendActDyn.c = 0;
  sendAbpDyn.c = 0;
  sendMbDyn.c = 0;
}

// Initialize data directory and record files required during the process.
// If the data directory doesn't exist, create it.
// If the data directory already exists, and if TGL_DELETE_FILE = 1, 
// remove all files in "fileL".
void InitRecFiles(void) {
  int n, cntFileL;
  char fileL[][30] = {"AbpLongFor", "AbpTurnover", "AbpUnb", "AbpDyn", 
		"AbpInstLenS", "AbpInstLenD", "AbpTraj", "AccuConf", 
		"AllChainList", "AllInstFor", "AllInstOri", "ActInstLenS", 
		"ActInstLenD", "ActFilaLen", "ActSev", "ActTraj", 
		"BeadBind", "BeadForce", "BeadPreforce", "BeadPrePos", "BeadPos", 
		"BndActMat", "BndLoc", "BndDiv", "BndUnbReb", "BunFilaEnd", 
		"ConfVmd.dat", "ConfVmd.psf", "ConfVmd.pdb", "Connect", 
		"CrossAng", "CrossDist", "ElasViscStre", 
		"InfoIndv", "InfoFila", "InfoFilaSeg", 
		"MechEall", "MechEperc", "MechEsupp", 
		"MembAdhes", "MembCen", "MembDim", "MembInstLenD", 
		"MotFilaLen", "MotUnbRate", "MotWalkRate", 
		"Parameter", "PercFila", "PercLen",	"PoreSize", 
		"Prestress", "Progress", "Stress", "Variable", OUTPUT_FILE};

  cntFileL = (int)sizeof(fileL) / (int)sizeof(fileL[0]);
  if (access(dataFold, 0) == -1) {
	mkdir(dataFold, 0755);
  }
  else if (DELETE_FILE != 0) {
	for(n = 0; n < cntFileL; n++) {
		CONT(!(access(GenFileName((const char *)fileL[n]), 0) == 0));
		remove(GenFileName((const char *)fileL[n]));
	}	
  }
}

void InitFileCheck(void) {
  int lenDataFold;
  FILE *fIn;

  // Check 'condition' and 'parallel'
  if (access("condition", 0) != 0) {
	Printf0("Error: the file 'condition' is missing!\n");
	exit(-1);
  }
  fIn = fopen("condition", "r");
  fscanf(fIn, "Directory for saving data = %s\n", dataFold);
  lenDataFold = (int)strlen(dataFold);
  if (dataFold[lenDataFold - 1] == '/' || dataFold[lenDataFold - 1] == '\\') {
	dataFold[lenDataFold - 1] = dataFold[lenDataFold];
  }
  fclose(fIn);

  // Initialize files for recording
  if (rank == 0) { InitRecFiles(); }

  if (access("parallel", 0) != 0) {
	Printf0("Error: the file 'parallel' is missing!\n");
	exit(-1);
  }
  if (gTglLoadNetData != 0 && access("Config", 0) != 0) {
	Printf0("Error: the file 'Config' is missing although the network data "
			"need to be loaded from the file!\n");
	exit(-1);
  }
}

// Feed the random seed from the timer information for generating 
// really random number sequences.
void InitRandSeed(void) {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  srand(tval.tv_usec / 1000 + rank);
  seed = (int)rand();
  // Initialize the seed for RNG
  init_genrand(seed); 
}

// The stochastic model is adopted from T. Erdmann and U.S. Schwarz, PRL, 2012.
// The rates of mechanochemical processes are defined in the 'condition' file.
// This function calculates the unbinding and walking rates of motors at 
// various forces.
void InitMotorWalkUnbindRates(int mode) {
  int m, n, k, nTotMotHead;
  double d, km, F0, Epp, Eel, f, gr, fac, vij, T10, *fac2, sumFac, *rate[2];
  long double *pji, *rij, *xij, *pi_inf, *g, *r, vb, unbRate, sum;
  FILE *fOut;
  MotUnbWalk *pL;

  pL = (mode == 0) ? &motWalk : &motUnb;
  nTotMotHead = motMC.nHead;
  if (mode == 0 && motSA.gTgl != 0) { 
	nTotMotHead *= motSA.nMotPerTF / 2; 
  }
  d = L_SCALE_IN_M / (double)nChAcX;
  km = STF_MOT_SPR;
  F0 = MOT_STALLF;
  Epp = -60e-21;

  MALLOC(fac2, double, nTotMotHead + 1);
  MALLOC(g, long double, nTotMotHead + 1);
  MALLOC(r, long double, nTotMotHead + 1);
  MALLOC(pi_inf, long double, nTotMotHead + 1);
  MALLOC(xij, long double, SQR(nTotMotHead + 1));
  MALLOC(pji, long double, SQR(nTotMotHead + 1));
  MALLOC(rij, long double, SQR(nTotMotHead + 1));

  for(m = 0; m < 2; m++) {
	MALLOC(rate[m], double, (int)F_N2S(1e-8));
	pL->maxF[m] = 0;
	while(1) {
		f = F_S2N((double)pL->maxF[m]) * ((m == 0) ? -1. : 1.);
		if (mode == 0 && motSA.gTgl != 0) { 
			f *= (double)(motSA.nMotPerTF / 2); 
		}
		for(n = 0; n < SQR(nTotMotHead + 1); n++) {
			xij[n] = 0;
			pji[n] = 0;
			rij[n] = 0;
		}
		for(n = 0; n < nTotMotHead + 1; n++) {
			r[n] = 0.;
			g[n] = 0.;
			pi_inf[n] = 0.;
		}
		for(n = 1; n <= nTotMotHead; n++) {
			for(k = 0; k <= n; k++) {
				P2A(xij,n,k,nTotMotHead + 1) = (f - k * km * d) / (n * km);
				Eel = km * ((n - k) * SQR(P2A(xij,n,k,nTotMotHead + 1)) 
						+ k * SQR(P2A(xij,n,k,nTotMotHead + 1) + d)) / 2.;
				fac2[k] = -1 / KT_IN_J * (Eel + Epp * k);
			}
			sumFac = 0;
			for(k = 0; k <= n; k++) {
				sumFac += fac2[k];
			}
			sumFac /= (double)(n + 1);
			for(k = 0; k <= n; k++) {
				fac2[k] -= sumFac;
				P2A(pji,n,k,nTotMotHead + 1) = expl((long double)fac2[k]);
			}
			sum = 0;
			for(k = 0; k < nTotMotHead + 1; k++) {
				sum += P2A(pji,n,k,nTotMotHead + 1);
			}
			if (sum > 0) {
				for(k = 0; k < nTotMotHead + 1; k++) {
					P2A(pji,n,k,nTotMotHead + 1) /= sum;
				}
			}
		}
       
		for(n = 0; n <= nTotMotHead; n++) {
			g[n] = (nTotMotHead - n) * motMC.k01;
		}
		for(n = 1; n <= nTotMotHead; n++) {
			for(k = 0; k <= n; k++) {
				P2A(rij,n,k,nTotMotHead + 1) = (n - k) * motMC.k10 
						+ k * motMC.k20	* expl(-km * (P2A(xij,n,k,nTotMotHead 
						+ 1) + d) / F0);
				r[n] += P2A(pji,n,k,nTotMotHead + 1) 
						* P2A(rij,n,k,nTotMotHead + 1);
			}
		}
		if (mode == 1) {
			T10 = 0;
			for(n = 1; n <= nTotMotHead; n++) {
				gr = 1;
				for(k = 1; k <= n - 1; k++) {
					if (r[k] > 0) {
						gr *= g[k] / r[k];
					}	
				}
				if (r[n] > 0) {
					T10 += 1. / r[n] * gr;
				}
			}
			unbRate = 1. / T10;
			rate[m][pL->maxF[m]] = K2P(pL->facK0 * unbRate);
		}
		else {
			sum = 1.;
			for(n = 1; n <= nTotMotHead; n++) {
				fac = 1;
				for(k = 0; k <= n - 1; k++) {
					if (r[k + 1] > 0) {
						fac *= g[k] / r[k + 1];
					}
				}
				pi_inf[n] = (long double)fac;
				sum += (long double)fac;
			}
			for(n = 1; n <= nTotMotHead; n++) {
				pi_inf[n] /= sum;
			}
			vb = 0;        
			for(n = 1; n <= nTotMotHead; n++) {
				for(k = 0; k <= n; k++) {
					vij = -g[n] * P2A(xij,n,k,nTotMotHead + 1) / (n + 1);
					if (n == 1 && k <= 1) {
						vij -= ((k == 0) ? motMC.k10 : motMC.k20) 
								* P2A(xij,n,k,nTotMotHead + 1);
					}		
						vb += vij * P2A(pji,n,k,nTotMotHead + 1) * pi_inf[n];
				}
			}
			rate[m][pL->maxF[m]] = K2P(pL->facK0 * vb / d);
		}
		(pL->maxF[m])++;
		if (mode == 1 && m == 0) {
			BREAK(rate[m][pL->maxF[m] - 1] > 0.999999);
		}
		else {
			BREAK(rate[m][pL->maxF[m] - 1] < 1e-16);
		}
	
	}
  }
  MALLOC(pL->p,double,pL->maxF[0] + pL->maxF[1] - 1);
  for(m = 0; m < 2; m++) {
	for(n = 0; n < pL->maxF[m]; n++) {
		pL->p[pL->maxF[0] + n * ((m == 0) ? -1 : 1) - 1] = rate[m][n];
	}
  }
  if (rank == 0) {
	fOut = fopen(GenFileName((mode == 0) ? "MotWalkRate" : "MotUnbRate"), "w");
	for(n = 0; n < pL->maxF[0] + pL->maxF[1] - 1; n++) {
	    fprintf(fOut, "%g\t%g\t%g\n", F_S2PN((double)n - (pL->maxF[0] - 1)), 
				pL->p[n], P2K(pL->p[n]));
	}
	fclose(fOut);
  }
  free(g);
  free(r);
  free(pi_inf);
  free(xij);
  free(pji);
  free(rij);
  free(fac2);
  for(n = 0; n < 2; n++) { free(rate[n]); }
}

