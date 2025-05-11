// This file contains functions related to bulk and segment rheology

/*---------------------------- Stress and strain -----------------------------*/

// Calculate and apply stress/strain.
void ApplyStressStrain(void) {
  int n, CS;
  double streErr;

  if (stra.acc > stra.lim[1] || stra.acc < stra.lim[0]) 
  { exit(0); }
  // If strain-controlled
  if (bulkRheoWay == 1) {
	SELECT(pres.dur + netForm.dur + motActiv.dur, currTimeStep, stra.dsp, 
			pres.magDsp, sinuStr.magDsp * cos(PI * 0.5 / ((double)sinuStr.prd 
			* 0.25)	* (double)((currTimeStep - (pres.dur + netForm.dur 
			+ motActiv.dur)) % sinuStr.prd)));
  }
  // If stress-controlled
  else {
	stre.goal += (pres.tgl != 0) ? pres.magDsp : (sinuStr.magDsp 
			* cos(PI * 0.5 / ((double)sinuStr.prd * 0.25)
			* (double)((currTimeStep - (pres.dur + netForm.dur 
			+ motActiv.dur)) % sinuStr.prd)));
	streErr = stre.goal - (stre.curr * (double)signStr);
	stre.accShErr += streErr;	
	// PI feedback
	// KPI_STRESS_P and KPI_STRESS_I are determined by trial-and-error
	stra.dsp = KPI_STRESS_P * streErr + KPI_STRESS_I * stre.accShErr;
  }
  stra.dsp *= (double)signStr;
  // Apply strain to actin segments fixed on the top boundary.
  for (n = 0; n < appStraParMe.c; n++) {
	P2(act.r,iAct[appStraParMe.l[n]],dirStr) += stra.dsp;
  }
  stra.accDsp += stra.dsp;
  if (bulkRheoType == 0) {
	stra.acc = stra.accDsp / dimDom[dirNoPBC];
  }
  else {
	stra.acc = stra.accDsp / (P2A(rGridInit,1,dirStr,NDIM) 
			- P2A(rGridInit,0,dirStr,NDIM));
	rGrid[dirStr][nGrid[dirStr] - 1] += stra.dsp;
	CS = 0;
	if (iCell[dirStr] == nCell[dirStr] - 1) {
		P2A(bnd.r,1,dirStr,NDIM) = rGrid[dirStr][nGrid[dirStr] - 1];
		CS = 1;
	}
	dimDom[dirStr] = rGrid[dirStr][nGrid[dirStr] - 1] - rGrid[dirStr][0];
	dimDomH[dirStr] = dimDom[dirStr] * 0.5;
 
	if (CS != 0)  {
	    minDimDomC = POS_LARGE_VALUE;
	    FOR_NDIM(n) {
			CONT(n == dir2D);
	        if (P2A(bnd.r,1,n,NDIM) - P2A(bnd.r,0,n,NDIM) < minDimDomC) {
	            minDimDomC = P2A(bnd.r,1,n,NDIM) - P2A(bnd.r,0,n,NDIM);
	        }
	    }
	}
  }
}

/*---------------------------- Stress and strain -----------------------------*/

/*----------- Procedures necessary with or without bulk rheology  ------------*/

// Perform initial procedures for bulk rheology.
void PrepareBulkRheology(void) {
  int direc;
  char dirCh[] = "xyz";

  if (rheoWay == 2) {
	if (recTraj.gTglCho == 0 && recTraj.act.c == 0 && recTraj.abp.c == 0) {
		Printf0("Warning: there is no information of actins or ABPs whose "
				"trajectories are recorded in Config. Nothing will be "
				"traced!\n\n");
	}
  }
  SeverLongFilament();
  direc = (bulkRheoType == 0) ? dirNoPBC : dirStr;
  if (confPbc[direc] != 0) {
	// In one direction, filaments crossing boundaries are severed and clamped,
	// and the periodic boundary condition is deactivated.
	SeverFilaOnPlane(direc, 1);
  }
  else {
	Printf0("Warning: the loaded network has no periodic boundary "
			"condition in the %c direction!\n\n", dirCh[direc]);
	if (gTglLoadNetDataFix == 0 && bndReb.gTgl == 0) {
		Printf0("Error: in this case, the information of clamped actins should "
				"be loaded from Config for bulk rheology, or the binding "
				"of actins on boundaries should be allowed!\n\n");
		exit(-1);
	}
	else if (meaStrePar.c == 0 && appStraPar.c == 0 
			&& gTglLoadNetDataFix != 0 && bndReb.gTgl == 0) {
		Printf0("Error: there is no information of clamped actins in Config "
				"for bulk rheology! In this case, for bulk rheology, the "
				"binding of actins on boundaries should be allowed!\n\n");
		exit(-1);
	}
  }
}

void PrepareNotBulkRheology(void) {
  int n, k;
  double *pR;
  char dirCh[] = "xyz";

  if (rheoWay == 0) {
	if (recTraj.gTglCho == 0 && recTraj.act.c == 0 && recTraj.abp.c == 0) {
		Printf0("Warning: there is no information of actins or ABPs whose "
				"trajectories are recorded in Config. Nothing will be "
				"traced!\n\n");
	}
  }
  FOR_NDIM(k) {
	if (confPbc[k] != 0) {
		if (pbc[k] == 0) {
			for(n = 0; n < nAct + nAbp; n++) {
				pR = (n < nAct) ? &P2(rAct,n,k) : &P2(rAbp,n - nAct,k);
				if (*pR >= rGrid[k][nGrid[k] - 1]) { *pR -= dimDom[k]; }
				else if (*pR < rGrid[k][0]) { *pR += dimDom[k]; }
			}
			SeverFilaOnPlane(k, 1);
		}
	}
	else {
		Printf0("Warning: the loaded network doesn't have a periodic boundary "
				"condition in the %c direction!\n\n", dirCh[k]);
	}
  }
}

/*----------- Procedures necessary with or without bulk rheology  ------------*/

/*---------------------- Related to segment rheology -------------------------*/

int ChooseTrajectoryListSubroutine(int n, int ind) {
  int CS, *pArr;

  CS = 1;
  pArr = (n == 0) ? &P2A(act.ch,ind,0,nChAc) : &P2A(abp.ch,ind,0,nChAb);
  if (n == 1 && motSA.gTgl != 0) {
	if (K_ABP(ind) == 2) {
		if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) 
		{ CS = 0; }
	}
	else { 
		if (pArr[0] < 0 && pArr[1] < 0) { CS = 0; }
	}
  }
  else { 
	if (pArr[0] < 0 && pArr[1] < 0) { CS = 0; }
  }
  if (n == 0) {
	if (act.fix[ind] > -1) { CS = 0; }
  }
  return CS;
}

// Select actin segments for segment-rheology. Trajectory of these actin 
// segments is traced during the simulation to calculate viscoelastic moduli.
void ChooseTrajectoryList(void) {
  int n, k, ind, cnt, nGoal, nTrajL, nTrajLme, nMe, sum, CS;
  int *chkCho, *id, *cntAll;
  ListInt *pL;

  MALLOC(cntAll, int, nCpu);
  for(n = 0; n < 2; n++) {
	if (n == 0) {
		CONT(recTraj.tgl == 0);
		nTrajL = recTraj.nActL;
		nMe = nActMe;
		pL = &recTraj.actMe;
		id = act.id;
	}
	else {
		CONT(recTraj.tgl2 == 0);
		nTrajL = recTraj.nAbpL;
		nMe = nAbpMe;
		pL = &recTraj.abpMe;
		id = abp.id;
	}
	pL->c = 0;
	MALLOC(chkCho,int,nMe);
	memset(chkCho, -1, sizeof(int) * nMe);

	cnt = 0;	
	for(k = 0; k < nMe; k++) {
		cnt += ChooseTrajectoryListSubroutine(n, k);
	}
	MPI_Allgather(&cnt, 1, MPI_INT, cntAll, 1, MPI_INT, MPI_COMM_WORLD);
	sum = SumArrInt(cntAll, nCpu);
	nTrajLme  = (int)((double)nTrajL * cnt / sum);

	nGoal = (cnt < nTrajLme) ? cnt : nTrajLme;
	CONT(nGoal == 0);
	while(1) {
        ind = GenRandIntIndex(nMe);
		CS = ChooseTrajectoryListSubroutine(n, ind);
		CONT(CS == 0);
		// Check whether the actin segment was previously chosen
        CONT(chkCho[ind] > -1);
		pL->l[pL->c] = id[ind];
		chkCho[ind] = 1;
		(pL->c)++;
		BREAK(pL->c == nGoal)
	}
	free(chkCho);
  }
  free(cntAll);
}

// mode = 0: actin, 1: ABP
void ReplaceElementInTrajList(int ind, int mode) {
  int n, CS, ind2, nMe, nCh, *chkL, *iL, *ch, *id;
  ListInt *trajL, parL;

  if (mode == 0) {
	trajL = &recTraj.actMe;
	iL = iAct;
	ch = act.ch;
	id = act.id;
	nMe = nActMe;
	nCh = nChAc;
  }
  else {
	trajL = &recTraj.abpMe;
	iL = iAbp;
	ch = abp.ch;
	id = abp.id;
	nMe = nAbpMe;
	nCh = nChAb;
  }

  CS = FindElementArray(trajL->l, trajL->c, ind, 0, 1);
  if (CS > -1) {
	MALLOC(chkL,int,nMe);
	memset(chkL, -1, sizeof(int) * nMe);
	for(n = 0; n < trajL->c; n++) {
		chkL[iL[trajL->l[n]]] = 1;
	}
	MALLOC(parL.l,int,nMe);
	parL.c = 0;
	for(n = 0; n < nMe; n++) {
		CONT((P2A(ch,n,0,nCh) < 0 && P2A(ch,n,1,nCh) < 0) 
				|| chkL[n] != -1);
		if (mode == 0) {
			CONT(act.fix[n] > -1);
			if (gTglBead != 0 && beadBind.gTgl != 0) { 
				CONT(act.bdFix[n] > -1);
			}
			if (gTglMb != 0) { 
				CONT(act.mbReb[n] > -1);
			}
		}	
		InsertElement1dArrayWoChk(parL.l, &parL.c, n);
	}
	if (parL.c > 0) {
		ind2 = GenRandIntIndex(parL.c);
		trajL->l[CS] = id[parL.l[ind2]];
	}
	else {
		DeleteElementArrayByIndex(trajL->l, &trajL->c, CS, 1);
	}
	free(parL.l);
	free(chkL);
  }
}

/*---------------------- Related to segment rheology -------------------------*/

/*-------------------------------- Recording ---------------------------------*/

// Record shear stress/strain
// During the application of prestrain/prestress, the information is stored
// in "Prestress". After that, the information is stored in "Stress".
void RecordStress(int period) {
  int n, k, oftTimeStep, direc, locActInd, abpInd;
  double fSum[NDIM], flowStre, area;
  double *streAll, streSum[NDIM], *accShStreUSall;
  FILE *fOut;

  // Shear
  if (bulkRheoType == 0) {
	area = L_S2M(dimDom[dirStr]) * L_S2M(dimDom[dirOther]);
  }
  // Normal
  else {
	area = L_S2M(dimDom[(dirStr + 1) % NDIM]) 
			* L_S2M(dimDom[(dirStr + 2) % NDIM]);
  }
  V3ZERO(fSum);

  for(n = 0; n < meaStreParMe.c; n++) {
	locActInd = iAct[meaStreParMe.l[n]];
	VV3ADD(fSum, &P2(act.f,locActInd,0));
  }
  // Count repulsive forces
  direc = (bulkRheoType == 0) ? dirNoPBC : dirStr;
  fSum[direc] += P2(bnd.f,2 * direc + 1,direc);

  VV3ADD(stre.acc, fSum);
  stra.accDspFl += stra.dsp;

  // bulkRheoWay = 0: stress-controlled, 1: strain-controlled
  // If stress-controlled
  oftTimeStep = netForm.dur + motActiv.dur;
  if (bulkRheoWay == 0) {
	stre.accShUS += fSum[dirStr];
	stra.accDspUS += stra.dsp;

	// For the stress-controlled mechanism, the total shear stress of a whole 
	// network has to be calculated more often to maintain the shear stress at
	// desired level.
	if ((currTimeStep - oftTimeStep) % prdUpdSinuStre == 0) {
		if (rank == rankMeaStre.l[0]) 
		{ MALLOC(accShStreUSall,double,rankMeaStre.c); }
		// Convert the force to stress
		stre.accShUS = -1. * F_S2N(stre.accShUS) 
				/ (double)prdUpdSinuStre / area;
		// Collect them
		CollectArrayDblFromSubdomainList(&stre.accShUS, accShStreUSall, 1, 
				&rankMeaStre, gotMeaStre);
		// Gather the values of stre.accShUS to the main node.
		if (rank == rankMeaStre.l[0]) { 
			stre.curr = 0.;
			for(n = 0; n < rankMeaStre.c; n++) 
			{ stre.curr += accShStreUSall[n]; }
			if (bulkRheoType == 0) {
				flowStre = VISCOSITY * stra.accDspUS  
						/ ((double)prdUpdSinuStre * dt) / dimDom[dirNoPBC] 
						* KT_IN_J / actF.dragR / SQR(L_SCALE_IN_M);
				stre.curr += flowStre;
			}
			free(accShStreUSall);
		}
		// Broadcast the calculated shear stress to all
		MPI_Bcast(&stre.curr, 1, MPI_DOUBLE, rankMeaStre.l[0], 
				MPI_COMM_WORLD);
		// If the shear stress exceeds the aimed value, the application of 
		// prestress is ended.
		if (pres.tgl != 0 && stre.curr * (double)signStr >= pres.mag) { 
			stre.goal = pres.mag;
			pres.tgl = 0; 
			pres.dur = currTimeStep - netForm.dur - motActiv.dur;
		    V3ZERO(stre.acc);
		    stra.accDspFl = 0.;
		}
		stre.accShUS = 0.;
		stra.accDspUS = 0.;
	}
	oftTimeStep += (pres.tgl != 0) ? 0 : pres.dur;
  }
  else {
	oftTimeStep += (currTimeStep < netForm.dur + pres.dur + motActiv.dur) 
			? 0 : pres.dur;
  }
  // Shear stress and strain are recorded at the interval of "period".
  if ((currTimeStep - oftTimeStep) % period == 0 
		&& currTimeStep != oftTimeStep) {
	if (rank == rankMeaStre.l[0]) 
	{ MALLOC(streAll,double,NDIM*rankMeaStre.c); }
	for(k = 0; k < NDIM; k++) {
	    stre.acc[k] = -1. * F_S2N(stre.acc[k]) / (double)period / area;
	}
	// Gather local stress.
	CollectArrayDblFromSubdomainList(stre.acc, streAll, NDIM, 
			&rankMeaStre, gotMeaStre);
	if (rank == rankMeaStre.l[0]) {
		V3ZERO(streSum);
		for(n = 0; n < rankMeaStre.c; n++) 
		{ VV3ADD(streSum, &streAll[3 * n]); }
		if (bulkRheoType == 0) {
			flowStre = VISCOSITY * stra.accDspFl / ((double)period * dt)
		            / dimDom[dirNoPBC] * KT_IN_J / actF.dragR 
					/ SQR(L_SCALE_IN_M);
		    streSum[dirStr] += flowStre;
		}
		if (currTimeStep < netForm.dur + pres.dur + motActiv.dur) 
		{ fOut = fopen(GenFileName("Prestress"), "a");	} 
		else { fOut = fopen(GenFileName("Stress"), "a"); }
 		if (bulkRheoType == 0) {
		    fprintf(fOut, "%lld\t%g\t%g\t%g\t%g\t%g\n", currTimeStep,
  					stra.acc, streSum[dirStr], 
					streSum[dirNoPBC], streSum[dirOther], flowStre);
		}
		else {
		    fprintf(fOut, "%lld\t%g\t%g\t%g\t%g\n", currTimeStep, 
					stra.acc, streSum[dirStr], 
					streSum[(dirStr + 1) % NDIM], streSum[(dirStr + 2) % NDIM]);
		}
	    fclose(fOut);
		free(streAll);
	}
    V3ZERO(stre.acc);
    stra.accDspFl = 0.;
  }
}

void RecordTrajectorySubroutine(double *arr, int mode) {
  int n, ch, nD;
  ListInt *pL;

  pL = (mode == 0) ? &recTraj.actMe : &recTraj.abpMe;
  nD = NDIM + 2;
  for(n = 0; n < pL->c; n++) {
	P2A(arr,n,0,nD) = (double)pL->l[n];
	if (mode == 0) {
		V3COPY(&P2A(arr,n,1,nD), &P2(act.r,iAct[pL->l[n]],0));
		ch = BinaryPackActinChainArray(pL->l[n]);
	}
	else {
		V3COPY(&P2A(arr,n,1,nD), &P2(abp.r,iAbp[pL->l[n]],0));
		ch = BinaryPackAbpChainArray(pL->l[n]);
	}
	P2A(arr,n,nD - 1,nD) = (double)ch;
  }
}

// Record the trajectory of selected components.
// mode = 0: actin, 1: ABP
void RecordTrajectory(int period) {
  int m, n, k, nD, ind, nSeq, cntAct, cntAbp, tag = 0, trajMeC[2], *trajMeCAll;
  double *sendArr, *recvArr;
  ListDbl actTrajAll, abpTrajAll, *trajAll;
  char fn[80];
  FILE *fOut;

  nD = NDIM + 2;
  MALLOC(trajMeCAll,int,nCpu * 2); 
  V2SET(trajMeC, recTraj.actMe.c, recTraj.abpMe.c);
  MPI_Gather(trajMeC, 2, MPI_INT, trajMeCAll, 2, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) { 
	actTrajAll.c = 0;
	abpTrajAll.c = 0;
	for(n = 0; n < nCpu; n++) {
		actTrajAll.c += trajMeCAll[2 * n];
		abpTrajAll.c += trajMeCAll[2 * n + 1];
	}
	MALLOC(recvArr,double,(actTrajAll.c + abpTrajAll.c) * nD);
	MALLOC(actTrajAll.l,double,actTrajAll.c * nD);
	MALLOC(abpTrajAll.l,double,abpTrajAll.c * nD);
	RecordTrajectorySubroutine(actTrajAll.l, 0);
	RecordTrajectorySubroutine(abpTrajAll.l, 1);
	cntAct = trajMeCAll[0];
	cntAbp = trajMeCAll[1];
	for(n = 1; n < nCpu; n++) {
		CONT(trajMeCAll[2 * n] + trajMeCAll[2 * n + 1] == 0);
		MPI_Recv(recvArr, (trajMeCAll[2 * n] + trajMeCAll[2 * n + 1]) * nD, 
				MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &status);
		for(k = 0; k < trajMeCAll[2 * n]; k++) {
			V5COPY(&P2A(actTrajAll.l,cntAct + k,0,nD), &P2A(recvArr,k,0,nD));
		}
		for(k = 0; k < trajMeCAll[2 * n + 1]; k++) {
			V5COPY(&P2A(abpTrajAll.l,cntAbp + k,0,nD), 
					&P2A(recvArr,trajMeCAll[2 * n] + k,0,nD));
		}
		cntAct += trajMeCAll[2 * n];
		cntAbp += trajMeCAll[2 * n + 1];
	}
	qsort(actTrajAll.l, actTrajAll.c, nD * sizeof(double), CompDbl);
	qsort(abpTrajAll.l, abpTrajAll.c, nD * sizeof(double), CompDbl);

	nSeq = (currTimeStep - netForm.dur) / period;
	for(n = 0; n < 2; n++) {	
		if (n == 0) {
			trajAll = &actTrajAll;
			strcpy(fn, "ActTraj");
		}
		else {
			trajAll = &abpTrajAll;
			strcpy(fn, "AbpTraj");
		}
	    fOut = fopen(GenFileName(fn), "a");
		if (nSeq == 0) {
			fprintf(fOut, "%d\t", ((n == 0) ? actTrajAll.c + 2
					: abpTrajAll.c + 2));
			Fprintf1dFillerInt(fOut, 0, NDIM + 2, 0);
		}
		for(m = 0; m < 2; m++) {
			fprintf(fOut, "%d\t", m);
			FOR_NDIM(k) {
				ind = (m == 0) ? 0 : nGrid[k] - 1;
				fprintf(fOut, "%g\t", rGrid[k][ind]);
			}
			Fprintf1dFillerInt(fOut, 0, 5 - NDIM, 0);
		}
		Fprintf2dArrayDouble(fOut, trajAll->l, trajAll->c, nD, 2, 0);
		fclose(fOut);
	}
	free(recvArr);
	free(actTrajAll.l);
	free(abpTrajAll.l);
  }
  else if (rank != 0 && recTraj.actMe.c + recTraj.abpMe.c > 0) {
	MALLOC(sendArr,double,(recTraj.actMe.c + recTraj.abpMe.c) * nD);
	RecordTrajectorySubroutine(sendArr, 0);
	RecordTrajectorySubroutine(&P2A(sendArr,recTraj.actMe.c,0,nD), 1);
	MPI_Send(sendArr, nD * (recTraj.actMe.c + recTraj.abpMe.c), 
			MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	free(sendArr); 
  }
  free(trajMeCAll); 
}

/*-------------------------------- Recording ---------------------------------*/
