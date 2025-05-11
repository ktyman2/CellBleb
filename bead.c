// This file contains functions related to microbeads.

/*--------------------------- Calculating force ------------------------------*/

double CalcBeadRepulsiveForcesSubroutine(double len, double dist) {
	double f;

	if (len < (1. - bead.thkRep) * dist) {
		f = bead.stfRep * bead.thkRep * dist;
	}
	else {
		f = bead.stfRep * (dist - len);
	}
	return f / len;
}

void CalcBeadRepulsiveForces(void) {
  int m, n, k, l, ind[3];
  double rPnt[3][NDIM], *f_p[3], ratio[2], dr[NDIM], fi[NDIM];
  double len, dist, f, *pR;

  for(n = 0; n < act.cyl.c; n++) {
	ind[0] = iAct[P2A(act.cyl.l,n,0,2)];
	ind[1] = iAct[P2A(act.cyl.l,n,1,2)];
	CONT(ind[0] < 0 || ind[1] < 0);
	CONT(ind[0] >= nActMe && ind[1] >= nActMe);
	for(m = 0; m < bead.n; m++) {
		for(k = 0; k < 2; k++) {
			pR = (k == 0) ? &P2(bead.r,m,0) : rPnt[0];
			V3COPY(rPnt[k],&P2(act.r,ind[k],0));
			FOR_NDIM(l) {
				CONT(pbc[l] != 1);
				CONT(!(fabs(rPnt[k][l] - pR[l]) > 1.5 * dimDomH[l]));
				rPnt[k][l] += dimDom[l] 
						* ((rPnt[k][l] < pR[l]) ? 1. : -1.);
			}
			if (dir2D > -1) {
				rPnt[k][dir2D] = dimDomH[dir2D];
			}
		}
		dist = bead.rad[m] + 0.5 * actF.dia;
		len = CalcSegPntDist(&rPnt[0], &P2(bead.r,m,0), dr, ratio);
		CONT(!(len < dist));
		f = CalcBeadRepulsiveForcesSubroutine(len, dist);
	    VS3COPY(fi, dr, f);
		VVS3ADD(&P2(act.f,ind[0],0), fi, REVSIGN(1. - ratio[0]));
		VVS3ADD(&P2(act.f,ind[1],0), fi, REVSIGN(ratio[0]));
		CONT(ind[0] >= nActMe);
		VV3ADD(&P2(bead.f,m,0), fi);
	}
  }
  for(n = 0; n < nAbpMe; n++) {
	for(m = 0; m < bead.n; m++) {
		V3COPY(rPnt[0],&P2(abp.r,n,0));
		FOR_NDIM(k) {
			CONT(pbc[k] != 1);
			CONT(!(fabs(rPnt[0][k] - P2(bead.r,m,k)) > 1.5 * dimDomH[k]));
			rPnt[0][k] += dimDom[k] 
					* ((rPnt[0][k] < P2(bead.r,m,k)) ? 1. : -1.);
		}
		if (dir2D > -1) {
			rPnt[0][dir2D] = dimDomH[dir2D];
		}
		dist = bead.rad[m] + 0.5 * abpF.repDia[K_ABP(n)];
		if (dir2D < 0) {
			len = CalcVecDist(dr, &P2(bead.r,m,0), rPnt[0], 0);
		}
		else {
			CalcVec(dr, &P2(bead.r,m,0), rPnt[0]);
			dr[dir2D] = 0.;
			len = V3LEN(dr);
		}
		CONT(!(len < dist));
		f = CalcBeadRepulsiveForcesSubroutine(len, dist);
	    VS3COPY(fi, dr, f);
		VV3SUB(&P2(abp.f,n,0), fi);
		VV3ADD(&P2(bead.f,m,0), fi);
	}
  }
  if (gTglMb != 0) { 
	for(n = 0; n < memb.unit.c; n++) {
		for(k = 0; k < dimMbNuc; k++) {
			ind[k] = iMb[P2A(memb.unit.l,n,k,dimMbNuc)];
			BREAK(ind[k] < 0);
		}
		CONT(k != dimMbNuc);
		for(m = 0; m < bead.n; m++) {
			for(k = 0; k < dimMbNuc; k++) {
				V3COPY(rPnt[k],&P2(memb.r,ind[k],0));
				f_p[k] = &P2(memb.f,ind[k],0);
				FOR_NDIM(l) {
					CONT(pbc[l] != 1);
					CONT(!(fabs(rPnt[k][l] - P2(bead.r,m,l)) 
							> 1.5 * dimDomH[l]));
					rPnt[k][l] += dimDom[l] 
							* ((rPnt[k][l] < P2(bead.r,m,l)) ? 1. : -1.);
				}
				if (dir2D > -1) {
					rPnt[k][dir2D] = dimDomH[dir2D];
				}
			}
			dist = bead.rad[m] + 0.5 * ((ISNUC(memb.idx[ind[0]]))
					? memb.nucThk : memb.thk);
			if (dimMbNuc == 2) {
				len = CalcSegPntDist(&rPnt[0], &P2(bead.r,m,0), dr, ratio);
			}
			else {
				len = CalcTriaPntDist(&rPnt[0], &P2(bead.r,m,0), dr, ratio);
			}
			CONT(!(len < dist));
			f = CalcBeadRepulsiveForcesSubroutine(len, dist);
		    VS3COPY(fi, dr, f);
			if (ind[0] < nMbMe) {	
				VV3ADD(&P2(bead.f,m,0), fi);
			}
			if (dimMbNuc == 2) {
				VVS3ADD(&P2(memb.f,ind[0],0), fi, REVSIGN(1. - ratio[0]));
				VVS3ADD(&P2(memb.f,ind[1],0), fi, REVSIGN(ratio[0]));
			}
			else {
				V3REVSIGN(fi);
				DistributeForceOnTria(&rPnt[0], &f_p[0], ratio, fi);
			}
		}
	}
  }
}

void CalcBeadActinBindingForces(void) {
  int n, k;
  double dr[NDIM], f, fi, len, stfBeadBind;
	
  stfBeadBind = KS_NPM2S(1e-3);

  FOR_ACTME(n) {
	CONT(!(act.bdFix[n] > -1));
	CalcVec(dr, &P2(act.r,n,0), &P2(act.bdRFix,n,0));
	len = V3LEN(dr);
	CONT(len == 0.);
	f = SPRING(stfBeadBind, len, 0.);
	f /= len;
	FOR_NDIM(k) {
		fi = f * dr[k];
		P2(act.f,n,k) += fi;
		P2(bead.f,act.bdFix[n],k) -= fi;
	}
  }
}

/*--------------------------- Calculating force ------------------------------*/

/*-------------------------- updating information ----------------------------*/

void UpdateBeadNewLocation(void) {
  int n, k, ind;
  double rand[2], fBr[NDIM], *beadFall, *prevBeadR;

  MALLOC(beadFall,double,bead.n * NDIM * nCpu);
  MPI_Gather(bead.f, bead.n * NDIM, MPI_DOUBLE, beadFall, bead.n * NDIM, 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (gTglBead != 0 && beadBind.gTgl != 0) { 
	MALLOC(prevBeadR,double,bead.n * NDIM);
	Copy1dArrayDouble(prevBeadR, bead.r, bead.n * NDIM);
  }
  if (rank == 0) {
	for(n = 1; n < nCpu; n++) {
		for(k = 0; k < bead.n * NDIM; k++) {
			bead.f[k] += beadFall[(n * bead.n * NDIM) + k];
		}
	}
	for(n = 0; n < bead.n; n++) {
		CalcBrownForcesSubroutine(fBr, bead.maxFbr[n], n, 0, rand);
		FOR_NDIM(k) {
			CONT(k == dir2D);
			P2(bead.r,n,k) += (P2(bead.f,n,k) + fBr[k]) / bead.drag[n] * dt;
		}
		ApplyBoundCondVector(&P2(bead.r,n,0), 6 + n, 0);
	}
  }
  MPI_Bcast(bead.r, bead.n * NDIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (gTglBead != 0 && beadBind.gTgl != 0) { 
	FOR_ACTME(n) {
		CONT(!(act.bdFix[n] > -1));
		ind = act.bdFix[n];
		FOR_NDIM(k) {
			P2(act.bdRFix,n,k) += P2(bead.r,ind,k) - P2(prevBeadR,ind,k);
		}
	}
	free(prevBeadR);
  }
  free(beadFall);
}

void ApplyBeadSinuSignal(void) {
  int n, k;
  double mag, err, *dsp;

  if (bead.rheoSig == 0 && beadBind.gTgl != 0) {
	MALLOC(dsp,double,bead.n);
  }
  if (bead.rheoSig == 0) {
	mag = (bead.rheoPreTgl != 0 ) ? bead.rheoPreMagDsp 
			: bead.rheoAmp2 * cos(PI * 0.5 / ((double)bead.rheoPrd * 0.25) 
			* (double)((currTimeStep - (pres.dur + netForm.dur + motActiv.dur)) 
			% bead.rheoPrd));
	bead.rheoFGoal += mag;
	for(n = 0; n < bead.n; n++) {	
		err = bead.rheoFGoal - bead.rheoFCurr[n];
		bead.rheoFAccErr[n] += err;
		// PI feedback
		mag = KPI_FORCE_P * err + KPI_FORCE_I * bead.rheoFAccErr[n];
		P2(bead.r,n,bead.rheoDir) += mag;
		if (beadBind.gTgl != 0) {
			dsp[n] = mag;
		}
	}
  }
  else {
	SELECT(pres.dur + netForm.dur + motActiv.dur, currTimeStep, mag, 
			bead.rheoPreMagDsp, bead.rheoAmp2 * cos(PI * 0.5 
			/ ((double)bead.rheoPrd * 0.25)	* (double)((currTimeStep 
			- (pres.dur + netForm.dur + motActiv.dur)) % bead.rheoPrd)));
	for(n = 0; n < bead.n; n++) {
		P2(bead.r,n,bead.rheoDir) += mag;
	}
  }
  if (beadBind.gTgl != 0) {
	FOR_ACTME(n) {
		CONT(!(act.bdFix[n] > -1));
		P2(act.bdRFix,n,bead.rheoDir) += (bead.rheoSig == 0) 
				? dsp[act.bdFix[n]] : mag;
	}
  }
  if (bead.rheoSig == 0 && beadBind.gTgl != 0) {
	free(dsp);
  }
}

void UpdateBeadActinBinding(void) {
  int m, n, *pArr;
  double dist, dr[NDIM], len;

  for(m = 0; m < bead.n; m++) {
	dist = bead.rad[m] + 0.5 * actF.dia;
	FOR_ACTME(n) {
		CONT(act.bdFix[n] > -1);
		pArr = &P2A(act.ch,n,0,nChAc);
		CONT((pArr[0] > -1 && pArr[1] > -1) || (pArr[0] < 0 && pArr[1] < 0));
		CalcVec(dr, &P2(act.r,n,0), &P2(bead.r,m,0));
		if (dir2D > -1) {
			dr[dir2D] = 0.;
		}
		len = V3LEN(dr);

		CONT(!(len <= dist));
		V3SCALE(dr, dist / len);
		act.bdFix[n] = m;
		V3ADD(&P2(act.bdRFix,n,0), &P2(bead.r,m,0), dr);
		if (dir2D > -1) {
			P2(act.bdRFix,n,dir2D) = P2(act.r,n,dir2D);
		}
		beadBind.cntMe[m]++;
	}
  }
}

/*-------------------------- updating information ----------------------------*/

/*-------------------------- Recording information ---------------------------*/

void RecordBeadForceLocation(void) {
  int n, k, l, oftTimeStep;
  double *beadFall, *pD;
  FILE *fOut;

  oftTimeStep = currTimeStep - netForm.dur;
  // Record position  
  if (rank == 0 && oftTimeStep >= 0 && oftTimeStep % recBeadLoc.prd == 0) {
	if (currTimeStep < netForm.dur + pres.dur + motActiv.dur)
	{ fOut = fopen(GenFileName("BeadPrepos"), "a");  }
	else { fOut = fopen(GenFileName("BeadPos"), "a"); }
	if ((netForm.dur > 0 && oftTimeStep ==  0)
			|| (netForm.dur == 0 && currTimeStep == recBeadLoc.prd)) {
		for(n = 0; n < bead.n; n++) { 
			fprintf(fOut, "%g\t", L_S2M(bead.rad[n]));
		}
		Fprintf1dFillerInt(fOut, 0, bead.n * (NDIM - 1) + 1, 0);
		fprintf(fOut, "%g\t%g\t%g\t", L_SCALE_IN_M, dtReal, VISCOSITY);
		Fprintf1dFillerInt(fOut, 0, bead.n * NDIM - 2, 0);
		if (bead.rheoWay == 1) {
			fprintf(fOut, "%g\t%d\t", bead.rheoAmp, bead.rheoPrd);
			Fprintf1dFillerInt(fOut, 0, bead.n * NDIM - 1, 0);
		}
	}
	fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayDouble(fOut, bead.r, bead.n * NDIM, 0);
	fclose(fOut);
  }

  if (bead.rheoWay == 1) {
	for(n = 0; n < 2; n++) {
		for(k = 0; k < bead.n; k++) {
			P2A(bead.rheoFAcc,n,k,bead.n) += P2(bead.f,k,bead.rheoDir);
		}
	}
	if ((oftTimeStep % bead.rheoUpdPrd == 0 && bead.rheoSig == 0) 
			|| (oftTimeStep % recBeadLoc.prd == 0 && recBeadLoc.tgl != 0)) {
		MALLOC(beadFall,double,bead.n * nCpu);
		for(n = 0; n < 2; n++) {
			CONT(n == 0 && (oftTimeStep % recBeadLoc.prd != 0 
					|| recBeadLoc.tgl == 0));
			CONT(n == 1 && (oftTimeStep % bead.rheoUpdPrd != 0 
					|| bead.rheoSig != 0));
			MPI_Gather(&P2A(bead.rheoFAcc,n,0,bead.n), bead.n, MPI_DOUBLE, 
					beadFall, bead.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (rank == 0) {
				pD = (n == 1) ? bead.rheoFCurr : bead.rheoFAcc;
				for(k = 0; k < bead.n; k++) {
					pD[k] = 0.;
					for(l = 0; l < nCpu; l++) {
						pD[k] += beadFall[l * bead.n + k];
					}
					pD[k] = -1. * ((n == 1) ? pD[k] / (double)bead.rheoUpdPrd 
							: F_S2PN(pD[k]) / (double)recBeadLoc.prd);
				}
			}
			if (n == 1) {
				MPI_Bcast(bead.rheoFCurr, bead.n, MPI_DOUBLE, 0, 
						MPI_COMM_WORLD);
				if (bead.rheoPreTgl != 0 
						&& bead.rheoFCurr[0] >= bead.rheoPreMag) {
					bead.rheoFGoal = bead.rheoPreMag;
					bead.rheoPreTgl = 0; 
					pres.dur = currTimeStep - netForm.dur - motActiv.dur;
				}
			}
			else {
				if (rank == 0) {
					if (currTimeStep < netForm.dur + pres.dur + motActiv.dur)
					{ fOut = fopen(GenFileName("BeadPreforce"), "a");  }
					else { fOut = fopen(GenFileName("BeadForce"), "a"); }
					fprintf(fOut, "%lld\t", currTimeStep);
					Fprintf1dArrayDouble(fOut, bead.rheoFAcc, bead.n, 0);
					fclose(fOut);
				}
			}
			SetAllValue1dArrayDouble(&P2A(bead.rheoFAcc,n,0,bead.n), 
					bead.n, 0.);
		}
		free(beadFall);
	}
  }
}

void RecordBeadActinBinding(void) {
  int n, k, *cntAll, *sumCnt;
  FILE *fOut; 

  fOut = fopen(GenFileName("BeadBind"), "a");
  MALLOC(cntAll,int,bead.n * 2 * nCpu);
  MPI_Gather(beadBind.cntMe, bead.n * 2, MPI_INT, cntAll, 
		bead.n * 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	MALLOC(sumCnt,int,bead.n * 2);
	memset(sumCnt, 0, sizeof(int) * bead.n * 2);
	for(n = 0; n < nCpu; n++) {
		for(k = 0; k < bead.n * 2; k++) {
			sumCnt[k] += cntAll[n * bead.n * 2 + k];
		}
	}
	fprintf(fOut, "%lld\t", currTimeStep);
	for(n = 0; n < bead.n * 2; n++) {
		fprintf(fOut, "%d\t", sumCnt[n]);
	}
	fprintf(fOut, "\n");
	free(sumCnt);
  }
  free(cntAll);
  fclose(fOut);
}

/*-------------------------- Recording information ---------------------------*/

/*------------------------------ Initialization ------------------------------*/

void InitBeadInformation(void) {
  int n, k, dir, ind[NDIM], ind2[NDIM], sft[4] = {0, 1, 2, 0};
  double dev, devMag[NDIM], sumCos;

  if (dir2D < 0) {
	if (bead.n != 1 && bead.n != 2 && bead.n != 4 && bead.n != 8) {
		Printf0("Error: the number of beads should be either 1, 2, 4, or 8!\n");
		exit(0);
	}
  }
  else {
	if (bead.n != 1 && bead.n != 2 && bead.n != 4) {
		Printf0("Error: the number of beads should be either 1, 2, or 4!\n");
		exit(0);
	}
  }
  bead.stfRep = KS_NPM2S(STF_REP_BEAD) / bead.thkRep * bead.facStfRep;
  CheckPeriodToggleParameter(&recBeadLoc);

  MALLOC(bead.r,double,bead.n * NDIM);
  MALLOC(bead.f,double,bead.n * NDIM);
  MALLOC(bead.rad,double,bead.n);
  MALLOC(bead.drag,double,bead.n);
  MALLOC(bead.maxFbr,double,bead.n);

  if (gTglBead != 0 && beadBind.gTgl != 0) {
	MALLOC(beadBind.cntMe,int,bead.n * 2);	
	memset(beadBind.cntMe, 0, sizeof(int) * bead.n * 2);
  }
  if (rank == 0) {
	VS3COPY(devMag, dimDom, INV(40.));
	V3IND_ASSIGN_CONST_INT(bead.n - 1, 2, ind2);
	for(n = 0; n < bead.n; n++) {
		V3IND_ASSIGN_CONST_INT(n, 2, ind);
		for(k = 0; k < NDIM; k++) {
			dir = (NDIM + k - sft[(int)log2(bead.n)]) % NDIM;
			dev = genrand_real3();
			dev = (2. * dev - 1.) * devMag[k];
			P2(bead.r,n,k) = (0.5 - 0.25 * ind2[dir] + 0.5 * ind[dir]) 
					* dimDom[k] + dev;
		}
	}
	if (!(dir2D < 0)) {
		for(n = 0; n < bead.n; n++) {
			P2(bead.r,n,dir2D) = dimDomH[dir2D];
		}
	}
  }
  MPI_Bcast(bead.r, bead.n * NDIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for(n = 0; n < bead.n; n++) {
	bead.rad[n] = L_UM2S(bead.rad2);
	if (dir2D < 0) { 
		bead.drag[n] = 6. * PI * VISCOSITY * L_S2M(bead.rad[n]) / actF.dragR;
		bead.maxFbr[n] = act.maxFbr[0] * sqrt(bead.drag[n]);
	}
	else {
		bead.drag[n] = 6. * PI * VISCOSITY * L_S2M(bead.rad[n]) / actF.dragR;
		bead.maxFbr[n] = act.maxFbr[0] * sqrt(bead.drag[n]);
	}
  }
  if (bead.rheoWay == 1) {
	MALLOC(bead.rheoFAcc,double,bead.n * 2);
	for(n = 0; n < bead.n * 2; n++) {
		bead.rheoFAcc[n] = 0.;
	}
	bead.rheoPrd = T_SEC2TS(bead.rheoPrdR);

	if (bead.rheoPreMag == 0.) { 
		bead.rheoPreMagDsp = 0.; 
		pres.dur = 0; 
		if (bead.rheoSig == 0) { bead.rheoPreTgl = 0; }
	}
	else { 
		if (bead.rheoSig == 1) {
			bead.rheoPreMagDsp = L_NM2S(bead.rheoPreRate) * dtReal; 
			pres.dur = T_SEC2TS(bead.rheoPreMag / bead.rheoPreRate);
			bead.rheoPreMag = L_NM2S(bead.rheoPreMag);	
		}
		else {
			bead.rheoPreMagDsp = F_PN2S(bead.rheoPreRate) * dtReal; 
			pres.dur = (int)POS_LARGE_VALUE; 
			bead.rheoPreMag = F_PN2S(bead.rheoPreMag);
			bead.rheoPreTgl = 1;
		}
	}

	sumCos = 0.;
	for(n = 0; n < bead.rheoPrd / 4; n++) {
		sumCos += cos(PI * 0.5 / ((double)bead.rheoPrd * 0.25) * (double)n);
	}
	bead.rheoAmp2 = ((bead.rheoSig == 0) 
			? F_PN2S(bead.rheoAmp) : L_NM2S(bead.rheoAmp)) / sumCos;

    if (bead.rheoSig == 0) { 
		bead.rheoUpdPrd = T_SEC2TS(PRD_UPD_SINU_FORCE);
		if (bead.rheoUpdPrd == 0) { bead.rheoUpdPrd = 1; }
		bead.rheoFGoal = 0.;
		MALLOC(bead.rheoFCurr,double,bead.n);
		MALLOC(bead.rheoFAccErr,double,bead.n);
		for(n = 0; n < bead.n; n++) {
			bead.rheoFCurr[n] = 0.;
			bead.rheoFAccErr[n] = 0.;
		}
	}
	else {
		bead.rheoUpdPrd = 1;
	}
  }
}

/*------------------------------ Initialization ------------------------------*/

int CheckActinAbpOverlapBead(double *r1, double *r2, int mode) {
  int n, k, l, CS;
  double len, dist, ratio, rPnt[2][NDIM], dr[NDIM], *pR;

  CS = 1;
  for(n = 0; n < bead.n; n++) {
	V3COPY(rPnt[0], r1);
	V3COPY(rPnt[1], r2);
	for(k = 0; k < 2; k++) {
		pR = (k == 0) ? &P2(bead.r,n,0) : rPnt[0];
		FOR_NDIM(l) {
			CONT(pbc[l] != 1);
			CONT(!(fabs(rPnt[k][l] - pR[l]) > 1.5 * dimDomH[l]));
			rPnt[k][l] += dimDom[l] 
					* ((rPnt[k][l] < pR[l]) ? -1. : 1.);
		}
		if (dir2D > -1) {
			rPnt[k][dir2D] = dimDomH[dir2D];
		}
	}
  	dist = bead.rad[n] + 0.5 * ((mode == 0) ? actF.dia : abpF.repDia[mode - 1]);
	len = CalcSegPntDist(&rPnt[0], &P2(bead.r,n,0), dr, &ratio);
	if (len < dist) {
		CS = 0;
		break;
	}
  }
  return CS;
}

