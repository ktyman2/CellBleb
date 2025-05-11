// This file contains functions related to membranes.

/*----------------------------- Calculate force ------------------------------*/

// Calculate volume-exclusion effects related to membrane.
void CalcMembraneRepulsiveForces(void) {
  int n, k, l, *nlPnt, ind[6], ind2, mode, CS;
  double rPnt[6][NDIM], *r_p[6], *f_p[6], dist[2];

  // Initialization of the list of elements interacting with membrane segments
  if (mbSld.act.gTgl != 0) {
	memset(mbSld.act.l, -1, sizeof(int) * (nActMe + nActCp) * dimMbNuc);
	FOR_ACTMECP(n) {
		P2A(mbSld.act.info,n,0,7) = POS_LARGE_VALUE;
	}
  } 
  if (mbSld.abp.gTgl != 0) {
	memset(mbSld.abp.l, -1, sizeof(int) * (nAbpMe + nAbpCp) * dimMbNuc);
	FOR_ABPMECP(n) {
		P2A(mbSld.abp.info,n,0,7) = POS_LARGE_VALUE;
	}
  } 

  nlPnt = neighMb.l;
  for(n = 0; n < neighMb.c; n++) {
	if (n > 0) { nlPnt += 2; }
	// From the neighboring list, required information is prepared
	mode = 2;
	SetAllValue1dArrayInt(ind, dimMbNuc * 2, -2);
	for (k = 0; k < 2; k++) {
		if (nlPnt[k] < nAct) {
			for(l = 0; l < 2; l++) {
				ind2 = dimMbNuc * k + l;
				ind[ind2] = iAct[P2A(act.cyl.l,nlPnt[k],l,2)];
				if (ind[ind2] > -1) {
					r_p[ind2] = &P2(act.r,ind[ind2],0);
					f_p[ind2] = &P2(act.f,ind[ind2],0);
				}
			}
			dist[k] = actF.dia;
			mode = 1;
		}
		else if (nlPnt[k] >= nAct && nlPnt[k] < nAct + nAbp) {
			ind2 = dimMbNuc * k;
			ind[ind2] = iAbp[nlPnt[k] - nAct];
			if (ind[ind2] > -1) {
				r_p[ind2] = &P2(abp.r,ind[ind2],0);
				f_p[ind2] = &P2(abp.f,ind[ind2],0);
				dist[k] = abpF.repDia[K_ABP(ind[ind2])];
			}
			mode = 0;
		}
		else {
			for(l = 0; l < dimMbNuc; l++) {
				ind2 = dimMbNuc * k + l;
				if (nlPnt[k] < nAct + nAbp + nUnitMb) {
					ind[ind2] = iMb[P2A(memb.unit.l,
							nlPnt[k] - nAct - nAbp,l,dimMbNuc)];
				}
				else {
					ind[ind2] = iMb[P2A(memb.unitCp.l, 
							nlPnt[k] - nAct - nAbp - nUnitMb,l,dimMbNuc)];
				}
				if (ind[ind2] > -1) {
					r_p[ind2] = &P2(memb.r,ind[ind2],0);
					f_p[ind2] = &P2(memb.f,ind[ind2],0);
				}
			}
			if (ind[dimMbNuc * k] > -1) {
				dist[k] = (ISNUC(memb.idx[ind[dimMbNuc * k]])) 
						? memb.nucThk : memb.thk;
			}
		}
	}
	CS = FindElementArray(ind, dimMbNuc * 2, -1, 0, 1);
	CONT(CS > -1);
	// Adjust segment positions
	for(k = 0; k < dimMbNuc * 2; k++) {
		CONT(ind[k] < 0);
		V3COPY(rPnt[k], r_p[k]);
		if (dimMbNuc == 2) {
			rPnt[k][dirNormMbNuc] = dimDomH[dirNormMbNuc]; 
		}
	}
	CalcRepulsiveForcesSubroutine(ind, nlPnt, rPnt, f_p, dist, mode + 3);
  }
}

void CalcMembraneSlideForces(void) {
  int m, n, k, mbInd[NDIM], locMbInd[NDIM], nMe, nCp, *pArr, cnt;
  double f, len, fac, fi[NDIM], dr[NDIM], ratio[3]; //dr2[NDIM], ratio[3];
  double *f_p[3], rPnt[3][NDIM];
  MembSldList *pM;

  for(m = 0; m < 2; m++) {
	if (m == 0) {
		CONT(mbSld.act.gTgl == 0);
		pM = &mbSld.act;
		nMe = nActMe;
		nCp = nActCp;
	}
	else {
		CONT(mbSld.abp.gTgl == 0);
		pM = &mbSld.abp;
		nMe = nAbpMe;
		nCp = nAbpCp;
	}
	for(n = 0; n < nMe + nCp; n++) {
		CONT(P2A(pM->l,n,0,dimMbNuc) == -1);
		if (m == 0) { CONT(ISACTM(n)); }
		else { CONT(ISABPIM(n)); }			
		len = P2A(pM->info,n,0,7);
		V3COPY(dr, &P2A(pM->info,n,1,7));
		V3COPY(ratio, &P2A(pM->info,n,4,7));
		cnt = 0;
		for(k = 0; k < dimMbNuc; k++) {
			mbInd[k] = P2A(pM->l,n,k,dimMbNuc);
			locMbInd[k] = iMb[mbInd[k]];
			CONT(locMbInd[k] >= nMbMe);
			cnt++;
		}

		CONT(n >= nMe && cnt == 0);

		for(k = 0; k < dimMbNuc; k++) {
			V3COPY(rPnt[k], &P2(memb.r,locMbInd[k],0));
			f_p[k] = &P2(memb.f,locMbInd[k],0);
		}
		f = -1. * mbSld.stf * (len - mbSld.eqDist) / len;
		VS3COPY(fi, dr, f);
		if (m == 0) {
			fac = (mbSld.act.side[n] == 0) ? 1. - ratio[dimMbNuc - 1] 
					: ratio[dimMbNuc - 1];
			V3SCALE(fi, fac);
			VV3ADD(&P2(act.f,n,0), fi);
		}
		else {
			VV3ADD(&P2(abp.f,n,0), fi);
		}
		// Apply the calculate force to a membrane segment
		if (cnt > 0) {
			if (dimMbNuc == 2) {
				VVS3ADD(f_p[0], fi, REVSIGN(1. - ratio[0]));
				VVS3ADD(f_p[1], fi, REVSIGN(ratio[0]));				
			}
			else {
				V3REVSIGN(fi);
				DistributeForceOnTria(&rPnt[0], &f_p[0], ratio, fi);
			}
		}
	}
  }
}

void CalcMembraneSpringForces(void) {
  int m, n, k, l, ind, side, kind, locInd, nBd, *chkBondForce;
  int nextInd, locNextInd, chkErr;
  double stf, len, eqLen, dr[NDIM], f, *membSprStf, conc;
  MembMyoCyto *pM; 
  Force *pS;

  nBd = nChMb - nMbAct;
  chkBondForce = allIntL;
  memset(chkBondForce, -1, sizeof(int) * nMbMe * nBd);

  if (gTglMbCont != 0) {
	ind = nObjMbNuc / (gTglNuc != 0 ? 2 : 1);
	MALLOC(membSprStf,double,ind);
	for(n = 0; n < ind; n++) {
		membSprStf[n] = (memb.area.val[n] <= memb.spr.hi) ? 0. : 
				memb.spr.stf * (memb.area.val[n] - memb.spr.hi);
	}
  }
  // Spring forces between membrane points
  FOR_MBME(m) {
	kind = (ISNUC(memb.idx[m])) ? 1 : 0;
	for(n = 0; n < nBd; n++) {
		ind = P2A(memb.ch,m,n,nChMb);
		CONT(ind < 0);
		locInd = iMb[ind];
		CONT(P2A(chkBondForce,m,n,nBd) != -1);
		len = CalcVecDist(dr, &P2(memb.r,m,0), &P2(memb.r,locInd,0), 0);
		if (gTglMbCont == 0 || kind == 1) {
			// Find the equilibrium length of a chain
			if (dimMbNuc == 2) {
				eqLen = (kind == 0) ? memb.len : memb.nucLen;
			}
			else {
				eqLen = P2A(memb.eqLen,m,n,nBd);
			}
			// Normalize the chain length using the equilibrium one
			len /= eqLen;
			pS = (kind == 0) ? &memb.spr : &memb.nucSpr;
			if (len < pS->lo || len > pS->hi) {
				f = SPRING(pS->stf, len, ((len < pS->lo) ? pS->lo : pS->hi));
				chkErr = CheckLargeForce(f, 14);
				if (chkErr != 1) {
					RecordErrorSpringForce(memb.id[m], ind, f, len, 3);
				}
				AddSpringForce(f, 1., dr, &P2(memb.f,m,0),
						&P2(memb.f,locInd,0));
				recMb.sprF[m] += REVSIGN(f); 
				if (locInd < nMbMe) {
					recMb.sprF[locInd] += REVSIGN(f); 
				}
			}
		}
		else {
			for(l = 0; l < 2; l++) { 
				pM = (l == 0) ? &memb.myo : &memb.cyto;
				if (mbPro.gTgl != 0) {
					conc = (pM->conc[m] < pM->conc[locInd]) 
							? pM->conc[m] : pM->conc[locInd];
				}
				else { conc = 1.; }
				if (l == 0) { 
					// From myosin contractility
					f = SPRING(membSprStf[memb.idx[m]] + memb.myo.stf * conc, 
							len, 0.);
				}
				else {
					// From actin cortex
					f = SPRING(memb.cyto.stf * conc, len, (dimMbNuc == 2 
							? memb.len : P2A(memb.eqLen,m,n,nBd)));
				}
				AddSpringForce(f, len, dr, &P2(memb.f,m,0), 
						&P2(memb.f,locInd,0));
			}
		}
		CONT(!(locInd < nMbMe));
		side = FindElementArray(&P2A(memb.ch,locInd,0,nChMb),
				nBd, memb.id[m], 0, 1);
		P2A(chkBondForce,locInd,side,nBd) = 1;
	}
  }
  // Spring forces between membrane points and the bound actins
  FOR_MBMECP(m) {
	pS = ISNUC(memb.idx[m]) ? &memb.nucSpr2 : &memb.spr2;
	for(n = 0; n < nMbAct; n++) {
		ind = P2A(memb.ch,m,n + nBd,nChMb);
		CONT(ind < 0);
		locInd = iAct[ind];
		// If both of the membrane point and actin belong to other subdomains,
		// the calculation doesn't need to be done here.
		CONT(m >= nMbMe && (locInd < 0 || locInd >= nActMe));
		CalcVec(dr, &P2(memb.r,m,0), &P2(act.r,locInd,0)); 
		if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
		len = V3LEN(dr);
		f = SPRING(pS->stf, len, pS->eq);
		chkErr = CheckLargeForce(f, 15);
		if (chkErr != 1) {
			RecordErrorSpringForce(memb.id[m], ind, f, len, 4);
		}
		AddSpringForce(f, len, dr, &P2(memb.f,m,0), &P2(act.f,locInd,0));
	}
  }
  if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
	FOR_MBME(n) {
		CONT(memb.fix[n] < 0);
		CONT(ISNUC(memb.idx[n]));
		CalcVec(dr, &P2(memb.r,n,0), &P2(memb.rFix,n,0));
		if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
		len = V3LEN(dr);
		CONT(len == 0.);
		f = SPRING(memb.spr3.stf, len, memb.spr3.eq);
		chkErr = CheckLargeForce(f, 16);
		if (chkErr != 1) {
			RecordErrorSpringForce(memb.id[n], -1, f, len, 5);
		}
		f /= len;
		FOR_NDIM(k) {
	  		P2(memb.f,n,k) += f * dr[k];
		}
	}
  }
  if (gTglMbCont != 0) {
	free(membSprStf);
  }
}

// Calculate the bending forces of a 2-D membrane
void Calc2dMembraneBendForces(void) {
  int n, ind1, ind2, chkErr, *pArr;
  double dr1[NDIM], dr2[NDIM], f1[NDIM], f2[NDIM], eqAng, stf, f;

  FOR_MBMECP(n) {
	pArr = &P2A(memb.ch,n,0,nChMb);
	ind1 = iMb[pArr[0]];
	ind2 = iMb[pArr[1]];
	CONT(n >= nMbMe && !((ind1 < nMbMe || ind2 < nMbMe) 
			&& (ind1 > -1 && ind2 > -1)));
	CONT(ind1 < 0 || ind2 < 0);
	// If nucleus
	if (ISNUC(memb.idx[n])) {
		eqAng = memb.nucBend.eq;
		stf = memb.nucBend.stf;
	}
	// If membrane
	else {
		// If the membrane underwent protrusion and is under recovery,
		// bending force shouldn't be applied.
		if (mbPro.gTgl != 0 && gTglMbCont != 0) {
			CONT(memb.cyto.conc[n] < 1. || memb.cyto.conc[ind1] < 1. 
					|| memb.cyto.conc[ind2] < 1.);
		}
		eqAng = memb.bend.eq;
		stf = memb.bend.stf * ((mbPro.gTgl != 0 && gTglMbCont != 0) 
				? memb.cyto.conc[n] : 1.);
	}
	CalcVec(dr1, &P2(memb.r,n,0), &P2(memb.r,ind1,0));	
	CalcVec(dr2, &P2(memb.r,ind2,0), &P2(memb.r,n,0));
	f = CalcBendForce(dr1, dr2, &P2(memb.f,ind1,0), &P2(memb.f,n,0), 
			&P2(memb.f,ind2,0), stf, eqAng, f1, f2);

	chkErr = CheckLargeForce(f, 17);
	if (chkErr != 1) { 
		RecordErrorBendingForce(memb.id[n], pArr[0], pArr[1], f, 5);
	}
  }
}

void Calc3dMembraneBendForces(void) {
  int n, k, l, loc, CS, *pArr, ind[4], listInd[2];
  double dr[6][NDIM], crs[4][NDIM], mag[3], e[NDIM], fi[4][NDIM], tempDbl[NDIM];
  double f, c, b11, b12, b22;
  ListInt list;

  MALLOC(list.l, int, memb.unit.c * 2 * 2);
  list.c = 0;
  for(n = 0; n < memb.unit.c; n++) {
	pArr = &P2A(memb.unit.l,n,0,3);
	for(k = 0; k < 3; k++) {
		CalcVec(dr[k], &P2(memb.r,iMb[pArr[(k + 1) % NDIM]],0), 
				&P2(memb.r,iMb[pArr[k]],0));
	}
	for(k = 0; k < 3; k++) {
		ind[0] = pArr[k];
		ind[1] = pArr[(k + 1) % NDIM];
		ind[2] = pArr[(k + 2) % NDIM];
		loc = FindElementArray(&P2A(memb.ch,iMb[ind[1]],0,nChMb), 
				nChMb - nMbAct, ind[0], 0, 1);
		loc = (loc + 2) % (P2A(memb.ch,iMb[ind[1]],5,nChMb) == -1 ? 5 : 6);
		ind[3] = P2A(memb.ch,iMb[ind[1]],loc,nChMb);

		if (ind[0] < ind[3]) { V2SET(listInd, ind[0], ind[3]); }
		else { V2SET(listInd, ind[3], ind[0]); }
		CS = Find2ElementArray(list.l, list.c, listInd[0], listInd[1], 0, 2);
		CONT(CS != -1);
		VS3COPY(dr[3], dr[(k + 2) % NDIM], -1.);
		V3CROSS(crs[0], dr[k], dr[3]);
		mag[0] = V3LEN_SQ(crs[0]);		

		CalcVec(dr[4], &P2(memb.r,iMb[ind[2]],0), &P2(memb.r,iMb[ind[3]],0));
		CalcVec(dr[5], &P2(memb.r,iMb[ind[1]],0), &P2(memb.r,iMb[ind[3]],0));
		V3CROSS(crs[1], dr[4], dr[5]);
		mag[1] = V3LEN_SQ(crs[1]);		
		mag[2] = sqrt(mag[0] * mag[1]);

		CONT(mag[0] < POS_SMALL_VALUE || mag[1] < POS_SMALL_VALUE
				|| mag[2] < POS_SMALL_VALUE);

		c = V3DOT(crs[0], crs[1]) / mag[2];
		b11 = -1. * c / mag[0];
		b22 = -1. * c / mag[1];
		b12 = 1. / mag[2];
		V3COPY(e, dr[(k + 1) % NDIM]);
		//
		V3CROSS(crs[2], crs[0], e);
		V3CROSS(crs[3], crs[1], e);
		VSS3ADD(fi[0], crs[2], crs[3], b11, b12);
		VSS3ADD(fi[3], crs[2], crs[3], -1. * b12, -1. * b22);
		//
		V3CROSS(crs[2], dr[3], crs[0]);
		V3CROSS(crs[3], crs[0], dr[4]);
		VSS3ADD(fi[1], crs[2], crs[3], b11, b12);
		V3CROSS(crs[2], dr[3], crs[1]);
		V3CROSS(crs[3], crs[1], dr[4]);
		VSS3ADD(tempDbl, crs[2], crs[3], b12, b22);
		VV3ADD(fi[1], tempDbl);
		//
		V3CROSS(crs[2], crs[0], dr[k]);
		V3CROSS(crs[3], dr[5], crs[0]);
		VSS3ADD(fi[2], crs[2], crs[3], b11, b12);
		V3CROSS(crs[2], crs[1], dr[k]);
		V3CROSS(crs[3], dr[5], crs[1]);
		VSS3ADD(tempDbl, crs[2], crs[3], b12, b22);
		VV3ADD(fi[2], tempDbl);
		//
		for(l = 0; l < 4; l++) {
			VVS3ADD(&P2(memb.f,iMb[ind[l]],0), fi[l], memb.bend.stf);
			if (iMb[ind[l]] < nMbMe) {
				recMb.bendF[iMb[ind[l]]] += V3LEN(fi[l]) * memb.bend.stf;
			}
		}
		V2COPY(&P2A(list.l,list.c,0,2), listInd);
		(list.c)++;
	}
  }
  free(list.l);
}

void CalcMembraneAreaForces(void) {
  int n, k, locInd, *pArr, kind;
  double stf, fac, dr[4][NDIM], xi[NDIM], alpha, lenXi;

  for(n = 0; n < memb.unit.c; n++) {
	pArr = &P2A(memb.unit.l,n,0,3);
	locInd = iMb[pArr[0]];
	kind = ISNUC(memb.idx[locInd]) ? 1 : 0;
	CONT((mbAre.gTgl == 0 && kind == 0) || (nucAre.gTgl == 0 && kind == 1));

	stf = (kind == 1) ? memb.nucArea.stf : memb.area.stf;
	for(k = 0; k < 3; k++) {
		CalcVec(dr[k], &P2(memb.r,iMb[pArr[(k + 1) % NDIM]],0), 
				&P2(memb.r,iMb[pArr[k]],0));
	}
	V3CROSS(xi, dr[2], dr[0]);
	lenXi = V3LEN(xi);
	if ((kind == 0 && mbAre.gTgl == 1) || (kind == 1 && nucAre.gTgl == 1)) { 
		alpha = (0.5 * lenXi - memb.unitEqAr[n]) / memb.unitEqAr[n];
	}
	else {
		alpha = memb.area.val[memb.idx[locInd]];
	}
	fac = -1. * stf * alpha / (2. * lenXi);

	for(k = 0; k < 3; k++) {
		V3CROSS(dr[3], xi, dr[(k + 1) % NDIM]);
		VVS3ADD(&P2(memb.f,iMb[pArr[k]],0), dr[3], fac);
	}
  }
}

void CalcMembraneVolumeForces(void) {
  int n, k, locInd, *pArr, kind;
  double stf, fac, dr[4][NDIM], xi[NDIM], t[NDIM];
 
  for(n = 0; n < memb.unit.c; n++) {
	pArr = &P2A(memb.unit.l,n,0,3);
	locInd = iMb[pArr[0]];
	kind = ISNUC(memb.idx[locInd]) ? 1 : 0;
	CONT((mbVol.gTgl == 0 && kind == 0) || (nucVol.gTgl == 0 && kind == 1));

	stf = (kind == 1) ? memb.nucVol.stf2 : memb.vol.stf2;
	for(k = 0; k < 3; k++) {
		CalcVec(dr[k], &P2(memb.r,iMb[pArr[(k + 1) % NDIM]],0), 
				&P2(memb.r,iMb[pArr[k]],0));
	}
	FOR_NDIM(k) {
		t[k] = (P2(memb.r,iMb[pArr[0]],k) + P2(memb.r,iMb[pArr[1]],k)
				+ P2(memb.r,iMb[pArr[2]],k)) / 3.;
	}
	VV3SUB(t, &P2(cenMb,memb.idx[locInd],0));

	V3CROSS(xi, dr[0], dr[2]);
	fac = -1. * stf * memb.vol.val[memb.idx[locInd]] / 6.;

	for(k = 0; k < 3; k++) {
		V3CROSS(dr[3], dr[(k + 1) % NDIM], t);
		VVS3ADD(&P2(memb.f,iMb[pArr[k]],0), xi, fac / 3.);
		VVS3ADD(&P2(memb.f,iMb[pArr[k]],0), dr[3], fac);
	}
  }
}

void CalcMembraneProtrusiveForces(void) {
  int n;
  double dr[NDIM];

  FOR_MBME(n) {
	// If a nucleus segment, skip it.
	CONT(ISNUC(memb.idx[n]));
	CONT(P2A(memb.pro,n,0,2) != 1);
	CalcUnitVec(dr, &P2(memb.r,n,0), &P2(cenMb,memb.idx[n],0));
	if (dimMbNuc == 2) {
		dr[dirNormMbNuc] = 0.;
	}
	VVS3ADD(&P2(memb.f,n,0), dr, mbPro.f);
  }
}

void CalcMembraneSubstrateProtrusiveForces(void) {
  int n, dir;
  double dr[NDIM], thres;

  dir = 2;
  thres = 1.;

  FOR_MBME(n) {
	// If a nucleus segment, skip it.
	CONT(ISNUC(memb.idx[n]));
	CONT(memb.fix[n] > -1);
	CONT(P2(memb.r,n,dir) >= rGrid[dir][0] + thres);
	CalcUnitVec(dr, &P2(memb.r,n,0), &P2(cenMb,memb.idx[n],0));
	dr[dir] = 0.;
	NormVec(dr);

	VVS3ADD(&P2(memb.f,n,0), dr, mbPro.f);
  }
}

/*----------------------------- Calculate force ------------------------------*/

/*---------------------------- Update information ----------------------------*/

// Update the fixation of membrane points in a space.
void UpdateMembraneFixation(void) {
  int n, CS;

  int dir = 2;
  double thres = 0.75;

  for(n = 0; n < nMbMe; n++) { 
	CONT(ISNUC(memb.idx[n]));
	// There should be an actin filament bound to the membrane point
	// If it is not fixed already
	CONT(memb.fix[n] > -1);
	//CONT(P2(memb.r,n,dir) >= rGrid[dir][0] +  thres && 
	//P2(memb.r,n,dir) < rGrid[dir][nGrid[dir] - 1] - thres);
	CONT(P2(memb.r,n,dir) >= rGrid[dir][0] +  thres);

	// Check a probability
	CONT(!(genrand_real3() < mbFix.p));
	CS = FindElementArray(noMbDyn.l, noMbDyn.c, memb.id[n], 0, 3);
	CONT(CS > -1);
	memb.fix[n] = 1;
	V3COPY(&P2(memb.rFix,n,0), &P2(memb.r,n,0));
	P2(memb.rFix,n,dir) = rGrid[dir][0];
	//P2(memb.rFix,n,dir) = (P2(memb.r,n,dir) < rGrid[dir][0] + thres) 
	//? rGrid[dir][0] : rGrid[dir][nGrid[dir] - 1];
	(mbFix.cntMe)++;
	V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), memb.id[n], -1, currTimeStep);
	(noMbDyn.c)++;
  }
}

// Update the release of membrane points fixed in a space.
void UpdateMembraneRelease(void) {
  int n, fMag, sumNFA, CS;
  double pDef, dr[NDIM], len;

  FOR_MBME(n) {
	CONT(ISNUC(memb.idx[n]));
	CONT(memb.fix[n] < 0);
	CalcVec(dr, &P2(memb.rFix,n,0), &P2(memb.r,n,0));
	if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }    
	len = V3LEN(dr);
	fMag = (int)(memb.spr3.stf * (len - memb.spr3.eq));
	fMag = TrimIntVal(fMag, 0, mbDef.maxF - 1);
	if (mbMat.gTgl != 0) {
		sumNFA = SumArrInt(&P2A(memb.nFA,n,0,nMbAct), nMbAct);
		fMag = (int)(fMag / (sumNFA + 1));
	}
	pDef = mbDef.p[fMag];
	CONT(!(genrand_real3() < pDef));
	CS = FindElementArray(noMbDyn.l, noMbDyn.c, memb.id[n], 0, 3);
	CONT(CS > -1);
	memb.fix[n] = -1;
	(mbDef.cntMe)++;
	V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), memb.id[n], -1, currTimeStep);
	(noMbDyn.c)++;
  }
}
void UpdateMembraneBindingSubroutine(int actInd, int mbInd, int side, 
		int mode) {
  int locActInd, locMbInd;

  locActInd = iAct[actInd];
  locMbInd = iMb[mbInd];
  if (locMbInd > -1) {
	if (mbMat.gTgl != 0) {
		P2A(memb.nFA,locMbInd,side,nMbAct) = 1;
	}
	P2A(memb.ch,locMbInd,side + nChMb - nMbAct,nChMb) = actInd;
  }
  if (locActInd > -1) { act.mbReb[locActInd] = mbInd; }
  if (locMbInd > -1 && locMbInd < nMbMe) { (mbReb.cntMe)++; };
  if (((locActInd > -1 && locActInd < nActMe) 
		|| (locMbInd > -1 && locMbInd < nMbMe)) && mpiMethod == 0) { 
	InsertLongChain(mbInd + nAct + nAbp, actInd, minDimDomC * 0.9);
  }
  V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), mbInd, actInd, currTimeStep);
  (noMbDyn.c)++;
}
// Update the binding of actin filaments on membranes.
void UpdateMembraneBinding(void) {
  int m, n, k, side, locActInd, locMbInd, CS;
  int *nlPnt, *actCyl, *mbCyl, *pArr;
  double dr[NDIM], len, pReb;
  ListInt chkPair;
  Force *pS;

  nlPnt = neighMb.l;
  chkPair.l = allIntL;
  chkPair.c = 0;
  for(n = 0; n < neighMb.c; n++) {
    if (n > 0) { nlPnt += 2; }
	// Consider the possible bond formation between membrane points and actins
	CONT(!((nlPnt[0] < nAct && nlPnt[1] >= nAct + nAbp) 
			|| (nlPnt[1] < nAct && nlPnt[0] >= nAct + nAbp)));
	side = (nlPnt[0] < nAct) ? 0 : 1;
	CONT(nlPnt[1 - side] >= nAct + nAbp + nUnitMb);
	actCyl = &P2A(act.cyl.l,nlPnt[side],0,2);
	mbCyl = &P2A(memb.unit.l,nlPnt[1 - side] - nAct - nAbp,0,dimMbNuc);
	for(m = 0; m < 2; m++) {
		locActInd = iAct[actCyl[m]];
		CONT(locActInd < 0);
		pArr = &P2A(act.ch,locActInd,0,nChAc);
		// Skip actin monomers
		CONT(pArr[0] < 0 && pArr[1] < 0);
		// Depending on the condition file, the middle part of actin filaments
		// cannot be bound.
		if (mbReb.gTglPa == 0) {
			CONT(pArr[0] > -1 && pArr[1] > -1);
		}
		// Allow the binding only on barbed ends
		CONT(!(pArr[0] < 0 && pArr[1] > -1));
		// If the actin is already bound to membrane, skip it.
		CONT(act.mbReb[locActInd] > -1);
		CS = CheckActinAvailability(actCyl[m], 2);
		CONT(CS > -1);
		for(k = 0; k < dimMbNuc; k++) {
			// If the actin is already bound to membrane, skip it.
			BREAK(act.mbReb[locActInd] > -1);
			locMbInd = iMb[mbCyl[k]];
			// If the information of the membrane point doesn't exist, skip it.
			CONT(locMbInd < 0 || locMbInd >= nMbMe);
			// Check whether the membrane is allowed to have bound actins
			CONT(memb.reb[locMbInd] < 0);
			// Check whether the membrane point has an available binding site
			side = FindElementArray(&P2A(memb.ch,locMbInd,nChMb - nMbAct,nChMb),
					nMbAct, -1, 0, 1);
			CONT(side == -1);
			pS = (ISNUC(memb.idx[locMbInd])) ? &memb.nucSpr2 : &memb.spr2;
			pReb = (ISNUC(memb.idx[locMbInd])) ? nucReb.p : mbReb.p;
			CalcVec(dr, &P2(act.r,locActInd,0), &P2(memb.r,locMbInd,0));
			if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
			len = V3LEN(dr);
			// If length between actin and membrane point doesn't lie within
			// a certain range, binding cannot occur.
			CONT(len < pS->lo || len >= pS->hi);
			// Check whether this paper was considered before.
			CS = Find2ElementArray(chkPair.l, chkPair.c, actCyl[m], 
					mbCyl[k], 0, 2);
			CONT(CS > -1);
			V2SET(&P2A(chkPair.l,chkPair.c,0,2), actCyl[m], mbCyl[k]);
			(chkPair.c)++;
			// Check a probability
        	CONT(!(genrand_real3() < pReb));
			CS = FindElementArray(noMbDyn.l, noMbDyn.c, mbCyl[k], 0, 3);
			CONT(CS > -1);
			UpdateMembraneBindingSubroutine(actCyl[m], mbCyl[k], side, 0);
			V3SET(&P2A(sendMbDyn.l,sendMbDyn.c,0,3), mbCyl[k], actCyl[m], 
					side + 1);
			(sendMbDyn.c)++;
			break;
		}
	}
  }
}

int UpdateMembraneUnbindMatureSubroutine(int actInd, int mbInd, int side,
		int mode) {
  int locActInd, locMbInd;
  
  locActInd = iAct[actInd];
  locMbInd = iMb[mbInd];
  if (locMbInd > -1) {
	if (side < 0) {
		side = FindElementArray(&P2A(memb.ch,locMbInd,nChMb - nMbAct,nChMb),
				nMbAct, actInd, 0, 1);
	}
	if (mbMat.gTgl != 0) {
		P2A(memb.nFA,locMbInd,side,nMbAct) = 0;
	}
	P2A(memb.ch,locMbInd,side + nChMb - nMbAct,nChMb) = -1;
  }
  if (locActInd > -1) { act.mbReb[locActInd] = -1; }
  // The pair should be deleted in longCh.l
  if (((locActInd > -1 && locActInd < nActMe) 
			|| (locMbInd > -1 && locMbInd < nMbMe)) && mpiMethod == 0) { 
	DeleteLongChain(mbInd + nAct + nAbp, actInd);
  }
  if (locMbInd > -1 && locMbInd < nMbMe) { (mbUnb.cntMe)++; }
  V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), mbInd, actInd, currTimeStep);
  (noMbDyn.c)++;
  return side;
}

// Update the unbinding of actin filaments from membranes.
void UpdateMembraneUnbindMature(void) {
  int n, k, fMag, CS, actInd, locActInd;
  double dr[NDIM], pUnb, len, f;
  FOR_MBME(n) {
	CS = FindElementArray(noMbDyn.l, noMbDyn.c, memb.id[n], 0, 3);
	CONT(CS > -1);
	for(k = 0; k < nMbAct; k++) {
		actInd = P2A(memb.ch,n,k + nChMb - nMbAct,nChMb);
		CONT(actInd < 0);
		CS = CheckActinAvailability(actInd, 2);
		CONT(CS > -1);
		locActInd = iAct[actInd];
		// Calculate the unit vector along the chain between membrane point 
	 	// and actin.
	  	CalcVec(dr, &P2(memb.r,n,0), &P2(act.r,locActInd,0));
	 	if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
		len = V3LEN(dr);
		f = memb.spr2.stf * (len - memb.spr2.eq);
		// Maturation
		if (mbMat.gTgl != 0) {
			CS = 1;
			// If actins exist inside the membrane, the maturation can occur
			// only at membrane points fixed in space.
			if (sideMb == 0) {
				if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
					if (memb.fix[n] < 0) {
						CS = 0;
					}
				}
			}
			if (mbMat.gTglPa == 0 && CS == 1) {
				if (P2A(act.ch,locActInd,0,nChAc) > -1 
						&& P2A(act.ch,locActInd,1,nChAc) > -1) {
					CS = 0;
				}
			}
			if (CS == 1) {
				if (P2A(memb.nFA,n,k,nMbAct) < mbMat.maxNFA) {
					fMag = TrimIntVal((int)f, 0, mbMat.maxF - 1);
					if (genrand_real3() < mbMat.p[fMag]) {
						P2A(memb.nFA,n,k,nMbAct)++;
						mbMat.cntMe2[0]++;
						continue;
					}
				}
			}
		}
		// Unbinding
		if (mbUnb.gTgl != 0) {
			fMag = TrimIntVal((int)f, 0, mbUnb.maxF - 1);
			if (mbMat.gTgl != 0) {
				fMag = (int)(fMag / (P2A(memb.nFA,n,k,nMbAct) + 1));
			}
		    pUnb = mbUnb.p[fMag];
			// Check a probability.
			CONT(!(genrand_real3() < pUnb));
			if (mbMat.gTgl != 0) {
		  		P2A(memb.nFA,n,k,nMbAct)--;
				mbMat.cntMe2[1]++;
				CONT(P2A(memb.nFA,n,k,nMbAct) > 0);
			}
			UpdateMembraneUnbindMatureSubroutine(actInd, memb.id[n], k, 0);
			V3SET(&P2A(sendMbDyn.l,sendMbDyn.c,0,3), memb.id[n], actInd, 
					-1 * (k + 1));
			(sendMbDyn.c)++;
		}
	}
  } 
}

void UpdateMembraneProtrusion(void) {
  int n, k, ind, CS;
  double pProt;

  int dir;
  double *distCen, radCT, distMb, dr[NDIM], dot;

  if (mbPro.gTglCT != 0) {
	MALLOC(distCen,double,nObjMbNuc);
	for(n = 0 ; n < nObjMbNuc; n++) {
		distCen[n] = 0.;
		for(k = 0; k < dimMbNuc; k++) {
			dir = (dimMbNuc == 2) ? dirMbNuc[k] : k;
			distCen[n] += SQR(dimDom[dir] * mbPro.rCT[dir] - P2(cenMb,n,dir));
		}
		distCen[n] = L_S2UM(sqrt(distCen[n]));
	}
	radCT = 1.2;
  }
  FOR_MBME(n) {
	CONT(ISNUC(memb.idx[n]));
	CS = FindElementArray(noMbDyn.l, noMbDyn.c, memb.id[n], 0, 3);
	CONT(CS > -1);
	if (P2A(memb.pro,n,0,2) == -1) {
		if (gTglMbCont != 0) {
			CONT(!(memb.myo.conc[n] > 0.8 && memb.cyto.conc[n] > 0.8));
		}
		CS = 1;
		for(k = 0; k < nChMb - nMbAct; k++) {
			ind = P2A(memb.ch,n,k,nChMb);
			CONT(ind < 0);
			if (iMb[ind] < 0) {
				CS = 0;
				break;
			}
			if (P2A(memb.pro,iMb[ind],0,2) != -1) {
				CS = 0;
				break;
			}
		}
		CONT(!(CS == 1));
		pProt = mbPro.p;
		// Chemotaxis
		if (mbPro.gTglCT != 0) {
			CONT(distCen[memb.idx[n]] < radCT);
			for(k = 0; k < dimMbNuc; k++) {
				dir = (dimMbNuc == 2) ? dirMbNuc[k] : k;
				dr[dir] = P2(memb.r,n,dir) - dimDom[dir] * mbPro.rCT[dir];
			}
			if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
			NormVec(dr);
			dot = V3DOT(&P2(memb.nDir,n,0), dr);
			if (dot > 0) {
				pProt *= SQR(dot);
			}
			else {
				pProt = 0.;
			}
		}
		// Reflect maturation of adhesion points
		CONT(!(genrand_real3() < pProt));
		V2SET(&P2A(memb.pro,n,0,2), 1, currTimeStep);
		if (gTglMbCont != 0) {
			memb.myo.conc[n] = 0.1;
			memb.cyto.conc[n] = 0.1;
		}
		mbPro.cntMe++;
	}
	else {
		if (currTimeStep - P2A(memb.pro,n,1,2) >= mbPro.dur) {
			V2SET_ALL(&P2A(memb.pro,n,0,2), -1);
		}
		if (mbPro.gTglCT != 0) {
			if (distCen[memb.idx[n]] < radCT) {
				V2SET_ALL(&P2A(memb.pro,n,0,2), -1);
			}
		}
	}
  }
  if (mbPro.gTglCT != 0) { 
	free(distCen); 
  }
}

void UpdateMembraneVolume(void) {
  int n, k, *pArr, idx;
  double *sum, *sumAll, vol, dr[2][NDIM], fac, eq; 
 
  MALLOC(sum,double,nObjMbNuc);
  MALLOC(sumAll,double,nCpu*nObjMbNuc);
  SetAllValue1dArrayDouble(sum, nObjMbNuc, 0.);
  for(n = 0; n < memb.unit.c; n++) {
    pArr = &P2A(memb.unit.l,n,0,dimMbNuc);
    CONT(!(iMb[pArr[0]] < nMbMe && iMb[pArr[0]] > -1));
	idx = memb.idx[iMb[pArr[0]]];
    if (dimMbNuc == 2) {
        vol = CalcTriaArea(&P2(cenMb,idx,0), &P2(memb.r,iMb[pArr[0]],0),
                &P2(memb.r,iMb[pArr[1]],0));
    }
    else {
        vol = CalcTetraVol(&P2(cenMb,idx,0), &P2(memb.r,iMb[pArr[0]],0),
                &P2(memb.r,iMb[pArr[1]],0), &P2(memb.r,iMb[pArr[2]],0));
    }
    CalcVec(dr[0], &P2(cenMb,idx,0), &P2(memb.r,iMb[pArr[0]],0));
    fac = (V3DOT(&P2(memb.unitNDir,n,0), dr[0]) < 0.) ? -1. : 1.;
    vol *= fac;
    sum[idx] += vol;
  }
  MPI_Gather(sum, nObjMbNuc, MPI_DOUBLE, sumAll, nObjMbNuc, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 1; n < nCpu; n++) {
		for(k = 0; k < nObjMbNuc; k++) {
		    sum[k] += P2A(sumAll,n,k,nObjMbNuc);
		}
	}
	if (gTglNuc != 0) {
		for(n = 0; n < nObjMbNuc / 2; n++) {
			sum[n] -= sum[n + nObjMbNuc / 2];
		}
	}
  }
  MPI_Bcast(sum, nObjMbNuc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (currTimeStep > 1) { 
  	for(n = 0; n < nObjMbNuc; n++) {
		eq = ISNUC(n) ? memb.nucVol.eq : memb.vol.eq;
		memb.vol.val[n] = (sum[n] - eq) / eq;
	}
  }
  else { 
	memb.vol.eq = sum[0]; 
	if (gTglNuc != 0) {
		memb.nucVol.eq = sum[nObjMbNuc / 2]; 
	}
	SetAllValue1dArrayDouble(memb.vol.val, nObjMbNuc, 0.);
  }
  free(sum);
  free(sumAll);
}

void UpdateMembraneUnitArea(void) {
  int n, k, nCh, *pArr;
  double area;

  FOR_MBME(n) {
	pArr = &P2A(memb.ch,n,0,nChMb);
	if (dimMbNuc == 2) { nCh = 2; }
	else { nCh = (pArr[5] < 0) ? 5 : 6; }
	memb.areaL[n] = 0.;
	for(k = 0; k < nCh; k++) {
		if (dimMbNuc == 2) { 
			area = CalcDist(&P2(memb.r,n,0), &P2(memb.r,iMb[pArr[k]],0), 0);
			area *= dimDom[dirNormMbNuc];
			memb.areaL[n] += area;
		}
		else {
			area = CalcTriaArea(&P2(memb.r,n,0), &P2(memb.r,iMb[pArr[k]],0), 
					&P2(memb.r,iMb[pArr[(k + 1) % nCh]],0));
			memb.areaL[n] += area;
		}
	}
	memb.areaL[n] /= (dimMbNuc == 2) ? 2. : 3.;
  }
}

void UpdateMembraneTotalArea(void) {
  int n, k, *pArr;
  double area, *sum, *sumAll, eq;
 
  MALLOC(sum,double,nObjMbNuc);
  MALLOC(sumAll,double,nCpu*nObjMbNuc);
  SetAllValue1dArrayDouble(sum, nObjMbNuc, 0.);

  for(n = 0; n < nMbMe; n++) {
	sum[memb.idx[n]] += memb.areaL[n];
  }
  MPI_Gather(sum, nObjMbNuc, MPI_DOUBLE, sumAll, nObjMbNuc, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 1; n < nCpu; n++) {
		for(k = 0; k < nObjMbNuc; k++) {
		    sum[k] += P2A(sumAll,n,k,nObjMbNuc);
		}
	}
  }
  MPI_Bcast(sum, nObjMbNuc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (currTimeStep > 1) { 
  	for(n = 0; n < nObjMbNuc; n++) {
		eq = ISNUC(n) ? memb.nucArea.eq : memb.area.eq;
		memb.area.val[n] = (sum[n] - eq) / eq;
	}
  }
  else { 
	memb.area.eq = sum[0]; 
	if (gTglNuc != 0) {
		memb.nucArea.eq = sum[nObjMbNuc / 2]; 
	}
	SetAllValue1dArrayDouble(memb.area.val, nObjMbNuc, 0.);
  }
  free(sum);
  free(sumAll);
}

// Update the center position of the membrane regularly.
void UpdateMembraneCenterPosition(void) {
  int n, k, idx;
  double *sum, *sumAll, sft;

  MALLOC(sum,double,NDIM*nObjMbNuc);
  MALLOC(sumAll,double,nCpu*NDIM*nObjMbNuc);
  SetAllValue1dArrayDouble(sum, NDIM * nObjMbNuc, 0.);
  FOR_MBME(n) {
	idx = memb.idx[n];
	FOR_NDIM(k) {
		sft = 0.;
		if (pbc[k] == 1) {
			if (P2(cenMb,idx,k) - P2(memb.r,n,k) > dimDomH[k]) { sft = 1.; }
			else if (P2(cenMb,idx,k) - P2(memb.r,n,k) < -1. * dimDomH[k]) 
			{ sft = -1.; }
		}
		P2(sum,idx,k) += P2(memb.r,n,k) + sft * dimDom[k];
	}
  }
  MPI_Gather(sum, NDIM * nObjMbNuc, MPI_DOUBLE, sumAll, NDIM * nObjMbNuc, 
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 1; n < nCpu; n++) {
		for(k = 0; k < NDIM * nObjMbNuc; k++) {
			sum[k] += P2A(sumAll,n,k,NDIM * nObjMbNuc);
		}
	}
	for(n = 0; n < nObjMbNuc; n++) {
		VS3COPY(&P2(cenMb,n,0), &P2(sum,n,0), 
				INV(ISNUC(n) ? nNucPerObj : nMbPerObj));
		ApplyBoundCondVector(&P2(cenMb,n,0), -1, 0);
	}
  }
  MPI_Bcast(cenMb, NDIM * nObjMbNuc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  free(sum);
  free(sumAll);
}

void UpdateMembraneProteinRecover(void) {
  int n;

  FOR_MBME(n) {
	CONT(ISNUC(memb.idx[n]));
	if (memb.myo.conc[n] < 1.) {
		memb.myo.conc[n] = memb.myo.conc[n]
				 + memb.myo.diff * (1. - memb.myo.conc[n]) * dt;
		if (memb.myo.conc[n] > 0.999) { memb.myo.conc[n] = 1.; }
	}
	if (memb.cyto.conc[n] < 1.) {
		memb.cyto.conc[n] = memb.cyto.conc[n]
				 + memb.cyto.diff * (1. - memb.cyto.conc[n]) * dt;
		if (memb.cyto.conc[n] > 0.999) { memb.cyto.conc[n] = 1.; }
	}
  }
}

void AssignMembraneSlideList(void) {
  int n, ind, *rebAct, *rebAbp;
  
  MALLOC(rebAct, int, nAct);
  MALLOC(rebAbp, int, nAbp);
  if (rank == 0) {
	SetAllValue1dArrayInt(rebAct, nAct, -1);
	SetAllValue1dArrayInt(rebAbp, nAbp, -1);

	for(n = 0; n < mbSld.act.nReb; n++) {
		while(1) {
			ind = GenRandIntIndex(nAct);
			CONT(rebAct[ind] != -1);
			rebAct[ind] = 1;
			break;
		}
	}

	for(n = 0; n < mbSld.abp.nReb; n++) {
		while(1) {
			ind = GenRandIntIndex(nAbp);
			CONT(rebAbp[ind] != -1);
			rebAbp[ind] = 1;
			break;
		}
	}
  }
  MPI_Bcast(rebAct, nAct, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rebAbp, nAbp, MPI_INT, 0, MPI_COMM_WORLD);

  FOR_ACT(n) {
  	CONT(rebAct[n] == -1);
	CONT(iAct[n] < 0);
	mbSld.act.reb[iAct[n]] = 1;
  }
  FOR_ABP(n) {
  	CONT(rebAbp[n] == -1);
	CONT(iAbp[n] < 0);
	mbSld.abp.reb[iAbp[n]] = 1;
  }

  if (rank == 0) {
	free(rebAct);
	free(rebAbp);  
  }
}

void UpdateMembraneSlideList(void) {
  int mm, m, n, ind, cnt, min, *pArr, nMe, nCp;
  ListInt reb;
  MembSldList *pM;
  // Calculate the value of cntMe
  if (currTimeStep % recProg.prd == 0 || currTimeStep == 1) {
	for(m = 0; m < 2; m++) {
		CONT(!((m == 0 && mbSld.act.gTgl != 0) 
				|| (m == 1 && mbSld.abp.gTgl != 0)));
		nMe = (m == 0) ? nActMe : nAbpMe;
		pM = (m == 0) ? &mbSld.act : &mbSld.abp;
		pM->cntMe = 0;
		for(n = 0; n < nMe; n++) {
			CONT(P2A(pM->l,n,0,dimMbNuc) == -1);
			if (m == 0) { CONT(ISACTM(n)); }
			else { CONT(ISABPIM(n)); }
			(pM->cntMe)++;
		}	
	}
  }
}

// Divide the value of mbSld.nReb into mbSld.nRebMe in each subdomain
void UpdateMembraneSlidePoint(void) {
  int m, n, sum, rem, nRebMe, *nMbAll, *nRebAll;
  MembSldList *pM;

  MALLOC(nMbAll,int,nCpu);
  MALLOC(nRebAll,int,nCpu * 2);
  MPI_Gather(&nMbMe, 1, MPI_INT, nMbAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	for(m = 0; m < 2; m++) {
		CONT(!((m == 0 && mbSld.act.gTgl != 0) 
				|| (m == 1 && mbSld.abp.gTgl != 0)));
		pM = (m == 0) ? &mbSld.act : &mbSld.abp;
		rem = pM->nReb;
		for(n = 0; n < nCpu; n++) {
			nRebMe = (int)round((double)nMbAll[n] / (double)nMb * pM->nReb);
			if (rem < nRebMe) { nRebMe = rem; }
			P2A(nRebAll,n,m,2) = nRebMe;
			rem -= P2A(nRebAll,n,m,2);
		}
		sum = 0;
		for(n = 0; n < nCpu; n++) {
			sum += P2A(nRebAll,n,m,2);
		}
		if (pM->nReb > sum) {
			for(n = nCpu - 1; n >= 0; n--) {
				if (P2A(nRebAll,n,m,2) > 0) {
					P2A(nRebAll,n,m,2) += pM->nReb - sum;
					break;
				}
			}
		}
	}
  }
  MPI_Bcast(nRebAll, nCpu * 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (mbSld.act.gTgl != 0) {
	  mbSld.act.nRebMe = P2A(nRebAll,rank,0,2);
  }
  if (mbSld.abp.gTgl != 0) {
	  mbSld.abp.nRebMe = P2A(nRebAll,rank,1,2);
  }
  free(nMbAll);
  free(nRebAll);
}

// Eliminate expired elements in noMbDyn.l
void UpdateNoMembraneDynamicsList(void) {
  int n;

  CheckArraySize(&noMbDyn, &noMbDyn.siz, 3, 1);
  for (n = noMbDyn.c - 1; n >= 0; n--) {
    CONT(!(currTimeStep - P2A(noMbDyn.l,n,2,3) > durNoMbDyn));
	DeleteElementArrayByIndex(noMbDyn.l,&noMbDyn.c,n,3);
  }
}

/*---------------------------- Update information ----------------------------*/

/*--------------------------- Recording information --------------------------*/

void RecordMembraneDimension(int period) {
  int n, idx, *nMbAll;
  double len, eq, *avgLen, *stdLen, *maxLen, *minLen;
  FILE *fOut;
 
  MALLOC(nMbAll,int,nCpu);
  if (rank == 0) { 
	MALLOC(rMb,double,nMb * NDIM); 
	MALLOC(avgLen,double,nObjMbNuc);
	MALLOC(stdLen,double,nObjMbNuc);
	MALLOC(maxLen,double,nObjMbNuc);
	MALLOC(minLen,double,nObjMbNuc);
	MALLOC(idxMb,int,nMb);
  }
  MPI_Gather(&nMbMe, 1, MPI_INT, nMbAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  Gather2dArrayDouble(nMbMe, nMbAll, NDIM, memb.id, memb.r, rMb);
  Gather2dArrayInt(nMbMe, nMbAll, 1, memb.id, memb.idx, idxMb);
  // Calculation of global area
  if (!(mbAre.gTgl == 2 || nucAre.gTgl == 2)) {
	UpdateMembraneTotalArea();
  }
  // Calculation of volume
  if (!(mbVol.gTgl != 0 || nucVol.gTgl != 0 
		|| (mbPro.gTgl != 0 && mbPro.gTglCT != 0))) {
	UpdateMembraneVolume();
  }
  if (rank == 0) {
	// Radius
	for(n = 0; n < nObjMbNuc; n++) {
	    avgLen[n] = 0.;
	    stdLen[n] = 0.;
	    maxLen[n] = NEG_LARGE_VALUE;
	    minLen[n] = POS_LARGE_VALUE;
	}
    FOR_MB(n) {
		idx = idxMb[n];
        len = CalcDist(&P2(rMb,n,0), &P2(cenMb,idx,0), 0);
        avgLen[idx] += len;
        stdLen[idx] += SQR(len);
        if (len > maxLen[idx]) { maxLen[idx] = len; }
        if (len < minLen[idx]) { minLen[idx] = len; }
    }
    fOut = fopen(GenFileName("MembRad"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t", nObjMbNuc);
		Fprintf1dFillerInt(fOut, 0, 5, 0);
	}
	for(n = 0; n < nObjMbNuc; n++) {
    	avgLen[n] /= (double)nMbPerObj;
	    stdLen[n] = sqrt(stdLen[n] / (double)nMbPerObj - SQR(avgLen[n]));
    	fprintf(fOut, "%lld\t%d\t%g\t%g\t%g\t%g\n", currTimeStep, n, 
				L_S2UM(avgLen[n]), L_S2UM(stdLen[n]), L_S2UM(maxLen[n]), 
				L_S2UM(minLen[n]));
	}
    fclose(fOut);
	// Volume
    fOut = fopen(GenFileName("MembVol"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t", nObjMbNuc);
		Fprintf1dFillerInt(fOut, 0, 3, 0);
	}
	for(n = 0; n < nObjMbNuc; n++) {
		eq = ISNUC(n) ? memb.nucVol.eq : memb.vol.eq;
    	fprintf(fOut, "%lld\t%d\t%g\t%g\n", currTimeStep, n, 
				(memb.vol.val[n] + 1.) * eq * CUBE(L_SCALE_IN_UM), 
				memb.vol.val[n]);
	}
    fclose(fOut);
	// Area
    fOut = fopen(GenFileName("MembArea"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t", nObjMbNuc);
		Fprintf1dFillerInt(fOut, 0, 3, 0);
	}
	for(n = 0; n < nObjMbNuc; n++) {
		eq = ISNUC(n) ? memb.nucArea.eq : memb.area.eq;
    	fprintf(fOut, "%lld\t%d\t%g\t%g\n", currTimeStep, n, 
				(memb.area.val[n] + 1.) * eq * SQR(L_SCALE_IN_UM), 
				memb.area.val[n]);
	}
    fclose(fOut);
    free(rMb);
	free(avgLen);
	free(stdLen);
	free(maxLen);
	free(minLen);
	free(idxMb);
  }
  free(nMbAll);
}

void RecordMembraneCenter(int period) {
  int n;
  FILE *fOut;
	
  if (rank == 0) {
	fOut = fopen(GenFileName("MembCen"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t", nObjMbNuc);
		Fprintf1dFillerInt(fOut, 0, 3, 0);
	}
	for(n = 0; n < nObjMbNuc; n++) {
    	fprintf(fOut, "%d\t%g\t%g\t%g\n", n, P2(cenMb,n,0), P2(cenMb,n,1),
				P2(cenMb,n,2));
	}
	fclose(fOut);
  }
}

void RecordMembraneAdhesion(void) {
  int n, *cntAll, *pCnt, sumCnt[2];
  double area, *areaAll;
  FILE *fOut;

  if (rank == 0) {
	  MALLOC(cntAll,int,nCpu);
	  MALLOC(areaAll,double,nCpu);
  }
  for(n = 0; n < 2; n++) {
	pCnt = (n == 0) ? &mbFix.cntMe : &mbDef.cntMe;
	MPI_Gather(pCnt, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt[n] = SumArrInt(cntAll, nCpu);
	}
  }
 area = 0.;
  for(n = 0; n < nMbMe; n++) {
	CONT(!(memb.fix[n] > -1));
	area += memb.areaL[n]; 
  }
  MPI_Gather(&area, 1, MPI_DOUBLE, areaAll, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	area = SumArrDbl(areaAll, nCpu);
	fOut = fopen(GenFileName("MembAdhes"), "a");
	fprintf(fOut, "%lld\t%d\t%d\t%d\t%g\n", currTimeStep, sumCnt[0], sumCnt[1], 
			sumCnt[0] - sumCnt[1], area * SQR(L_SCALE_IN_UM));
	fclose(fOut);

	free(cntAll);
	free(areaAll);
  }
}

/*--------------------------- Recording information --------------------------*/

/*------------------------------ Initialization ------------------------------*/

void InitMembraneInformation(void) {
  int m, n, k, cnt, ind, ind2, CS, CS2, nObj, nPerObj, *pArr;
  double r[NDIM], dr[NDIM], v[NDIM], rBnd[NDIM*2], cen[NDIM];
  double ang, por, rad, cNeiEdge;

  if (mbReb.gTgl != 0 && gTglLoadMbNucData == 0) {
	if (rank == 0) {
		memset(rebMb, -1, sizeof(int) * nMb);
		for(n = 0; n < nObjMbNuc; n++) {
			if (ISNUC(n)) {
				ind2 = nObjMbNuc / 2 * nMbPerObj 
						+ (n - nObjMbNuc / 2) * nNucPerObj;
				nPerObj = nNucPerObj;
				por = nucReb.por;
			}
			else {
				ind2 = n * nMbPerObj;
				nPerObj = nMbPerObj;
				por = mbReb.por;
			}
			cnt = 0;
			while(cnt < (int)(por * nPerObj)) {
				ind = GenRandIntIndex(nPerObj) + ind2;
				if (rebMb[ind] == -1) {
					rebMb[ind] = 1;
					cnt++;
				}
			}
		}
	}
	MPI_Bcast(rebMb, nMb, MPI_INT, 0, MPI_COMM_WORLD);
  }

  // Initialize the position and chain information of membrane points.
  // For 2D membrane
  if (dimMbNuc == 2) {
	if (gTglLoadMbNucData == 0) { 
		nObj = nObjMbNuc / ((gTglNuc != 0) ? 2 : 1);
		for(m = 0; m < nObjMbNuc; m++) {
			ind = m % nObj;
			// Nucleus
			if (ISNUC(m)) {
				ind2 = nObjMbNuc / 2 * nMbPerObj 
						+ (m - nObjMbNuc / 2) * nNucPerObj;
				nPerObj = nNucPerObj;
				rad = radNuc;
			}
			// Membrane
			else {
				ind2 = m * nMbPerObj;
				nPerObj = nMbPerObj;
				rad = radMb;
			}
			ang = 2. * PI / (double)nPerObj;

			if (nObj == 1) {
				V3COPY(&P2(cenMb,m,0), dimDomH);
			}
			else if (nObj == 2) {
				P2(cenMb,m,dirMbNuc[0]) = dimDom[dirMbNuc[0]] 
						* ((ind == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirMbNuc[1]) = dimDom[dirMbNuc[1]] 
						* ((ind == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirNormMbNuc) = dimDomH[dirNormMbNuc];
			}
			else if (nObj == 4) {
				P2(cenMb,m,dirMbNuc[0]) = dimDom[dirMbNuc[0]] 
						* ((ind / 2 == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirMbNuc[1]) = dimDom[dirMbNuc[1]] 
						* ((ind % 2 == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirNormMbNuc) = dimDomH[dirNormMbNuc];
			}
			for(n = 0; n < nPerObj; n++) {
				P2(rMb,ind2 + n,dirMbNuc[0]) = P2(cenMb,m,dirMbNuc[0])
						+ rad * cos(ang * (double)n);
				P2(rMb,ind2 + n,dirMbNuc[1]) = P2(cenMb,m,dirMbNuc[1])
						+ rad * sin(ang * (double)n);
				P2(rMb,ind2 + n,dirNormMbNuc) = P2(cenMb,m,dirNormMbNuc);
				V2SET(&P2A(chMb,ind2 + n,0,nChMb), ind2 + n + 1, ind2 + n - 1);
				memset(&P2A(chMb,ind2 + n,2,nChMb), -1, sizeof(int) * nMbAct);
				CalcUnitVec(&P2(nDirMb,ind2 + n,0), &P2(cenMb,m,0), 
						&P2(rMb,ind2 + n,0));
				idxMb[ind2 + n] = m;
			}
			P2A(chMb,ind2,1,nChMb) = ind2 + nPerObj - 1;
			P2A(chMb,ind2 + nPerObj - 1,0,nChMb) = ind2;
			for(n = 0; n < nPerObj; n++) {
				V2SET(&P2A(unitMb,ind2 + n,0,2), ind2 + n, 
						ind2 + (n + 1) % nPerObj);
				V3AVG(&P2(nDirUnitMb,ind2 + n,0), &P2(nDirMb,ind2 + n,0), 
						&P2(nDirMb,ind2 + (n + 1) % nPerObj,0));
			}
		}
		nUnitMb = nMb;
	}
  }
  // For 3D membrane
  else {
	Init3dMembraneInformation();
  }
  // Distribute the information to subdomains.
  memset(iMb, -1, sizeof(int) * nMb);
  cNeiEdge = (motWalk.gTgl != 0) ? neiEdge + 1. : neiEdge;
  V6COPY(rBnd,bnd.r);
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 1));
	if (iCell[k] == 0 && nCell[k] > 1)
	{ P2A(rBnd,0,k,NDIM) = rGrid[k][nGrid[k] - 1]; }
	else if (iCell[k] == nCell[k] - 1 && nCell[k] > 1)
	{ P2A(rBnd,1,k,NDIM) = rGrid[k][0]; }
  }
  nMbMe = 0;
  nMbCp = 0;
  for(m = 0; m < 2; m++) {
	FOR_MB(n) {
		CS = 0;
		CS2 = 0;
		V3COPY(r, &P2(rMb,n,0));
		ExtractConfigSubroutine(r, rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) ||
				(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			ind = (m == 0) ? nMbMe : nMbMe + nMbCp;
			V3COPY(&P2(memb.r,ind,0), &P2(rMb,n,0));
			V3COPY(&P2(memb.nDir,ind,0), &P2(nDirMb,n,0));
			memb.idx[ind] = idxMb[n];
			memb.id[ind] = n;
			iMb[n] = ind;
			if (mbReb.gTgl != 0 || mbUnb.gTgl != 0 
					|| nucReb.gTgl != 0 || nucUnb.gTgl != 0) {
				memb.reb[ind] = rebMb[n]; 
			}
			if (m == 0) { nMbMe++; }
			else { nMbCp++; }
		}
	}
  } 
  FOR_MB(n) {
    CONT(iMb[n] < 0);
    Copy1dArrayInt(&P2A(memb.ch,iMb[n],0,nChMb), &P2A(chMb,n,0,nChMb), nChMb);
	if (dimMbNuc == 3) {
		Copy1dArrayDouble(&P2A(memb.eqLen,iMb[n],0,nChMb - nMbAct), 
				&P2A(eqLenMb,n,0,nChMb - nMbAct), nChMb - nMbAct);
	}
  }
  memb.unit.c = 0;
  memb.unitCp.c = 0;
  for(n = 0; n < nUnitMb; n++) {
 	pArr = &P2A(unitMb,n,0,dimMbNuc);
	for(k = 0; k < dimMbNuc; k++) {
		BREAK(iMb[pArr[k]] > -1 && iMb[pArr[k]] < nMbMe);
	}
	CONT(k == dimMbNuc);
	Copy1dArrayInt(&P2A(memb.unit.l,memb.unit.c,0,dimMbNuc), pArr, dimMbNuc);
	V3COPY(&P2(memb.unitNDir,memb.unit.c,0), &P2(nDirUnitMb,n,0));
	if (dimMbNuc == 3 && (mbAre.gTgl == 1 || nucAre.gTgl == 1)) { 
		memb.unitEqAr[memb.unit.c] = eqArUnitMb[n];
	}
	(memb.unit.c)++;
  }

  radMbInit = radMb;
  radNucInit = radNuc;

if (tglBias != 0 && typeBias == 2) {
	FOR_MBME(n) {
		memb.bias[n] = (P2(memb.r,n,0) > rGrid[0][0]
            + dimDom[0] * 0.5 + radMb * (1. - thkBias)) ? 1 : 0;
	}
}

  free(rMb);
  free(nDirMb);
  free(chMb);
  free(unitMb);
  free(nDirUnitMb);
  free(idxMb);
  if (dimMbNuc == 3 && (mbAre.gTgl == 1 || nucAre.gTgl == 1)) { 
	free(eqArUnitMb);
  }
  if (mbReb.gTgl != 0) { free(rebMb); }
  if (dimMbNuc == 3) {
	free(eqLenMb);
  }
  UpdateSubdomSectionLocation(1);
  volMe = HowManyInOutMembrane(bnd.r);
  if (actNuc.gTgl != 0) {
	gTglMbActNuc = (volMe > 0.) ? 1 : 0;
  }
  UpdateMembraneSlidePoint();
}

void Init3dMembraneInformation(void){
  // Array used for initializing the vertices of the icosahedron
  int initTriaArr[60] = {0,2,1,0,3,2,0,4,3,0,5,4,0,1,5,1,2,7,2,3,8,3,4,9,4,5,
		10,5,1,6,1,7,6,2,8,7,3,9,8,4,10,9,5,6,10,6,7,11,7,8,11,8,9,11,9,10,11,
		10,6,11};
  int mm, m, n, k, l, CS, nMb2, nUnitMb2, *pArr, *pArr2, nCh, nObj, nPerObj;
  int nUnitMbPerObj, nUnitNucPerObj, nUnitPerObj;
  int ind[2], vert[6], order[12] = {0,3,5,3,1,4,5,4,2,5,3,4};
  int *chMb2, *unitMb2;	
  double theta, sinThe, cosThe, phi, len, lev, rad, min;
  double dr[4][NDIM], r[NDIM], rIcos[12][NDIM], ang[2];
  double *rMb2, *nDirMb2, *eqLenMb2, *nDirUnitMb2, *eqArUnitMb2;
  ListInt newVert;

  nUnitMbPerObj = 20 * (int)pow(4., levMb);
  nUnitNucPerObj = 20 * (int)pow(4., levNuc);
  nObj = nObjMbNuc / ((gTglNuc != 0) ? 2 : 1);
  for(mm = 0; mm < ((gTglNuc != 0) ? 2 : 1); mm++) {
	if (mm == 0) {
		lev = levMb;
		nPerObj = nMbPerObj;
		rad = radMb;
		nUnitPerObj = nUnitMbPerObj;
	}
	else {
		lev = levNuc;
		nPerObj = nNucPerObj;
		rad = radNuc;
		nUnitPerObj = nUnitNucPerObj;
	}
	MALLOC(rMb2,double,NDIM*nPerObj);
	MALLOC(nDirMb2,double,NDIM*nPerObj);
	MALLOC(chMb2,int,6*nPerObj);
	MALLOC(eqLenMb2,double,6*nPerObj);
	MALLOC(unitMb2,int,3*nUnitPerObj);
	MALLOC(nDirUnitMb2,double,NDIM*nUnitPerObj);
	if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
		MALLOC(eqArUnitMb2,double,nUnitPerObj);
	}

	memset(chMb2, -1, (sizeof(int) * 6 * nPerObj));

	theta = DEG2RAD(26.56505117707799);
	sinThe = sin(theta);
	cosThe = cos(theta);
	// The lower and upper vertex of icosahedron
	V3SET(rIcos[0], 0, 0, -1.); 
	V3SET(rIcos[11], 0, 0, 1.); 	
	// Lower pentagon of icosahedron
	phi = PI / 5;
	for (n = 1; n < 6; n++) {
		V3SET(rIcos[n], cosThe * cos(phi), cosThe * sin(phi), -1 * sinThe);
		phi += 0.4 * PI;
	}
	// Upper pentagon of icosahedron
	phi = 0;
	for (n = 6; n < 11; n++) {
	  	V3SET(rIcos[n], cosThe * cos(phi), cosThe * sin(phi), sinThe);
	  	phi += 0.4 * PI;
	}
	for (n = 0; n < 12; n++) {
		V3SCALE(rIcos[n], rad);
		V3COPY(&P2(rMb2,n,0), rIcos[n]);
	}
	Copy1dArrayInt(unitMb2, initTriaArr, 60);

	nMb2 = 12;
	nUnitMb2 = 20;

	for(m = 0; m < lev; m++) {
		MALLOC(newVert.l,int,(int)(30 * pow(4, m + 1)) * 3);
		newVert.c = 0;
		for(n = 0; n < nUnitMb2; n++) { 
			V3COPY(vert, &P2A(unitMb2,n,0,3));
			// Determines the next 3 nodes for the the current triangle's 
			// subdivision
			for(k = 0; k < 3; k++){
				V2SET(ind, k, (k + 1) % 3);
				for(l = 0; l < newVert.c; l++) {
					pArr = &P2A(newVert.l,l,0,3);
					BREAK((pArr[0] == vert[ind[0]] && pArr[1] == vert[ind[1]])
							|| (pArr[0] == vert[ind[1]] 
							&& pArr[1] == vert[ind[0]]));
				}
				if (l < newVert.c) {
					vert[k + 3] = P2A(newVert.l,l,2,3);
				}
				else {
					vert[k + 3] = nMb2;
					V3ADD(r, &P2(rMb2,vert[ind[0]],0), 
							&P2(rMb2,vert[ind[1]],0));
					NormVec(r);
					V3SCALE(r, rad);
					V3COPY(&P2(rMb2,nMb2,0), r);
					V3SET(&P2A(newVert.l,newVert.c,0,3), vert[ind[0]], 
							vert[ind[1]], nMb2);
					(newVert.c)++;
					nMb2++;
				}
			}
			for(k = nUnitMb2 - 1; k >= n; k--){
				V3COPY(&P2A(unitMb2,k + 3,0,3), &P2A(unitMb2,k,0,3));
			}		
			for(k = 0; k < 12; k++){
				unitMb2[n * 3 + k] = vert[order[k]];
			}		
			nUnitMb2 += 3;
			n += 3;
		}
		free(newVert.l);
	}
	// Error check
	if (nPerObj != nMb2 || nUnitPerObj != nUnitMb2) {
		printf("Error: %d vs %d\n", nMb, nPerObj);
	}

	// Offset the membrane position to the center of a domain
	for(n = 0; n < nMb2; n++) {
		VS3COPY(&P2(nDirMb2,n,0), &P2(rMb2,n,0), -1.);
		NormVec(&P2(nDirMb2,n,0));
	}

	for(m = 0; m < nUnitMb2; m++) {
		qsort(&P2A(unitMb2,m,0,3), 3, sizeof(int), CompInt);
	    pArr = &P2A(unitMb2,m,0,3);
		CalcVec(dr[0], &P2(rMb2,pArr[1],0), &P2(rMb2,pArr[0],0));
		CalcVec(dr[1], &P2(rMb2,pArr[2],0), &P2(rMb2,pArr[0],0));
		V3CROSS(dr[2], dr[0], dr[1]);	
		if (V3COS(dr[2], &P2(nDirMb2,pArr[0],0)) < 0.) {
			SwapInt(&pArr[1], &pArr[2]);
			V3REVSIGN(dr[2]);
		}
		if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
			eqArUnitMb2[m] = 0.5 * V3LEN(dr[2]);
		}
		NormVec(dr[2]);
		V3COPY(&P2(nDirUnitMb2,m,0), dr[2]);
		
	    for(n = 0; n < 3; n++) {
	        ind[0] = pArr[n];
	        ind[1] = pArr[(n + 1) % 3];
	        len = CalcDist(&P2(rMb2,ind[0],0), &P2(rMb2,ind[1],0), 0);
	        for(k = 0; k < 2; k++) {
	            pArr2 = &P2A(chMb2,ind[k],0,6);
	            CS = -1;
	            for(l = 0; l < 6; l++) {
	                if (pArr2[l] == ind[1 - k]) {
	                    CS = 1;
	                    break;
	                }
	                BREAK(pArr2[l] == -1);
	            }
	            if (CS == -1 && l < 6) {
	                pArr2[l] = ind[1 - k];
	                P2A(eqLenMb2,ind[k],l,6) = len;
	            }
	        }
	    }
	}
 
	for(m = 0; m < nMb2; m++) {
		nCh = (P2A(chMb2,m,5,6) < 0) ? 5 : 6;
		CalcVec(dr[3], &P2(rMb2,P2A(chMb2,m,0,6),0), &P2(rMb2,m,0));
		V3CROSS(dr[0], dr[3], &P2(nDirMb2,m,0));
		for(n = 1; n < nCh; n++) {
			for(k = n + 1; k < nCh; k++) {
				for(l = 0; l < 2; l++) {
					CalcVec(dr[3], &P2(rMb2,P2A(chMb2,m,(l == 0) ? n : k,6),0), 
							&P2(rMb2,m,0));
					V3CROSS(dr[l + 1], dr[3], &P2(nDirMb2,m,0));
					ang[l] = V3ANG(dr[0], dr[l + 1]);
					V3CROSS(dr[3], dr[0], dr[l + 1]);
					if (V3DOT(&P2(nDirMb2,m,0), dr[3]) > 0.) {
						ang[l] = 2 * PI - ang[l];
					}
				}
				if (ang[0] > ang[1]) {
				 	SwapInt(&P2A(chMb2,m,n,6), &P2A(chMb2,m,k,6));
				 	SwapDbl(&P2A(eqLenMb2,m,n,6), &P2A(eqLenMb2,m,k,6));
				}	
			}
		}
	}
	for(m = 0; m < nObj; m++) {
		if (nObj == 1) {
			V3COPY(cenMb, dimDomH);
		}
		else if (nObj == 2) {
			VS3COPY(&P2(cenMb,m,0), dimDom, (m % nObj == 0) ? 0.25 : 0.75);
		}
		ind[0] = mm * nMbPerObj * nObjMbNuc / 2 + m * nPerObj;
		for(n = 0; n < nPerObj; n++) {
			V3ADD(&P2(rMb,ind[0] + n,0), &P2(rMb2,n,0), &P2(cenMb,m,0));
			V3COPY(&P2(nDirMb,ind[0] + n,0), &P2(nDirMb2,n,0));
			for(k = 0; k < 6; k++) {
				P2A(eqLenMb,ind[0] + n,k,6) = P2A(eqLenMb2,n,k,6);
				P2A(chMb,ind[0] + n,k,nChMb) = P2A(chMb2,n,k,6) 
						+ (P2A(chMb2,n,k,6) > -1 ? ind[0] : 0);
			}
	  		memset(&P2A(chMb,ind[0] + n,6,nChMb), -1, 
					sizeof(int) * (nChMb - 6));
			idxMb[ind[0] + n] = m + mm * nObj;
		}
		ind[1] = mm * nUnitMbPerObj * nObjMbNuc / 2 + m * nUnitPerObj;
		for(n = 0; n < nUnitPerObj; n++) {
			for(k = 0; k < 3; k++) {
				P2A(unitMb,ind[1] + n,k,3) = P2A(unitMb2,n,k,3) + ind[0];
			}
			V3COPY(&P2(nDirUnitMb,ind[1] + n,0), &P2(nDirUnitMb2,n,0));
			if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
				eqArUnitMb[ind[1] + n] = eqArUnitMb2[n];
			}
		}
	}

	min = POS_LARGE_VALUE;
	for(n = 0; n < nMb; n++) {
		for(k = 0; k < 6; k++) {
			CONT(P2A(chMb,n,k,nChMb) < 0);
			CONT(P2A(eqLenMb,n,k,6) > min);
			min = P2A(eqLenMb,n,k,6);
		}
	}
	free(rMb2);
	free(nDirMb2);
	free(eqLenMb2);
	free(chMb2);
	free(unitMb2);
	free(nDirUnitMb2);
	if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
		free(eqArUnitMb2);
	}
  }
}

/*------------------------------ Initialization ------------------------------*/

/*----------------------- Process transfered messages ------------------------*/

void UpdateMembraneDynamicsEventsSubroutine(ListInt *all) {
  int n, side, CS, CS2, mode;
  int mbInd, actInd, locMbInd, locActInd, mbInd2;
  int actRankMol, mbRankMol, mbRankMol2; 
  ListInt confMb;

  confMb.l = allIntL;
  confMb.c = 0;
  for(n = sendMbDyn.c; n < all->c; n++) {
	// mode = 0: unbind, 2: bind
	mode = (P2A(all->l,n,2,3) < 0) ? 0 : 1;
	mbInd = P2A(all->l,n,0,3);
	actInd = P2A(all->l,n,1,3);
	side = abs(P2A(all->l,n,2,3)) - 1;	
	locMbInd = iMb[mbInd];
	locActInd = iAct[actInd];
	if (locMbInd >= nMbMe || locMbInd < 0) {
		CS = 1;
		if (locActInd > -1) {
			actRankMol = CalcRankMolecule(&P2(act.r,locActInd,0));
			mbInd2 = act.mbReb[locActInd];
			if (mbInd2 > -1 && mbInd2 != mbInd) {
				CS2 = Find2ElementArray(confMb.l, confMb.c, locActInd, 
						side, 0, 2);
				if (CS2 == -1) {
					V2SET(&P2A(confMb.l,confMb.c,0,2), locActInd, side);
					(confMb.c)++;
				}
				if (iMb[mbInd2] > -1) {
					mbRankMol2 = CalcRankMolecule(&P2(memb.r,iMb[mbInd2],0));
					if (mbRankMol2 != actRankMol) {
						UpdateMembraneUnbindMatureSubroutine(actInd, 
								mbInd2, side, 1);
						if (locMbInd > -1) {
							mbRankMol 
									= CalcRankMolecule(&P2(memb.r,locMbInd,0));
							CS = (actRankMol == mbRankMol) ? 1 : 0;
						}
						else { CS = 1; }
					}
					else { CS = 0; }
				}
			}
			else if (mbInd2 < 0) {
				CS2 = Find2ElementArray(confMb.l, confMb.c, locActInd, 
						side, 0, 2);
				if (CS2 > -1) {
					if (locMbInd > -1) {
						mbRankMol = CalcRankMolecule(&P2(memb.r,locMbInd,0));
						CS = (actRankMol == mbRankMol) ? 1 : 0;
					}
					else { CS = 1; }
				}
				else { CS = 1; }
			}
			// Modify act.ch
			if (CS == 1) {
				if (mode == 0) {
					UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, 
							side, 2);
				}
				else {
					UpdateMembraneBindingSubroutine(actInd, mbInd, side, 1);
				}
				if (locMbInd >= nMbMe) { CS = 2; }
			}
		}
		if (locMbInd >= nMbMe && CS == 1) {
			if (mode == 0) {
				UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, side, 3);
			}
			else {
				UpdateMembraneBindingSubroutine(actInd, mbInd, side, 2);
			}
		}
	}
	if (locMbInd > -1 && locMbInd < nMbMe) {
		if (mode == 0) {
			UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, side, 10);
		}
	}
  }
  sendMbDyn.c = 0;
}
void UpdateMembraneDynamicsEvents(void) {
  int n, sizeArr;
  ListInt all;

  sizeArr = sendMbDyn.siz;
  sizeArr *= (mpiMethod == 0) ? cntAdjRank[0] : 2;
  MALLOC(all.l,int,sizeArr);

  all.c = sendMbDyn.c;
  for(n = 0; n < sendMbDyn.c; n++) {
	V3COPY(&P2A(all.l,n,0,3), &P2A(sendMbDyn.l,n,0,3));
  }
  CollectArrayIntFromAdjacentSubdomain(&all, 3);
  UpdateMembraneDynamicsEventsSubroutine(&all);

  free(all.l);
}

/*----------------------- Process transfered messages ------------------------*/

/*----------------- Interactions between membrane and others -----------------*/

int TrimNetworkForMembraneSubroutine(double len) {
  int CS = 1;

  if (sideMb == 0) {
	if (len < radMb - 0.6 * memb.thk) {
		if (gTglNuc != 0) {
			if (len > radNuc + 0.6 * memb.nucThk) {
				CS = 0;
			}
		}
		else { CS = 0; }
	}
  }
  else {
	if (len > radMb + 0.6 * memb.thk) { 
		CS = 0;
	}
  }
  return CS;
}

void TrimNetworkForMembrane(void) {
  int n, k, ind, CS;
  double len, dr[NDIM];
 
  FOR_ACT(n) {
	CalcVec(dr, &P2(rAct,n,0), dimDomH);
	if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
	len = V3LEN(dr);
	CS = TrimNetworkForMembraneSubroutine(len);
	CONT(CS != 1);
	DetachAbpOnActin(n);
	for(k = 0; k < 2; k++) {
		ind = P2A(chAct,n,k,nChAc);
		CONT(!(ind > -1));
		DetachAbpOnActin(ind);
		P2A(chAct,ind,1 - k,nChAc) = -1;
	}
	SetAllValue1dArrayInt(&P2A(chAct,n,0,nChAc), nChAc, 1);
  }
  FOR_ABP(n) {
	CalcVec(dr, &P2(rAbp,n,0), dimDomH);
	if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
	len = V3LEN(dr);
	CS = TrimNetworkForMembraneSubroutine(len);
	CONT(CS != 1);
	DetachActinOnAbp(n);
	V2SET_ALL(&P2A(chAbp,n,0,nChAb), -1);
	V2SET_ALL(&P2A(chAbp,n,3,nChAb), -1);
  }
}

int CheckActinAbpOverlapMembrane(double *r1, double *r2, int mode) { 
  int m, n, k, *cylInd, CS, ind, kind;
  double rPnt[5][NDIM], ratio[3];
  double len, distMb, distNuc, minDist, dr[NDIM], dr2[NDIM], dot;

  // Calculate the central position and length of the actin cylinder to put
  V3COPY(rPnt[0], r1);
  V3COPY(rPnt[1], r2);
  if (dimMbNuc == 2) { 
	for(n = 0; n < 2; n++) {
		rPnt[n][dirNormMbNuc] = dimDomH[dirNormMbNuc];
	}
  }
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	CONT(!(fabs(r1[k] - r2[k]) > 1.5 * dimDomH[k]));
	if (iCell[k] == 0) { 
		ind = (r1[k] > r2[k]) ? 0 : 1;
		rPnt[ind][k] -= dimDom[k];
	}
	else {
		ind = (r1[k] < r2[k]) ? 0 : 1;
		rPnt[ind][k] += dimDom[k];
	}
  }
  // To figure out which is the inner side, a reference point is required.
  // During network formation, the center of the entire domain is the refernece
  // point. 
  
  len = (mode == 0) ? actF.dia : abpF.repDia[mode - 1];
  distMb = AVG2(len, memb.thk);
  distNuc = AVG2(len, memb.nucThk);
  CS = 1;
  if (currTimeStep <= netForm.dur) {
	for(m = 0; m < nObjMbNuc; m++) {
		if (sideMb == 0) { CS = 1; }
		for(n = 0; n < 2; n++) { 
			len = CalcDist(rPnt[n], &P2(cenMb,m,0), 0);
			if (sideMb == 0) {
//				if (len < radMb - distMb && len > (radMb - distMb) 
				if (len < radMb - 0.05 && len > (radMb - distMb) 
						* (1. - thkMbActNuc)) {
					if (gTglNuc != 0) { 
						CONT(len > radNuc + distNuc);
					}
					else { continue; }
				}
			}
			else {
				CONT(len > radMb + distMb);
			}
			CS = 0;
			break;
		}
		if (sideMb == 0) { BREAK(CS == 1); }
		else { BREAK(CS == 0); }
	}
  }
  else {
	if (memb.unit.c + memb.unitCp.c > 0) {
		minDist = POS_LARGE_VALUE;
		for(n = 0; n < memb.unit.c + memb.unitCp.c ; n++) {
			cylInd = (n < memb.unit.c) ? &P2A(memb.unit.l,n,0,dimMbNuc)
					: &P2A(memb.unitCp.l,n - memb.unit.c,0,dimMbNuc);
			CONT(iMb[cylInd[0]] < 0 || iMb[cylInd[1]] < 0);
			if (dimMbNuc == 3) { CONT(iMb[cylInd[2]] < 0);	}
			for(m = 0; m < dimMbNuc; m++) {
			    V3COPY(rPnt[m + 2],&P2(memb.r,iMb[cylInd[m]],0));
			    FOR_NDIM(k) {
					CONT(pbc[k] != 1);
			        CONT(fabs(rPnt[m + 2][k] - rPnt[0][k]) <= 1.5 * dimDomH[k]);
		            rPnt[m + 2][k] += dimDom[k]
		                    * ((rPnt[m + 2][k] < rPnt[0][k]) ? 1. : -1.);
				}
			}
			if (dimMbNuc == 2) {
				len = CalcSegSegDist(&rPnt[2], &rPnt[0], dr, ratio, 0);
			}
			else {
				len = CalcTriaSegDist(&rPnt[2], &rPnt[0], dr, ratio);
			}
			// If the distance between two segments is smaller than the sum of
			// their radii (overlapping), repulsive forces should act.
			if (len < (ISNUC(memb.idx[iMb[cylInd[0]]]) ? distNuc : distMb)) {
				CS = 0;
				break;
			}
			// Find the nearest membrane unit
			if (len < minDist) { 
				ind = n;
				kind = ISNUC(memb.idx[iMb[cylInd[0]]]) ? 1 : 0;
				minDist = len;
				V3COPY(dr2, dr);
			}
		}
		if (CS == 1) {
			if (ind < memb.unit.c) {
				dot = V3DOT(&P2(memb.unitNDir,ind,0), dr2);
			}
			else {
				CalcMembUnitNormalDirec(&P2A(memb.unitCp.l,
						ind - memb.unit.c,0,dimMbNuc), dr);
				dot = V3DOT(dr, dr2);
			}
			if (kind == 1) {
				if (dot < 0.) { CS = 0; }
			}
			else {
				if ((sideMb == 0 && dot > 0.) || (sideMb == 1 && dot < 0.))	{ 
					CS = 0; 
				}
			}
		}
    }
	else {
		if (nActMe == 0) { CS = 0; }
	}
  }
  return CS;
}

double HowManyInOutMembrane(double *rBnd) {
  int k, cnt, nObj;
  int ind[NDIM], begin[NDIM], end[NDIM];
  double dist, thres[2], r[NDIM], resol, min, max, bound, cenDom[NDIM];

  resol = 0.5;
  if (bnd.gTglRnd != 0) {
	FOR_NDIM(k) {
    	cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
	}
  }
  if (sideMb == 0) {
	thres[0] = (radMb - 0.5 * memb.thk) * (1. - thkMbActNuc);
	thres[1] = radMb - 0.5 * memb.thk;
	if (gTglNuc != 0) {
		if (thres[0] < radNuc + 0.5 * memb.nucThk) {
			thres[0] = radNuc + 0.5 * memb.nucThk;
		}
	}
  }
  else {
	thres[0] = radMb + 0.5 * memb.thk;
	thres[1] = POS_LARGE_VALUE;
  }
  nObj = nObjMbNuc / (gTglNuc != 0 ? 2 : 1);

  if (sideMb == 0) {
  min = POS_LARGE_VALUE;
  max = NEG_LARGE_VALUE;
  for(k = 0; k < nObj; k++) {
	min = (min > P2(cenMb,k,0) - radMb) ? P2(cenMb,k,0) - radMb : min;
	max = (max < P2(cenMb,k,0) + radMb) ? P2(cenMb,k,0) + radMb : max;
  }

  for(k = 0; k < dimMbNuc; k++) {
	bound = (min > P2A(rBnd,0,dirMbNuc[k],NDIM)) 
			? min : P2A(rBnd,0,dirMbNuc[k],NDIM);
	begin[k] = (int)ceil(bound / resol);
	bound = (max < P2A(rBnd,1,dirMbNuc[k],NDIM)) 
			? max : P2A(rBnd,1,dirMbNuc[k],NDIM);
	end[k] = (int)floor(bound / resol);
  }
  }
  else {
	for(k = 0; k < dimMbNuc; k++) {
		begin[k] = (int)ceil(P2A(rBnd,0,dirMbNuc[k],NDIM) / resol);
		end[k] = (int)floor(P2A(rBnd,1,dirMbNuc[k],NDIM) / resol);
	}
  }

  cnt = 0;
  if (dimMbNuc == 3) {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			for(ind[2] = begin[2]; ind[2] < end[2]; ind[2]++) {	
				VS3COPY(r, ind, resol);
				if (bnd.gTglRnd != 0) {
					dist = CalcDist(cenDom, r, 0);
					CONT(dist > bnd.radRnd);
				}
				for(k = 0; k < nObj; k++) {
					dist = CalcDist(&P2(cenMb,k,0), r, 0);
					if (dist >= thres[0] && dist <= thres[1]) { 
						if (sideMb == 0) {
							cnt++; 
							break;
						}
						else { continue; }
					}
					if (sideMb == 1) { 
						cnt++;
						break;
					}
				}
			}
		}
	}
	return (cnt * CUBE(resol));
  }
  else {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			r[dirMbNuc[0]] = ind[0] * resol;
			r[dirMbNuc[1]] = ind[1] * resol;
			r[dirNormMbNuc] = cenMb[dirNormMbNuc];
			if (bnd.gTglRnd != 0) {
				dist = CalcDist(cenDom, r, 0);
				CONT(dist > bnd.radRnd);
			}
			for(k = 0; k < nObj; k++) {
				dist = CalcDist(&P2(cenMb,k,0), r, 0);
				if (dist >= thres[0] && dist <= thres[1]) { 
					if (sideMb == 0) {
						cnt++; 
						break;
					}
					else { continue; }
				}
				if (sideMb == 1) { 
					break;
				}
			}
			if (sideMb == 1 && k == nObj) {
				cnt++;
			}
		}
	}
	return (cnt * SQR(resol));
  }
}

/*----------------- Interactions between membrane and others -----------------*/

