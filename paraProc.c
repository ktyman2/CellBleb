// This file contains functions for parallel processing.

/*----------------------- Related to handling long chains --------------------*/

// Subroutine for UpdateLongChainNormalSubroutine()
// This finds the ranks where a particle is possibly located at a next time 
// step based on its location. The ranks found here are stored in pRank.
void UpdateLongChainNormalSubSubroutine1(double *rMol, ListInt *pRank, 
		double limit) {
  int k, n, oftInd[][2] = {{-1, 0}, {0}, {0, 1}}, cntOftInd[] = {2, 1, 2};
  int oft[NDIM], oft2[NDIM], offset[NDIM], idx[NDIM];
  double r[NDIM], bound;

  pRank->c = 0;
  V3COPY(r, rMol);
  if (rheoWay > 0 && bulkRheoType == 0) { ConvertRectDomainVector(r, 0); }
  V3SET_ALL(offset, 1);
  V3SET_ALL(oft, 1);

  // oft[] indicates a current rank, and oft2 represents possible ranks.
  FOR_NDIM(k) {
	if (iCell[k] > 0 && iCell[k] < nCell[k] - 1) {
		if (r[k] < P2A(bnd.r,0,k,NDIM)) { oft[k] = 0; }
		else if (r[k] > P2A(bnd.r,1,k,NDIM)) { oft[k] = 2; }
		for(n = 0; n < 2; n++) {
			CONT(!(fabs(r[k] - P2A(bnd.r,n,k,NDIM)) < limit));
			if (r[k] < P2A(bnd.r,n,k,NDIM)) { offset[k] = 2; }
			else { offset[k] = 0; }
		}
	}
	else if ((iCell[k] == 0 || iCell[k] == nCell[k] - 1) && nCell[k] > 1) {
		if (iCell[k] == 0) {
			if (r[k] >= rGrid[k][1] && r[k] < rGrid[k][2]) { oft[k] = 2; }
			else if (r[k] >= rGrid[k][2]) { oft[k] = 0; }
		}
		else {
			if (r[k] >= rGrid[k][nGrid[k] - 3] && r[k] < rGrid[k][nGrid[k] - 2])
			{ oft[k] = 0; }
			else if (r[k] < rGrid[k][nGrid[k] - 3]) { oft[k] = 2; }
		}
        for(n = 0; n < 2; n++) {
			if (n == 0) {
				bound = (iCell[k] == 0) ? rGrid[k][nGrid[k] - 1]
						: rGrid[k][nGrid[k] - 2];
			}
			else {
				bound = (iCell[k] == 0) ? rGrid[k][1] : rGrid[k][0];
			}
			if (r[k] > rGrid[k][nGrid[k] - 1] - limit) { offset[k] = 2; }
			else if (r[k] < rGrid[k][0] + limit) { offset[k] = 0; }
			else if (fabs(r[k] - bound) < limit) {
				if (r[k] < bound) { offset[k] = 2; }
				else { offset[k] = 0; }
			}
		}
	}
  }
  for(idx[0] = 0; idx[0] < cntOftInd[offset[0]]; idx[0]++) {
	for(idx[1] = 0; idx[1] < cntOftInd[offset[1]]; idx[1]++) {
		for(idx[2] = 0; idx[2] < cntOftInd[offset[2]]; idx[2]++) {
			FOR_NDIM(k) { 
				oft2[k] = oft[k] + oftInd[offset[k]][idx[k]]; 
				if (oft2[k] == 3) { oft2[k] = 1; }
				else if (oft2[k] == -1) { oft2[k] = 1; }
			}
			V3IND_BACK_CONST_INT(pRank->l[pRank->c], oft2, NDIM);
			(pRank->c)++;
		}
	}
  }
}

// Subroutine for UpdateLongChainNormalSubroutine()
void UpdateLongChainNormalSubSubroutine2(int rankMol, int rankL, int ind1, 
		int ind2, int ind3) {
  int CS;

  if (rankMol != rank) {
	// "longChIntMsg.l" has the list of actins which a 
	// current subdomain has to transfer to other subdomains 
	// where ABP can belong at next step
	if (rankL == 13) {
		CS = Find2ElementArray(longChIntMsg.l,
				longChIntMsg.c, ind1, ind2, 0, 2);
		if (CS < 0) {
			V2SET(&P2A(longChIntMsg.l,longChIntMsg.c,0,2), ind1, ind2);
			(longChIntMsg.c)++;
		}
	}
	// "longChExtMsg.l" has the list of actins which the 
	// other subdomains have to transfer to subdomains where 
	// ABP can belong at next step
	else {
		if (ind2 != ind3) {
			CS = Find3ElementArray(longChExtMsg.l,
					longChExtMsg.c, ind1, ind2, ind3, 0, 3);
			if (CS < 0) {
				V3SET(&P2A(longChExtMsg.l,longChExtMsg.c,
						0,3), ind1, ind2, ind3);
				(longChExtMsg.c)++;
			}
		}
	}
  }
}

// Subroutine for UpdateLongChainNormal()
void UpdateLongChainNormalSubroutine(int ind1, int ind2, double len) {
  int n, k, CS, CS2; 
  int abpRankInd, actRankInd, abpRankMol, actRankMol, ind;
  double r[NDIM], *pR;
  ListInt abpRank, actRank;

  MALLOC(abpRank.l, int, nCpu);
  MALLOC(actRank.l, int, nCpu);

  CS = 0; 
  for(n = 0; n < 2; n++) { 
	ind = (n == 0) ? ind1 : ind2;
	if (ind < nAct) { V3COPY(r, &P2(act.r,iAct[ind],0)); }
	else if (ind >= nAct && ind < nAct + nAbp) 
	{ V3COPY(r, &P2(abp.r,iAbp[ind - nAct],0)); }
	else { V3COPY(r, &P2(memb.r,iMb[ind - nAct - nAbp],0)); }
	if (rheoWay > 0 && bulkRheoType == 0) 
	{ ConvertRectDomainVector(r, 0); }
	// either of two elements should be placed near boundary
	FOR_NDIM(k) {
		CONT(r[k] >= P2A(bnd.r,0,k,NDIM) + P2A(edge,0,k,NDIM) 
				&& r[k] < P2A(bnd.r,1,k,NDIM) - P2A(edge,1,k,NDIM));
		CS++;	
		break; 
	}
	BREAK(CS > 0);
  }

  // If the current chain length is long enough, they should be in longCh.l
  if (CS > 0 && len >= neiEdge * 0.9) {
	CS2 = InsertLongChain(ind1, ind2, len * 1.1);
	// If it is in longCh.l already, the chain length is just updated
	if (CS2 > -1) { longChDist[CS2] = len * 1.1; }
  }
  // If chain length is not long enough, or if both of them exist at outside
  // of synchronized regions, they are removed from longCh.l
  else { DeleteLongChain(ind1, ind2); }

  // If this is an ABP-actin chain, another actin located in a barbed direction
  // should be transferred for calculation of bending in actin-ABP-actin 
  // or right angles between ABP arms and axes of actin filaments.
  if (CS > 0 && ind1 >= nAct && ind1 < nAct + nAbp && ind2 < nAct) {
	// Find rank to which ABP or actin belongs 
	abpRankMol = CalcRankMolecule(&P2(abp.r,iAbp[ind1 - nAct],0));
	abpRank.c = 0;
	// Find all possible ranks to which ABP or actin can belong at next step
	UpdateLongChainNormalSubSubroutine1(&P2(abp.r,iAbp[ind1 - nAct],0), 	
			&abpRank, maxDisp);
	ind = P2A(act.ch,iAct[ind2],0,nChAc);
	actRank.c = 0;
	if (iAct[ind] < 0) {
		actRankMol = -1;
		// If the current subdomain doesn't have information of actin,
		// it scans all possible subdomains using other actin. 
		UpdateLongChainNormalSubSubroutine1(&P2(act.r,iAct[ind2],0), 
				&actRank, maxActCh);
	}
	else {
		actRankMol = CalcRankMolecule(&P2(act.r,iAct[ind],0)); 
		// Find all possible ranks to which actin can belong at next step
		UpdateLongChainNormalSubSubroutine1(&P2(act.r,iAct[ind],0),
				&actRank, maxDisp);
	}
	for(n = 0; n < actRank.c; n++) {
		actRankInd = P2A(adjRank,actRank.l[n],0,2);
		// If the subdomain does not exist, skip this procedure
		CONT(actRank.l[n] != 13 && actRankInd < 0);
		for(k = 0; k < abpRank.c; k++) {
			abpRankInd = P2A(adjRank,abpRank.l[k],0,2);
			// If the subdomain does not exist, skip this procedure
			CONT(abpRank.l[k] != 13 && abpRankInd < 0);
			// If ABP is not in the current subdomain
			UpdateLongChainNormalSubSubroutine2(abpRankMol, actRank.l[n],
					ind, abpRankInd, actRankInd);
			// If actin belongs to other subdomains
			UpdateLongChainNormalSubSubroutine2(actRankMol, abpRank.l[k],
					ind1, actRankInd, abpRankInd);
		}
	}
  }
  free(actRank.l);
  free(abpRank.l);
}

// normal method.
void UpdateLongChainNormal(void) {
  int n, k, nBd, side;
  int abpInd, actInd, locActInd, locAbpInd, locMbInd;
  int *pArr, *chkL;
  double len;

  chkL = allIntL;
  memset(chkL, -1, sizeof(int) * nActMe * 2);

  CheckArraySize(&longChExtMsg, &longChExtMsg.siz, 3, 0);
  CheckArraySize(&longChIntMsg, &longChIntMsg.siz, 2, 0);

  // Initialize messgages
  longChExtMsg.c = 0;
  longChIntMsg.c = 0;

  FOR_ABPME(n) {
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	// If motors don't walk, inactive motors cannot have long chains, so
	// they don't need to be considered here
	for (k = 0; k < 2; k++) {
		actInd = pArr[k];
		CONT(actInd < 0);
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
		len = CalcDist(&P2(abp.r,n,0), &P2(act.r,locActInd,0), 0);
		UpdateLongChainNormalSubroutine(abp.id[n] + nAct, actInd, len);
	}
	if (ISMTF(K_ABP(n)) && motSA.cenDist * 1.5 > neiEdge) {
		if (abp.mId[n] == 0) {
			for(k = 0; k < 2; k++) {
				abpInd = pArr[k + 3];
				CONT(abpInd < 0);
				locAbpInd = iAbp[abpInd];
				CONT(locAbpInd < 0);
				CONT(abp.mId[locAbpInd] != 0);
				len = CalcDist(&P2(abp.r,locAbpInd,0), &P2(abp.r,n,0), 0);
				UpdateLongChainNormalSubroutine(abp.id[n] + nAct, 
						abpInd + nAct, len);
			}	
		}
	}
  }

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	// Actins usually have very high extensional stiffness. However, If a chain
	// between actins can be extended a lot due to motors or external forces. 
	// The following part takes care of such long actin chains.
	if (maxActCh > neiEdge * 1.25) {
		for(k = 0; k < 2; k++) {
			CONT(P2A(chkL,n,k,2) == 1);
			actInd = P2A(act.ch,n,k,nChAc);
			CONT(actInd < 0);
			locActInd = iAct[actInd];
			CONT(locActInd < 0);
			len = CalcDist(&P2(act.r,locActInd,0), &P2(act.r,n,0), 0);
			UpdateLongChainNormalSubroutine(actInd, act.id[n], len);
			if (locActInd < nActMe) {
				P2A(chkL,locActInd,1 - k,2) = 1;
			}
		}
	}

	// This is for a case with actin in a current subdomain and 
	// ABP in other subdomain
	for(k = 0; k < nChAc - 2; k++) {
		abpInd = P2A(act.ch,n,k + 2,nChAc);
		CONT(abpInd < 0);
		locAbpInd = iAbp[abpInd];
		CONT(locAbpInd < nAbpMe);
		len = CalcDist(&P2(abp.r,locAbpInd,0), &P2(act.r,n,0), 0);
		UpdateLongChainNormalSubroutine(abpInd + nAct, act.id[n], len);
	}
  }
  if (gTglMb != 0) {
	nBd = nChMb - nMbAct;
	memset(chkL, -1, sizeof(int) * nMbMe * nBd);
	FOR_MBME(n) {
		pArr = &P2A(memb.ch,n,0,nChMb);
		for(k = 0; k < nBd; k++) {
			CONT(pArr[k] < 0);
			locMbInd = iMb[pArr[k]];
			CONT(locMbInd < 0);
			CONT(P2A(chkL,n,k,nBd) != -1);		
			len = CalcDist(&P2(memb.r,n,0), &P2(memb.r,locMbInd,0), 0);
			UpdateLongChainNormalSubroutine(memb.id[n] + nAct + nAbp, 
					pArr[k] + nAct + nAbp, len);
			if (locMbInd < nMbMe) {
				side = FindElementArray(&P2A(memb.ch,locMbInd,0,nChMb),
						nBd, memb.id[n], 0, 1);
				P2A(chkL,locMbInd,side,nBd) = 1;
			}
		}
		for(k = nChMb - nMbAct; k < nChMb; k++) {
			actInd = pArr[k];
			CONT(actInd < 0);
			locActInd = iAct[actInd];		
			CONT(locActInd < 0);
			len = CalcDist(&P2(memb.r,n,0), &P2(act.r,locActInd,0), 0);
			UpdateLongChainNormalSubroutine(memb.id[n] + nAct + nAbp, 
					actInd, len);
		}
	}
  }
}
	
// Update "longChIntMsg" with Plympton method.
void UpdateLongChainPlympton(void) {
  int n, k, actRank, abpRank;
  int actInd, actInd2, abpInd, locActInd, locActInd2, locAbpInd;

  // Initialize messgages
  longChIntMsg.c = 0;
  for (n = 0; n < nAbpMe; n++) {
	for (k = 0; k < 2; k++) {
		actInd = P2A(abp.ch,n,k,nChAb);
		CONT(actInd < 0);
		locActInd = iAct[actInd];
		CONT(locActInd < 0);
		actInd2 = P2A(act.ch,locActInd,0,nChAc);
		CONT(actInd2 < 0);
		locActInd2 = iAct[actInd2];
		CONT(locActInd2 > -1 && locActInd2 < nActMe);
		InsertElement1dArrayWChk(longChIntMsg.l, &longChIntMsg.c, actInd2);
	}
  }
  FOR_ACTME(n) {
	actInd = P2A(act.ch,n,0,nChAc);
	CONT(actInd < 0);
	locActInd = iAct[actInd];
	CONT(locActInd < 0);
	actRank = CalcRankMolecule(&P2(act.r,locActInd,0));
	for(k = 0; k < nChAc - 2; k++) {
		abpInd = P2A(act.ch,n,k + 2,nChAc);
		CONT(abpInd < 0);
		locAbpInd = iAbp[abpInd];
		CONT(locAbpInd < nAbpMe);
		abpRank = CalcRankMolecule(&P2(abp.r,locAbpInd,0));
		CONT(actRank == abpRank);
		InsertElement1dArrayWChk(longChIntMsg.l, &longChIntMsg.c, actInd);
		break;
	}
  }
}

// Delete an element in "longCh.l"
void DeleteLongChain(int ind1, int ind2) {
  int n, k, *pArr;

  for(n = 0; n < longCh.c; n++) {
	pArr = &P2A(longCh.l,n,0,2);
	CONT(!((pArr[0] == ind1 && pArr[1] == ind2) 
			|| (pArr[0] == ind2 && pArr[1] == ind1)));
	DeleteElementArrayByIndex(longCh.l, &longCh.c, n, 2);
	for(k = n; k < longCh.c; k++) {
		longChDist[k] = longChDist[k + 1];
	}
	break;			
  }
}

// Insert an element in "longCh.l"
int InsertLongChain(int ind1, int ind2, double len) {
  int CS, ind[2];

  if (ind1 > ind2) { V2SET(ind, ind1, ind2); }
  else { V2SET(ind, ind2, ind1);  }
  CS = Find2ElementArray(longCh.l, longCh.c, ind[0], ind[1], 0, 2);
  if (CS < 0) {
	V2SET(&P2A(longCh.l,longCh.c,0,2), ind[0], ind[1]);
	longChDist[longCh.c] = len;
	(longCh.c)++;
  }
  return CS;
}

/*----------------------- Related to handling long chains --------------------*/

/*-------------------------------- Subdomains --------------------------------*/

// This function adjusts the lines which divide a whole domain into
// several subdomains. The aim is to make each subdomain take relatively 
// the same number of particles to maintain computational efficiency.
void UpdateSubdomSectionLocation(int mode) {
  int m, n, k, CS, CS2, end, *nActAll, *nAbpAll, *nMbAll;
  int *tagAct, *tagAbp, *tagActMe, *tagAbpMe;
  int iGrid[NDIM], nGridInd[NDIM], dire[NDIM];
  double *arrGrid, dimDomC[NDIM], rBnd[NDIM * 2];
  double outThick[NDIM], inThick[NDIM], lThick, rThick, disp, nPar[2], mbMaxSiz;
  double **rGridPrev;

  // Displacement of grid lines at each trial
  disp = 0.5;
  if (mode == 0) {
	// Gather the information of positions required for adjustment.
	MALLOC(nActAll,int,nCpu);
	MALLOC(nAbpAll,int,nCpu);
	MALLOC(tagActMe,int,nActMe);
	MALLOC(tagAbpMe,int,nAbpMe);
	MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (gTglMb != 0) { 
		MALLOC(nMbAll,int,nCpu);
		MPI_Gather(&nMbMe, 1, MPI_INT, nMbAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	if (rank == 0) {
		MALLOC(rAct,double,nAct*NDIM);
		MALLOC(rAbp,double,nAbp*NDIM);
	    MALLOC(tagAct,int,nAct);
	    MALLOC(tagAbp,int,nAbp);
		if (gTglMb != 0) { MALLOC(rMb,double,nMb*NDIM); }
	}
	Gather2dArrayDouble(nActMe, nActAll, NDIM, act.id, act.r, rAct);
	Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM, abp.id, abp.r, rAbp);
	if (gTglMb != 0) {
		Gather2dArrayDouble(nMbMe, nMbAll, NDIM, memb.id, memb.r, rMb);
	}
	FOR_ACTME(n) {
	    tagActMe[n] = (!ISACTM(n)) ? 1 : 0;
	}
	FOR_ABPME(n) {
		tagAbpMe[n] = (!ISABPIM(n)) ? 1 : 0;
	}
	Gather2dArrayInt(nActMe, nActAll, 1, act.id, tagActMe, tagAct);
	Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, tagAbpMe, tagAbp);
  }
  end = 1;
  FOR_NDIM(k) { 
	if (nGrid[k] > 2) { 
		end *= nGrid[k] - 2; 
		nGridInd[k] = nGrid[k] - 2;
	}
	else { nGridInd[k] = 1; }
  }
  MALLOC(arrGrid,double,V3SUM(nGrid) - 6);

  // Only main CPU calculates the new boundaries.
  if (rank == 0) {
	MALLOC2(rGridPrev,double,NDIM);
	FOR_NDIM(k) {
		MALLOC(rGridPrev[k],double,nGrid[k]);
		for(n = 0; n < nGrid[k]; n++) {
			rGridPrev[k][n] = rGrid[k][n];
		}
	}
	FOR_NDIM(k) {
		outThick[k] = neiEdge + ((pbc[k] == 1 && nGrid[k] > 2) ? neiEdge : 0.);
		inThick[k] = neiEdge * 2.;
	}
	for(m = 0; m < end; m++) {
		V3IND_ASSIGN_INT(m, nGridInd, 1, iGrid);
		FOR_NDIM(k) { iGrid[k] += 1; }
		CS2 = 1;
		do {
			CS = 0;
			FOR_NDIM(k) {
				if (nGrid[k] <= 2) {
					CS++;
					continue;
				}
				lThick = (iGrid[k] == 1) ? outThick[k] : inThick[k];
				rThick = (iGrid[k] == nGrid[k] - 2) ? outThick[k] : inThick[k];
                if (rGrid[k][iGrid[k]] + rThick >= rGridPrev[k][iGrid[k] + 1]) {
					if (CS2 == 0 && dire[k] == 1) { CS++; }
					else { rGrid[k][iGrid[k]] -= disp; }
				}
				else if (rGrid[k][iGrid[k]]
						<= rGridPrev[k][iGrid[k] - 1] + lThick) {
					if (CS2 == 0 && dire[k] == -1) { CS++; }
					else { rGrid[k][iGrid[k]] += disp; }
				}
				// Count the number of particles belonging to two spaces 
				// divided by a line corresponding to "rGrid[k][iGrid[k]]".
				V2ZERO(nPar);
				if (mode == 0) {
					FOR_ACT(n) {
						CONT(tagAct[n] != 1);
						nPar[(P2(rAct,n,k) < rGrid[k][iGrid[k]]) ? 0 : 1] += 1.;
					}
					FOR_ABP(n) {
						CONT(tagAbp[n] != 1);
						nPar[(P2(rAbp,n,k) < rGrid[k][iGrid[k]]) ? 0 : 1] += 1.;
					}
					if (gTglMb != 0) {
						FOR_MB(n) {
							nPar[(P2(rMb,n,k) < rGrid[k][iGrid[k]]) ? 0 : 1] 
									+= 1.;
						}
					}
				}
				else {
					for(n = 0; n < NDIM; n++) {
						P2A(rBnd,0,n,NDIM) = 0;
						P2A(rBnd,1,n,NDIM) = dimDom[n];
					}
					P2A(rBnd,1,k,NDIM) = rGrid[k][iGrid[k]];
					nPar[0] = HowManyInOutMembrane(rBnd);
					P2A(rBnd,0,k,NDIM) = rGrid[k][iGrid[k]];
					P2A(rBnd,1,k,NDIM) = dimDom[k];
					nPar[1] = HowManyInOutMembrane(rBnd);
				}
				nPar[0] *= (double)(nGrid[k] - iGrid[k] - 1);
				nPar[1] *= (double)iGrid[k];
				// If the ratio of two numbers is not ideal, move the line
				// in either of two directions.
				if (CS2 == 1) { 
					dire[k] = (nPar[0] < nPar[1]) ? 1 : -1;
				}
				if (nPar[0] < nPar[1] && dire[k] == 1) {
					if (rGrid[k][iGrid[k]] + rThick + disp
							< rGrid[k][iGrid[k] + 1]) { 
						rGrid[k][iGrid[k]] += disp; 
					}
					else { CS++; }
				}
				else if (nPar[0] > nPar[1] && dire[k] == -1) { 
					if (rGrid[k][iGrid[k]] 
							> rGrid[k][iGrid[k] - 1] + lThick + disp) { 
						rGrid[k][iGrid[k]] -= disp; 
					}
					else { CS++; }
				}
				else { CS++; }
			}
			CS2 = 0;
		} while(CS < 3);
	}
	// Pack the calculated new boundaries
	for(n = 0; n < nGrid[0] - 2; n++) { arrGrid[n] = rGrid[0][n + 1]; }
	for(n = 0; n < nGrid[1] - 2; n++) 
	{ arrGrid[n + nGrid[0] - 2] = rGrid[1][n + 1]; }
	for(n = 0; n < nGrid[2] - 2; n++) 
	{ arrGrid[n + nGrid[0] + nGrid[1] - 4] = rGrid[2][n + 1]; }
  }
  MPI_Bcast(arrGrid, V3SUM(nGrid) - 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Unpack the calculated new boundaries
  if (rank != 0) {
	for(n = 0; n < nGrid[0] - 2; n++) { rGrid[0][n + 1] = arrGrid[n]; }
	for(n = 0; n < nGrid[1] - 2; n++) 
	{ rGrid[1][n + 1] = arrGrid[n + nGrid[0] - 2]; }
	for(n = 0; n < nGrid[2] - 2; n++) 
	{ rGrid[2][n + 1] = arrGrid[n + nGrid[0] + nGrid[1] - 4]; }
  }
  // Post-process related to neighboring list
  minDimDomC = POS_LARGE_VALUE;
  if (gTglMb != 0) {
	mbMaxSiz = (abpF.maxSiz > memb.len) ? abpF.maxSiz : memb.len;
  }
  if (gTglNuc != 0) {
	if (mbMaxSiz < memb.nucLen) { mbMaxSiz = memb.nucLen; }
  }
  FOR_NDIM(k) {
	if (nCell[k] > 1) {
		P2A(bnd.r,0,k,NDIM) = rGrid[k][iCell[k]];
		P2A(bnd.r,1,k,NDIM) = rGrid[k][iCell[k] + 1];
		dimDomC[k] = P2A(bnd.r,1,k,NDIM) - P2A(bnd.r,0,k,NDIM);

		cell.n[k] = (int)((dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM))
	            / ((abpF.maxSiz + DT_NL_UPDATE) * 1.2));
	    cell.wid[k] = (dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM))
	            / (double)cell.n[k];

	    cell.base[k] = P2A(bnd.r,0,k,NDIM) - P2A(edge,0,k,NDIM);
		if (gTglMb != 0) {
			cellMb.n[k] = (k == dirNormMbNuc && dimMbNuc == 2) ? 1 : 
					(int)((dimDomC[k] + P2A(edge,0,k,NDIM) + P2A(edge,1,k,NDIM))
					/ ((mbMaxSiz + DT_NL_UPDATE) * 1.2));
			cellMb.wid[k] = (dimDomC[k] + P2A(edge,0,k,NDIM) 
					+ P2A(edge,1,k,NDIM)) / (double)cellMb.n[k];
			cellMb.base[k] = cell.base[k];
		}
	}
	if (dimDomC[k] < minDimDomC && k != dir2D) { minDimDomC = dimDomC[k]; }
  }
  maxNeiLenSq = POS_LARGE_VALUE;
  free(arrGrid);  
  if (rank == 0) {
    FOR_NDIM(k) { free(rGridPrev[k]); }
    free(rGridPrev);
  }
  if (mode == 0) {
	free(nActAll);
	free(nAbpAll);
	free(tagActMe);
	free(tagAbpMe);
	if (gTglMb != 0) { free(nMbAll); }
	if (rank == 0) {
		free(rAct);
		free(rAbp);
	    free(tagAct);
	    free(tagAbp);
		if (gTglMb != 0) { free(rMb); }
	}
  }
}

void RecordSubdomSectionLocation (void) {
  int n, k;
  FILE *fOut;

  fOut = fopen(GenFileName("BndDiv"), "a");
  fprintf(fOut, "%lld\t", currTimeStep);
  FOR_NDIM(k) {
	for(n = 0; n < nGrid[k]; n++) {
		fprintf(fOut, "%g\t", rGrid[k][n]);
	}
  }
  fprintf(fOut, "\n");
  fclose(fOut);	
}

/*-------------------------------- Subdomains --------------------------------*/

/*--------------------- Move particles between subdomains --------------------*/

// Subroutine for MovePartilces()
void MoveParticlesSubroutine3(ListInt *unitMbMv, double *unitMbInfoMv, 
		int ind, int nInfo) {
  int n, k, loc, *pArr; 

  for(n = memb.unit.c - 1; n >= 0; n--) {
	pArr = &P2A(memb.unit.l,n,0,dimMbNuc);
	for(k = 0; k < dimMbNuc; k++) {
		BREAK(pArr[k] == ind);
	}
	CONT(k == dimMbNuc);
	loc = k;
	Copy1dArrayInt(&P2A(unitMbMv->l,unitMbMv->c,0,dimMbNuc), pArr, dimMbNuc);
	V3COPY(&P2A(unitMbInfoMv,unitMbMv->c,0,nInfo), &P2(memb.unitNDir,n,0));
    if (tglMbNucLocAre != 0) { 
		P2A(unitMbInfoMv,unitMbMv->c,3,nInfo) = memb.unitEqAr[n];
	}
	(unitMbMv->c)++;
	for(k = 1; k < dimMbNuc; k++) {
		BREAK(iMb[pArr[(loc + k) % dimMbNuc]] > -1 
				&& iMb[pArr[(loc + k) % dimMbNuc]] < nMbMe);
	}
	CONT(k < dimMbNuc);
	DeleteElementInNeighborList(n, 2);
	for(k = n; k < memb.unit.c - 1; k++) {
		V3COPY(&P2(memb.unitNDir,k,0), &P2(memb.unitNDir,k + 1,0));
    	if (tglMbNucLocAre != 0) { 
			memb.unitEqAr[k] = memb.unitEqAr[k + 1];
		}
	}
	Copy1dArrayInt(&P2A(memb.unitCp.l,memb.unitCp.c,0,dimMbNuc), pArr, 
			dimMbNuc);
	DeleteElementArrayByIndex(memb.unit.l, &(memb.unit.c), n, dimMbNuc);
	for(k = 0; k < insNeiPar.c; k++) {
		CONT(!(insNeiPar.l[k] > nAct + nAbp + n 
				&& insNeiPar.l[k] < nAct + nAbp + nUnitMb));
		insNeiPar.l[k]--;
	}
	InsertElement1dArrayWoChk(insNeiPar.l, &insNeiPar.c,
			nAct + nAbp + nUnitMb + memb.unitCp.c);
	(memb.unitCp.c)++;
  }
}
// Subroutine for MovePartilces()
void MoveParticlesSubroutine1(ListInt *longChMv, double *longChDistMv, 
		int ind) {
  int n, k, oppInd; 

  for(n = longCh.c - 1; n >= 0; n--) {
	for(k = 0; k < 2; k++) {
		BREAK(P2A(longCh.l,n,k,2) == ind);
	}
	CONT(k == 2);
	// copy information from longCh.l to longChMv.l
	V2COPY(&P2A(longChMv->l,longChMv->c,0,2), &P2A(longCh.l,n,0,2));
	longChDistMv[longChMv->c] = longChDist[n];
	(longChMv->c)++;
	// delete them in longCh.l
	oppInd = P2A(longCh.l,n,1 - k,2);
	if (oppInd < nAct) { CONT(iAct[oppInd] > -1); }
	else if (oppInd >= nAct && oppInd < nAct + nAbp) 
	{ CONT(iAbp[oppInd - nAct] > -1); }
	else { CONT(iMb[oppInd - nAct - nAbp] > -1); }
	DeleteElementArrayByIndex(longCh.l, &longCh.c, n, 2);
	for(k = n; k < longCh.c; k++) {
		longChDist[k] = longChDist[k + 1];
	}
  }
}

// mode = 0: actin, 1: ABP, 2: Membrane
void MoveParticlesSubroutine2(int locInd, int mode) {
  int m, n, k;
  
  if (mode == 0) {
	for(n = locInd; n < nActMe - 1; n++) {
	    act.id[n] = act.id[n + 1];
	    V3COPY(&P2(act.r,n,0), &P2(act.r,n + 1,0));
	    V3COPY(&P2(act.f,n,0), &P2(act.f,n + 1,0));
	    Copy1dArrayInt(&P2A(act.ch,n,0,nChAc), 
				&P2A(act.ch,n + 1,0,nChAc), nChAc);
	    V3COPY(&P2(act.rPrev,n,0), &P2(act.rPrev,n + 1,0));
	    act.fix[n] = act.fix[n + 1];
		if (bndReb.gTgl != 0 && bndReb.gTglSpr == 1) {
	    	V3COPY(&P2(act.bndRFix,n,0), &P2(act.bndRFix,n + 1,0));
		}
		if (gTglBead != 0 && beadBind.gTgl != 0) { 
			act.bdFix[n] = act.bdFix[n + 1]; 
	    	V3COPY(&P2(act.bdRFix,n,0), &P2(act.bdRFix,n + 1,0));
		}
		if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) 
		{ act.nFA[n] = act.nFA[n + 1]; }
		if (gTglActCapAll != 0) { act.cap[n] = act.cap[n + 1]; }
		if (actAge.gTgl != 0) { act.age[n] = act.age[n + 1]; }
		act.len[n] = act.len[n + 1];
		if (gTglMb != 0) { 
			act.mbReb[n] = act.mbReb[n + 1]; 
			if (mbSld.act.gTgl != 0) {
				mbSld.act.reb[n] = mbSld.act.reb[n + 1];
				Copy1dArrayInt(&P2A(mbSld.act.l,n,0,dimMbNuc), 
						&P2A(mbSld.act.l,n + 1,0,dimMbNuc), dimMbNuc);
			}
		}
	    act.iF[n] = act.iF[n + 1];
  		if (confVmdInfo.tgl != 0) {
	        recAct.len[n] = recAct.len[n + 1];
	        recAct.sprF[n] = recAct.sprF[n + 1];
	        recAct.bendF[n] = recAct.bendF[n + 1];
	        V4COPY(&P2A(recAct.allF,n,0,NDIM + 1), 
					&P2A(recAct.allF,n + 1,0,NDIM + 1));
	        recAct.cnt[n] = recAct.cnt[n + 1];
	    }
	    iAct[act.id[n]] = n;
	}
	nActMe--;
  }
  else if (mode == 1) {
	for(m = locInd; m < nAbpMe - 1; m++) {
		abp.id[m] = abp.id[m + 1];
		V3COPY(&P2(abp.r,m,0), &P2(abp.r,m + 1,0));
		V3COPY(&P2(abp.f,m,0), &P2(abp.f,m + 1,0));
		V4COPY(&P2A(recInstSprFabp,m,0,4), &P2A(recInstSprFabp,m + 1,0,4));
		Copy1dArrayInt(&P2A(abp.ch,m,0,nChAb), 
				&P2A(abp.ch,m + 1,0,nChAb), nChAb);
		if (abpAge.gTgl != 0) { abp.age[m] = abp.age[m + 1]; }
		if (motSA.gTgl != 0) {
			abp.mId[m] = abp.mId[m + 1];
		}
		if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
			mbSld.abp.reb[m] = mbSld.abp.reb[m + 1];
			Copy1dArrayInt(&P2A(mbSld.abp.l,m,0,dimMbNuc), 
					&P2A(mbSld.abp.l,m + 1,0,dimMbNuc), dimMbNuc);
		}
		if (recAbpDyn.tgl != 0) {
			abpDyn.cntU[m] = abpDyn.cntU[m + 1];
			abpDyn.cntB[m] = abpDyn.cntB[m + 1];
			abpDyn.cntW[m] = abpDyn.cntW[m + 1];
		}
  		if (confVmdInfo.tgl != 0) {
			Copy1dArrayDouble(&P2A(recAbp.len,m,0,recAbp.nL), 
					&P2A(recAbp.len,m + 1,0,recAbp.nL), recAbp.nL);
			recAbp.sprF[m] = recAbp.sprF[m + 1];
			recAbp.bendF[m] = recAbp.bendF[m + 1];
	        V4COPY(&P2A(recAbp.allF,m,0,NDIM + 1), 
					&P2A(recAbp.allF,m + 1,0,NDIM + 1));
			recAbp.cnt[m] = recAbp.cnt[m + 1];
		}
		if (recLongF.tgl != 0) {
			V4COPY(&P2A(recLongSprFabp,m,0,4), &P2A(recLongSprFabp,m + 1,0,4));
		}
		if (tglRecAbpTurn != 0) {
			Copy1dArrayDouble(&P2A(abpTurn,m,0,7), &P2A(abpTurn,m + 1,0,7), 7);
		}
		iAbp[abp.id[m]] = m;
	}
	nAbpMe--;
  }
  else {
	for(n = locInd; n < nMbMe - 1; n++) {
	    memb.id[n] = memb.id[n + 1];
	    V3COPY(&P2(memb.r,n,0), &P2(memb.r,n + 1,0));
	    V3COPY(&P2(memb.f,n,0), &P2(memb.f,n + 1,0));
	    V3COPY(&P2(memb.nDir,n,0), &P2(memb.nDir,n + 1,0));
		Copy1dArrayInt(&P2A(memb.ch,n,0,nChMb), 
				&P2A(memb.ch,n + 1,0,nChMb), nChMb);
		V3COPY(&P2(memb.rPrev,n,0), &P2(memb.rPrev,n + 1,0));
		memb.idx[n] = memb.idx[n + 1];
		if (dimMbNuc == 3) {
			Copy1dArrayDouble(&P2A(memb.eqLen,n,0,nChMb - nMbAct), 
					&P2A(memb.eqLen,n + 1,0,nChMb - nMbAct), nChMb - nMbAct);
		}
		if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		    memb.fix[n] = memb.fix[n + 1];
	    	V3COPY(&P2(memb.rFix,n,0), &P2(memb.rFix,n + 1,0));
		}
		if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) { 
			memb.reb[n] = memb.reb[n + 1];
		}
		if (mbMat.gTgl != 0) {
			Copy1dArrayInt(&P2A(memb.nFA,n,0,nMbAct), 
					&P2A(memb.nFA,n + 1,0,nMbAct), nMbAct);
		}
		if (mbPro.gTgl != 0) {
			V2COPY(&P2A(memb.pro,n,0,2), &P2A(memb.pro,n + 1,0,2));
			if (gTglMbCont != 0) {
				memb.myo.conc[n] = memb.myo.conc[n + 1];
				memb.cyto.conc[n] = memb.cyto.conc[n + 1];
			}
		}
if (tglBias != 0 && typeBias == 2) {
	memb.bias[n] = memb.bias[n + 1];
}
	    iMb[memb.id[n]] = n;
	}
	nMbMe--;
  }
}
// Transfer particles which move out of the current subdomain to adjacent
// subdomains.
void MoveParticles(void) {
  int m, n, k, l, cnt, cntRecvMsg, cntElem, CS, *pArr, kind;
  int ind[3], ind2, actInd, abpInd, mbInd, indAdjCellMv; 
  int *chkL, oftMv[NDIM];
  int nActMvSub, nAbpMvSub, nMbMvSub, nChAbMv;
  int begin, end, rep, posi, tag = 0, nInfo;
  double r[NDIM], *longChDistMv, *unitMbInfoMv;
  ListInt longChMv, actCh, unitMbMv;


  nInfo = NDIM + tglMbNucLocAre;  

  nChAbMv = (motSA.gTgl != 0) ? nChAb : nChAb - 2;
  MALLOC(chkL,int,(mpiMethod == 0) ? cntAdjRank[0] : 2);
  if (mpiMethod == 0) {
	longChMv.l = allIntL;
	longChDistMv = allDblL;
  }
  if (gTglMb != 0) {
	MALLOC(unitMbMv.l,int,dimMbNuc*nMbC);
	MALLOC(unitMbInfoMv,double,nInfo*nMbC);
  }
  if (modeActCh != 0) {
	MALLOC(actCh.l,int,2 * nChAc - 1);
  }

  FOR_ACTCP(n) { iAct[act.id[n + nActMe]] = -1; }
  FOR_ABPCP(n) { iAbp[abp.id[n + nAbpMe]] = -1; }
  FOR_MBCP(n) { iMb[memb.id[n + nMbMe]] = -1; }

  nActCp = 0;
  nAbpCp = 0;
  nMbCp = 0;

  insNeiPar.c = 0; 

  rep = (mpiMethod == 0) ? 1 : NDIM;
  for(m = 0; m < rep; m++) {
	if (mpiMethod == 0) { 
		begin = 0;   
		end = NDIM; 
	}
	else { 
		begin = m;   
		end = m + 1; 
	}
	if (cntAdjRank[m] > 0) {
		// Move particles
		memset(cntMvPar, 0, sizeof(int) * cntAdjRank[m] * 3);
		for(n = 0; n < cntAdjRank[m]; n++) { 
			mvPar[n].c = 0; 
		}
		// Find particles to be transferred
		for(n = 0; n < nActMe + nAbpMe + nMbMe; n++) {
			kind = SetKind(n, nActMe, nAbpMe);
			if (kind == 0) { 
				CONT(ISACTM(n));
				V3COPY(r, &P2(act.r,n,0)); 
			}
			else if (kind == 1) { 
				CONT(ISABPIM(n - nActMe));
				V3COPY(r, &P2(abp.r,n - nActMe,0)); 
			}
			else {
				V3COPY(r, &P2(memb.r,n - nActMe - nAbpMe,0)); 
			}
			// if bulk rheology is used, PBC should be applied.
			if (rheoWay > 0 && bulkRheoType == 0) {
				ConvertRectDomainVector(r, 0);
			}
			CS = 1;
			V3SET_ALL(oftMv, 1);

			for(k = begin; k < end; k++) {
				// If there is no periodic boundary condition, particles 
				// crossing the outer boundary of the domain are ignored.
				if (pbc[k] == 0 && ((iCell[k] == 0 && r[k] < rGrid[k][0])
						|| (iCell[k] == nCell[k] - 1 
						&& r[k] >= rGrid[k][nGrid[k] - 1]))) { 
					CS = 0;   
					continue;   
				}
				// Treat particles crossing the outer boundary.
				else if (nCell[k] > 1) {
					// For middle cells
					if (iCell[k] > 0 && iCell[k] < nCell[k] - 1) {
						if (r[k] < P2A(bnd.r,0,k,NDIM)) { oftMv[k] -= 1; }
						else if (r[k] >= P2A(bnd.r,1,k,NDIM))
						{ oftMv[k] += 1; }
					}
					// For cells located at edges
					else if (iCell[k] == 0 && r[k] >= P2A(bnd.r,1,k,NDIM)) {
						if (r[k] >= rGrid[k][2]) { oftMv[k] -= 1; }
						else { oftMv[k] += 1; }
					}
					else if (iCell[k] == nCell[k] - 1 
							&& r[k] < P2A(bnd.r,0,k,NDIM)) {
						if (r[k] >= rGrid[k][nGrid[k] - 3]) { oftMv[k] -= 1; }
						else { oftMv[k] += 1; }
					}
				}
			}
			CONT(CS == 0 && mpiMethod == 1);
			CS = 0;
			// Normal method
			if (mpiMethod == 0) {
				V3IND_BACK_CONST_INT(indAdjCellMv, oftMv, NDIM);
			    ind2 = P2A(adjRank,indAdjCellMv,1,2);
				// "indAdjCellMv == 13" corrsponds to the current cell. (1,1,1)
			    if (indAdjCellMv != 13 && ind2 > -1) { CS = 1; }
			}
			// Plympton method
			else {
				ind2 = P2A(adjRank,m * 2 + oftMv[m] / 2,1,2);
				if (oftMv[m] != 1 && ind2 > -1) { CS = 1; }
			}
			// Add particles in list
			if (CS == 1) {
	            InsertElement1dArrayWoChk(mvPar[ind2].l, &mvPar[ind2].c, n);
				P2A(cntMvPar,ind2,kind,3)++;
			}
		}
	}
	// Send chosen particles
	// The information of the particles is packed into a message
	for(n = 0; n < cntAdjRank[m]; n++) {
		posi = 0;
		V3COPY(ind, &P2A(cntMvPar,n,0,3));
		if (mpiMethod == 0) { longChMv.c = 0; }
		if (gTglMb != 0) { unitMbMv.c = 0; }

		MPI_PACK_INT(&P2A(cntMvPar,n,0,3), 3, n);
		// actin
		for(k = 0; k < ind[0]; k++) {
			actInd = mvPar[n].l[k];
			MPI_PACK_INT(&act.id[actInd], 1, n);
			MPI_PACK_DBL(&P2(act.r,actInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(act.f,actInd,0), NDIM, n);
			if (modeActCh == 0) { 
				MPI_PACK_INT(&P2A(act.ch,actInd,0,nChAc), nChAc, n);
			}
			else {		
				V2COPY(&actCh.l[1], &P2A(act.ch,actInd,0,nChAc));
				actCh.c = 2;
				for(l = 2; l < nChAc; l++) {
					if (P2A(act.ch,actInd,l,nChAc) > -1) {
						V2SET(&actCh.l[actCh.c + 1], l, 
								P2A(act.ch,actInd,l,nChAc));
						actCh.c += 2;
					}
				}
				actCh.l[0] = actCh.c;
				MPI_PACK_INT(actCh.l, actCh.c + 1, n);
			}
			MPI_PACK_DBL(&P2(act.rPrev,actInd,0), NDIM, n);
			MPI_PACK_INT(&act.fix[actInd], 1, n);
			if (bndReb.gTgl != 0 && bndReb.gTglSpr == 1) {
				MPI_PACK_DBL(&P2(act.bndRFix,actInd,0), NDIM, n);
			}
			if (gTglBead != 0 && beadBind.gTgl != 0) { 
				MPI_PACK_INT(&act.bdFix[actInd], 1, n);
				MPI_PACK_DBL(&P2(act.bdRFix,actInd,0), NDIM, n);
			}
			if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) 
			{ MPI_PACK_INT(&act.nFA[actInd], 1, n); }
			if (gTglActCapAll != 0) { MPI_PACK_INT(&act.cap[actInd], 1, n); }
			if (actAge.gTgl != 0) { MPI_PACK_INT(&act.age[actInd], 1, n); }
			MPI_PACK_INT(&act.len[actInd], 1, n);
			if (gTglMb != 0) { 
				MPI_PACK_INT(&act.mbReb[actInd], 1, n); 
				if (mbSld.act.gTgl != 0) {
					MPI_PACK_INT(&mbSld.act.reb[actInd], 1, n);
					MPI_PACK_INT(&P2A(mbSld.act.l,actInd,0,dimMbNuc), 
							dimMbNuc, n);
				}
			}
			MPI_PACK_INT(&act.iF[actInd], 1, n);
			if (recTraj.tgl != 0) {
				CS = FindElementArray(recTraj.actMe.l, recTraj.actMe.c,	
						act.id[actInd], 0, 1);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					DeleteElement1dArray(recTraj.actMe.l, &recTraj.actMe.c, 
							act.id[actInd]);
				}
			}
			if (rheoWay > 0) {
				if (act.fix[actInd] > 10) {
					DeleteElement1dArray(meaStreParMe.l, 
							&meaStreParMe.c, act.id[actInd]);
					DeleteElement1dArray(appStraParMe.l, 
							&appStraParMe.c, act.id[actInd]);
			
				}
			}
  			if (confVmdInfo.tgl != 0) {
				MPI_PACK_DBL(&recAct.len[actInd], 1, n);
				MPI_PACK_DBL(&recAct.sprF[actInd], 1, n);
				MPI_PACK_DBL(&recAct.bendF[actInd], 1, n);
				MPI_PACK_DBL(&P2A(recAct.allF,actInd,0,NDIM + 1), NDIM + 1, n);
				MPI_PACK_INT(&recAct.cnt[actInd], 1, n);
			}
			if (tglActDyn != 0) {
				CS = FindElementArray(noActDyn.l, noActDyn.c, 
						act.id[actInd], 0, 2);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noActDyn.l,CS,0,2), 2, n);
					DeleteElementArrayByIndex(noActDyn.l, &noActDyn.c, CS, 2);
				}
			}
			if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 || motMoBind.tgl != 0 
					|| tglActDisBstSevAbp != 0) { 
				CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, 
						act.id[actInd], 1, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noAbpDyn.l,CS,0,3), 3, n);
		        }
			}
			if (gTglMb != 0 && (mbReb.gTgl != 0 || mbUnb.gTgl != 0)) {
				CS = FindElementArray(noMbDyn.l, noMbDyn.c, 
						act.id[actInd], 1, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noMbDyn.l,CS,0,3), 3, n);
		        }
			}
			if (mpiMethod == 0) {
				MoveParticlesSubroutine1(&longChMv, longChDistMv, 
						act.id[actInd]);
			}
			iAct[act.id[actInd]] = -1;
			act.id[actInd] = -1;
		}
		// ABP
		for(k = 0; k < ind[1]; k++) {
			abpInd = mvPar[n].l[k + ind[0]] - nActMe;
			pArr = &P2A(abp.ch,abpInd,0,nChAb);
			MPI_PACK_INT(&abp.id[abpInd], 1, n);
			MPI_PACK_DBL(&P2(abp.r,abpInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(abp.f,abpInd,0), NDIM, n);
			MPI_PACK_DBL(&P2A(recInstSprFabp,abpInd,0,4), 4, n);
			MPI_PACK_INT(pArr, nChAbMv, n);
			if (abpAge.gTgl != 0) { MPI_PACK_INT(&abp.age[abpInd], 1, n); }
			if (motSA.gTgl != 0) {
				MPI_PACK_INT(&abp.mId[abpInd], 1, n);
			}
			if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
				MPI_PACK_INT(&mbSld.abp.reb[abpInd], 1, n);
				MPI_PACK_INT(&P2A(mbSld.abp.l,abpInd,0,dimMbNuc), dimMbNuc, n);
			}
			if (recAbpDyn.tgl != 0) {
				MPI_PACK_INT(&abpDyn.cntU[abpInd], 1, n);
				MPI_PACK_INT(&abpDyn.cntB[abpInd], 1, n);
				MPI_PACK_INT(&abpDyn.cntW[abpInd], 1, n);
			}
			if (recTraj.tgl2 != 0) {
				CS = FindElementArray(recTraj.abpMe.l, recTraj.abpMe.c,	
						abp.id[abpInd], 0, 1);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					DeleteElement1dArray(recTraj.abpMe.l, &recTraj.abpMe.c, 
							abp.id[abpInd]);
				}
			}
  			if (confVmdInfo.tgl != 0) {
				MPI_PACK_DBL(&P2A(recAbp.len,abpInd,0,recAbp.nL), recAbp.nL, n);
				MPI_PACK_DBL(&recAbp.sprF[abpInd], 1, n);
				MPI_PACK_DBL(&recAbp.bendF[abpInd], 1, n);
				MPI_PACK_DBL(&P2A(recAbp.allF,abpInd,0,NDIM + 1), NDIM + 1, n);
				MPI_PACK_INT(&recAbp.cnt[abpInd], 1, n);
			}
			if (recLongF.tgl != 0) {
				MPI_PACK_DBL(&P2A(recLongSprFabp,abpInd,0,4), 4, n);
			}
			if (mpiMethod == 0) {
				MoveParticlesSubroutine1(&longChMv, longChDistMv,
						abp.id[abpInd] + nAct);
			}
			if (pArr[0] > -1 && pArr[1] > -1) {
				if (tglNeiAbpDC != 0) {
					DeleteElementInNeighborList(abp.id[abpInd], 1);
				}
			}
			else {
				if (tglNeiAbpSC != 0) {
					DeleteElementInNeighborList(abp.id[abpInd], 1);
				}
				if (K_ABP(abpInd) == 2) { 
					if (pArr[0] > -1) { nMotInaMe--; }
					else { nMotMme--; }
				}
				else {
					if (pArr[0] > -1) { nAcpInaMe--; }
					else { nAcpMme--; }
				}
			}
			if ((tglAcpAcInaDyn != 0 && K_ABP(abpInd) != 2) 
					|| (tglMotAcInaDyn != 0 && K_ABP(abpInd) == 2)
					|| tglActDisBstSevAbp != 0) {
				CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, 
						abp.id[abpInd], 0, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noAbpDyn.l,CS,0,3), 3, n);
		        }
			}
			if (tglRecAbpTurn != 0) {
				MPI_PACK_DBL(&P2A(abpTurn,abpInd,0,7), 7, n);
			}
			iAbp[abp.id[abpInd]] = -1;
			abp.id[abpInd] = -1;
		}
		// Membrane
		for(k = 0; k < ind[2]; k++) {
			mbInd = mvPar[n].l[k + ind[0] + ind[1]] - nActMe - nAbpMe;
			MPI_PACK_INT(&memb.id[mbInd], 1, n);
			MPI_PACK_DBL(&P2(memb.r,mbInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(memb.f,mbInd,0), NDIM, n);
			MPI_PACK_DBL(&P2(memb.nDir,mbInd,0), NDIM, n);
			MPI_PACK_INT(&P2A(memb.ch,mbInd,0,nChMb), nChMb, n);
			MPI_PACK_DBL(&P2(memb.rPrev,mbInd,0), NDIM, n);
			MPI_PACK_INT(&memb.idx[mbInd], 1, n);
			if (dimMbNuc == 3) {
				MPI_PACK_DBL(&P2A(memb.eqLen,mbInd,0,nChMb - nMbAct), 
						nChMb - nMbAct, n);
			}
			if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
				MPI_PACK_INT(&memb.fix[mbInd], 1, n);
				MPI_PACK_DBL(&P2(memb.rFix,mbInd,0), NDIM, n);
			}
			if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
				MPI_PACK_INT(&memb.reb[mbInd], 1, n);
			}
			if (mbMat.gTgl != 0) {
				MPI_PACK_INT(&P2A(memb.nFA,mbInd,0,nMbAct), nMbAct, n);
			}
			if (mbPro.gTgl != 0) {
				MPI_PACK_INT(&P2A(memb.pro,mbInd,0,2), 2, n);
				if (gTglMbCont != 0) {
					MPI_PACK_DBL(&memb.myo.conc[mbInd], 1, n);
					MPI_PACK_DBL(&memb.cyto.conc[mbInd], 1, n);
				}
			}
			if (tglBias != 0 && typeBias == 2) {
				MPI_PACK_INT(&memb.bias[mbInd], 1, n);
			}
			if (tglMbNucDyn != 0) {
				CS = FindElementArray(noMbDyn.l, noMbDyn.c, 
						memb.id[mbInd], 0, 3);
				MPI_PACK_INT(&CS, 1, n);
				if (CS > -1) {
					MPI_PACK_INT(&P2A(noMbDyn.l,CS,0,3), 3, n);
		        }
			}
			if (mpiMethod == 0) {
				MoveParticlesSubroutine1(&longChMv, longChDistMv,
						memb.id[mbInd] + nAct + nAbp);
			}
			MoveParticlesSubroutine3(&unitMbMv, unitMbInfoMv, 
					memb.id[mbInd], nInfo);

			iMb[memb.id[mbInd]] = -1;
			memb.id[mbInd] = -1;
			
		}
		// Additional information to be packed
		if (mpiMethod == 0) {
			MPI_PACK_INT(&longChMv.c, 1, n);
			MPI_PACK_INT(longChMv.l, longChMv.c * 2, n);
			MPI_PACK_DBL(longChDistMv, longChMv.c, n);
			// If it's the subdomain where actin in longChExtMsg.l can exist
			// at next step, the message is transferred. 
			longChMv.c = 0;
			for(k = 0; k < longChExtMsg.c; k++) {
				if (iRank[n] == P2A(longChExtMsg.l,k,2,3)) { 
					V2COPY(&P2A(longChMv.l,longChMv.c,0,2),
							&P2A(longChExtMsg.l,k,0,3));
					(longChMv.c)++;
				}
			}
			MPI_PACK_INT(&longChMv.c, 1, n);
			MPI_PACK_INT(longChMv.l, longChMv.c * 2, n);
		}
		else {
			MPI_PACK_INT(&longChIntMsg.c, 1, n);
			MPI_PACK_INT(longChIntMsg.l, longChIntMsg.c, n);
		}
	 	if (gTglMb != 0) {
			MPI_PACK_INT(&unitMbMv.c, 1, n);
			MPI_PACK_INT(unitMbMv.l, unitMbMv.c * dimMbNuc, n);
			MPI_PACK_DBL(unitMbInfoMv, unitMbMv.c * nInfo, n);
		}	
		// Send and receive the messages to adjacent subdomains
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, iRank[(mpiMethod == 0) 
				? n : 2 * m + n], tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, 
				iRank[(mpiMethod == 0) ? n : 2 * m + n], 
				tag, MPI_COMM_WORLD, &rReq[n]);
	}
	if (cntAdjRank[m] > 0) {
		// Delete moved particles
		// actin
		for(n = nActMe - 1; n >= 0; n--) {  
			CONT(!(act.id[n] < 0));
			MoveParticlesSubroutine2(n, 0);
		}
		// ABP
		for(n = nAbpMe - 1; n >= 0; n--) {
			CONT(!(abp.id[n] < 0));
			MoveParticlesSubroutine2(n, 1);
		}
		// Membrane
		for(n = nMbMe - 1; n >= 0; n--) {  
			CONT(!(memb.id[n] < 0));
			MoveParticlesSubroutine2(n, 2);
		}

		// Receive information of particles transferred from adjacent domains
		cnt = 0;
		cntRecvMsg = 0;
		memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(chkL, 0, sizeof(int) * cntAdjRank[m]);
		while(cntRecvMsg < cntAdjRank[m]) {
			if (mpiTestSendFlag[cnt] == 0) {
				MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
			}
			if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0 
					&& chkL[cnt] == 0) {
				MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
				if (mpiTestRecvFlag[cnt] != 0) { 
					posi = 0;
			  		MPI_UNPACK_INT(&nActMvSub, 1, cnt);
			  		MPI_UNPACK_INT(&nAbpMvSub, 1, cnt);
			  		MPI_UNPACK_INT(&nMbMvSub, 1, cnt);
					// actin
					for(n = 0; n < nActMvSub; n++) {
						actInd = nActMe + n;
				  		pArr = &P2A(act.ch,actInd,0,nChAc);
				  		MPI_UNPACK_INT(&act.id[actInd], 1, cnt);
				  		MPI_UNPACK_DBL(&P2(act.r,actInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(act.f,actInd,0), NDIM, cnt);
						if (modeActCh == 0) {
				  			MPI_UNPACK_INT(pArr, nChAc, cnt);
						}
						else {
							memset(&pArr[2], -1, sizeof(int) * (nChAc - 2));
				  			MPI_UNPACK_INT(&actCh.c, 1, cnt);
				  			MPI_UNPACK_INT(actCh.l, actCh.c, cnt);
							V2COPY(pArr, actCh.l);
							for(k = 2; k < actCh.c; k += 2) {
								pArr[actCh.l[k]] =  actCh.l[k + 1];
							}
						}
						iAct[act.id[actInd]] = actInd;
				  		MPI_UNPACK_DBL(&P2(act.rPrev,actInd,0), NDIM, cnt);
					  	MPI_UNPACK_INT(&act.fix[actInd], 1, cnt);
						if (bndReb.gTgl != 0 && bndReb.gTglSpr == 1) {
				  			MPI_UNPACK_DBL(&P2(act.bndRFix,actInd,0), 
									NDIM, cnt);
						}
						if (gTglBead != 0 && beadBind.gTgl != 0) { 
					  		MPI_UNPACK_INT(&act.bdFix[actInd], 1, cnt);
				  			MPI_UNPACK_DBL(&P2(act.bdRFix,actInd,0), NDIM, cnt);
						}
						if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) 
						{ MPI_UNPACK_INT(&act.nFA[actInd], 1, cnt); }
						if (gTglActCapAll != 0) {
						  	MPI_UNPACK_INT(&act.cap[actInd], 1, cnt);
						}
						if (actAge.gTgl != 0) {
						  	MPI_UNPACK_INT(&act.age[actInd], 1, cnt);
						}
						MPI_UNPACK_INT(&act.len[actInd], 1, cnt);
						if (gTglMb != 0) { 
					  		MPI_UNPACK_INT(&act.mbReb[actInd], 1, cnt);
							if (mbSld.act.gTgl != 0) {
								MPI_UNPACK_INT(&mbSld.act.reb[actInd], 1, cnt);
								MPI_UNPACK_INT(&P2A(mbSld.act.l,actInd,0,
										dimMbNuc), dimMbNuc, cnt);
							}
						}
					  	MPI_UNPACK_INT(&act.iF[actInd], 1, cnt);
						if (recTraj.tgl != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
								InsertElement1dArrayWoChk(recTraj.actMe.l, 
										&recTraj.actMe.c, act.id[actInd]);
							}
						}
						if (rheoWay > 0) {
							if (act.fix[actInd] > 10) {
								InsertElement1dArrayWoChk(meaStreParMe.l,
										&meaStreParMe.c, act.id[actInd]);
								InsertElement1dArrayWoChk(appStraParMe.l,
										&appStraParMe.c, act.id[actInd]);
							}
						}
			  			if (confVmdInfo.tgl != 0) {
					  		MPI_UNPACK_DBL(&recAct.len[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAct.sprF[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAct.bendF[actInd], 1, cnt);
					  		MPI_UNPACK_DBL(&P2A(recAct.allF,actInd,0,NDIM + 1), 
									NDIM + 1, cnt);
					  		MPI_UNPACK_INT(&recAct.cnt[actInd], 1, cnt);
						}
						if (tglActDyn != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noActDyn.l, 
										noActDyn.c,0,2), 2, cnt);
								(noActDyn.c)++;
							}
						}
						if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 
								|| motMoBind.tgl != 0 
								|| tglActDisBstSevAbp != 0) { 
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noAbpDyn.l, 
									noAbpDyn.c,0,3), 3, cnt);
								(noAbpDyn.c)++;
							}
						}
						if (gTglMb != 0 && (mbReb.gTgl != 0 
								|| mbUnb.gTgl != 0)) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noMbDyn.l, 
									noMbDyn.c,0,3), 3, cnt);
								(noMbDyn.c)++;
							}
						}
					}
					// ABP
					for(n = 0; n < nAbpMvSub; n++) {
						abpInd = nAbpMe + n;
						pArr = &P2A(abp.ch,abpInd,0,nChAb);
				  		MPI_UNPACK_INT(&abp.id[abpInd], 1, cnt);
				  		MPI_UNPACK_DBL(&P2(abp.r,abpInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(abp.f,abpInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2A(recInstSprFabp,abpInd,0,4), 4, cnt);
				  		MPI_UNPACK_INT(pArr, nChAbMv, cnt);
						if (abpAge.gTgl != 0) {
						  	MPI_UNPACK_INT(&abp.age[abpInd], 1, cnt);
						}
						if (motSA.gTgl != 0) {
				  			MPI_UNPACK_INT(&abp.mId[abpInd], 1, cnt);
						}
						else {
							V2SET_ALL(&pArr[3], -1);
						}
						if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
							MPI_UNPACK_INT(&mbSld.abp.reb[abpInd], 1, cnt);
							MPI_UNPACK_INT(&P2A(mbSld.abp.l,abpInd,0,dimMbNuc), 
									dimMbNuc, cnt);
						}
						if (recAbpDyn.tgl != 0) {
							MPI_UNPACK_INT(&abpDyn.cntU[abpInd], 1, cnt);
							MPI_UNPACK_INT(&abpDyn.cntB[abpInd], 1, cnt);
							MPI_UNPACK_INT(&abpDyn.cntW[abpInd], 1, cnt);	
						}
						if (recTraj.tgl2 != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
								InsertElement1dArrayWoChk(recTraj.abpMe.l, 
										&recTraj.abpMe.c, abp.id[abpInd]);
							}
						}
			  			if (confVmdInfo.tgl != 0) {
					  		MPI_UNPACK_DBL(&P2A(recAbp.len,abpInd,0,recAbp.nL), 
									recAbp.nL, cnt);
					  		MPI_UNPACK_DBL(&recAbp.sprF[abpInd], 1, cnt);
					  		MPI_UNPACK_DBL(&recAbp.bendF[abpInd], 1, cnt);
					  		MPI_UNPACK_DBL(&P2A(recAbp.allF,abpInd,0,NDIM + 1), 
									NDIM + 1, cnt);
					  		MPI_UNPACK_INT(&recAbp.cnt[abpInd], 1, cnt);
						}
						if (recLongF.tgl != 0) {
				  			MPI_UNPACK_DBL(&P2A(recLongSprFabp,
									abpInd,0,4), 4, cnt);
						}
						iAbp[abp.id[abpInd]] = abpInd;
		
						if (pArr[0] > -1 && pArr[1] > -1) {
							if (tglNeiAbpDC != 0) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
										&insNeiPar.c, nAct + abp.id[abpInd]);
							}
						}
						else {
							if (tglNeiAbpSC != 0) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
									&insNeiPar.c, nAct + abp.id[abpInd]);
							}
							if (pArr[2] == 2) { 
								if (pArr[0] > -1) { nMotInaMe++; }
								else { nMotMme++; }
							}
							else {
								if (pArr[0] > -1) { nAcpInaMe++; }
								else { nAcpMme++; }
							}
						}
						if ((tglAcpAcInaDyn != 0 && pArr[2] != 2) 
								|| (tglMotAcInaDyn != 0 && pArr[2] == 2)
								|| tglActDisBstSevAbp != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noAbpDyn.l, 
									noAbpDyn.c,0,3), 3, cnt);
								(noAbpDyn.c)++;
							}
						}
						if (tglRecAbpTurn != 0) {
							MPI_UNPACK_DBL(&P2A(abpTurn,abpInd,0,7), 7, cnt);
						}
					}
					// Membrane
					for(n = 0; n < nMbMvSub; n++) {
						mbInd = nMbMe + n;
				  		MPI_UNPACK_INT(&memb.id[mbInd], 1, cnt);
				  		MPI_UNPACK_DBL(&P2(memb.r,mbInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(memb.f,mbInd,0), NDIM, cnt);
				  		MPI_UNPACK_DBL(&P2(memb.nDir,mbInd,0), NDIM, cnt);
				  		MPI_UNPACK_INT(&P2A(memb.ch,mbInd,0,nChMb), nChMb, cnt);
				  		MPI_UNPACK_DBL(&P2(memb.rPrev,mbInd,0), NDIM, cnt);
				  		MPI_UNPACK_INT(&memb.idx[mbInd], 1, cnt);
						if (dimMbNuc == 3) {
					  		MPI_UNPACK_DBL(&P2A(memb.eqLen,mbInd,0,nChMb 
									- nMbAct), nChMb - nMbAct, cnt);
						}
						if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
					  		MPI_UNPACK_INT(&memb.fix[mbInd], 1, cnt);
					  		MPI_UNPACK_DBL(&P2(memb.rFix,mbInd,0), NDIM, cnt);
						}
						if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
				  			MPI_UNPACK_INT(&memb.reb[mbInd], 1, cnt);
						}
						if (mbMat.gTgl != 0) {
				  			MPI_UNPACK_INT(&P2A(memb.nFA,mbInd,0,nMbAct), 
									nMbAct, cnt);
						}
						if (mbPro.gTgl != 0) {
				  			MPI_UNPACK_INT(&P2A(memb.pro,mbInd,0,2), 2, cnt);
							if (gTglMbCont != 0) {
					  			MPI_UNPACK_DBL(&memb.myo.conc[mbInd], 1, cnt);
						  		MPI_UNPACK_DBL(&memb.cyto.conc[mbInd], 1, cnt);
							}
						}

						if (tglBias != 0 && typeBias == 2) {
							MPI_UNPACK_INT(&memb.bias[mbInd], 1, cnt);
						}

						if (tglMbNucDyn != 0) {
					  		MPI_UNPACK_INT(&CS, 1, cnt);
							if (CS > -1) {
					  			MPI_UNPACK_INT(&P2A(noMbDyn.l, 
									noMbDyn.c,0,3), 3, cnt);
								(noMbDyn.c)++;
							}
						}
						iMb[memb.id[mbInd]] = mbInd;
					}
					// Additional information to be received
					if (mpiMethod == 0) {
						MPI_UNPACK_INT(&cntElem, 1, cnt);
						MPI_UNPACK_INT(longChMv.l, cntElem * 2, cnt);
						MPI_UNPACK_DBL(longChDistMv, cntElem, cnt);
						for(n = 0; n < cntElem; n++) {
							pArr = &P2A(longChMv.l,n,0,2);
							InsertLongChain(pArr[0], pArr[1], longChDistMv[n]);
						}
						// The transferred longChExtMsg is stored in 
						// longChIntMsg.l.
						MPI_UNPACK_INT(&cntElem, 1, cnt);
						MPI_UNPACK_INT(&P2A(longChIntMsg.l,longChIntMsg.c
								,0,2) , cntElem * 2, cnt);
						longChIntMsg.c += cntElem;					
					}
					else {
						MPI_UNPACK_INT(&cntElem, 1, cnt);
						MPI_UNPACK_INT(&longChIntMsg.l[longChIntMsg.c], 
								cntElem, cnt);
						longChIntMsg.c += cntElem;
					}
	 				if (gTglMb != 0) {
						MPI_UNPACK_INT(&cntElem, 1, cnt);
						MPI_UNPACK_INT(unitMbMv.l, cntElem * dimMbNuc, cnt);
						MPI_UNPACK_DBL(unitMbInfoMv, cntElem * nInfo, cnt);
						for(n = 0; n < cntElem; n++) {
							pArr = &P2A(unitMbMv.l,n,0,dimMbNuc);
							if (dimMbNuc == 2) {
								CS = Find2ElementArray(memb.unit.l, memb.unit.c,
									pArr[0], pArr[1], 0, 2);
							}
							else {
								CS = Find3ElementArray(memb.unit.l, memb.unit.c,
									pArr[0], pArr[1], pArr[2], 0, 3);
							}
							CONT(CS > -1);
							if (dimMbNuc == 2) {
								CS = Find2ElementArray(memb.unitCp.l,
										memb.unitCp.c,pArr[0], pArr[1], 0, 2);
							}
							else {
								CS = Find3ElementArray(memb.unitCp.l, 
										memb.unitCp.c, pArr[0], pArr[1], 
										pArr[2], 0, 3);
							}
							if (CS > -1) {
								DeleteElementInNeighborList(CS, 3);
								for(k = 0; k < insNeiPar.c; k++) {
									ind2 = nAct + nAbp + nUnitMb;
									CONT(!(insNeiPar.l[k] > ind2 + CS
										&& insNeiPar.l[k] < ind2 + nUnitMb));
									insNeiPar.l[k]--;
								}
								DeleteElementArrayByIndex(memb.unitCp.l, 
										&memb.unitCp.c, CS, dimMbNuc);
							}
							Copy1dArrayInt(&P2A(memb.unit.l,memb.unit.c,0,
									dimMbNuc), pArr, dimMbNuc);
							V3COPY(&P2(memb.unitNDir,memb.unit.c,0), 
									&P2A(unitMbInfoMv,n,0,nInfo));
    						if (tglMbNucLocAre != 0) { 
								memb.unitEqAr[memb.unit.c] 
										= P2A(unitMbInfoMv,n,3,nInfo);
							}
							InsertElement1dArrayWoChk(insNeiPar.l, &insNeiPar.c,
									nAct + nAbp + memb.unit.c);
							(memb.unit.c)++;
						}
					}	
					nActMe += nActMvSub;
					nAbpMe += nAbpMvSub;
					nMbMe += nMbMvSub;
					cntRecvMsg++;
					chkL[cnt] = 1;
				}

			}
			cnt++;
			if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
  } 
  free(chkL);
  if (gTglMb != 0) {
	free(unitMbMv.l);
	free(unitMbInfoMv);
  }
  if (modeActCh != 0) {
	free(actCh.l);
  }
}

/*--------------------- Move particles between subdomains --------------------*/

/*--------------------- Copy particles between subdomains --------------------*/

// Subroutine for CopyPartilcesSubroutine()
void CopyParticlesSubSubroutine(int id, int *idx, int ind) {
  int k, CS, end;

  // If actin, two actins are transferred
  if (id < nActMe + nActCp) {
	for(k = 0; k < 2; k++) {
		CONT(!(idx[k] > -1));
		CS = InsertElement1dArrayWChk(cpPar[ind].l, &cpPar[ind].c, idx[k]);
		if (CS == -1) { P2A(cntCpPar,ind,0,3)++; }
	}
  }
  // If ABP, only one actin is transferred
  else if (id >= nActMe + nActCp
		&& id < nActMe + nActCp + nAbpMe + nAbpCp) {
	for(k = 0; k < 2; k++) {
		CONT(!(idx[k] > -1));
		CS = FindElementArray(cpPar[ind].l, 
				P2A(cntCpPar,ind,0,3), idx[k], 0, 1);
		CONT(!(CS == -1));
		// Squeeze it in the list
		InsertElementArrayByIndex(cpPar[ind].l,
				&cpPar[ind].c, &idx[k], P2A(cntCpPar,ind,0,3), 1);
		P2A(cntCpPar,ind,0,3)++;
	}
  }
  else {
	end = (dimMbNuc == 2) ? 2 : 5;
	for(k = 0; k < end; k++) {
		CONT(!(idx[k] > -1));
		CS = FindElementArray(&cpPar[ind].l[P2A(cntCpPar,ind,0,3) 
				+ P2A(cntCpPar,ind,1,3)], P2A(cntCpPar,ind,2,3), idx[k], 0, 1);
		CONT(!(CS == -1));
		InsertElement1dArrayWoChk(cpPar[ind].l, &cpPar[ind].c, idx[k]);
		P2A(cntCpPar,ind,2,3)++;
	}
  }
}

// Subroutine for CopyParticles() with normal method
// mode: 0 - normal, 1 - long
void CopyParticlesNormalSubroutine(int id, int *oftCp, int *idx, int mode) { 
  int m, n, k, l, ind, kind, cntCpuList, *cpuList;
  int *cpuListTemp, CS, indAdjCellCp;
  int oftInd[][2] = {{0}, {1}, {2}, {0, 2}}, cntOftInd[4] = {1, 1, 1, 2};
  int offCpuList[][7] = OFFSET_CPU_LIST, offCpuLen[] = OFFSET_CPU_LEN;

  if (mode == 0 && oftCp[0] != 3 && oftCp[1] != 3 && oftCp[2] != 3) {
	indAdjCellCp = ((oftCp[0] * NDIM) + oftCp[1]) * NDIM + oftCp[2];
	if (indAdjCellCp != 13) {
		cntCpuList = offCpuLen[indAdjCellCp];
		cpuList = offCpuList[indAdjCellCp];
	}
	else { return; }
  }
  else {
	MALLOC(cpuList,int,27);
	MALLOC(cpuListTemp,int,27);
	memset(cpuListTemp, -1, sizeof(int) * 27);
	for(m = 0; m < cntOftInd[oftCp[0]]; m++) {
		for(n = 0; n < cntOftInd[oftCp[1]]; n++) {
			for(k = 0; k < cntOftInd[oftCp[2]]; k++) {
	  			indAdjCellCp = ((oftInd[oftCp[0]][m] * NDIM) 
						+ oftInd[oftCp[1]][n]) * NDIM + oftInd[oftCp[2]][k];
				if (indAdjCellCp != 13) {
					for(l = 0; l < offCpuLen[indAdjCellCp]; l++) {
						cpuListTemp[offCpuList[indAdjCellCp][l]] = 1;
					}
				}
			}
		}
	}
	cntCpuList = 0;
	for(n = 0; n < 27; n++) {
		CONT(!(cpuListTemp[n] == 1));
		cpuList[cntCpuList] = n;
		cntCpuList++;
	}
	free(cpuListTemp);
  }
  for(n = 0; n < cntCpuList; n++) {
	ind = P2A(adjRank,cpuList[n],1,2);
	CONT(!(ind > -1));
	CS = InsertElement1dArrayWChk(cpPar[ind].l, &cpPar[ind].c, id);
	if (CS == -1) {
  		kind = SetKind(id, nActMe + nActCp, nAbpMe + nAbpCp);
		P2A(cntCpPar,ind,kind,3)++; 
	}
	CopyParticlesSubSubroutine(id, idx, ind);
  }
  if (!(mode == 0 && oftCp[0] != 3 && oftCp[1] != 3 && oftCp[2] != 3)) {
	free(cpuList); 
  }
}

// Subroutine for CopyParticles() with Plympton method
int CopyParticlesPlymptonSubroutine2(int id) { 
  int k, ind, kind, CS, nBd; 

  kind = SetKind(id, nActMe + nActCp, nAbpMe + nAbpCp);
  CS = -1;
  if (kind == 0) {
	for(k = 0; k < nChAc; k++) {
		ind = P2A(act.ch,id,k,nChAc);
		CONT(ind < 0);
		if (k < 2) { CONT(iAct[ind] > -1); }
		else { CONT(iAbp[ind] > -1); }
		CS = 1;   
		break; 
	}
  }
  else if (kind == 1) {
	id -= nActMe + nActCp;
	for(k = 0; k < nChAb; k++) {
		if (k == 2) {
			if (ISMTF(K_ABP(id))) { continue; }
			else { break; }
		}
		ind = P2A(abp.ch,id,k,nChAb);
		CONT(ind < 0);
		if (k < 2) { CONT(iAct[ind] > -1); }
		else { CONT(iAbp[ind] > -1); }
		CS = 1;   
		break; 
	}
  }
  else {
	id -= nActMe + nActCp + nAbpMe + nAbpCp;
	nBd = (dimMbNuc == 2) ? 2 : 6;	
	for(k = 0; k < nChMb; k++) {
		ind = P2A(memb.ch,id,k,nChMb);
		CONT(ind < 0);
		if (k < nBd) { CONT(iMb[ind] > -1); }
		else { CONT(iAct[ind] > -1); }
		CS = 1;
		break;
	}	
  }
  return CS;
}

// Subroutine for CopyParticles() with Plympton method
void CopyParticlesPlymptonSubroutine(int id, int dire, int *idx, int m) { 
  int k, ind, begin, end, CS, kind;
  double *r;

  kind = SetKind(id, nActMe + nActCp, nAbpMe + nAbpCp);
  if (dire < 3) {
	begin = dire / 2;
	end = dire / 2 + 1;
  }
  else {
	begin = 0;
	end = 2;
  }

  for(k = begin; k < end; k++) {
	ind = P2A(adjRank,m * 2 + k,1,2);
	CONT(!(ind > -1));
	CS = InsertElement1dArrayWChk(cpPar[ind].l, &cpPar[ind].c, id);
	if (CS == -1) {
		P2A(cntCpPar,ind,kind,3)++; 
	}
	CopyParticlesSubSubroutine(id, idx, ind);
  }
}

void CopyParticlesSubroutine(int begin, int end, int *oftCp,
	double *r, double dist) {
  int k;

  V3SET_ALL(oftCp, 1);
  for(k = begin; k < end; k++) {
	if (!(pbc[k] == 0 && iCell[k] == 0) && r[k] >= P2A(bnd.r,0,k,NDIM)
			&& r[k] < P2A(bnd.r,0,k,NDIM) + dist) {
		oftCp[k]--;
	}
	if (!(pbc[k] == 0 && iCell[k] == nCell[k] - 1)
			&& r[k] >= P2A(bnd.r,1,k,NDIM) - dist
			&& r[k] < P2A(bnd.r,1,k,NDIM)) {
		oftCp[k] = (oftCp[k] == 1) ? 2 : 3;
	}
  }
}

// Send information of particles located near boundaries to adjacent boundaries
// for synchronization.
void CopyParticles(void) {
  int m, n, k, l, cnt, cntRecvMsg, cntElem, CS, oftCp[NDIM], cntMb, kind;
  int ind[3], ind2, ind3, idx[5], *chkL, *pArr;
  int actInd, abpInd, mbInd, absInd, oppInd;
  int begin, end, rep,  posi, maxPosi, tag = 0, nChAbCp;
  int nActCpSub, nAbpCpSub, nActCpPre, nAbpCpPre, nMbCpSub, nMbCpPre;
  double r[NDIM], dist;
  ListInt actBurst, actCh;

  maxPosi = NEG_LARGE_VALUE; 
  nChAbCp = (motSA.gTgl != 0) ? nChAb : nChAb - 2;
  if (modeActCh != 0) {
	MALLOC(actCh.l,int,2 * nChAc - 1);
  }
  MALLOC(chkL,int,(mpiMethod == 0) ? cntAdjRank[0] : 2);
  if (actBst.tgl != 0) { MALLOC(actBurst.l, int, nActMe * 3); }
  rep = (mpiMethod == 0) ? 1 : NDIM;

  for(m = 0; m < rep; m++) {
	if (mpiMethod == 0) { 
		begin = 0;	
		end = NDIM; 
	}
	else { 
		begin = m;	
		end = m + 1; 
	}
	if (cntAdjRank[m] > 0) {
		memset(cntCpPar, 0, sizeof(int) * cntAdjRank[m] * 3); 
		for(n = 0; n < cntAdjRank[m]; n++) { 
			cpPar[n].c = 0; 
		}
		// Find particles to be synchronized
		for(n = 0; n < nActMe + nActCp + nAbpMe + nAbpCp + nMbMe + nMbCp; n++) {
			kind = SetKind(n, nActMe + nActCp, nAbpMe + nAbpCp);
			if (n < nActMe) { 
				CONT(ISACTM(n));
			}
			else if (n >= nActMe + nActCp && n < nActMe + nActCp + nAbpMe) {
				ind2 = n - nActMe - nActCp;
				CONT(ISABPIM(ind2));
			}
			if (kind == 0) { 
				V3COPY(r, &P2(act.r,n,0)); 
				absInd = act.id[n]; 
			}
			else if (kind == 1) {
				V3COPY(r, &P2(abp.r,n - nActMe - nActCp,0)); 
				absInd = abp.id[n - nActMe - nActCp] + nAct;
			}
			else {
				V3COPY(r, &P2(memb.r,n - nActMe - nActCp - nAbpMe - nAbpCp,0)); 
				absInd = memb.id[n - nActMe - nActCp - nAbpMe - nAbpCp] 
						+ nAct + nAbp;
			}
			if (rheoWay > 0 && bulkRheoType == 0) {
				ConvertRectDomainVector(r, 0);
			}
			CopyParticlesSubroutine(begin, end, oftCp, r, neiEdge);
			V2SET_ALL(idx, -1);
			// actin
			if (kind == 0) {
				pArr = &P2A(act.ch,n,0,nChAc);
				// If the actin is not the end of filaments,
				if (pArr[0] > -1) {
					// If no ABP is attached on the actin,
					// another actin connected to the actin is transferred.
					// This is for calculation of bending of actin filament.
					if (pArr[1] > -1) {
						for(k = 0; k < 2; k++) {
							CONT(!(iAct[pArr[k]] > -1 
									&& iAct[pArr[1 - k]] < 0));
							idx[k] = iAct[pArr[k]];
							break;
						}
					}
					// If ABP is attached, but if the ABP doesn't exist
					// in the subdomain, both actins connected to 
					// the actin are transferred.
					// This is for calculation of bending of right angle.
					if (idx[0] == -1) {
						for(k = 2; k < nChAc; k++) {
							CONT(pArr[k] < 0);
							CONT(iAbp[pArr[k]] > -1);
							idx[0] = iAct[pArr[0]];
							break;
						}
					}
				}
			}
			// ABP
			else if (kind == 1) {
				ind2 = n - nActMe - nActCp;
				if (abpF.bend[K_ABP(ind2)].facStf > 0. 
						&& !(ISMTF(K_ABP(ind2)))) {
					pArr = &P2A(abp.ch,ind2,0,nChAb);
					// If the ABP is connected to two act.filaments
					if (pArr[0] > -1 && pArr[1] > -1) {
						// If one actin belongs to a current domain,
						// and if the other actin belongs to the other domain,
						// the actin in the current domain should be 
						// transferred. This is for calculating bending between
						// two arms of ABPs.
						for(k = 0; k < 2; k++) {
							CONT(!(iAct[pArr[k]] > -1
									&& iAct[pArr[1 - k]] < 0));
							idx[0] = iAct[pArr[k]];
							idx[1] = iAct[P2A(act.ch,iAct[pArr[k]],0,nChAc)];
							break;
						}
					}
				}
			}
			// Membrane	
			else {
				ind2 = n - nActMe - nActCp - nAbpMe - nAbpCp;
				pArr = &P2A(memb.ch,ind2,0,nChMb);
				// 2D membrane
				if (dimMbNuc == 2) {
					for(k = 0; k < 2; k++) {
						CONT(!(iMb[pArr[k]] > -1 && iMb[pArr[1 - k]] < 0));
						idx[k] = iMb[pArr[k]] + nActMe 
								+ nActCp + nAbpMe + nAbpCp;
						break;
					}
				}
				// 3D membrane
				else {
					V3SET_ALL(&idx[2], -1);
					for(k = 0; k < 6; k++) {
						CONT(pArr[k] < 0);
						BREAK(iMb[pArr[k]] < 0);
					}
					if (k < 6) {
						cntMb = 0;
						for(k = 0; k < 6; k++) {
							if (iMb[pArr[k]] > -1) {
								idx[cntMb] = iMb[pArr[k]] + nActMe 
									+ nActCp + nAbpMe + nAbpCp;
								cntMb++;
							}
						}
					}
				}
			}
			// Very long chains are handled here..
			if (mpiMethod == 0) {
				CopyParticlesNormalSubroutine(n, oftCp, idx, 0); 
				// Check longCh.l and transfer to more adjacent subdomains
				// if necessary
				for(k = 0; k < longCh.c; k++) {
					pArr = &P2A(longCh.l,k,0,2);
					CS = -1;
					for(l = 0; l < 2; l++) {
						CONT(absInd != pArr[l]);
						oppInd = pArr[1 - l];
						if (oppInd < nAct) {
							if (iAct[oppInd] < 0) { CS = 1; }
						}
						else if (oppInd >= nAct && oppInd < nAct + nAbp) {
							if (iAbp[oppInd - nAct] < 0) { CS = 1; }
						}
						else {
							if (iMb[oppInd - nAct - nAbp] < 0) { CS = 1; }
						}
						break;
					}
					if (CS == 1) {
						dist = longChDist[k] + 1.5;
						CopyParticlesSubroutine(0, NDIM, oftCp, r, dist);
						CopyParticlesNormalSubroutine(n, oftCp, idx, 1); 
					}
				}
			}
			else {
				if (oftCp[m] != 1) {
					CopyParticlesPlymptonSubroutine(n, oftCp[m], 
							idx, m); 
				}
				else {
					CS = CopyParticlesPlymptonSubroutine2(n);
					if (CS == 1) {
						CopyParticlesPlymptonSubroutine(n, 3, idx, m); 
					}
				}
			}
			if (mpiMethod == 0) {
				// If actin, check longChIntMsg.l and transfer actin 
				// if it's in the list.
				for(k = 0; k < longChIntMsg.c; k++) {
					pArr = &P2A(longChIntMsg.l,k,0,2);
					CONT(pArr[0] != absInd);
					// Find the destination
					ind2 = -1;
					for(l = 0; l < 27; l++) {
						CONT(!(pArr[1] == P2A(adjRank,l,0,2)));
						ind2 = P2A(adjRank,l,1,2);
						break;
					}
					CONT(!(ind2 > -1));
					CS = InsertElement1dArrayWChk(cpPar[ind2].l, 
							&cpPar[ind2].c, n);
					CONT(!(CS == -1));
					P2A(cntCpPar,ind2,kind,3)++; 
				}
			}
			if (mpiMethod != 0 && kind == 0) {
				CS = FindElementArray(longChIntMsg.l, longChIntMsg.c, 
						absInd, 0, 1);
				if (CS != -1) {
					CopyParticlesPlymptonSubroutine(n, 3, idx, m); 
				}
			}
		}
	}
	// Pack information of the chosen particles into a message
	for(n = 0; n < cntAdjRank[m]; n++) {
		posi = 0;
		V3COPY(ind, &P2A(cntCpPar,n,0,3));
		MPI_PACK_INT(&P2A(cntCpPar,n,0,3), 3, n);
		if (tglActFormDyn != 0) {
			MPI_PACK_INT(&cntNucAssAnn, 1, n);
		}
		// actin
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_INT(&act.id[cpPar[n].l[k]], 1, n);
		}
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_DBL(&P2(act.r,cpPar[n].l[k],0), NDIM, n);
		}
		// Send the whole or filled chain information, depending on which
		// is a more efficient way..
		if (modeActCh == 0) { 
			for(k = 0; k < ind[0]; k++) {
				MPI_PACK_INT(&P2A(act.ch,cpPar[n].l[k],0,nChAc), nChAc, n);
			}
		}
		else {		
			for(k = 0; k < ind[0]; k++) {
				actInd = cpPar[n].l[k];
				V2COPY(&actCh.l[1], &P2A(act.ch,actInd,0,nChAc));
				actCh.c = 2;
				for(l = 2; l < nChAc; l++) {
					CONT(!(P2A(act.ch,actInd,l,nChAc) > -1));
					V2SET(&actCh.l[actCh.c + 1], l, 
							P2A(act.ch,actInd,l,nChAc));
					actCh.c += 2;
				}
				actCh.l[0] = actCh.c;
				MPI_PACK_INT(actCh.l, actCh.c + 1, n);
			}
		}
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_INT(&act.iF[cpPar[n].l[k]], 1, n);
		}
		for(k = 0; k < ind[0]; k++) {
			MPI_PACK_INT(&act.len[cpPar[n].l[k]], 1, n);
		}
		if (gTglMb != 0) { 
			for(k = 0; k < ind[0]; k++) {
				MPI_PACK_INT(&act.mbReb[cpPar[n].l[k]], 1, n);
			}
			if (mbSld.act.gTgl != 0) {
				for(k = 0; k < ind[0]; k++) {
					MPI_PACK_INT(&mbSld.act.reb[cpPar[n].l[k]], 1, n);
				}
			}
		}
		if (actBst.tgl != 0) {
			actBurst.c = 0;
			for(k = 0; k < actBst.fil.c; k++) { 
				for(l = 0; l < ind[0]; l++) {
					CONT(!(act.iF[cpPar[n].l[l]] == P2A(actBst.fil.l,k,0,3)));
					V3COPY(&P2A(actBurst.l,actBurst.c,0,3), 
							&P2A(actBst.fil.l,k,0,3));
					(actBurst.c)++;
					break;
				}
			}
			MPI_PACK_INT(&actBurst.c, 1, n);
			if (actBurst.c > 0) { 
				MPI_PACK_INT(actBurst.l, actBurst.c * 3, n);
			}
		}
		// ABP
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_INT(&abp.id[cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp], 1, n);
		}
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_DBL(&P2(abp.r,cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp,0), NDIM, n);
		}
		for(k = 0; k < ind[1]; k++) {
			MPI_PACK_INT(&P2A(abp.ch,cpPar[n].l[k + ind[0]] 
					- nActMe - nActCp,0,nChAb), nChAbCp, n);
		}
		if (motSA.gTgl != 0) {
			for(k = 0; k < ind[1]; k++) {
				MPI_PACK_INT(&abp.mId[cpPar[n].l[k + ind[0]] 
						- nActMe - nActCp], 1, n);
			}
		}
		if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
			for(k = 0; k < ind[1]; k++) {
				MPI_PACK_INT(&mbSld.abp.reb[cpPar[n].l[k + ind[0]] 
						- nActMe - nActCp], 1, n);
			}
		}
		if (gTglMb != 0) {
			// Membrane
			for(k = 0; k < ind[2]; k++) {
				mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
						- nActMe - nActCp - nAbpMe - nAbpCp;
				MPI_PACK_INT(&memb.id[mbInd], 1, n);
			}
			for(k = 0; k < ind[2]; k++) {
				mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
						- nActMe - nActCp - nAbpMe - nAbpCp;
				MPI_PACK_DBL(&P2(memb.r,mbInd,0), NDIM, n);
			}
			for(k = 0; k < ind[2]; k++) {
				mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
						- nActMe - nActCp - nAbpMe - nAbpCp;
				MPI_PACK_DBL(&P2(memb.nDir,mbInd,0), NDIM, n);
			}
			for(k = 0; k < ind[2]; k++) {
				mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
						- nActMe - nActCp - nAbpMe - nAbpCp;
				MPI_PACK_INT(&P2A(memb.ch,mbInd,0,nChMb), nChMb, n);
			}
			for(k = 0; k < ind[2]; k++) {
				mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
						- nActMe - nActCp - nAbpMe - nAbpCp;
				MPI_PACK_INT(&memb.idx[mbInd], 1, n);
			}
			if (mbMat.gTgl != 0) {
				for(k = 0; k < ind[2]; k++) {
					mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
							- nActMe - nActCp - nAbpMe - nAbpCp;
					MPI_PACK_INT(&P2A(memb.nFA,mbInd,0,nMbAct), nMbAct, n);
				}
			}
			if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
				for(k = 0; k < ind[2]; k++) {
					mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
							- nActMe - nActCp - nAbpMe - nAbpCp;
					MPI_PACK_INT(&memb.reb[mbInd], 1, n);
				}
			}
			if (mbPro.gTgl != 0) {
				for(k = 0; k < ind[2]; k++) {
					mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
							- nActMe - nActCp - nAbpMe - nAbpCp;
					MPI_PACK_INT(&P2A(memb.pro,mbInd,0,2), 2, n);
				}
				if (gTglMbCont != 0) {
					for(k = 0; k < ind[2]; k++) {
						mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
								- nActMe - nActCp - nAbpMe - nAbpCp;
						MPI_PACK_DBL(&memb.myo.conc[mbInd], 1, n);
					}
					for(k = 0; k < ind[2]; k++) {
						mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
								- nActMe - nActCp - nAbpMe - nAbpCp;
						MPI_PACK_DBL(&memb.cyto.conc[mbInd], 1, n);
					}
				}
			}
			if (tglBias != 0 && typeBias == 2) {
				for(k = 0; k < ind[2]; k++) {
					mbInd = cpPar[n].l[k + ind[0] + ind[1]] 
							- nActMe - nActCp - nAbpMe - nAbpCp;
					MPI_PACK_INT(&memb.bias[mbInd], 1, n);
				}
			}
		}	
		// Send and receive the messages to adjacent subdomains
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, iRank[(mpiMethod == 0) 
				? n : 2 * m + n], tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, 
				iRank[(mpiMethod == 0) ? n : 2 * m + n], 
				tag, MPI_COMM_WORLD, &rReq[n]);
		if (posi > maxPosi) { maxPosi = posi; }
	}
	if (cntAdjRank[m] > 0) {
		// Receive messages and unpack the transferred information.
		cnt = 0;
		cntRecvMsg = 0;
		memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
		memset(chkL, 0, sizeof(int) * cntAdjRank[m]);

		nActCpPre = nActCp;
		nAbpCpPre = nAbpCp;
		nMbCpPre = nMbCp;
		while(cntRecvMsg < cntAdjRank[m]) {
			if (mpiTestSendFlag[cnt] == 0) {
				MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
			}
			if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0 
					&& chkL[cnt] == 0) {
				MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
				if (mpiTestRecvFlag[cnt] != 0) { 
					posi = 0;
			  		MPI_UNPACK_INT(&nActCpSub, 1, cnt); 
			  		MPI_UNPACK_INT(&nAbpCpSub, 1, cnt);
			  		MPI_UNPACK_INT(&nMbCpSub, 1, cnt);
					if (tglActFormDyn != 0) { 
						MPI_UNPACK_INT(&cntNucAssAnn, 1, cnt);
					}
					// actin
			  		MPI_UNPACK_INT(&act.id[nActMe + nActCp], 
							nActCpSub, cnt);
			  		MPI_UNPACK_DBL(&P2(act.r,nActMe + nActCp,0), 
							nActCpSub * NDIM, cnt);
					if (modeActCh == 0) {
			  			MPI_UNPACK_INT(&P2A(act.ch,nActMe + nActCp,0,nChAc), 
								nActCpSub * nChAc, cnt);
					}
					else {
						memset(&P2A(act.ch,nActMe + nActCp,0,nChAc), -1, 
								sizeof(int) * nActCpSub * nChAc);
						for(k = 0; k < nActCpSub; k++) {
							pArr = &P2A(act.ch,nActMe + nActCp + k,0,nChAc);
				  			MPI_UNPACK_INT(&actCh.c, 1, cnt);
				  			MPI_UNPACK_INT(actCh.l, actCh.c, cnt);
							V2COPY(pArr, actCh.l);
							for(l = 2; l < actCh.c; l += 2) {
								pArr[actCh.l[l]] = actCh.l[l + 1];
							}
						}
					}
			  		MPI_UNPACK_INT(&act.iF[nActMe + nActCp], nActCpSub, cnt);
					MPI_UNPACK_INT(&act.len[nActMe + nActCp], nActCpSub, cnt);
					if (gTglMb != 0) { 
			  			MPI_UNPACK_INT(&act.mbReb[nActMe + nActCp], 
								nActCpSub, cnt);
						if (mbSld.act.gTgl != 0) {
							MPI_UNPACK_INT(&mbSld.act.reb[nActMe + nActCp], 
									nActCpSub, cnt);
						}
					}
					if (actBst.tgl != 0) {
			  			MPI_UNPACK_INT(&actBurst.c, 1, cnt);
						if (actBurst.c > 0) {
			  				MPI_UNPACK_INT(actBurst.l, actBurst.c * 3, cnt);
						}
						for(k = 0; k < actBurst.c; k++) {
							CS = FindElementArray(actBst.fil.l, actBst.fil.c, 
									P2A(actBurst.l,k,0,3), 0, 3);
							CONT(CS != -1);
							V3COPY(&P2A(actBst.fil.l,actBst.fil.c,0,3), 
									&P2A(actBurst.l,k,0,3));
							(actBst.fil.c)++;
						}
					}
					// ABP
			  		MPI_UNPACK_INT(&abp.id[nAbpMe + nAbpCp], nAbpCpSub, cnt);
			  		MPI_UNPACK_DBL(&P2(abp.r,nAbpMe + nAbpCp,0), 
							nAbpCpSub * NDIM, cnt);
					if (motSA.gTgl == 0) {
						for(k = 0; k < nAbpCpSub; k++) {
				  			MPI_UNPACK_INT(&P2A(abp.ch,nAbpMe + nAbpCp + k,0,
								nChAb),	nChAb - 2, cnt);
							V2SET_ALL(&P2A(abp.ch,nAbpMe + nAbpCp + k,3,nChAb), 
									-1);
						}
					}
					else {
			  			MPI_UNPACK_INT(&P2A(abp.ch,nAbpMe + nAbpCp,0,nChAb), 
								nAbpCpSub * nChAb, cnt);
		  				MPI_UNPACK_INT(&abp.mId[nAbpMe + nAbpCp], nAbpCpSub, 
								cnt);
					}
					if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
						MPI_UNPACK_INT(&mbSld.abp.reb[nAbpMe + nAbpCp], 
								nAbpCpSub, cnt);
					}
					if (gTglMb != 0) {
						// Membrane
						mbInd = nMbMe + nMbCp;
				  		MPI_UNPACK_INT(&memb.id[mbInd], nMbCpSub, cnt);
				  		MPI_UNPACK_DBL(&P2(memb.r,mbInd,0), nMbCpSub * NDIM, 
								cnt);
				  		MPI_UNPACK_DBL(&P2(memb.nDir,mbInd,0), 
								nMbCpSub * NDIM, cnt);
				  		MPI_UNPACK_INT(&P2A(memb.ch,mbInd,0,nChMb), 
								nMbCpSub * nChMb, cnt);
				  		MPI_UNPACK_INT(&memb.idx[mbInd], nMbCpSub, cnt);
						if (mbMat.gTgl != 0) {
					  		MPI_UNPACK_INT(&P2A(memb.nFA,mbInd,0,nMbAct), 
									nMbCpSub * nMbAct, cnt);
						}
						if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
					  		MPI_UNPACK_INT(&memb.reb[mbInd], nMbCpSub, cnt);
						}
						if (mbPro.gTgl != 0) {
					  		MPI_UNPACK_INT(&P2A(memb.pro,mbInd,0,2), 
									nMbCpSub * 2, cnt);
							if (gTglMbCont != 0) {
						  		MPI_UNPACK_DBL(&memb.myo.conc[mbInd], 
										nMbCpSub, cnt);
						  		MPI_UNPACK_DBL(&memb.cyto.conc[mbInd], 
										nMbCpSub, cnt);
							}
						}
						if (tglBias != 0 && typeBias == 2) {
							MPI_UNPACK_INT(&memb.bias[mbInd], nMbCpSub, cnt);
						}
					}
					if (tglActFormDyn != 0) {
						if (cntNucAssAnn > 0) {
							for(k = 0; k < nActCpSub; k++) {
								InsertElement1dArrayWoChk(insNeiPar.l, 
									&insNeiPar.c, act.id[nActMe + nActCp + k]);
							}
						}
					}
				
					nActCp += nActCpSub;
					nAbpCp += nAbpCpSub;
					if (gTglMb != 0) { nMbCp += nMbCpSub; }
					cntRecvMsg++;
					chkL[cnt] = 1;
				}
			}
			cnt++;
			if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
		for(n = nActCpPre; n < nActCp; n++) {
			iAct[act.id[n + nActMe]] = n + nActMe;
		}
		for(n = nAbpCpPre; n < nAbpCp; n++) {
			iAbp[abp.id[n + nAbpMe]] = n + nAbpMe;
		}
		for(n = nMbCpPre; n < nMbCp; n++) {
			iMb[memb.id[n + nMbMe]] = n + nMbMe;
		}
	}
  }
  free(chkL);
  if (actBst.tgl != 0) { free(actBurst.l); } 
  if (modeActCh != 0) { free(actCh.l); }
  
  if (maxPosi > sizeBufMsg / 2) {
	end = (mpiMethod == 0) ? cntAdjRank[0] : 2;
	sizeBufMsg *= 2;
	for(n = 0; n < end; n++) {
		free(bufSendMsg[n]);
		free(bufRecvMsg[n]);
		MALLOC(bufSendMsg[n],char,sizeBufMsg);
		MALLOC(bufRecvMsg[n],char,sizeBufMsg);
	}	
  }
  for(n = 0; n < insNeiPar.c; n++) {
	ind[0] = insNeiPar.l[n];
	if (ind[0] < nAct) {
		ind2 = iAct[ind[0]];
		pArr = &P2A(act.ch,ind2,0,nChAc);
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			CONT(iAct[pArr[k]] < 0);
			CS = FindElementArray(act.cyl.l, act.cyl.c, ind[0], k, 2);
			CONT(CS > -1);
			InsertElementInNeighborList(ind[0], pArr[k], k);
		}
	}
	else {
  		kind = SetKind(ind[0], nAct + nAbp, nUnitMb);
		if (kind == 0) { ind[0] -= nAct; }
		else if (kind == 1) { ind[0] -= nAct + nAbp; }
		else if (kind == 2) { ind[0] -= nAct + nAbp + nUnitMb; }
		InsertElementInNeighborList(ind[0], -1, kind + 2);
	}	
  }
}

/*--------------------- Copy particles between subdomains --------------------*/

/*-------------------- Process conflicts and dynamic events ------------------*/

void UpdateAbpUnbRebLists(int abpInd, int actInd, int side, int mode) {
  V3SET(&P2A(sendAbpDyn.l,sendAbpDyn.c,0,3), abpInd, actInd, side);
  (sendAbpDyn.c)++;
  V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, currTimeStep);
  (noAbpDyn.c)++;
}

void UpdateActinSeverEvents(void) {
  int n, CS, actInd, iFila, curr, sizeArr;
  ListInt all, sendActSevL;

  sendActSevL.c = 0;
  sizeArr = sendActDyn.siz;
  MALLOC(sendActSevL.l,int,sizeArr); 
  sizeArr *= (mpiMethod == 0) ? cntAdjRank[0] : 2;
  MALLOC(all.l,int,sizeArr);
  all.c = sendActDyn.c;
  // Gather sendActDyn.l from adjacent subdomains
  for(n = 0; n < sendActDyn.c; n++) {
	V3COPY(&P2A(all.l,n,0,3), &P2A(sendActDyn.l,n,0,3));
  }
  CollectArrayIntFromAdjacentSubdomain(&all, 3);
  for(n = 0; n < all.c; n++) {
	actInd = P2A(all.l,n,0,3);
	iFila = P2A(all.l,n,1,3);
	CONT(!(iAct[actInd] > -1));
	CS = (iAct[actInd] < nActMe) ? 1 : 0;
	curr = actInd;
	while (curr > -1 && curr < nAct) {
		if ((iAct[curr] >= nActMe && CS == 1) || iAct[curr] < 0) {
			if (n >= sendActDyn.c) {
				V2SET(&P2A(sendActSevL.l,sendActSevL.c,0,2), curr, iFila);
				(sendActSevL.c)++;
			}
			if (iAct[curr] < 0) { break; }
			CS = 0;
		}
		else if (iAct[curr] > -1 && iAct[curr] < nActMe && CS == 0) {
			CS = 1;
		}	
		if  (act.iF[iAct[curr]] == iFila) { break; }
		else { act.iF[iAct[curr]] = iFila; }
		curr = P2A(act.ch,iAct[curr],0,nChAc);
	}
  }
  CheckArraySize(&sendActDyn, &sendActDyn.siz, 3, 0);
  for(n = 0; n < sendActSevL.c; n++) {
	V2COPY(&P2A(sendActDyn.l,n,0,3), &P2A(sendActSevL.l,n,0,2));
	P2A(sendActDyn.l,n,2,3) = -4;
  }
  sendActDyn.c = sendActSevL.c;
  free(all.l);
  free(sendActSevL.l);
}

// Subroutine for UpdateActinAbpDynamicsEventsSubroutine()
void UpdateActinAbpDynamicsEventsSubSubroutine(ListInt *all, int abpInd, 
		int actInd, int side, ListInt *confAct, int *CS) {
  int locActInd, actInd2, locActInd2, locAbpInd, abpInd2, locAbpInd2, loc[2];
  int abpRankMol, abpRankMol2, actRankMol, abpSide, *pArr, *pArr2, CS2;

  locActInd = iAct[actInd];
  locAbpInd = iAbp[abpInd];
  abpInd2 = P2A(act.ch,locActInd,side,nChAc);
  actRankMol = CalcRankMolecule(&P2(act.r,locActInd,0));
  // It is possible that multiple ABPs belonging to different subdomains 
  // try to bind to the same actin. Check the possibility here.
  // If so, only ABP in the same subdomain as that of actin can bind to actin.
  if (abpInd2 > -1 && abpInd2 != abpInd) {
	CS2 = Find2ElementArray(confAct->l, confAct->c, locActInd, side, 0, 2);
	if (CS2 == -1) {
		V2SET(&P2A(confAct->l,confAct->c,0,2), locActInd, side);
		(confAct->c)++;
	}
	//P2A(confAct,locActInd,side - 2,2) = 1;
	locAbpInd2 = iAbp[abpInd2];
	if (locAbpInd2 > -1) {
		pArr2 = &P2A(abp.ch,locAbpInd2,0,nChAb); 
		abpRankMol2 = CalcRankMolecule(&P2(abp.r,locAbpInd2,0));
		// If they belong to different subdomains
		if (abpRankMol2 != actRankMol) {
			CS2 = 1;
			if (pArr2[0] == actInd) { abpSide = 0; }
			else if (pArr2[1] == actInd) { abpSide = 1; }
			else { 
				CS2 = -1;
			}
			if (CS2 > -1) {
				// Sever a link between actin and the pre-existing ABP
				UpdateActinDisassemblySubroutine2(abpInd2, actInd);
				// Check whether the pre-existing ABP originates from walking
				CS2 = FindElementArray(all->l, all->c, abpInd2, 0, 3);
			}
			if (CS2 > -1 && CS2 < all->c - 1) {
				pArr = &P2A(all->l,CS2,0,3);
				// If it originates from walking, the pre-existing ABP should
				// return to the previous actin.
				loc[0] = (int)((P2A(pArr,0,2,3) - 2) / nChAcY);
				loc[1] = (int)((P2A(pArr,1,2,3) - 2) / nChAcY);
				if (abpInd2 == P2A(pArr,1,0,3) 
						&& ((P2A(pArr,0,1,3) != P2A(pArr,1,1,3) 
						&& loc[0] == nChAcX - 1 && loc[1] == 0) 
						|| (P2A(pArr,0,1,3) == P2A(pArr,1,1,3) 
						&& loc[1] - loc[0] == 1))) { 
					// actInd2 is the previous actin to return.
					actInd2 = P2A(pArr,0,1,3);
					locActInd2 = iAct[actInd2];
					if (locActInd2 > -1) {
						P2A(act.ch,locActInd2,P2A(pArr,0,2,3),nChAc) = abpInd2;
					}
					if (abpSide == 0) {
						pArr2[1] = pArr2[0];
					}
					pArr2[abpSide] = actInd2;
					if (locAbpInd2 < nAbpMe) {
						if (pArr2[0] > -1 && pArr2[1] > -1) {
							if (tglNeiAbpSC != 0 && tglNeiAbpDC == 0) {
								DeleteElementInNeighborList(abpInd2, 1);
							}
						}
						(motWalk.cntMe)--;
						if (pArr2[1] < 0) { (motInaUnb.cntMe)--; } 
						else { (motUnb.cntMe)--; }
						if (recAbpDyn.tgl != 0) {
							abpDyn.cntW[locAbpInd2]--;
							abpDyn.cntU[locAbpInd2]--;
						}
					}
					V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), 
							abpInd2, actInd2, currTimeStep);
					(noAbpDyn.c)++;
					if (mpiMethod == 0 && actInd != actInd2
							&& ((locActInd2 > -1 && locActInd2 < nActMe) 
							|| (locAbpInd2 > -1 && locAbpInd2 < nAbpMe))) {
						InsertLongChain(abpInd2 + nAct, actInd2, 
								minDimDomC * 0.9);
					}
				}
			}
			// if (abpRankMol2 != actRankMol)
			// and if information of the current ABP is available
			if (locAbpInd > -1) {
				abpRankMol = CalcRankMolecule(&P2(abp.r,locAbpInd,0));
				*CS = (abpRankMol == actRankMol) ? 1 : 0;
			}
			// if information of the current ABP is not available
			else { *CS = 1; }
		}
		// if (abpRankMol2 == actRankMol)
		else { *CS = 0; }
	}
  }
  // If there was a conflict on the binding spot, but if it is currently
  // avalable
  else if (abpInd2 < 0) {
	CS2 = Find2ElementArray(confAct->l, confAct->c, locActInd, side, 0, 2);
	if (CS2 > -1) {
		if (locAbpInd > -1) {
			abpRankMol = CalcRankMolecule(&P2(abp.r,locAbpInd,0));
			*CS = (actRankMol == abpRankMol) ? 1 : 0;
		}
		else { *CS = 1; }
	}
	else { *CS = 1; }
  }
}

// Subroutine for UpdateActinAbpDynamicsEvents()
void UpdateActinAbpDynamicsEventsSubroutine(ListInt *all, int mode) {
  int n, k, side, walk, begin, curr, CS, CS2, *pArr, *pArr2, loc[2];
  int abpInd, actInd, actInd2, locAbpInd, locActInd, nextActInd, locNextActInd;
  ListInt sendActSevL, confAct;

  if (actSev.tgl != 0) { 
	MALLOC(sendActSevL.l,int,nActMin*2); 
	sendActSevL.c = 0;
  }
  confAct.l = allIntL;
  confAct.c = 0;
  begin = sendAbpDyn.c + ((mode == 0) ? sendActDyn.c : 0);
  for(n = begin; n < all->c; n++) {
	// ABP unbinding/binding/walking
	if (P2A(all->l,n,2,3) >= 0) {
		abpInd = P2A(all->l,n,0,3);
		actInd = P2A(all->l,n,1,3);
		locAbpInd = iAbp[abpInd];
		locActInd = iAct[actInd];
		// Check whether it is a walking action or not.
		// If the first element is the same in two consecutive rows,
		// it indicates the walking.
		walk = 0;
		CS2 = 1;
		if (n < (all->c) - 1 && motWalk.tgl != 0) {
			if (abpInd == P2A(all->l,n + 1,0,3) && P2A(all->l,n + 1,2,3) >= 0) {
				loc[0] = (int)((P2A(all->l,n,2,3) - 2) / nChAcY);
				loc[1] = (int)((P2A(all->l,n + 1,2,3) - 2) / nChAcY);
				if ((actInd != P2A(all->l,n + 1,1,3) && loc[0] == nChAcX - 1
						&& loc[1] == 0) || (actInd == P2A(all->l,n + 1,1,3)
						&& loc[1] - loc[0] == 1)) {
					nextActInd = P2A(all->l,n + 1,1,3);
					walk = 1;
					if (iAct[nextActInd] > -1) {
						side = P2A(all->l,n + 1,2,3);
						// Check whether a conflict exists on the binding spot 
						// to walk.
						UpdateActinAbpDynamicsEventsSubSubroutine(all, abpInd,
								nextActInd, side, &confAct, &CS2);
					}
				}
			}
		}
		if (locAbpInd >= nAbpMe || locAbpInd < 0) {
			CS = 1;
			if (locActInd > -1 && !(CS2 == 0 && walk == 1)) {
				pArr = &P2A(act.ch,locActInd,0,nChAc);
				side = P2A(all->l,n,2,3);
				UpdateActinAbpDynamicsEventsSubSubroutine(all, abpInd, 
						actInd, side, &confAct, &CS);
				// Modify act.ch
				if (CS == 1) {
					if (mpiMethod == 0 && locActInd < nActMe) {
						if (pArr[side] == abpInd) {
							DeleteLongChain(abpInd + nAct, actInd);
						}
						else {
							InsertLongChain(abpInd + nAct, actInd, 
									minDimDomC * 0.9);
						}
					}
					pArr[side] = (pArr[side] == abpInd) ? -1 : abpInd;
				}
			}	
			if (locAbpInd >= nAbpMe) {
				// Modify abp.ch
				pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
				// If it is a successful walking motion
				if (CS2 == 1 && walk == 1) {
					pArr2[(pArr2[0] == actInd) ? 0 : 1] = nextActInd;
				}
				// If it is a simple binding or unbinding
				else if (CS == 1 && walk == 0) {
					if (pArr2[0] == actInd) {
				       	pArr2[0] = pArr2[1];
				       	pArr2[1] = -1;	
					}
					else if (pArr2[1] == actInd) {
						pArr2[1] = -1; 
					}
					else { 
						pArr2[(pArr2[0] < 0) ? 0 : 1] = actInd; 
					}
				}
			}
			V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), abpInd, actInd, 
					currTimeStep);
			(noAbpDyn.c)++;
		}
		// From actin dissembly or severing
		// (Among ABP dynamics, unbinding of ABPs by actin disassembly and
		// severing is only a thing which can occur at the side of actin.
		if (locAbpInd > -1 && locAbpInd < nAbpMe && locActInd > -1) {
			pArr2 = &P2A(abp.ch,locAbpInd,0,nChAb);
			// By chance, ABP might have been unbound already by its unbinding 
			// or walking at the same time step, so it needs to check the bond.
			if (pArr2[0] == actInd || pArr2[1] == actInd) {
				// Sever the chain
				UpdateActinDisassemblySubroutine2(abpInd, actInd);
			}
		}
		// If walking event
		if (walk == 1) {
			if (CS2 == 1) {
				locNextActInd = iAct[nextActInd];
				if (locNextActInd > -1) {
					side = P2A(all->l,n + 1,2,3);
					P2A(act.ch,locNextActInd,side,nChAc) = abpInd;
					if (locNextActInd < nActMe) {
						if (mpiMethod == 0 && actInd != nextActInd) {
							InsertLongChain(abpInd + nAct, nextActInd, 
									minDimDomC * 0.9);
						}
					}
					V3SET(&P2A(noAbpDyn.l,noAbpDyn.c,0,3), 
							abpInd, nextActInd, -1 * currTimeStep);
					(noAbpDyn.c)++;
				}
			}	
			n++; 
		}
	}
	// Actin disassembly
	else if (P2A(all->l,n,2,3) == -1 || P2A(all->l,n,2,3) == -2) {
		actInd = P2A(all->l,n,0,3);
		actInd2 = P2A(all->l,n,1,3);
		side = -1 * P2A(all->l,n,2,3) - 1;
		locActInd = iAct[actInd];
		pArr = &P2A(act.ch,locActInd,0,nChAc);
		// This must be among copied particles
		if (locActInd > -1) { 
	        if (pArr[0] > -1 || pArr[1] > -1) {
				// Check whether or not ABPs are bound on it
				if (side == 1) {
					for(k = 2; k < nChAc; k++) {
						abpInd = pArr[k];
						CONT(!(abpInd > -1));
						UpdateActinDisassemblySubroutine2(abpInd, actInd);
					}
				}
		        V3SET_ALL(&P2(act.r,locActInd,0), 0.);
		        SetAllValue1dArrayInt(pArr, nChAc, -1);
				act.iF[locActInd] = -1;
			}
		  	V2SET(&P2A(noActDyn.l,noActDyn.c,0,2), actInd, currTimeStep);
		  	(noActDyn.c)++;
		}
		// This can be among copied particles or among current subdomain
		if (iAct[actInd2] > -1) { 
			pArr2 = &P2A(act.ch,iAct[actInd2],0,nChAc);
	        if (pArr2[0] > -1 || pArr2[1] > -1) {
				// Check whether or not ABPs are bound on it
				if (side == 0) {
					for(k = 2; k < nChAc; k++) {
						abpInd = pArr2[k];
						CONT(!(abpInd > -1));
						UpdateActinDisassemblySubroutine2(abpInd, actInd2);
					}
				}
				UpdateActinDisassemblySubroutine(actInd2, side);
			}
		}
		// Update the neighboring list
		// Find and delete an actin cylinder corresponding the disassembled one
		//DeleteElementInNeighborList(P2A(all->l,n,0,3), 0);
		DeleteElementInNeighborList(actInd, 0);
	}
	// Adjust chain information for actin severing
	else if (P2A(all->l,n,2,3) == -3) {
		UpdateActinSeveringSubroutine(&P2A(all->l,n,0,3));
	}
	// Adjust chain information for actin severing
	else if (P2A(all->l,n,2,3) == -5) {
		UpdateActinAnnealingSubroutine(&P2A(all->l,n,0,3));
	}
	// Change the index of actin filament for actin severing
	else if (P2A(all->l,n,2,3) == -4) {
		actInd = P2A(all->l,n,0,3);
		locActInd = iAct[actInd];
		// If it is in the current subdomain or is one of the copied particles
		if (locActInd > -1) {
			// If it has a different filament index from one in the message
			// Change toward pointed end
			curr = P2A(act.ch,locActInd,1,nChAc);
			while (curr > -1 && curr < nAct) {
				if (iAct[curr] > -1) {
					if  (act.iF[iAct[curr]] == P2A(all->l,n,1,3)) { break; }
					else { act.iF[iAct[curr]] = P2A(all->l,n,1,3); }
				}
				else { break; }
				curr = P2A(act.ch,iAct[curr],1,nChAc);
			}
			// Change toward barbed end
			UpdateActinSeverAnnealSubroutine(actInd, P2A(all->l,n,1,3), 
					&sendActSevL, 1);
		}
	}
  }
  CheckArraySize(&sendAbpDyn, &sendAbpDyn.siz, 3, 0);
  sendAbpDyn.c = 0;
  if (mode == 0) {
	CheckArraySize(&sendActDyn, &sendActDyn.siz, 3, 0);
  }
  if (actSev.tgl != 0) { 
	for(n = 0; n < sendActSevL.c; n++) {
		V2COPY(&P2A(sendActDyn.l,n,0,3), &P2A(sendActSevL.l,n,0,2));
		P2A(sendActDyn.l,n,2,3) = -4;
	}
	sendActDyn.c = sendActSevL.c;
	free(sendActSevL.l);
  }
  else { sendActDyn.c = 0; }
}

// This function let other subdomains know the happening of the dynamic 
// behaviors of ABPs (e.g. unbinding, binding, and walking). 
void UpdateActinAbpDynamicsEvents(int mode) {
  int n, sizeArr;
  ListInt all;

  sizeArr = sendAbpDyn.siz + ((mode == 0) ? sendActDyn.siz : 0.);
  sizeArr *= (mpiMethod == 0) ? cntAdjRank[0] : 2;
  MALLOC(all.l,int,sizeArr);

  all.c = sendAbpDyn.c;
  for(n = 0; n < sendAbpDyn.c; n++) {
	V3COPY(&P2A(all.l,n,0,3), &P2A(sendAbpDyn.l,n,0,3));
  }
  if (mode == 0) {
	(all.c) += sendActDyn.c;
	for(n = 0; n < sendActDyn.c; n++) {
		V3COPY(&P2A(all.l,n + sendAbpDyn.c,0,3), &P2A(sendActDyn.l,n,0,3));
	}
  }
  CollectArrayIntFromAdjacentSubdomain(&all, 3);
  UpdateActinAbpDynamicsEventsSubroutine(&all, mode);
  free(all.l);
}

/*-------------------- Process conflicts and dynamic events ------------------*/

/*-------------------- Balancing counters between subdomains -----------------*/

// It is possible that some subdomain has no actin for motor nucleation. Then,
// the counter for possible nucleation, "motSA.nNucMe", is balanced between
// subdomains in this function.
// This function is not called after all the nucleation processes occur.
void UpdateMotorNucleCounter(void) {
  int *motSAnNucAll, sum, base, rem, ind;

  MALLOC(motSAnNucAll,int,nCpu);
  MPI_Gather(&motSA.nNucMe, 1, MPI_INT, motSAnNucAll, 1, MPI_INT, 
		0, MPI_COMM_WORLD);
  if (rank == 0) {
	sum = SumArrInt(motSAnNucAll, nCpu);
	if (sum > 0) {
		base = (int)(sum / nCpu);
		rem = sum % nCpu;
		SetAllValue1dArrayInt(motSAnNucAll, nCpu, base);
		while(rem > 0) {
			ind = GenRandIntIndex(nCpu);
			CONT(motSAnNucAll[ind] > base);
			motSAnNucAll[ind]++;
			rem--;
		}	
	}
	else {
		SetAllValue1dArrayInt(motSAnNucAll, nCpu, -1);
	}
  }
  MPI_Scatter(motSAnNucAll, 1, MPI_INT, &motSA.nNucMe, 1, MPI_INT, 
			0, MPI_COMM_WORLD);
  if (motSA.nNucMe == -1) {
	updMotN.tgl = 0;
	motSA.nNucMe = 0;
  }
  free(motSAnNucAll);
}

// Fill gap by information of transferred elements
// mode = 0: actin, 1: ACP, 2: motor
void UpdateActinAbpMonomerListSubroutine2(int locInd, int ind, int mode) {
  int k;

  if (mode == 0) {
	V3SET_ALL(&P2(act.r,locInd,0), 0.);
	SetAllValue1dArrayInt(&P2A(act.ch,locInd,0,nChAc), nChAc, -1);
	act.id[locInd] = ind;
	act.iF[locInd] = -1;
	act.fix[locInd] = -1;
	if (gTglBead != 0 && beadBind.gTgl != 0) { act.bdFix[locInd] = -1; }
	if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) { act.nFA[locInd] = 0; }
	if (gTglActCapAll != 0) { act.cap[locInd] = -1; }
	if (actAge.gTgl != 0) { act.age[locInd] = 0; }
	act.len[locInd] = 0;
	if (gTglMb != 0) { 
		act.mbReb[locInd] = -1; 
		if (mbSld.act.gTgl != 0) {
			mbSld.act.reb[locInd] = -1;
			SetAllValue1dArrayInt(&P2A(mbSld.act.l,locInd,0,dimMbNuc), 
					dimMbNuc, -1);
		}
	}
	iAct[ind] = locInd;
	if (confVmdInfo.tgl != 0) {
		recAct.len[locInd] = 0.;
		recAct.sprF[locInd] = 0.;
		recAct.bendF[locInd] = 0.;
		V4SET_ALL(&P2A(recAct.allF,locInd,0,NDIM + 1), 0.);
		recAct.cnt[locInd] = 0;
	}
  }
  else {
	V3SET_ALL(&P2(abp.r,locInd,0), 0.);
	SetAllValue1dArrayInt(&P2A(abp.ch,locInd,0,nChAb), nChAb, -1);
	P2A(abp.ch,locInd,2,nChAb) = mode - 1;
	if (abpAge.gTgl != 0) { abp.age[locInd] = 0; }
	if (motSA.gTgl != 0) {
		abp.mId[locInd] = -1;
	}
	if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
		mbSld.abp.reb[locInd] = -1;
		SetAllValue1dArrayInt(&P2A(mbSld.abp.l,locInd,0,dimMbNuc), 
				dimMbNuc, -1);
	}
	if (recAbpDyn.tgl != 0) {
		abpDyn.cntU[locInd] = 0;
		abpDyn.cntB[locInd] = 0;
		abpDyn.cntW[locInd] = 0;
	}
	abp.id[locInd] = ind;
	iAbp[ind] = locInd;

	if (confVmdInfo.tgl != 0) {
		SetAllValue1dArrayDouble(&P2A(recAbp.len,locInd,0,recAbp.nL), 
				recAbp.nL, 0.);
		recAbp.sprF[locInd] = 0.;
		recAbp.bendF[locInd] = 0.;
		V4SET_ALL(&P2A(recAbp.allF,locInd,0,NDIM + 1), 0.);
		recAbp.cnt[locInd] = 0;
	}
	if (recLongF.tgl != 0) {	
		V4SET_ALL(&P2A(recLongSprFabp,locInd,0,4), 0.);
	}
  }
}

// Shift all the information of copied elements from locInd2 to locInd1
// mode = 0: actin, 1: ABP
void UpdateActinAbpMonomerListSubroutine(int locInd1, int locInd2, int mode) {
  int n, k;

  if (mode == 0) {
	V3COPY(&P2(act.r,locInd1,0), &P2(act.r,locInd2,0));
	Copy1dArrayInt(&P2A(act.ch,locInd1,0,nChAc), 
			&P2A(act.ch,locInd2,0,nChAc), nChAc);
	if (gTglMb != 0 && mbSld.act.gTgl != 0) {
		mbSld.act.reb[locInd1] = mbSld.act.reb[locInd2];
		Copy1dArrayInt(&P2A(mbSld.act.l,locInd1,0,dimMbNuc), 
				&P2A(mbSld.act.l,locInd2,0,dimMbNuc), dimMbNuc);
	}
	act.iF[locInd1] = act.iF[locInd2];
	act.id[locInd1] = act.id[locInd2];
	iAct[act.id[locInd2]] = locInd1;
  }
  else {
	V3COPY(&P2(abp.r,locInd1,0), &P2(abp.r,locInd2,0));
	Copy1dArrayInt(&P2A(abp.ch,locInd1,0,nChAb), 
			&P2A(abp.ch,locInd2,0,nChAb), nChAb);
	if (motSA.gTgl != 0) {
		abp.mId[locInd1] = abp.mId[locInd2];
	}
	if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
		mbSld.abp.reb[locInd1] = mbSld.abp.reb[locInd2];
		Copy1dArrayInt(&P2A(mbSld.abp.l,locInd1,0,dimMbNuc), 
				&P2A(mbSld.abp.l,locInd2,0,dimMbNuc), dimMbNuc);
	}
	if (recAbpDyn.tgl != 0) {
		abpDyn.cntU[locInd1] = abpDyn.cntU[locInd2];
		abpDyn.cntB[locInd1] = abpDyn.cntB[locInd2];
		abpDyn.cntW[locInd1] = abpDyn.cntW[locInd2];
	}
	abp.id[locInd1] = abp.id[locInd2];
	iAbp[abp.id[locInd2]] = locInd1;
  }
}

void UpdateActinAbpMonomerList(void) {
  int m, n, k, l, posi, rep, tag = 0, idRank, kind, ind, sft, CS;
  int gTglActSN, tglActAD, tglActAbp;
  int oriMC[3], subMC[3], nAdjRank, cnt, cntRecvMsg, cntSendM, cntRecvM;
  int *nMe, *nCp, *pArr, *pArr2, *actAbpMCall, *chkL, *abpKindL;
  int *iFilaPCall, oriIFilaPC; 
  double *volAdj, fac, dimDomC[NDIM];
  ListInt *pM;  

  nAdjRank = (mpiMethod == 0) ? cntAdjRank[0] : 2;
  tglActAD = (actAss.tgl != 0 || actDis.tgl != 0 || actBst.tgl != 0) ? 1 : 0;
  gTglActSN = (actSev.gTgl != 0 || actNuc.gTgl != 0) ? 1 : 0;
  tglActAbp = tglActAD + tglAbpInaMoDyn * 2;

  if (bnd.gTglRnd != 0) {
	volMe = (double)HowManyInOutBound();
  }
  else {
	V3SUB(dimDomC, &P2A(bnd.r,1,0,NDIM), &P2A(bnd.r,0,0,NDIM));
	volMe = V3PROD(dimDomC);
  }

  MALLOC(actAbpMCall,int,nAdjRank * tglActAbp);
  MALLOC(chkL,int,nAdjRank);
  if (tglAbpInaMoDyn != 0) 
  { abpKindL = allIntL; }
  if (gTglActSN != 0) 
  { MALLOC(iFilaPCall,int,nAdjRank); }
  if (updSubdSize.tgl != 0) 
  { MALLOC(volAdj,double,nAdjRank); }

  V3SET_ALL(subMC, 0);
  for(n = 0; n < 3; n++) {
	if (n == 0) { pM = &actM; }
	else if (n == 1) { 
		pM = &acpM; 
	}
	else { 
		pM = &motM;
	}
	for(k = pM->c - 1; k >= 0; k--) {
		if (n == 0) {
			CS = FindElementArray(noActDyn.l, noActDyn.c, pM->l[k], 0, 2);
			if (CS == -1) {
				CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, pM->l[k], 1, 3);
			}
		}
		else {
			CS = FindElementArray(noAbpDyn.l, noAbpDyn.c, pM->l[k], 0, 3);
		}
		CONT(CS == -1);
		subMC[n]++;
	}
  }

  rep = (mpiMethod == 0) ? 1 : NDIM;
  for(m = 0; m < rep; m++) {
	V3SET(oriMC, actM.c, acpM.c, motM.c);
	if (gTglActSN != 0) { oriIFilaPC = iFilaP.c; }
	// gather actM.c, acpM.c, or motM.c of adjacent nodes
	for(n = 0; n < cntAdjRank[m]; n++) {
		idRank = iRank[(mpiMethod == 0) ? n : 2 * m + n];
		posi = 0;
		if (tglActAD != 0) { MPI_PACK_INT(&actM.c, 1, n);	}
		if (tglAbpInaMoDyn != 0) { 
			MPI_PACK_INT(&acpM.c, 1, n); 
			MPI_PACK_INT(&motM.c, 1, n); 
		}
		if (gTglActSN != 0) { MPI_PACK_INT(&iFilaP.c, 1, n); }
		if (updSubdSize.tgl != 0) { MPI_PACK_DBL(&volMe, 1, n); }
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, idRank, 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, idRank,
				tag, MPI_COMM_WORLD, &rReq[n]);
	}
    if (cntAdjRank[m] > 0) {
        cnt = 0;
        cntRecvMsg = 0;
        memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(chkL, 0, sizeof(int) * cntAdjRank[m]);
        while(cntRecvMsg < cntAdjRank[m]) {
            if (mpiTestSendFlag[cnt] == 0) {
                MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
            }
            if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0
                    && chkL[cnt] == 0) {
                MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
                if (mpiTestRecvFlag[cnt] != 0) {
                    posi = 0;
					MPI_UNPACK_INT(&actAbpMCall[tglActAbp * cnt], 
							tglActAbp, cnt);
					if (gTglActSN != 0) 
					{ MPI_UNPACK_INT(&iFilaPCall[cnt], 1, cnt);  }
					if (updSubdSize.tgl != 0) 
					{ MPI_UNPACK_DBL(&volAdj[cnt], 1, cnt); }
                    cntRecvMsg++;
                    chkL[cnt] = 1;
				}
			}
            cnt++;
            if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
	sft = GenRandIntIndex(cntAdjRank[m]);
	// Pack elements from actM.l, acpM.l, or motM.l
	for(ind = 0; ind < cntAdjRank[m]; ind++) {
		n = (ind + sft) % cntAdjRank[m];
		idRank = iRank[(mpiMethod == 0) ? n : 2 * m + n];
		posi = 0;
		for(k = 0; k < 3; k++) {
			CONT((k == 0 && tglActAD == 0) || (k != 0 && tglAbpInaMoDyn == 0));
			if (k == 0) {
				pM = &actM;
				nMe = &nActMe;
				nCp = &nActCp;
				pArr = act.id;
				pArr2 = iAct;
			}
			else {
				if (k == 1) {
					CONT(gTglImpAcpM == 0);
					pM = &acpM;
				}
				else {
					CONT(gTglImpMotM == 0);
					pM = &motM;
				}
				nMe = &nAbpMe;
				nCp = &nAbpCp;
				pArr = abp.id;
				pArr2 = iAbp;
			}
			fac = (updSubdSize.tgl != 0) ? volAdj[n] / volMe : 1.;
			// cntSendM is the number of monomeric elements to send to 
			// adjacent subdomains
			cntSendM = (int)floor((double)(oriMC[k] * fac 
					- actAbpMCall[tglActAbp * n + k - (1 - tglActAD)]) 
					/ (double)nAdjRank);
			if (cntSendM == 0) { cntSendM = 1; }
			if (cntSendM > pM->c - subMC[k]) { cntSendM = pM->c - subMC[k]; }
			if (cntSendM < 0) { cntSendM = 0; }
			MPI_PACK_INT(&cntSendM, 1, n);
			MPI_PACK_INT(pM->l, cntSendM, n);
			if (k == 1) { 
				for(l = 0; l < cntSendM; l++) {
					abpKindL[l] = K_ABP(pArr2[pM->l[l]]);
				}
				MPI_PACK_INT(abpKindL, cntSendM, n);
			}
			// Delete the information of sent elements
			for(l = 0; l < cntSendM; l++) {
				pArr[pArr2[pM->l[l]]] = -1;
				pArr2[pM->l[l]] = -1;
			}
			for(l = *nMe - 1; l >= 0; l--) {  
				CONT(!(pArr[l] < 0));
				MoveParticlesSubroutine2(l, (int)((k + 1) / 2));
			}
			// Shift the information of copied elements to fill the gap
		    for (l = 0; l < *nCp; l++) {
				UpdateActinAbpMonomerListSubroutine(*nMe + l, 
						*nMe + l + cntSendM, k);
			}
			// Delete sent elements in actM.l
			(pM->c) -= cntSendM;
			if (cntSendM > 0) {
				for(l = 0; l < pM->c; l++) {
					pM->l[l] = pM->l[l + cntSendM];
				}
			}
		}
		if (gTglActSN != 0) {
			cntSendM = (int)((double)(oriIFilaPC - iFilaPCall[n]) 
					/ (double)nAdjRank);
			if (cntSendM < 0) { cntSendM = 0; }
			if (cntSendM > iFilaP.c) { cntSendM = iFilaP.c; }
			MPI_PACK_INT(&cntSendM, 1, n);
			MPI_PACK_INT(iFilaP.l, cntSendM, n);
			// Delete sent elements in actM.l
			iFilaP.c -= cntSendM;
			if (cntSendM > 0) {
				for(k = 0; k < iFilaP.c; k++) {
					iFilaP.l[k] = iFilaP.l[k + cntSendM];
				}
			}
		}
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, idRank, 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, idRank,
				tag, MPI_COMM_WORLD, &rReq[n]);
	}
    if (cntAdjRank[m] > 0) {
        cnt = 0;
        cntRecvMsg = 0;
        memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(chkL, 0, sizeof(int) * cntAdjRank[m]);
        while(cntRecvMsg < cntAdjRank[m]) {
            if (mpiTestSendFlag[cnt] == 0) {
                MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
            }
            if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0
                    && chkL[cnt] == 0) {
                MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
                if (mpiTestRecvFlag[cnt] != 0) {
                    posi = 0;
					for(n = 0; n < 3; n++) {
						CONT((n == 0 && tglActAD == 0) 
								|| (n != 0 && tglAbpInaMoDyn == 0));
						if (n == 0) {
							pM = &actM;
							nMe = &nActMe;
							nCp = &nActCp;
						}
						else {
							if (n == 1) {
								CONT(gTglImpAcpM == 0);
								pM = &acpM;
							}
							else {
								CONT(gTglImpMotM == 0);
								pM = &motM;
							}
							nMe = &nAbpMe;
							nCp = &nAbpCp;
						}
						MPI_UNPACK_INT(&cntRecvM, 1, cnt);
						MPI_UNPACK_INT(&pM->l[pM->c], cntRecvM, cnt);
						// Shift back the copied elements to make a room
						if (cntRecvM > 0) {
						    for (k = *nCp - 1; k >= 0; k--) {
								UpdateActinAbpMonomerListSubroutine(*nMe 
										+ cntRecvM + k, *nMe + k, n);
							}
					    }
						if (n == 1) {
							MPI_UNPACK_INT(abpKindL, cntRecvM, cnt);
						}
						// Fill gap by information of transferred elements
						for(k = 0; k < cntRecvM; k++) {
							if (n == 0) { kind = 0; }
							else if (n == 1) { kind = abpKindL[k] + 1; }
							else { kind = 2 + 1; }
 							UpdateActinAbpMonomerListSubroutine2(*nMe + k, 
									pM->l[pM->c + k], kind);
						}
						pM->c += cntRecvM;
						*nMe += cntRecvM;
					}
					if (gTglActSN != 0) {
						MPI_UNPACK_INT(&cntRecvM, 1, cnt);
						MPI_UNPACK_INT(&iFilaP.l[iFilaP.c], cntRecvM, cnt);
						(iFilaP.c) += cntRecvM;
					}
                    cntRecvMsg++;
                    chkL[cnt] = 1;
				}
			}
            cnt++;
            if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
  }
  free(chkL);
  free(actAbpMCall);
  if (gTglActSN != 0) 
  { free(iFilaPCall); }
  if (updSubdSize.tgl != 0) 
  { free(volAdj); }
}

/*-------------------- Balancing counters between subdomains -----------------*/

/*--------------- Collect information from adjacent subdomains ---------------*/
void CollectArrayIntFromAdjacentSubdomain(ListInt *all, int col) {
  int m, n, posi, rep, tag = 0, idRank;
  int cnt, cntRecvMsg, cntElem, *chkList;

  MALLOC(chkList,int,((mpiMethod == 0) ? cntAdjRank[0] : 2));
  rep = (mpiMethod == 0) ? 1 : NDIM;
  for(m = 0; m < rep; m++) {
	for(n = 0; n < cntAdjRank[m]; n++) {
		idRank = iRank[(mpiMethod == 0) ? n : 2 * m + n];
		posi = 0;
		MPI_PACK_INT(&(all->c), 1, n);
		MPI_PACK_INT(all->l, (all->c) * col, n);
		MPI_Isend(bufSendMsg[n], posi, MPI_PACKED, idRank, 
				tag, MPI_COMM_WORLD, &sReq[n]);
		MPI_Irecv(bufRecvMsg[n], sizeBufMsg, MPI_PACKED, idRank,
				tag, MPI_COMM_WORLD, &rReq[n]);
	}
    if (cntAdjRank[m] > 0) {
        cnt = 0;
        cntRecvMsg = 0;
        memset(mpiTestRecvFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(mpiTestSendFlag, 0, sizeof(int) * cntAdjRank[m]);
        memset(chkList, 0, sizeof(int) * cntAdjRank[m]);
        while(cntRecvMsg < cntAdjRank[m]) {
            if (mpiTestSendFlag[cnt] == 0) {
                MPI_Test(&sReq[cnt], &mpiTestSendFlag[cnt], &status);
            }
            if (mpiTestRecvFlag[cnt] == 0 && mpiTestSendFlag[cnt] != 0
                    && chkList[cnt] == 0) {
                MPI_Test(&rReq[cnt], &mpiTestRecvFlag[cnt], &status);
                if (mpiTestRecvFlag[cnt] != 0) {
                    posi = 0;
					MPI_UNPACK_INT(&cntElem, 1, cnt)
					MPI_UNPACK_INT(&P2A(all->l,all->c,0,col), 
							cntElem * col, cnt);
					(all->c) += cntElem;
                    cntRecvMsg++;
                    chkList[cnt] = 1;
				}
			}
            cnt++;
            if (cnt == cntAdjRank[m]) { cnt = 0; }
		}
	}
  }
  free(chkList);
}

void CollectArrayDblFromSubdomainList(double *sendData, double *recvData, 
		int size, ListInt *subdList, int gotData) {
  int n, posi, tag = 0, *flag, cnt, cntRecvMsg, cntMsg, sizBufMsg;
  char *bufMsg, **bufMsg2;
  MPI_Request *req;

  // Gather the local sums to main CPU
  sizBufMsg = sizeof(double) * size;
  if (rank != subdList->l[0] && gotData != 0) {
    MALLOC(bufMsg,char,sizBufMsg);
	posi = 0;
	MPI_Pack(sendData, size, MPI_DOUBLE, bufMsg, sizBufMsg,
			&posi, MPI_COMM_WORLD);
	MPI_Send(bufMsg, posi, MPI_PACKED, subdList->l[0], tag, MPI_COMM_WORLD);
    free(bufMsg);
  }
  else if (rank == subdList->l[0]) {
	cntMsg = (subdList->c) - 1;
	for(n = 0; n < size; n++) {
		recvData[n] = sendData[n];
	}
    MALLOC(flag,int,cntMsg);
    MALLOC(req,MPI_Request,cntMsg);
    MALLOC2(bufMsg2,char,cntMsg);
    for(n = 0; n < cntMsg; n++) {
        MALLOC(bufMsg2[n],char,sizBufMsg);
        MPI_Irecv(bufMsg2[n], sizBufMsg, MPI_PACKED, subdList->l[n + 1], tag,
                MPI_COMM_WORLD, &req[n]);
    }
    cnt = 0;
    cntRecvMsg = 0;
    memset(flag, 0, sizeof(int) * cntMsg);
    while(cntRecvMsg < cntMsg) {
        if (flag[cnt] == 0) {
            MPI_Test(&req[cnt], &flag[cnt], &status);
            if (flag[cnt] != 0) {
                posi = 0;
                MPI_Unpack(bufMsg2[cnt], sizBufMsg, &posi, 
						&recvData[size * (cntRecvMsg + 1)], 
						size, MPI_DOUBLE, MPI_COMM_WORLD);
                cntRecvMsg++;
            }
        }
        cnt++;
        if (cnt == cntMsg) { cnt = 0; }
    }
    for(n = 0; n < cntMsg; n++) {
		free(bufMsg2[n]);
	}
	free(bufMsg2);
	free(req);
	free(flag);
  }
}

/*--------------- Collect information from adjacent subdomains ---------------*/
