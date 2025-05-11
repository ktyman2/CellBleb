// This file contains functions related to boundaries or domain

/*----------------------- Related to force calculation -----------------------*/
void CalcBoundSpringForces(void) {
  int n, ind;
  double f, len, dr[NDIM];

  // Calculate the sum of local forces acting on the ends of clamped actin 
  // filaments
  FOR_ACTME(n) {
	CONT(act.fix[n] < 0);
	ind = (act.fix[n] % 10) - 1;
	if (bndReb.gTglSpr == 0) {
		VV3ADD(&P2(bnd.f,ind,0), &P2(act.f,n,0));
	}
	else {
		len = CalcVecDist(dr, &P2(act.r,n,0), &P2(act.bndRFix,n,0), 0);
		f = SPRING(bnd.stfSpr[ind], len, 0.);
		AddSpringForce(f, len, dr, &P2(act.f,n,0), &P2(bnd.f,ind,0));
	}
  }
}

void CalcBoundRepulsiveForces(void) {
  int n, k, locActInd, *pArr;
  double f, dist, cenDom[NDIM], dr[NDIM];

  FOR_NDIM(k) {
    cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
  }
  FOR_ACTME(n) {
	CONT(ISACTM(n));
    dist = CalcVecDist(dr, cenDom, &P2(act.r,n,0), 0);
    CONT(!(dist > bnd.radRnd));
    f = bnd.stfRep * (dist - bnd.radRnd);
    VVS3ADD(&P2(act.f,n,0), dr, f / dist);
  }
  for(n = 0; n < nAbpMe; n++) {
	CONT(ISABPIM(n));
    pArr = &P2A(abp.ch,n,0,nChAb);
    dist = CalcVecDist(dr, cenDom, &P2(abp.r,n,0), 0);
    CONT(!(dist > bnd.radRnd));
    f = bnd.stfRep * (dist - bnd.radRnd);
    VVS3ADD(&P2(abp.f,n,0), dr, f / dist);
  }
}

/*----------------------- Related to force calculation -----------------------*/

/*---------------------------- Update information ----------------------------*/

void UpdateBoundaryActinUnbindMature(void) {
  int n, ind, side, direc, CS;
  double fMag, fMag2, pUnb, len;

  FOR_ACTME(n) {
	CONT(act.fix[n] < 0);
	ind = (act.fix[n] % 10) - 1;
	direc = ind / 2;
	side = ind % 2;

	if (bndReb.gTglSpr == 1) {
		len = CalcDist(&P2(act.r,n,0), &P2(act.bndRFix,n,0), 0);
		fMag = bnd.stfSpr[ind] * len;
	}
	else {
		fMag = P2(act.f,n,direc) * (1 - side * 2);
	}
	if (bndMat.gTgl != 0) {
		CS = 1;
		if (bndMat.gTglPa == 0) {
			if (P2A(act.ch,n,0,nChAc) > -1 && P2A(act.ch,n,1,nChAc) > -1) {
				CS = 0;
			}
		}
		if (act.nFA[n] < bndMat.maxNFA && CS == 1) {
			fMag2 = TrimDblVal(fMag, 0., (double)(bndMat.maxF - 1));
		    pUnb = bndMat.p[(int)fMag2];
			if (genrand_real3() < pUnb) { 
				act.nFA[n]++; 
				bndMat.cntMe2[ind]++;
				continue;
			}
		}
		fMag /= (double)act.nFA[n];
	}
	fMag2 = TrimDblVal(fMag, 0., (double)(bndUnb.maxF - 1));
    pUnb = P2A(bndUnb.p,ind,(int)fMag2,bndUnb.maxF);
	CONT(!(genrand_real3() < pUnb));
	if (bndMat.gTgl != 0) { 
		act.nFA[n]--; 
		bndMat.cntMe2[ind + 6]++;
		CONT(act.nFA[n] > 0);
	}
	act.fix[n] = -1;
	if (bndMat.gTgl != 0) { act.nFA[n] = 0; }
	(bndUnb.cntMe[ind])++;

  	if (rheoWay > 0) { 
		if (direc == ((bulkRheoType == 0) ? dirNoPBC : dirStr)) {
			DeleteElement1dArray(meaStrePar.l, &meaStrePar.c, act.id[n]);
			DeleteElement1dArray(appStraPar.l, &appStraPar.c, act.id[n]);
		}
	}
  }
}

void UpdateBoundaryActinBinding(void) {
  int n, k, ind, side; 
  double bndLoc, r;

  FOR_ACTME(n) {
	CONT(act.fix[n] > -1);
	CONT(ISACTM(n));
	if (bndReb.gTglPa == 0) {
		CONT(P2A(act.ch,n,0,nChAc) > -1 && P2A(act.ch,n,1,nChAc) > -1);
	}
	FOR_NDIM(k) {
		CONT(bndReb.p[2 * k] == 0. && bndReb.p[2 * k + 1] == 0.);
		CONT(k == dir2D && bndReb.gTglSpr == 0);
		CONT(pbc[k] != 0);
		CONT(iCell[k] != 0 && iCell[k] != nCell[k] - 1);
		r = P2(act.r,n,k);
		if (nCell[k] > 1) {
			if (iCell[k] == 0) {
				CONT(!(r >= rGrid[k][0] && r < rGrid[k][0] + bndReb.dep));
				ind = 2 * k;
			}
			else if (iCell[k] == nCell[k] - 1) {
				CONT(!(r < rGrid[k][nGrid[k] - 1] 
						&& r >= rGrid[k][nGrid[k] - 1] - bndReb.dep));
				ind = 2 * k + 1;
			}
			else { continue; }
		}
		else {
			if (bndReb.dep > dimDomH[k]) {
				CONT(!(r >= rGrid[k][0] && r < rGrid[k][1]));
				if (bndReb.p[2 * k] > 0. && bndReb.p[2 * k + 1] > 0.) {
					ind = 2 * k + (r - rGrid[k][0] < dimDomH[k]) ? 0 : 1;
				}
				else {
					ind = 2 * k + (bndReb.p[2 * k] == 0.) ? 1 : 0;
				}
			}
			else {
				if (r >= rGrid[k][0] && r < rGrid[k][0] + bndReb.dep) {
					ind = 2 * k;
				}
				else if (r < rGrid[k][1] && r >= rGrid[k][1] - bndReb.dep) {
					ind = 2 * k + 1;
				}
				else { continue; }
			}
		}
		CONT(!(genrand_real3() < bndReb.p[ind]));
		act.fix[n] = ind + 1;
		if (bndReb.gTglSpr == 1) {
			V3COPY(&P2(act.bndRFix,n,0), &P2(act.r,n,0));
			P2(act.bndRFix,n,k) = (ind % 2 == 0) ? rGrid[k][0] 
					: rGrid[k][nGrid[k] - 1];
		}
		if (bndMat.gTgl != 0) {
			act.nFA[n] = 1;
		}
		(bndReb.cntMe[ind])++;

		if (rheoWay > 0) {
			if (k == ((bulkRheoType == 0) ? dirNoPBC : dirStr)) {
				InsertElement1dArrayWoChk(meaStreParMe.l, &meaStreParMe.c, 
						act.id[n]);
				InsertElement1dArrayWoChk(appStraParMe.l, &appStraParMe.c, 
						act.id[n]);
			}
		}
		break;
	}
  }
}

// Move boundaries 
void UpdateBoundaryLocation(void) {
  int n, k, direc, side, ind, CS;
  double bndFsum[NDIM * 2 * NDIM], *bndFall;
  double arrGrid[NDIM * 2], prevRGrid[NDIM * 2];
  double area, dR = 0.006, diff, vol, p, *pR;

  if (rank == bndMv.rankF.l[0]) {
	MALLOC(bndFall,double,(bndMv.rankF.c) * NDIM * 2 * NDIM);
  }
  CollectArrayDblFromSubdomainList(bnd.f, bndFall, NDIM * 2 * NDIM, 
		&bndMv.rankF, bnd.gotF);
  if (rank == bndMv.rankF.l[0]) {
	for(k = 0; k < NDIM * 2 * NDIM; k++) { bndFsum[k] = 0.; }
	for(n = 0; n < bndMv.rankF.c; n++) {
		for(k = 0; k < NDIM * 2 * NDIM; k++) {
			bndFsum[k] += P2A(bndFall,n,k,NDIM * 2 * NDIM); 
		}
	}

	if (bndVol.gTgl != 0) {
		vol = V3PROD(dimDom);
		p = -1. * bndVol.stf * (vol - bndVol.eq);
		FOR_NDIM(k) {
			area = dimDom[(k + 1) % NDIM] * dimDom[(k + 2) % NDIM];
			for(n = 0; n < 2; n++) {
				P2(bndFsum,2 * k + n,k) += p * area * (n == 0 ? -1. : 1.);
			}
		}
	}

	// Determine the displacement of each boundary toward an ideal location
	// determined by stress-strain relationship.
	FOR_NDIM(k) {
	    area = dimDom[(k + 1) % NDIM] * dimDom[(k + 2) % NDIM]; 
		for(n = 0; n < 2; n++) {
			ind = (n == 0) ? 0 : nGrid[k] - 1;
			P2A(prevRGrid,n,k,NDIM) = rGrid[k][ind];
			// Pa
			if (bndMv.stfUnit == 0) {
                diff = P2A(rGridInit,n,k,NDIM) + P2A(bndMv.thk,n,k,NDIM)
                        * P2(bndFsum,2 * k + n,k) / area 
						/ P2A(bndMv.stf,n,k,NDIM) - rGrid[k][ind];
			}
			// N/m
			else {
				diff = P2(bndFsum,2 * k + n,k) / P2A(bndMv.stf,n,k,NDIM) 
						- (rGrid[k][ind] - P2A(rGridInit,n,k,NDIM));
			}
			rGrid[k][ind] += (fabs(diff) > dR) ? SignDbl(diff) * dR : diff;
			P2A(arrGrid,n,k,NDIM) = rGrid[k][ind];
		}
	}
	free(bndFall);
  }
  MPI_Bcast(arrGrid, NDIM * 2, MPI_DOUBLE, bndMv.rankF.l[0], MPI_COMM_WORLD);

  // Post-process
  if (rank != bndMv.rankF.l[0]) {
	FOR_NDIM(k) {
		for(n = 0; n < 2; n++) { 
			ind = (n == 0) ? 0 : nGrid[k] - 1;
			P2A(prevRGrid,n,k,NDIM) = rGrid[k][ind];
			rGrid[k][ind] = P2A(arrGrid,n,k,NDIM);
		}
	}
  }
  CS = 0;
  FOR_NDIM(k) {
	if (iCell[k] == 0) { 
		P2A(bnd.r,0,k,NDIM) = rGrid[k][0]; 
		CS = 1;
	}
	if (iCell[k] == nCell[k] - 1) {
		P2A(bnd.r,1,k,NDIM) = rGrid[k][nGrid[k] - 1]; 
		CS = 1;
	}
 	dimDom[k] = rGrid[k][nGrid[k] - 1] - rGrid[k][0];
  	dimDomH[k] = dimDom[k] * 0.5;
  }
  if (CS != 0)  {
	minDimDomC = POS_LARGE_VALUE;
	FOR_NDIM(k) {
		CONT(k == dir2D);
		CONT(!(P2A(bnd.r,1,k,NDIM) - P2A(bnd.r,0,k,NDIM) < minDimDomC));
		minDimDomC = P2A(bnd.r,1,k,NDIM) - P2A(bnd.r,0,k,NDIM);
	}
  }
  if (bnd.gotF != 0) {
	// Clamped actin filaments also have to be displaced with boundaries.
	FOR_ACTME(n) {
		CONT(act.fix[n] < 0);
		pR = (bndReb.gTglSpr == 1) ? &P2(act.bndRFix,n,0) : &P2(act.r,n,0);
		direc = ((act.fix[n] % 10) - 1) / 2;
		side = ((act.fix[n] % 10) - 1) % 2;
		pR[direc] += rGrid[direc][(side == 0) ? 0 : 
					nGrid[direc] - 1] - P2A(prevRGrid,side,direc,NDIM);
		FOR_NDIM(k) {
			CONT(k == direc);
			pR[k] = (pR[k] - P2A(prevRGrid,0,k,NDIM)) 
					/ (P2A(prevRGrid,1,k,NDIM) - P2A(prevRGrid,0,k,NDIM)) 
					* dimDom[k] + rGrid[k][0];
		}
  	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		// Adjust the position of fixed points for membranes.
		FOR_MBME(n) {
			CONT(ISNUC(memb.idx[n]));
			CONT(memb.fix[n] < 0);
			side = (P2(memb.rFix,n,2) - rGrid[2][0] < dimDomH[2]) ? 0 : 1;
			P2(memb.rFix,n,2) += (rGrid[2][side * (nGrid[2] - 1)] 
					- P2A(prevRGrid,side,2,NDIM));
		}
	}
  }
}

/*---------------------------- Update information ----------------------------*/

/*------------------------ Handle boundary conditions ------------------------*/

// Check whether the chain crosses a boundary or not.
void CheckCrossBound(int *sft, double *dr, double *dr2) {
  int k;
 
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	CONT(!(fabs(dr2[k] - dr[k]) > 1.5 * dimDomH[k]));
	sft[k] += (dr2[k] > dr[k]) ? -1 : 1;
  }
}

// Apply periodic boundary condition (PBC) to a chain vector.
// If the absolute value of a component of the vector is greater than half of 
// domain width with PBC, it is added or subtracted by domain width.
void ApplyBoundCondVecDiff(double *dr) {
  int k;
  FOR_NDIM(k) {
	CONT(pbc[k] != 1);
	if (dr[k] >= dimDomH[k]) { dr[k] -= dimDom[k]; }
	else if (dr[k] <= REVSIGN(dimDomH[k])) { dr[k] += dimDom[k]; }
  }
}

// If mode is -1, neglect repulsive condition.
// Otherwise, if mode is 0, actin. If mode is 1-3, ABP.
// If 5, nucleus. If 4, membrane
void ApplyBoundCondVector(double *r, int mode, int mode2) {
  int k, indBndF;
  double disp, dispAll, repF, drag, sft;

  FOR_NDIM(k) {
	disp = (k == dirStr && rheoWay > 0 && bulkRheoType == 0) 
			? r[dirNoPBC] * stra.acc : 0.;
	// If PBC
    if (pbc[k] == 1) {
		if (r[k] < disp + rGrid[k][0])
		{ r[k] += dimDom[k]; }
		else if (r[k] >= disp + rGrid[k][nGrid[k] - 1])
		{ r[k] -= dimDom[k]; }
	}
	// If no PBC and mode >= 0, repulsive condition is applied.
	else if (pbc[k] != 1 && mode >= 0) {
		sft = 0.;
		if (mode == 0) { drag = 1.; }
		else if (mode >= 1 && mode <= 3) { drag = abpF.drag[mode - 1].inv; }
		else if (mode == 4) { drag = INV(memb.dragR); }
		else if (mode == 5) { drag = INV(memb.nucDragR); }
		else { 
			drag = INV(bead.drag[mode - 6]); 
			sft = bead.rad[mode - 6];
		}
		if ((r[k] < disp + rGrid[k][0] + sft) 
				|| (r[k] >= disp + rGrid[k][nGrid[k] - 1] - sft)) {
			if (r[k] < disp + rGrid[k][0] + sft) {
				dispAll = disp + rGrid[k][0] + sft;
				indBndF = 0;
			}
			else {
				dispAll = disp + rGrid[k][nGrid[k] - 1] - sft;
				indBndF = 1;
			}
			repF = bnd.stfRep * (r[k] - dispAll);
			// If shear-bulk-rheology (against a tilted boundary)
			if (k == dirStr && rheoWay > 0 && bulkRheoType == 0) {
				r[dirStr] -= repF / (SQR(stra.acc) + 1.) * drag * dt;
				r[dirNoPBC] -= (repF  * stra.acc / (SQR(stra.acc) + 1.)) 
						* drag * dt;
				if (mode2 != 0) {
					P2A(bnd.f,2 * dirStr + indBndF,dirStr,NDIM) += repF 
							/ (SQR(stra.acc) + 1.);
					P2A(bnd.f,2 * dirNoPBC + indBndF,dirNoPBC,NDIM) 
							+= repF * stra.acc / (SQR(stra.acc) + 1.);
				}
			}
			// If a flat boundary
			else {
				r[k] -= repF * drag * dt;
				if (bndMv.gTgl != 0 && mode2 != 0) {
					P2A(bnd.f,2 * k + indBndF,k,NDIM) += repF;
				}
			}
		}
	}
  }
}

// Apply the periodic boundary condition on the positions of whole things
void ApplyBoundCondAll(void) {
  int n;

  FOR_ACTME(n) {
	CONT(ISACTM(n));
	ApplyBoundCondVector(&P2(act.r,n,0), 0, (act.fix[n] > -1 ? 0 : 1));
  }
 
  FOR_ABPME(n) {
	CONT(ISABPIM(n));
	ApplyBoundCondVector(&P2(abp.r,n,0), K_ABP(n) + 1, 1);
  }
  if (gTglMb != 0) {
	FOR_MBME(n) {
		ApplyBoundCondVector(&P2(memb.r,n,0), 
				((ISNUC(memb.idx[n])) ? 5 : 4), 0);
	}
  }
}

// Regardless of shear deformation, offset positions of all particles into 
// retangular-solid domain
// If mode = 0, adjust the position only in shearing direction
// If mode > 0, apply PBC in other directions.
void ConvertRectDomainVector(double *r, int mode) {
  int k;

  if (mode == 0) {
	if (rheoWay > 0) {
		if (r[dirStr] < rGrid[dirStr][0]) { 
			while(r[dirStr] < rGrid[dirStr][0]) { 
				r[dirStr] += dimDom[dirStr];
			}
		}
		else if (r[dirStr] >= rGrid[dirStr][nGrid[dirStr] - 1]) {
			while(r[dirStr] >= rGrid[dirStr][nGrid[dirStr] - 1]) {
				r[dirStr] -= dimDom[dirStr];
			}
		}
	}
  }
  else {
	FOR_NDIM(k) {
		if (pbc[k] != 0 || (pbc[k] == 0 && confPbc[k] != 0)) {
			if (r[k] >= rGrid[k][nGrid[k] - 1]) { r[k] -= dimDom[k]; }
			else if (r[k] < rGrid[k][0]) { r[k] += dimDom[k]; }
		}
	}
  }
}

/*------------------------ Handle boundary conditions ------------------------*/

/*---------------- Interactions between boundaries and others ----------------*/

int CheckParticleInDomain(double *r) {
  int CS, k;
  double disp;
  CS = 1;
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 0));
	disp = (k == dirStr && rheoWay > 0 && bulkRheoType == 0) 
			? r[dirNoPBC] * stra.acc : 0.; 
	CONT(!(r[k] < disp + rGrid[k][0] || r[k] >= disp + rGrid[k][nGrid[k] - 1]));
	CS = 0;
	break;
  }
  return CS;
}

double HowManyInOutBound(void) {
  int k, cnt, dirX, dirY;
  int ind[NDIM], begin[NDIM], end[NDIM];
  double dist, resol, r[NDIM], cenDom[NDIM];

  resol = 0.5;
  FOR_NDIM(k) {
    cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
  }

  FOR_NDIM(k) {
    begin[k] = (int)ceil(P2A(bnd.r,0,k,NDIM) / resol);
    end[k] = (int)floor(P2A(bnd.r,1,k,NDIM) / resol);
  }

  if (dir2D == -1) {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			for(ind[2] = begin[2]; ind[2] < end[2]; ind[2]++) {
				VS3COPY(r, ind, resol);
				dist = CalcDist(cenDom, r, 0);
				if (dist <= bnd.radRnd) {
					cnt++;
				}
			}
		}
	}
	return (int)(cnt * CUBE(resol));
  }
  else {
	dirX = (dir2D + 1) % NDIM;
	dirY = (dir2D + 2) % NDIM;
	cnt = 0;
	for(ind[dirX] = begin[dirX]; ind[dirX] < end[dirX]; ind[dirX]++) {
		for(ind[dirY] = begin[dirY]; ind[dirY] < end[dirY]; ind[dirY]++) {
			r[dirX] = ind[dirX] * resol;
			r[dirY] = ind[dirY] * resol;
			r[dir2D] = cenDom[dir2D];
			dist = CalcDist(cenDom, r, 0);
			if (dist <= bnd.radRnd) {
				cnt++;
			}
		}
	}
	return (int)(cnt * SQR(resol));
  }
}

int CheckActinAbpOverlapBoundary(double *r) {
  int CS, k;
  double dist, cenDom[NDIM];

  FOR_NDIM(k) {
    cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
  }
  dist = CalcDist(cenDom, r, 0);
  CS = (dist > bnd.radRnd) ? 0 : 1;
  return CS;
}

/*---------------- Interactions between boundaries and others ----------------*/

/*-------------------------- Recording information ---------------------------*/

void RecordBoundaryLocation(void) {
  int k;
  FILE *fOut;

  // Record the locations of boundaries at the interval of "recBndLoc.prd".
  fOut = fopen(GenFileName("BndLoc"), "a");
  fprintf(fOut, "%lld\t", currTimeStep);
  FOR_NDIM(k) {
	fprintf(fOut, "%g\t%g\t", rGrid[k][0], rGrid[k][nGrid[k] - 1]);
  }
  fprintf(fOut, "\n");
  fclose(fOut);
}

void RecordBoundaryActinUnbindBind(void) {
  int n, k, *cntBndUnbRebAll, *cntBndUnbRebSum;
  FILE *fOut;

  MALLOC(cntBndUnbRebAll,int,nCpu*NDIM*4);
  MALLOC(cntBndUnbRebSum,int,NDIM*4);
  V6COPY(cntBndUnbRebSum, bndUnb.cntMe);
  V6COPY(&cntBndUnbRebSum[NDIM*2], bndReb.cntMe);

  MPI_Gather(cntBndUnbRebSum, NDIM * 4, MPI_INT, cntBndUnbRebAll, NDIM * 4, 
		MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 0; n < NDIM * 4; n++) {
		for(k = 1; k < nCpu; k++) {
			cntBndUnbRebSum[n] += P2A(cntBndUnbRebAll,k,n,NDIM*4);
		}
	}
	fOut = fopen(GenFileName("BndUnbReb"), "a");
	fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayInt(fOut, cntBndUnbRebSum, NDIM * 4, 0);
	fclose(fOut);
  }
  free(cntBndUnbRebSum);
  free(cntBndUnbRebAll);
}

void RecordBoundaryActinMature(void) {
  int n, k, *cntAll;
  FILE *fOut;

  MALLOC(cntAll,int,nCpu*NDIM*4);
  MPI_Gather(bndMat.cntMe2, NDIM * 4, MPI_INT, cntAll, NDIM * 4, 
		MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 0; n < NDIM * 4; n++) {
		for(k = 1; k < nCpu; k++) {
			cntAll[n] += P2A(cntAll,k,n,NDIM * 4);
		}
	}
	fOut = fopen(GenFileName("BndActMat"), "a");
	fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayInt(fOut, cntAll, NDIM * 4, 0);
	fclose(fOut);
  }
  free(cntAll);
}

void RecordBoundaryTractionForce(void) {
  int n, k, ind, *sendCall;
  double dr[NDIM], f[NDIM];
  FILE *fOut;
  ListDbl send, recv;

  MALLOC(send.l,double,nActMe*7);
  if (rank == 0) { MALLOC(recv.l,double,nAct*7); }
  send.c = 0;
  FOR_ACTME(n) {
	CONT(act.fix[n] < 0);
	ind = (act.fix[n] % 10) - 1;
	CalcVec(dr, &P2(act.r,n,0), &P2(act.bndRFix,n,0));
	FOR_NDIM(k) {
		f[k] = F_S2PN(dr[k] * bnd.stfSpr[ind]);
	}
	P2A(send.l,send.c,0,7) = ind;
	V3COPY(&P2A(send.l,send.c,1,7), &P2(act.bndRFix,n,0));
	V3COPY(&P2A(send.l,send.c,4,7), f);
	send.c++;
  }
  MALLOC(sendCall,int,nCpu);
  MPI_Gather(&send.c, 1, MPI_INT, sendCall, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDoubleWoIndWoSort(send.c, sendCall, send.l, recv.l, 7);
  recv.c = SumArrInt(sendCall, nCpu);
  if (rank == 0) {
	fOut = fopen(GenFileName("BndTracF"), "a");
	fprintf(fOut, "%d\t0\t0\t0\t0\t0\t0\t0\n", recv.c);
	Fprintf2dArrayDouble(fOut, recv.l, recv.c, 7, 0, 0);
	fclose(fOut);
  }

  free(sendCall);
  free(send.l);
  if (rank == 0) { free(recv.l); }
}

/*-------------------------- Recording information ---------------------------*/
