// This file contains functions checking and recording errors.

// Record the kind of an error in the file, "Error".
void RecordError(int kind) {
  FILE *fErr;
  char fn[80];

  sprintf(fn, "Error_%d", rank);
  fErr = fopen(GenFileName(fn), "a");
  fprintf(fErr, "\ncurrTimeStep : %lld\n", currTimeStep);
  fprintf(fErr, "rank : %d\n", rank);

  switch (kind) {
	case 1 :
		fprintf(fErr, "Repulsive force between actins\n");
		break;  
	case 2 :
		fprintf(fErr, "Repulsive between ABPs\n");
		break;  
	case 3 :
		fprintf(fErr, "Repulsive force between actin and ABP\n");
		break;
	case 4 :
		fprintf(fErr, "Spring force between actins\n");
		break;
	case 5 :
		fprintf(fErr, "Spring between actin and ABP\n");
		break;
	case 6 :
		fprintf(fErr, "Bending force for actin filaments\n");
		break;
	case 7 :
		fprintf(fErr, "Bending force for actin-ABP-actin\n");
		break;  
	case 8 :
		fprintf(fErr, "Bending force for right angle between actin and ABP\n");
		break;  
	case 9 :
		fprintf(fErr, "Force for cross-linking angle\n");
		break;  
	case 10 :
		fprintf(fErr, "Total force on actin\n");
		break;
	case 11 :
		fprintf(fErr, "Total force on ABP\n");
		break;
	case 12 :
		fprintf(fErr, "Spring force on the backbone of motors\n");
		break;
	case 13 :
		fprintf(fErr, "Bending force on the backbone of motors\n");
		break;
	case 14 :
		fprintf(fErr, "Spring force between membranes\n");
		break;
	case 15 :
		fprintf(fErr, "Spring force between membrane and actin\n");
		break;
	case 16 :
		fprintf(fErr, "Spring force between membrane and space\n");
		break;
	case 17 :
		fprintf(fErr, "Bending force for membranes\n");
		break;
	case 18 :
		fprintf(fErr, "Total force on membrane\n");
		break;
  }
  fclose(fErr);
}

// Subroutine for functions recording errors
void RecordErrorArrayInt(FILE *fOut, const char *name, int ind, 
		int *arr, int cnt) {
  fprintf(fOut, "%s of %d: ", name, ind);
  Fprintf1dArrayInt(fOut, arr, cnt, 0);
}

// Subroutine for functions recording errors
void RecordErrorArrayDouble(FILE *fOut, const char *name, int ind, double *arr, 
		int cnt) {
  fprintf(fOut, "%s of %d: ", name, ind);
  Fprintf1dArrayDouble(fOut, arr, cnt, 0);
}

void RecordErrorSubroutine1(FILE *fErr, int ind) {
  if (iAct[ind] > -1) {
	RecordErrorArrayInt(fErr, "act.ch", ind, &P2A(act.ch,iAct[ind],0,nChAc), 
			nChAc);
	RecordErrorArrayDouble(fErr, "act.r", ind, &P2(act.r,iAct[ind],0), NDIM);
  }
  else {
	fprintf(fErr, "Actin %d doesn't exist in the subdomain!\n", ind);
  }
  fprintf(fErr, "\n");
}

void RecordErrorSubroutine2(FILE *fErr, int ind) {
  if (iAbp[ind] > -1) {
	RecordErrorArrayInt(fErr, "abp.ch", ind, &P2A(abp.ch,iAbp[ind],0,nChAb), 
			nChAb);
	RecordErrorArrayDouble(fErr, "abp.r", ind, &P2(abp.r,iAbp[ind],0), NDIM);
  }
  else {
	fprintf(fErr, "ABP %d doesn't exist in the subdomain!\n", ind);
  }
  fprintf(fErr, "\n");
}

void RecordErrorSubroutine3(FILE *fErr, int ind) {
  if (iMb[ind] > -1) {
	RecordErrorArrayInt(fErr, "memb.ch", ind, &P2A(memb.ch,iMb[ind],0,nChMb), 
			nChMb);
	RecordErrorArrayDouble(fErr, "memb.r", ind, &P2(memb.r,iMb[ind],0), NDIM);
  }
  else {
	fprintf(fErr, "Membrane %d doesn't exist in the subdomain!\n", ind);
  }
  fprintf(fErr, "\n");
}

// Record detailed information related to an error caused by spring force.
void RecordErrorSpringForce(int ind1, int ind2, double f, double len, 
		int mode) {
  char fn[80];
  FILE *fErr;

  sprintf(fn, "Error_%d", rank);
  fErr = fopen(GenFileName(fn), "a");
  // A spring force between actins
  if (mode == 0) {
	fprintf(fErr, "Spring: actin %d (%d/%d) and actin %d (%d/%d)\n\n",
			ind2, iAct[ind2], nActMe, ind1, iAct[ind1], nActMe);
	RecordErrorSubroutine1(fErr, ind1);
	RecordErrorSubroutine1(fErr, ind2);
  }
  // A spring force between actin and ABP
  else if (mode == 1) {
	fprintf(fErr, "Spring: ABP %d (%d/%d) and actin %d (%d/%d)\n\n",
			ind2, iAbp[ind2], nAbpMe, ind1, iAct[ind1], nActMe);
	RecordErrorSubroutine2(fErr, ind2);
	RecordErrorSubroutine1(fErr, ind1);
	RecordErrorSubroutine1(fErr, P2A(act.ch,iAct[ind1],0,nChAc));
  } 
  // A spring force between ABPs
  else if (mode == 2) {
	fprintf(fErr, "Spring: ABP %d (%d/%d) and ABP %d (%d/%d)\n\n",
			ind2, iAbp[ind2], nAbpMe, ind1, iAbp[ind1], nAbpMe);
	RecordErrorSubroutine2(fErr, ind2);
	RecordErrorSubroutine2(fErr, ind1);
  }
  // A spring force between membranes
  else if (mode == 3) {
	fprintf(fErr, "Spring: membrane %d (%d/%d) and membrane %d (%d/%d)\n\n",
			ind1, iMb[ind1], nMbMe, ind2, iMb[ind2], nMbMe);
	RecordErrorSubroutine3(fErr, ind1);
	RecordErrorSubroutine3(fErr, ind2);
  }
  // A spring force between membrane and actin
  else if (mode == 4) {
	fprintf(fErr, "Spring: membrane %d (%d/%d) and actin %d (%d/%d)\n\n",
			ind1, iMb[ind1], nMbMe, ind2, iAct[ind2], nActMe);
	RecordErrorSubroutine3(fErr, ind1);
	RecordErrorSubroutine1(fErr, ind2);
  }
  // A spring force between membrane and space
  else {
	fprintf(fErr, "Spring: membrane %d (%d/%d) and space\n\n",
			ind1, iMb[ind1], nMbMe);
	RecordErrorSubroutine3(fErr, ind1);
  }
  fprintf(fErr, "Force = %g\n", f);
  fprintf(fErr, "Length = %g\n", len);
  fclose(fErr);
}

// Record detailed information related to an error caused by repulsive force.
void RecordErrorRepForce(int *ind, int *nlPnt, double (*rPnt)[NDIM], 
	  double f, double dist, double *ratio, int mode) {
  int n, k, neighInd[2], *pArr[2];
  char fn[80];
  FILE *fErr;

  sprintf(fn, "Error_%d", rank);
  fErr = fopen(GenFileName(fn), "a");
  if (mode < 3) {
	fprintf(fErr, "Repulsive: ");
	for(n = 0; n < 2; n++) {
		if (nlPnt[n] < nAct) {
			fprintf(fErr, "actin segment %d", nlPnt[n]);
			pArr[n] = &P2A(act.cyl.l,nlPnt[n],0,2);
		}
		else {
			fprintf(fErr, "ABP %d", nlPnt[n] - nAct);
		}
		if (n == 0) {
			fprintf(fErr, " and ");
		}
		else {
			fprintf(fErr, "\n\n");
		}
	}
	for(n = 0; n < 2; n++) {
		if (nlPnt[n] < nAct) {
			fprintf(fErr, "actin segment %d: %d (%d/%d) %d (%d/%d)\n\n", 
					nlPnt[n], ((ind[2 * n] > -1) ? act.id[ind[2 * n]] : -1), 
					ind[2 * n], nActMe, ((ind[2 * n + 1] > -1) ? 
					act.id[ind[2 * n + 1]] : -1), ind[2 * n + 1], nActMe);
			for(k = 0; k < 2; k++) {
				RecordErrorSubroutine1(fErr, act.id[ind[2 * n + k]]);
			}
			fprintf(fErr, "Ratio = %g\n", ratio[n]);
		}
		else {
			fprintf(fErr, "ABP: %d (%d/%d)\n", nlPnt[n] - nAct, 
					ind[dimMbNuc * n], nAbpMe);
			RecordErrorSubroutine2(fErr, nlPnt[n] - nAct);
		}
	}
  }
  else {
	if (nlPnt[0] < nAct + nAbp + nUnitMb) {
		neighInd[0] = nlPnt[0] - nAct - nAbp;
		fprintf(fErr, "Repulsive: membrane segment %d and ", neighInd[0]);
		pArr[0] = &P2A(memb.unit.l,neighInd[0],0,dimMbNuc);
	}
	else {
		neighInd[0] = nlPnt[0] - nAct - nAbp - nUnitMb;
		fprintf(fErr, "Repulsive: copied membrane segment %d and ", 
				neighInd[0]);
		pArr[0] = &P2A(memb.unitCp.l,neighInd[0],0,dimMbNuc);
	}
	if (mode == 4) {
		fprintf(fErr, "actin segment %d\n\n", nlPnt[1]);
		pArr[1] = &P2A(act.cyl.l,nlPnt[1],0,2);
	}
	else if (mode == 3) {
		neighInd[1] = nlPnt[1] - nAct;
		fprintf(fErr, "ABP %d\n\n", neighInd[1]);
	}
	else {
		if (nlPnt[1] < nAct + nAbp + nUnitMb) {
			neighInd[1] = nlPnt[1] - nAct - nAbp;
			fprintf(fErr, "membrane segment %d\n", neighInd[1]);
			pArr[1] = &P2A(memb.unit.l,neighInd[1],0,dimMbNuc);
		}
		else {
			neighInd[1] = nlPnt[1] - nAct - nAbp - nUnitMb;
			fprintf(fErr, "copied membrane segment %d\n", neighInd[1]);
			pArr[1] = &P2A(memb.unitCp.l,neighInd[1],0,dimMbNuc);
		}
	}
	fprintf(fErr, "membrane segment %d: ", neighInd[0]);
	for(k = 0; k < dimMbNuc; k++) {
		fprintf(fErr, "%d (%d/%d) ", pArr[0][k], ind[k], nMbMe);
	}
	fprintf(fErr, "\n");
	for(k = 0; k < dimMbNuc; k++) {
		RecordErrorSubroutine3(fErr, pArr[0][k]);
	}
	if (mode == 4) {
		fprintf(fErr, "actin segment %d: ", nlPnt[1]);
		for(k = 0; k < 2; k++) {
			fprintf(fErr, "%d (%d/%d) ", pArr[1][k], ind[dimMbNuc + k], nActMe);
		}
		fprintf(fErr, "\n");
		for(k = 0; k < 2; k++) {
			RecordErrorSubroutine1(fErr, pArr[1][k]);
		}
	}
	else if (mode == 3) {
		fprintf(fErr, "ABP: %d (%d/%d)\n", nlPnt[1] - nAct, ind[dimMbNuc], 
				nAbpMe);
		RecordErrorSubroutine2(fErr, nlPnt[1] - nAct);
	}
	else {
		fprintf(fErr, "membrane segment %d: ", neighInd[1]);
		for(k = 0; k < dimMbNuc; k++) {
			fprintf(fErr, "%d (%d/%d) ", pArr[1][k], ind[dimMbNuc + k], nMbMe);
		}
		fprintf(fErr, "\n");
		for(k = 0; k < dimMbNuc; k++) {
			RecordErrorSubroutine3(fErr, pArr[1][k]);
		}
	}
	if (mode == 3 && dimMbNuc == 2) {
		fprintf(fErr, "Ratio = %g\n", ratio[0]);
	}
	else if (((mode == 4 || mode == 5) && dimMbNuc == 2) 
			|| (mode == 3 && dimMbNuc == 3)) {
		fprintf(fErr, "Ratio = %g\t%g\n", ratio[0], ratio[1]);
	}
	else if (mode == 4 && dimMbNuc == 3) {
		fprintf(fErr, "Ratio = %g\t%g\t%g\n", ratio[0], ratio[1], ratio[2]);
	}
	else {
		fprintf(fErr, "Ratio = %g\t%g\t%g\t%g\n", ratio[0], ratio[1], ratio[2], 
				ratio[3]);
	}
  }
  fprintf(fErr, "Force = %g\n", f);
  fprintf(fErr, "Distance = %g\n", dist);
  fclose(fErr);
}

// Record detailed information related to an error caused by bending force.
void RecordErrorBendingForce(int ind1, int ind2, int ind3, double f, int mode) {
  char fn[80];
  FILE *fErr;

  sprintf(fn, "Error_%d", rank);
  fErr = fopen(GenFileName(fn), "a");

  switch (mode) {
    case 0:
		fprintf(fErr, "Filament bending: actin %d (%d/%d), actin "
				"%d (%d/%d), actin %d (%d/%d)\n\n", ind1, iAct[ind1], nActMe, 
				ind2, iAct[ind2], nActMe, ind3, iAct[ind3], nActMe);    
		break;
    case 1:
		fprintf(fErr, "ABP bending: actin %d (%d/%d), ABP %d (%d/%d), "
				"actin %d (%d/%d)\n\n",	ind2, iAct[ind2], nActMe, 
				ind1, iAbp[ind1], nAbpMe, ind3, iAct[ind3], nActMe);  
		break;
    case 2:
		fprintf(fErr, "Crosslink angle: ABP %d (%d/%d), actin %d (%d/%d), "
				"actin %d (%d/%d)\n\n", ind1, iAbp[ind1], nAbpMe, 
				ind2, iAct[ind2], nActMe, ind3, iAct[ind3], nActMe);  
		break;
    case 3:
		fprintf(fErr, "90-angle bending: ABP %d (%d/%d), actin %d (%d/%d), "
				"actin %d (%d/%d)\n\n", ind1, iAbp[ind1], nAbpMe, 
				ind2, iAct[ind2], nActMe, ind3, iAct[ind3], nActMe);  
		break;
	case 4:
		fprintf(fErr, "Motor backbone bending: ABP %d (%d/%d), ABP %d (%d/%d), "
				"ABP %d (%d/%d)\n\n", ind2, iAbp[ind2], nAbpMe, 
				ind1, iAbp[ind1], nAbpMe, ind3, iAbp[ind3], nAbpMe);  
		break;
	case 5:
		fprintf(fErr, "Membrane bending: membrane %d (%d/%d), membrane %d "
				"(%d/%d), membrane %d (%d/%d)\n\n", ind1, iMb[ind1], nMbMe, 
				ind2, iMb[ind2], nMbMe, ind3, iMb[ind3], nMbMe);  
		break;
  }
  if (mode < 4) {
	if (mode == 0) {
		RecordErrorSubroutine1(fErr, ind1);
	}
	else {
		RecordErrorSubroutine2(fErr, ind1);
	}
	RecordErrorSubroutine1(fErr, ind2);
	RecordErrorSubroutine1(fErr, ind3);
  }
  else if (mode == 4) {
	RecordErrorSubroutine2(fErr, ind1);
	RecordErrorSubroutine2(fErr, ind2);
	RecordErrorSubroutine2(fErr, ind3);
  }
  else {
	RecordErrorSubroutine3(fErr, ind1);
	RecordErrorSubroutine3(fErr, ind2);
	RecordErrorSubroutine3(fErr, ind3);
  }
  fprintf(fErr, "Force = %g\n", f);
  fclose(fErr);
}

void RecordErrorTotalForceSubroutine1(FILE *fErr, int ind1) {
  int ind2;

  if (ind1 > -1) {
	ind2 = iAct[ind1];
	if (ind2 > -1) {
		fprintf(fErr, "actin %d (%d/%d):\n", ind1, ind2, nActMe);	
    	RecordErrorArrayInt(fErr, "act.ch", ind1, &P2A(act.ch,ind2,0,nChAc), 
				nChAc);
		RecordErrorArrayDouble(fErr, "act.r", ind1, &P2(act.r,ind2,0), NDIM);
		RecordErrorArrayDouble(fErr, "act.f", ind1, &P2(act.f,ind2,0), NDIM);
		fprintf(fErr, "fixAct = %d\n", act.fix[ind2]);
	}
	else {
		fprintf(fErr, "Actin %d doesn't exist in the subdomain!\n", ind1);
	}
	fprintf(fErr, "\n");
  }
}

void RecordErrorTotalForceSubroutine2(FILE *fErr, int ind1) {
  int ind2;

  if (ind1 > -1) {
	ind2 = iAbp[ind1];
	fprintf(fErr, "ABP %d (%d/%d):\n", ind1, ind2, nAbpMe);
	if (ind2 > -1) {
    	RecordErrorArrayInt(fErr, "abp.ch", ind1, &P2A(abp.ch,ind2,0,nChAb), 
				nChAb);
		RecordErrorArrayDouble(fErr, "abp.r", ind1, &P2(abp.r,ind2,0), NDIM);
		RecordErrorArrayDouble(fErr, "abp.f", ind1, &P2(abp.f,ind2,0), NDIM);
	}
	else {
		fprintf(fErr, "ABP %d doesn't exist in the subdomain!\n", ind1);
	}
	fprintf(fErr, "\n");
  }
}

void RecordErrorTotalForceSubroutine3(FILE *fErr, int ind1) {
  int ind2;

  if (ind1 > -1) {
	ind2 = iMb[ind1];
	fprintf(fErr, "membrane %d (%d/%d):\n", ind1, ind2, nMbMe);
	if (ind2 > -1) {
    	RecordErrorArrayInt(fErr, "memb.ch", ind1, &P2A(memb.ch,ind2,0,nChMb), 
				nChMb);
		RecordErrorArrayDouble(fErr, "memb.r", ind1, &P2(memb.r,ind2,0), NDIM);
		RecordErrorArrayDouble(fErr, "memb.f", ind1, &P2(memb.f,ind2,0), NDIM);
	}
	else {
		fprintf(fErr, "Membrane %d doesn't exist in the subdomain!\n", ind1);
	}
	fprintf(fErr, "\n");
  }
}

// Record detailed information related to an error caused by total force.
void RecordErrorTotalForce(int ind, int mode) {
  int ind2, k; 
  FILE *fErr;
  char fn[80];

  sprintf(fn, "Error_%d", rank);
  fErr = fopen(GenFileName(fn), "a");

  fprintf(fErr, "Total force: ");
  if (mode == 0) {
	RecordErrorTotalForceSubroutine1(fErr, ind);
    fprintf(fErr, "\n---- Bound on the actin ----\n");
    for(k = 0; k < 2; k++) {
        ind2 = P2A(act.ch,iAct[ind],k,nChAc);
		RecordErrorTotalForceSubroutine1(fErr, ind2);
    }
    for(k = 2; k < nChAc; k++) {
        ind2 = P2A(act.ch,iAct[ind],k,nChAc);
		RecordErrorTotalForceSubroutine2(fErr, ind2);
    }
  }
  else if (mode == 1) {
	RecordErrorTotalForceSubroutine2(fErr, ind);
    fprintf(fErr, "\n---- Bound on the ABP ----\n");
    for(k = 0; k < nChAb; k++) {
		CONT(k == 2);
        ind2 = P2A(abp.ch,iAbp[ind],k,nChAb);
		if (k < 2) { 
			RecordErrorTotalForceSubroutine1(fErr, ind2);
		}
		else {
			RecordErrorTotalForceSubroutine2(fErr, ind2);
		}
    }
  }
  else {
	RecordErrorTotalForceSubroutine3(fErr, ind);
    fprintf(fErr, "\n---- Bound on the membrane ----\n");
    for(k = 0; k < nChMb; k++) {
        ind2 = P2A(memb.ch,iMb[ind],k,nChMb);
		if (k < nChMb - nMbAct) { 
			RecordErrorTotalForceSubroutine3(fErr, ind2);
		}
		else {
			RecordErrorTotalForceSubroutine1(fErr, ind2);
		}
    }
  }
  fclose(fErr);
}

void RecordErrorElement(int ind, int mode) {
  int k, locInd, ind2;

  if (mode == 0) {
	locInd = iAct[ind];
	printf("actin: %d (%d/%d)\n", ind, locInd, nActMe);
	if (locInd > -1 && locInd < nActMe + nActCp) {
		printf("act.r: %g\t%g\t%g\n", P2(act.r,locInd,0), P2(act.r,locInd,1),
				P2(act.r,locInd,2));
		printf("act.f: %g\t%g\t%g\n", P2(act.f,locInd,0), P2(act.f,locInd,1),
				P2(act.f,locInd,2));
	    for(k = 0; k < nChAc; k++) {
	        ind2 = P2A(act.ch,locInd,k,nChAc);
	        CONT(ind2 < 0);
	        if (k < 2) {
	            printf("%d: %d(%d/%d): ", k, ind2, iAct[ind2], nActMe);
	            if (iAct[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(act.r,iAct[ind2],0),
	                        P2(act.r,iAct[ind2],1), P2(act.r,iAct[ind2],2));
    	            printf("%g\t%g\t%g\t", P2(act.f,iAct[ind2],0),
	                        P2(act.f,iAct[ind2],1), P2(act.f,iAct[ind2],2));
	            }
	            printf("\n");
			}
			else {
	            printf("%d: %d(%d/%d): ", k, ind2, iAbp[ind2], nAbpMe);
	            if (iAbp[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(abp.r,iAbp[ind2],0),
	                        P2(abp.r,iAbp[ind2],1), P2(abp.r,iAbp[ind2],2));
    	            printf("%g\t%g\t%g\t", P2(abp.f,iAbp[ind2],0),
	                        P2(abp.f,iAbp[ind2],1), P2(abp.f,iAbp[ind2],2));
	            }
	            printf("\n");
			}
		}
	}
  }
  else if (mode == 1) {
	locInd = iAbp[ind];
	printf("ABP: %d (%d/%d)\n", ind, locInd, nAbpMe);
	if (locInd > -1 && locInd < nAbpMe + nAbpCp) {
		printf("abp.r: %g\t%g\t%g\n", P2(abp.r,locInd,0), P2(abp.r,locInd,1),
				P2(abp.r,locInd,2));
		printf("abp.f: %g\t%g\t%g\n", P2(abp.f,locInd,0), P2(abp.f,locInd,1),
				P2(abp.f,locInd,2));
	    for(k = 0; k < nChAb; k++) {
	        ind2 = P2A(abp.ch,locInd,k,nChAb);
	        CONT(ind2 < 0 || k == 2);
	        if (k < 2) {
	            printf("%d: %d(%d/%d): ", k, ind2, iAct[ind2], nActMe);
	            if (iAct[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(act.r,iAct[ind2],0),
	                        P2(act.r,iAct[ind2],1), P2(act.r,iAct[ind2],2));
    	            printf("%g\t%g\t%g\t", P2(act.f,iAct[ind2],0),
	                        P2(act.f,iAct[ind2],1), P2(act.f,iAct[ind2],2));
	            }
	            printf("\n");
			}
			else {
	            printf("%d: %d(%d/%d): ", k, ind2, iAbp[ind2], nAbpMe);
	            if (iAbp[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(abp.r,iAbp[ind2],0),
	                        P2(abp.r,iAbp[ind2],1), P2(abp.r,iAbp[ind2],2));
    	            printf("%g\t%g\t%g\t", P2(abp.f,iAbp[ind2],0),
	                        P2(abp.f,iAbp[ind2],1), P2(abp.f,iAbp[ind2],2));
	            }
	            printf("\n");
			}
		}
	}
  }
  else {
	locInd = iMb[ind];
	printf("membrane: %d (%d/%d)\n", ind, locInd, nMbMe);
	if (locInd > -1 && locInd < nMbMe + nMbCp) {
		printf("memb.idx: %d\n", memb.idx[locInd]);
		printf("memb.r: %g\t%g\t%g\n", P2(memb.r,locInd,0), P2(memb.r,locInd,1),
				P2(memb.r,locInd,2));
		printf("memb.f: %g\t%g\t%g\n", P2(memb.f,locInd,0), P2(memb.f,locInd,1),
				P2(memb.f,locInd,2));
	    for(k = 0; k < nChMb; k++) {
	        ind2 = P2A(memb.ch,locInd,k,nChMb);
	        CONT(ind2 < 0);
	        if (k < 2) {
	            printf("%d: %d(%d/%d): ", k, ind2, iMb[ind2], nMbMe);
	            if (iMb[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(memb.r,iMb[ind2],0),
	                        P2(memb.r,iMb[ind2],1), P2(memb.r,iMb[ind2],2));
	                printf("%g\t%g\t%g\t", P2(memb.f,iMb[ind2],0),
	                        P2(memb.f,iMb[ind2],1), P2(memb.f,iMb[ind2],2));
	            }
	            printf("\n");
			}
	        else {
    	        printf("%d: %d(%d/%d): ", k, ind2, iAct[ind2], nActMe);
	            if (iAct[ind2] > -1) {
    	            printf("%g\t%g\t%g\t", P2(act.r,iAct[ind2],0),
	                        P2(act.r,iAct[ind2],1), P2(act.r,iAct[ind2],2));
	                printf("%g\t%g\t%g\t", P2(act.f,iAct[ind2],0),
	                        P2(act.f,iAct[ind2],1), P2(act.f,iAct[ind2],2));
	            }
	            printf("\n");
	        }
	    }
	}
  }
}

// Check large forces beyond "magUnstF".
int CheckLargeForce(double f, int kind) {
  int chkErr = 1;

  if  (fabs(f) > magUnstF) {
    Printf("Error: too large force is produced.\n");
    Printf("currTimeStep: %d, rank %d, force: %lf\n", currTimeStep, rank, f);
    stopSig = 1;
    RecordError(kind);
    chkErr = -1;
  }
  if (isnan(f) != 0) {
    Printf("Error: some force has a 'nan' value.\n");
    Printf("currTimeStep: %d, rank %d, force: nan\n", currTimeStep, rank);
    stopSig = 1;
    chkErr = -1;
  } 
  return chkErr;
}

void CheckLargeTotalForceAll(int mode) {
  int n;
  for(n = 0; n < nActMe; n++) {
	CONT(!(fabs(P2(act.f,n,0)) > magUnstF));
	printf("Large force!\n");
	printf("mode = %d\n, rank = %d, currTimeStep = %lld\n", mode, rank, 
			currTimeStep);
	RecordErrorElement(act.id[n], 0);	
	exit(-1);
  }
  for(n = 0; n < nAbpMe; n++) {
	CONT(!(fabs(P2(abp.f,n,0)) > magUnstF));
	printf("Large force!\n");
	printf("mode = %d\n", mode);
	printf("rank = %d, currTimeStep = %lld\n", rank, currTimeStep);
	RecordErrorElement(abp.id[n], 1);	
	exit(-1);
  }
  for(n = 0; n < nMbMe; n++) {
	CONT(!(fabs(P2(memb.f,n,0)) > magUnstF));
	printf("Large force!\n");
	printf("mode = %d\n", mode);
	printf("rank = %d, currTimeStep = %lld\n", rank, currTimeStep);
	RecordErrorElement(memb.id[n], 2);	
	exit(-1);
  }
}

void CheckNanForce(int mode) {
  int n, k, ind;
  FOR_ACTME(n) {
	CONT(!(isnan(P2(act.r,n,0)) != 0 || isnan(P2(act.r,n,1)) != 0 
			|| isnan(P2(act.r,n,2)) != 0 || isnan(P2(act.f,n,0)) != 0 
			|| isnan(P2(act.f,n,1)) != 0 || isnan(P2(act.f,n,2)) != 0));
	printf("NaN force!\n");
	printf("mode = %d\n, rank = %d, currTimeStep = %lld\n", mode, rank, 
			currTimeStep);
	RecordErrorElement(act.id[n], 0);
	exit(-1);
  }

  FOR_ABPME(n) {
	CONT(!(isnan(P2(abp.r,n,0)) != 0 || isnan(P2(abp.r,n,1)) != 0 
			|| isnan(P2(abp.r,n,2)) != 0 || isnan(P2(abp.f,n,0)) != 0 
			|| isnan(P2(abp.f,n,1)) != 0 || isnan(P2(abp.f,n,2)) != 0));
	printf("NaN force!\n");
	printf("mode = %d\n, rank = %d, currTimeStep = %lld\n", mode, rank, 
			currTimeStep);
	RecordErrorElement(abp.id[n], 1);
	exit(-1);
  }

  FOR_MBME(n) {
	CONT(!(isnan(P2(memb.r,n,0)) != 0 || isnan(P2(memb.r,n,1)) != 0 
			|| isnan(P2(memb.r,n,2)) != 0 || isnan(P2(memb.f,n,0)) != 0 
			|| isnan(P2(memb.f,n,1)) != 0 || isnan(P2(memb.f,n,2)) != 0));
	printf("NaN force!\n");
	printf("mode = %d\n, rank = %d, currTimeStep = %lld\n", mode, rank, 
			currTimeStep);
	RecordErrorElement(memb.id[n], 2);
	exit(-1);
  }
}


