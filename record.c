// This file contains functions which record data.

/*------------------- Recording the progress of simulations ------------------*/

void RecordProgressSubroutine(int *cntMe, FILE *fOut) {
  int *cntAll, sumCnt;

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt);
	Printf("%5d ", sumCnt);
  }
  free(cntAll);
}

void RecordProgressSubroutine2(int *nPar1, int *nPar2, int nPar3, 
		FILE *fOut, int mode) {
  int n, *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  V2SET_ALL(sumCnt, 0);
  for(n = 0; n < 2; n++) {
	CONT(!((n == 0 && (int)(mode / 2) == 1) || (n == 1 && mode % 2 == 1)));
	MPI_Gather(((n == 0) ? nPar1 : nPar2), 1, MPI_INT, cntAll, 
			1, MPI_INT, 0, MPI_COMM_WORLD);
	CONT(rank != 0);
	sumCnt[n] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[n]);
	Printf("%5d/", sumCnt[n]);
  }
  if (rank == 0) {
	if (mode == 0) { 
		fprintf(fOut, "%d\t", nPar3);
		Printf("%5d ", nPar3);
	}
	else {
		fprintf(fOut, "%d\t%d\t", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
		Printf("%5d/%5d ", nPar3 - sumCnt[0] - sumCnt[1], nPar3);
	}
  }
  free(cntAll);
}

void RecordProgressSubroutine3(int *cntMe, int *cntMe2, FILE *fOut) {
  int *cntAll, sumCnt[2];

  MALLOC(cntAll,int,nCpu);
  MPI_Gather(cntMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[0] = SumArrInt(cntAll, nCpu);
  }
  MPI_Gather(cntMe2, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt[1] = SumArrInt(cntAll, nCpu);
	fprintf(fOut, "%d\t", sumCnt[0] + sumCnt[1]);
	Printf("%5d ", sumCnt[0] + sumCnt[1]);
  }
  free(cntAll);
}

// Record the progress of the simulation.
void RecordProgress(void) {
  int n, nAcp, mode, nMot2[3], *pArr;
  double vol, area;
  time_t now;
  FILE *fOut;

  if (rank == 0) {
	fOut = fopen(GenFileName("Progress"), "a");
	fprintf(fOut, "%lld\t%g\t", currTimeStep, (double)currTimeStep * dtReal);
	Printf("%12d ", currTimeStep);
  }
  nAcp = nAbp - nMot;
  // Gather information from all subdomains first.
  if (nAct > 0) {
	mode = 0;
	if (actAss.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&actM.c, &actM.c, nAct, fOut, mode);
	if (actNuc.gTgl != 0 || actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ RecordProgressSubroutine(&nActFilaMe, fOut); }
	if (actNuc.gTgl != 0) { RecordProgressSubroutine(&actNuc.cntMe, fOut); }
	if (actAss.gTgl != 0) { RecordProgressSubroutine(&actAss.cntMe, fOut); }
	if (actDis.gTgl != 0 || actBst.gTgl != 0) 
	{ RecordProgressSubroutine(&actDis.cntMe, fOut); }
	if (actBst.gTgl != 0) { RecordProgressSubroutine(&actBst.cntMe, fOut); }
	if (actSev.gTgl != 0) { RecordProgressSubroutine(&actSev.cntMe, fOut); }
	if (actAnn.gTgl != 0) { RecordProgressSubroutine(&actAnn.cntMe, fOut); }
	if (actCap.gTgl != 0) { RecordProgressSubroutine(&actCap.cntMe, fOut); }
	if (actUnc.gTgl != 0) { RecordProgressSubroutine(&actUnc.cntMe, fOut); }
  }
  if (nAcp > 0) { 
	mode = 0;
	if (acpInaUnb.gTgl != 0 || acpMoBind.gTgl != 0) { mode += 2; }
	if (acpUnb.gTgl != 0 || acpReb.gTgl != 0) { mode++; }
	RecordProgressSubroutine2(&nAcpMme, &nAcpInaMe, nAcp, fOut, mode); 
	if (acpInaUnb.gTgl != 0) 
	{ RecordProgressSubroutine(&acpInaUnb.cntMe, fOut); }
	if (acpMoBind.gTgl != 0) 
	{ RecordProgressSubroutine(&acpMoBind.cntMe, fOut); }
	if (acpUnb.gTgl != 0) { RecordProgressSubroutine(&acpUnb.cntMe, fOut); }
	if (acpReb.gTgl != 0) { RecordProgressSubroutine(&acpReb.cntMe, fOut); }
  }
  if (nMot > 0) { 
	mode = 0;
	if (motInaUnb.gTgl != 0 || motMoBind.gTgl != 0) { mode += 2; }
	if (motUnb.gTgl != 0 || motReb.gTgl != 0) { mode++; }
	if (motSA.gTgl == 0) { 
		RecordProgressSubroutine2(&nMotMme, &nMotInaMe, nMot, fOut, mode); 
		if (motInaUnb.gTgl != 0) 
		{ RecordProgressSubroutine(&motInaUnb.cntMe, fOut); }
		if (motMoBind.gTgl != 0) 
		{ RecordProgressSubroutine(&motMoBind.cntMe, fOut); }
		if (motUnb.gTgl != 0) { RecordProgressSubroutine(&motUnb.cntMe, fOut); }
		if (motReb.gTgl != 0) { RecordProgressSubroutine(&motReb.cntMe, fOut); }
	}
	else {
		nMot2[0] = 0;
		FOR_ABPME(n) {
			pArr = &P2A(abp.ch,n,0,nChAb);
			CONT(pArr[2] != 2);
			CONT(!(pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0));
			nMot2[0] += 2;
		}
		nMot2[1] = 2 * nMotMme - nMot2[0] + nMotInaMe;
		nMot2[2] = nMot * 2;
		RecordProgressSubroutine2(&nMot2[0], &nMot2[1], nMot2[2], fOut, mode); 
		if (motUnb.gTgl != 0) {
			RecordProgressSubroutine3(&motUnb.cntMe, &motInaUnb.cntMe, fOut); 
		}
		if (motReb.gTgl != 0) { 
			RecordProgressSubroutine3(&motReb.cntMe, &motMoBind.cntMe, fOut); 
		}
	}
	if (motWalk.gTgl != 0) { RecordProgressSubroutine(&motWalk.cntMe, fOut); }
    if (motSA.gTgl != 0) {
		RecordProgressSubroutine(&motSA.cntNucMe, fOut);
		RecordProgressSubroutine(&motSA.cntAssMe, fOut);
		if (motSA.to.gTgl != 0) {
			RecordProgressSubroutine(&motSA.cntTurnMe, fOut);
		}
    }
  }
  if (gTglMb != 0) {
	RecordProgressSubroutine(&nMbMe, fOut);

	if (mbSld.act.gTgl != 0) { 
		RecordProgressSubroutine(&mbSld.act.cntMe, fOut); 
	}
	if (mbSld.abp.gTgl != 0) { 
		RecordProgressSubroutine(&mbSld.abp.cntMe, fOut); 
	}
	if (mbUnb.gTgl != 0) { RecordProgressSubroutine(&mbUnb.cntMe, fOut); }
	if (mbReb.gTgl != 0) { RecordProgressSubroutine(&mbReb.cntMe, fOut); }
	if (mbDef.gTgl != 0) { RecordProgressSubroutine(&mbDef.cntMe, fOut); }
	if (mbFix.gTgl != 0) { RecordProgressSubroutine(&mbFix.cntMe, fOut); }
	if (mbPro.gTgl != 0) { RecordProgressSubroutine(&mbPro.cntMe, fOut); }
	if (mbMat.gTgl != 0) { 
		RecordProgressSubroutine(&mbMat.cntMe2[0], fOut); 
		RecordProgressSubroutine(&mbMat.cntMe2[1], fOut); 
	}
	if (gTglNuc != 0) {
		if (nucUnb.gTgl != 0) { RecordProgressSubroutine(&nucUnb.cntMe, fOut); }
		if (nucReb.gTgl != 0) { RecordProgressSubroutine(&nucReb.cntMe, fOut); }
	}
  }
  if (rank == 0) {
	now = time(NULL);
	fprintf(fOut, "%s", ctime(&now));
	fclose(fOut);
	Printf("\n");
  }

}

/*------------------- Recording the progress of simulations ------------------*/

/*------------------------- Recording parameter values -----------------------*/

void RecordInitParameterSubroutine(const char *str, int tgl1, int tgl2, 
		int tgl3, FILE *fOut) {
  fprintf(fOut, "%s = %s / %s / %s\n", str, 
		((tgl1 != 0 && tgl3 != 0) ? "on" : "off"),
		((tgl1 != 0 && tgl2 != 0) ? "on" : "off"), (tgl1 != 0 ? "on": "off"));
}

void RecordInitParameterSubroutine2(const char *str, 
		FuncCont *cont, FILE *fOut) {
  fprintf(fOut, "%s = %s\n", str, (cont->tgl != 0) ? "on" : "off");
  if (cont->tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", cont->prdR, cont->prd);
  }
}

// Record the parameter values which were initially loaded from "condition" 
// and "Config".
void RecordInitParameter(void) {
  int chk[NDIM][2], chk2, k, nAcp;
  double dimDomC[NDIM] , vol, ratio, ratio2, cAbp, rAbp;
  char direc[4] = "xyz", ch[80];
  FILE *fOut;  

  nAcp = nAbp - nMot;

  V3DIV(dimDomC, dimDom, nCell);
  vol = V3PROD(dimDom) * CUBE(L_SCALE_IN_M) * 1e3;
  cAbp = (double)nAbp * 2. / N_AVO / vol * 1e6;
  rAbp = (double)nAbp / ((double)nAct * nActPerSeg);
  fOut = fopen(GenFileName("Parameter"), "w");

  fprintf(fOut, "Duration of initial network formation = %g [s] "
		"(%lld timesteps)\n", netForm.durR, netForm.dur);
  fprintf(fOut, "Time step = %g [s] (%g)\n", dtReal, dt);
  fprintf(fOut, "Threshold of unstable force = %g [pN] (%g)\n",
		MAG_UNSTABLE_FORCE * 1e12, magUnstF);

  fprintf(fOut, "\n================== Number and concentration of actins and "
		"ABPs =================\n");
  fprintf(fOut, "Number of actin segments = %d (C_A = %g uM)\n", nAct, 
		(double)nAct * 2. * nActPerSeg / N_AVO / vol * 1e6);
  fprintf(fOut, "Length of actin segments = %g [nm]\n", L_SCALE_IN_NM);
  fprintf(fOut, "Number of ABPs = %d (C_ABP = %g uM, R_ABP = %g)\n", nAbp, 
		cAbp, rAbp);
  if (nAbp > 0) {
	ratio = (nAbp > 0) ? (double)nMot / (double)nAbp : 0.;
	ratio2 = (nAcp > 0) ? (double)nAbpDet[0] / (double)nAcp : 0.;
	fprintf(fOut, "(N_ACP^C = %d, C_ACP^C = %g uM, R_ACP^C = %g)\n", nAbpDet[0],
			cAbp * (1. - ratio) * ratio2, rAbp * (1. - ratio) * ratio2);
	fprintf(fOut, "(N_ACP^B = %d, C_ACP^B = %g uM, R_ACP^B = %g)\n", nAbpDet[1],
			cAbp * (1. - ratio) * (1. - ratio2), rAbp * (1. - ratio) 
			* (1. - ratio2));
	fprintf(fOut, "(N_M = %d, C_M = %g uM, R_M = %g)\n", nMot, 
			cAbp * ratio, rAbp * ratio);
  }

  fprintf(fOut, "\n===================== Conditions for domain and boundaries "
		"=====================\n");
  fprintf(fOut, "Domain width: x = %g [um] (%g), y = %g [um] (%g),"
		" z = %g [um] (%g)\n", L_S2UM(dimDom[0]), dimDom[0], 
		L_S2UM(dimDom[1]), dimDom[1], L_S2UM(dimDom[2]), dimDom[2]);
  fprintf(fOut, "Periodic boundary condition: x = %s, y = %s, z = %s\n", 
		(pbc[0] == 1 ? "on" : "off"), (pbc[1] == 1 ? "on" : "off"),
		(pbc[2] == 1 ? "on" : "off"));
  fprintf(fOut, "Geometry of a domain = ");
  if (dir2D > -1) {
	fprintf(fOut, "2-D with a normal direction in %c\n", direc[dir2D]);
  }
  else {
	fprintf(fOut, "3-D\n");
  }
  fprintf(fOut, "Movement of boundaries following stress-strain relation"
		" = %s\n", ((bndMv.gTgl != 0) ? "on" : "off"));
  if (bndMv.gTgl != 0) {
	if (bndMv.stfUnit == 0) {
		fprintf(fOut, "Young's modulus = %g on -x, %g on +x, %g on -y, %g on "
				"+y, %g on -z, %g on +z [Pa]\n", S2PA(bndMv.stf[0]), 
				S2PA(bndMv.stf[3]), S2PA(bndMv.stf[1]), S2PA(bndMv.stf[4]), 
				S2PA(bndMv.stf[2]), S2PA(bndMv.stf[5]));
	}
	else {
		fprintf(fOut, "Young's modulus = %g on -x, %g on +x, %g on -y, %g on "
				"+y, %g on -z, %g on +z [N/m]\n", KS_S2NPM(bndMv.stf[0]), 
				KS_S2NPM(bndMv.stf[3]), KS_S2NPM(bndMv.stf[1]), 
				KS_S2NPM(bndMv.stf[4]), KS_S2NPM(bndMv.stf[2]), 
				KS_S2NPM(bndMv.stf[5]));
	}
	fprintf(fOut, "Conservation of volume encapsulated by boundaries = %s\n",
			((bndVol.gTgl != 0) ? "on" : "off"));
	if (bndVol.gTgl != 0) {
		fprintf(fOut, "Stiffnesses of volume conservation = %g [N/m^5] (%g)\n",
				bndVol.stf * (KT_IN_J / pow(L_SCALE_IN_M, 6.)), bndVol.stf);
	}
  }
  if (recBndLoc.prdR > 0.) {
	fprintf(fOut, "Period of recording the boundary location = %g [s] "
			"(%d timesteps)\n", recBndLoc.prdR, recBndLoc.prd);
  }
  if (nAct > 0) {
	fprintf(fOut, "Unbinding of actin filaments from boundaries"
	 		" = %s\n", ((bndUnb.gTgl != 0) ? "on" : "off"));
	if (bndUnb.gTgl != 0) { 
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				K0_BND_UNB);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_BND_UNB);
		fprintf(fOut, "Factor that adjusts k0 = %g on -x, %g on "
				"+x, %g on -y, %g on +y, %g on -z, %g on +z\n", bndUnb.facK0[0],
				bndUnb.facK0[1], bndUnb.facK0[2], bndUnb.facK0[3], 
				bndUnb.facK0[4], bndUnb.facK0[5]);
		fprintf(fOut, "Factor that adjusts x = %g on -x, %g on "
				"+x, %g on -y, %g on +y, %g on -z, %g on +z\n", bndUnb.facX[0], 
				bndUnb.facX[1], bndUnb.facX[2], bndUnb.facX[3], bndUnb.facX[4], 
				bndUnb.facX[5]);
	}
	fprintf(fOut, "Binding of actin filaments on boundaries"
			" = %s\n", ((bndReb.gTgl != 0) ? "on" : "off"));
	if (bndReb.gTgl != 0) {
		fprintf(fOut, "Binding of entire parts of the filaments"
				" = %s\n", ((bndReb.gTglPa != 0) ? "on" : "off"));
		fprintf(fOut, "Thickness of space at which filaments can bind = %g "
				"[nm] (%g)\n", bndReb.depR, bndReb.dep);
		fprintf(fOut, "Binding rate(k) = %g [1/uM s]\n", K_BND_REB);
		fprintf(fOut, "Factor that adjusts k = %g on -x, %g on +x, %g on -y, "
				"%g on +y, %g on -z, %g on +z\n", bndReb.facK[0], 
				bndReb.facK[1], bndReb.facK[2], bndReb.facK[3], bndReb.facK[4],
				bndReb.facK[5]);
	}
	if (bndMat.gTgl != 0 && bndUnb.gTgl != 0) {
		fprintf(fOut, "Manuration rate of links between actin filaments and "
				"boundary(k0) = %g [1/s]\n", K0_ACT_MAT);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_ACT_MAT);
		fprintf(fOut, "Factor that adjusts k0 = %g\n", bndMat.facK0);
		fprintf(fOut, "Factor that adjusts x = %g\n", bndMat.facX);
		fprintf(fOut, "Maximum level of maturation = %d\n", bndMat.maxNFA);
	}
	fprintf(fOut, "Moving of actin filaments on boundaries"
			" = %s\n", ((bnd.gTglActMv != 0) ? "on" : "off"));
	if (bnd.gTglActMv != 0) {
		fprintf(fOut, "Drag coefficient for the movement (relative to actin "
				"segment) = %g on -x, %g on +x, %g on -y, %g on +y, %g on -z, "
				"%g on +z\n", bnd.drag[0], bnd.drag[1], bnd.drag[2], 
				bnd.drag[3], bnd.drag[4], bnd.drag[5]);
	}
	if (bndUnb.gTgl != 0 || bndReb.gTgl != 0) {
		if (recBndUnbReb.prdR > 0.) {
			fprintf(fOut, "Period of recording the unbinding and binding on "
					"boundaries = %g [s] (%d timesteps)\n", recBndUnbReb.prdR, 
					recBndUnbReb.prd);
		}
	}
  }

  fprintf(fOut, "\n======================== Information about measurement "
			"=========================\n");
  fprintf(fOut, "Net duration of simulation = %g [s] "
		"(%lld timesteps)\n", rheo.durR, rheo.dur);
  fprintf(fOut, "Duration of motor activation before rheological measurement "
		"= %g [s] (%lld timesteps)\n", motActiv.durR, motActiv.dur);

  strcpy(ch, "");
  if (rheoWay == -1) { strcpy(ch, "neither"); }
  else if (rheoWay == 0 || rheoWay == 2) { strcpy(ch, "segment-rheology"); }
  if (rheoWay > 0) {
	if (rheoWay == 2) { strcat(ch, " and "); }
	if (bulkRheoType == 0) { strcat(ch, "shear-"); }
	else { strcat(ch, "normal-"); }
	if (bulkRheoWay == 1) { strcat(ch, "strain-controlled bulk rheology"); }
	else { strcat(ch, "stress-controlled bulk rheology"); }
  }
  fprintf(fOut, "Measurement method = %s\n", ch);
  if (rheoWay > 0) {
	fprintf(fOut, "---------------------------------- Bulk rheology "
			"-------------------------------\n");
    fprintf(fOut, "Direction of applied %s = %c%c\n", ((bulkRheoType == 0) 
			? "shear" : "normal"), ((signStr > 0) ? '+' : '-'),  
			direc[dirStr]);
	if (bulkRheoType == 0) {
		fprintf(fOut, "Direction normal to the sheared surface = %c\n", 
				direc[dirNoPBC]);
	}
	if (bulkRheoWay == 1) {
		fprintf(fOut, "Amount of prestrain = %g [percent] (%g [um])\n", 
				pres.mag * 100., pres.mag * L_S2UM(dimDom[(bulkRheoType == 0) 
				? dirNoPBC : dirStr]));
	    fprintf(fOut, "Prestrain rate = %g [1/s] (%g [um/s])\n", 
				pres.rate, pres.rate * L_S2UM(dimDom[(bulkRheoType == 0) 
				? dirNoPBC : dirStr]));
		fprintf(fOut, "Amplitude of sinusoidal strain = %g "
				"[percent]\n", 100. * sinuStr.amp);
		fprintf(fOut, "Period of sinusoidal strain = %g [s] (%d timesteps)\n",
				sinuStr.prdR, sinuStr.prd);
	}
	else {
		fprintf(fOut, "Amount of prestress = %g [Pa]\n", pres.mag);
	    fprintf(fOut, "Prestress rate = %g [Pa/s]\n", pres.rate);
		fprintf(fOut, "Amplitude of sinusoidal stress = %g [Pa]\n",
				sinuStr.amp);
		fprintf(fOut, "Period of sinusoidal stress = %g [s] (%d timesteps)\n",
				sinuStr.prdR, sinuStr.prd);
	}
	fprintf(fOut, "Lower and uppder limits of strain beyond which a run is "
			"terminated = %g, %g\n", stra.lim[0], stra.lim[1]);
	fprintf(fOut, "Period of recording %s = %g [s] (%d timesteps)\n", 
			((bulkRheoWay == 0) ? "stress" : "strain"), 
			recStre.prdR, recStre.prd);
  }
  if (rheoWay == 0 || rheoWay == 2) {
	fprintf(fOut, "-------------------------------- Segment rheology "
			"------------------------------\n");
	fprintf(fOut, "Period of recording trajectory = %g [s] (%d timesteps)\n", 
			recTraj.prdR, recTraj.prd);
	fprintf(fOut, "The number of recorded actin segments = %d (%g [percent])\n",
			recTraj.nActL, (double)recTraj.nActL / (double)nActGoal * 100.);
	fprintf(fOut, "The number of recorded ABPs = %d (%g [percent])\n", 
			recTraj.nAbpL, (double)recTraj.nAbpL / (double)nAbpGoal * 100.);
	fprintf(fOut, "Choose new components for recording the trajectory = %s\n",
			((recTraj.gTglCho != 0) ? "on" : "off"));
  }
  fprintf(fOut, "\n============= Conditions for dynamic behaviors of actins "
		"and ABPs ==============\n");
  RecordInitParameterSubroutine2("Equilibration of counter for monomeric "
		"actins and ABPs", &updMono, fOut);

  if (nAct > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors of "
			"actins ==================\n");
	RecordInitParameterSubroutine("Nucleation of actins", 
			actNuc.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	fprintf(fOut, "Control of the number of actin filaments by the "
			"nucleation = %s\n", (actNuc.gTglFN != 0 ? "on" : "off"));
	RecordInitParameterSubroutine("Assembly of actins", 
			actAss.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Disassembly of actins", 
			actDis.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Bursting disassembly of actins", 
			actBst.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Severing of actins", 
			actSev.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Annealing of actins", 
			actAnn.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Capping of actins at ends", 
			actCap.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Capping of actins by severing", 
			gTglActSevCap, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Uncapping of actins", 
			actUnc.gTgl, gTglActDynPres, gTglActDynNF, fOut);
	RecordInitParameterSubroutine("Degradation of actins", 
			actDgd.gTgl, gTglActDynPres, gTglActDynNF, fOut);
  }
  if (nAcp > 0) {
	fprintf(fOut, "\n=================== Conditions for dynamic behaviors "
			"of ACPs ===================\n");
	RecordInitParameterSubroutine("Unbinding of inactive ACPs", 
			acpInaUnb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of monomeric ACPs", 
			acpMoBind.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Unbinding of active ACPs", 
			acpUnb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	RecordInitParameterSubroutine("Binding of inactive ACPs", 
			acpReb.gTgl, gTglAcpDynPres, gTglAcpDynNF, fOut);
	fprintf(fOut, "Check cross-linking angle at binding of inactive ACPs = "
			"%s\n", (acpReb.gTglCrsAng != 0 ? "on" : "off"));
	fprintf(fOut, "Implicit consideration of monomeric ACPs = %s\n", 
			(gTglImpAcpM != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "\n================== Conditions for dynamic behaviors "
			"of motors ==================\n");
	RecordInitParameterSubroutine("Walking of active motors", 
			motWalk.gTgl, gTglMotWalkPres, gTglMotWalkNF, fOut);
	fprintf(fOut, "Slide-off by the walking = %s\n",
			(gTglMotWalkSld != 0 ? "on" : "off"));
	fprintf(fOut, "Self-assembly of motors = %s\n",
			(motSA.gTgl != 0 ? "on" : "off"));
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Regulate relative constant size of the self-assembled "
				"structure = %s\n",	(motSA.gTglConSiz != 0 ? "on" : "off"));
		RecordInitParameterSubroutine("Unbinding of motors", 
				motUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of motors", 
				motReb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Check the alignment between actin and motor filaments "
				"for binding of motors = %s\n", (motReb.gTglCrsAng != 0 
				? "on" : "off"));
		fprintf(fOut, "Allow the binding of motors in a dimer only in opposite "
				"directions = %s\n", (motReb.gTglOppDir != 0 ? "on" : "off"));
		fprintf(fOut, "Turnover of motor filaments = ");
		if (motSA.to.gTgl == 0) {
			fprintf(fOut, "no\n");
		}
		else {
			if (motSA.to.mode == 0) {
				fprintf(fOut, "free\n");
			}
			else {
				fprintf(fOut, "free & force-dependent\n");
			}
		}
  		RecordInitParameterSubroutine2("Equilibration of counter for motor "
				"self-assembly", &updMotN, fOut);
	}
	else {
		RecordInitParameterSubroutine("Unbinding of inactive motors", 
				motInaUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of monomeric motors", 
				motMoBind.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Unbinding of active motors", 
				motUnb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		RecordInitParameterSubroutine("Binding of inactive motors", 
				motReb.gTgl, gTglMotUnbRebPres, gTglMotUnbRebNF, fOut);
		fprintf(fOut, "Check cross-linking angle for binding of inactive "
				"motors = %s\n", (motReb.gTglCrsAng != 0 ? "on" : "off"));
		fprintf(fOut, "Implicit consideration of monomeric motors = %s\n", 
				(gTglImpMotM != 0 ? "on" : "off"));
	}
  }
  fprintf(fOut, "\n======== Adjustment for rates of the dynamic behaviors of "
		"actin and ABP ========\n");
  fprintf(fOut, "During network formation(i*k) = %g\n", netForm.facK);
  fprintf(fOut, "During prestrain/prestress(i*k) = %g\n", pres.facK);
  fprintf(fOut, "During motor activation(i*k) = %g\n", motActiv.facK);

  if (nAct > 0) {
	if (actNuc.gTgl != 0 || actAss.gTgl != 0 || actDis.gTgl != 0 
			|| actBst.gTgl != 0 || actSev.gTgl != 0 || actAnn.gTgl != 0 
			|| actCap.gTgl != 0 || actUnc.gTgl != 0 || actDgd.gTgl != 0) {
		fprintf(fOut, "\n======================== Parameters of actin dynamics "
				"=========================\n");
	}
	if (actNuc.gTgl != 0) {
		fprintf(fOut, "Nucleation rate = %g [1/uM s]\n", actNuc.k);
	}
	if (actAss.gTgl != 0) {
		fprintf(fOut, "Assembly rate = %g [1/uM s] at barbed, %g [1/uM s] at "
				"pointed\n", actAss.k[0], actAss.k[1]);
	}
	if (actDis.gTgl != 0) {
		fprintf(fOut, "Disassembly rate = %g [1/s] at barbed, %g [1/s] at "
				"pointed\n", actDis.k[0], actDis.k[1]);
		fprintf(fOut, "Factor for varying disassembly rate of actin with "
				"ABPs = %g\n", actDis.facKWA);
	}
	if (actBst.gTgl != 0) {
		fprintf(fOut, "Bursting disassembly rate = %g [1/s] at barbed, %g "
				"[1/s] at pointed\n", actBst.k[0], actBst.k[1]);
		fprintf(fOut, "Factor for varying bursting disassembly rate of actin "
				"with ABPs = %g\n", actBst.facKWA);
	}
	if (actSev.gTgl != 0) {
		fprintf(fOut, "k0 for severing rate = %g [1/s]\n", actSev.k);
		fprintf(fOut, "Sensitivity of severing rate to angle = %g [deg]\n", 
				actSev.facX);
		fprintf(fOut, "Factor for varying severing rate of actin with ABPs = "
				"%g\n", actSev.facKWA);
	}
	if (actAnn.gTgl != 0) {
		fprintf(fOut, "Annealing rate = %g [1/uM s]\n", actAnn.k);
		fprintf(fOut, "Angular constraint of annealing = 0 +- %g [deg]\n", 
				RAD2DEG(actAnn.ang));	
	}
	if (actCap.gTgl != 0) {
		fprintf(fOut, "Capping rate = %g [1/uM s] at barbed, %g [1/uM s] at "
				"pointed\n", actCap.k[0], actCap.k[1]);
	}
	if (actUnc.gTgl != 0) {
		fprintf(fOut, "Uncapping rate = %g [1/uM s] at barbed, %g [1/uM s] at "
				"pointed\n", actUnc.k[0], actUnc.k[1]);
	}
	if (actDgd.gTgl != 0) {
		fprintf(fOut, "Degradation rate = %g [1/s]\n", actDgd.k);
		fprintf(fOut, "Sensitivity to distance = %g\n", actDgd.x1);
		fprintf(fOut, "Sensitivity to force = %g\n", actDgd.x2);
		fprintf(fOut, "Maximum distance of degradation = %g\n", actDgd.dist);
	}
  }

  if (nAcp > 0) {
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0 || acpInaUnb.gTgl != 0 
			|| acpUnb.gTgl != 0) {
		fprintf(fOut, "\n================= Parameters for unbinding and "
				"binding of ACPs =================\n");
	}
	if (acpMoBind.gTgl != 0 || acpReb.gTgl != 0) {
		fprintf(fOut, "Binding and binding rate(k) = %g [1/uM s]\n", 
				K_ACP_BIND);
		fprintf(fOut, "Factor that adjusts k = %g\n", acpReb.facK);
	}
  }

  if (nMot > 0) {
	if (motSA.gTgl != 0) { 
		fprintf(fOut, "\n====== Parameters for motor self-assembly, "
				"unbinding, binding, and walking =====\n");
	}
	else {
		fprintf(fOut, "\n============== Parameters for motor unbinding, "
				"binding, and walking ============\n");
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0) {  
		fprintf(fOut, "Binding rate(k) = %g [1/uM s]\n", 
				40. * motMC.nHead);
		fprintf(fOut, "Factor that adjusts k = %g\n", motReb.facK);
	}
	if (motInaUnb.gTgl != 0 || motUnb.gTgl != 0) {
		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]\n", 
				(motUnb.gTgl != 0 || motInaUnb.gTgl != 0) ? 
				log(1 - motUnb.p[motUnb.maxF[0] - 1]) / REVSIGN(dtReal) : 0.);
	}
	if (motWalk.gTgl != 0) {
		fprintf(fOut, "Walking rate without applied force(k0) = %g [1/s]\n", 
				(motWalk.gTgl != 0) ? log(1 - motWalk.p[motWalk.maxF[0] - 1]) 
				/ REVSIGN(dtReal) : 0.);
		fprintf(fOut, "Factor that adjusts k0 of unbinding and walking = %g\n", 
				motUnb.facK0);
		fprintf(fOut, "Each walking distance = %g [nm]\n", 
				L_SCALE_IN_NM / (double)nChAcX);
	}
	if (motMoBind.gTgl != 0 || motReb.gTgl != 0 || motInaUnb.gTgl != 0 
			|| motUnb.gTgl != 0 || motWalk.gTgl != 0) {
		fprintf(fOut, "Transition rates between mechanochemical states in each "
				"head(k01, k10, k12, k21, k20) = %g, %g, %g, %g, %g [1/s]\n", 
				motMC.k01, motMC.k10, motMC.k12, motMC.k21, motMC.k20);
		fprintf(fOut, "Number of heads which each motor arm represents = %d\n", 
				motMC.nHead);
	}
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Average number of motors per each self-assembled "
				"structure = %d\n", motSA.nMotPerTF);
		fprintf(fOut, "k for motor assembly = %g [1/s]\n", motSA.kAss);
		fprintf(fOut, "Turnover of free motor filaments = %s\n", 
				(motSA.to.gTgl != 0) ? "on" : "off");
		if (motSA.to.gTgl != 0) {
			fprintf(fOut, "k for motor free turnover = %g [1/s]\n", motSA.to.k);
			if (motSA.to.mode == 1) {
				fprintf(fOut, "k0 for motor force-dependent turnover = %g "
						"[1/s]\n", motSA.to.k0);
				fprintf(fOut, "Sensitivity of motor force-dependent turnover "
						"to force = %g [m]\n", motSA.to.x);
			}
		}
	}
  }

  if (nAbp > 0) {
	fprintf(fOut, "\n============================== Geometry of ABP"
			" =================================\n");
	fprintf(fOut, "Number of binding sites for ABPs on each actin segment in "
			"longitudinal and transverse directions = %d, %d\n", nChAcX, 
			nChAcY);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Length of ACP^C arm = %g [nm] (%g)\n", 
			L_ACPC_ARM * 1.0e9, L_M2S(L_ACPC_ARM));
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Length of ACP^B arm = %g [nm] (%g)\n", 
			L_ACPB_ARM * 1.0e9, L_M2S(L_ACPB_ARM));
  }
  if (nMot > 0) {
	fprintf(fOut, "Length of motor arm = %g [nm] (%g)\n", 
			L_MOT_ARM * 1.0e9, L_M2S(L_MOT_ARM));
	if (motSA.gTgl != 0) {
		fprintf(fOut, "Length of bare zone in self-assembled motors = %g [nm] "
				"(%g)\n", L_MOTBACK_CEN_DIST * 1.0e9, 
				L_M2S(L_MOTBACK_CEN_DIST));
		fprintf(fOut, "Spacing between motors in self-assembled motors = %g "
				"[nm] (%g)\n", L_MOTBACK_DIST * 1.0e9, L_M2S(L_MOTBACK_DIST));
	}
  }
  fprintf(fOut, "\n=============================== Drag coefficients "
		"==============================\n");
  fprintf(fOut, "Viscosity of medium = %g [Pa s] (%g)\n", VISCOSITY, 
		VISCOSITY * L_SCALE_IN_M / actF.dragR);
  if (nAct > 0) {
	fprintf(fOut, "Drag coefficient of actin segment = %g [kg/s] (%g)\n", 
			actF.dragR, 1.0);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^C = %g [kg/s] (%g)\n", 
			abpF.drag[0].n * actF.dragR, abpF.drag[0].n);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "Drag coefficient of ACP^B = %g [kg/s] (%g)\n", 
			abpF.drag[1].n * actF.dragR, abpF.drag[1].n);
  }
  if (nMot > 0) {
	fprintf(fOut, "Drag coefficient of motor = %g [kg/s] (%g)\n", 
			abpF.drag[2].n * actF.dragR, abpF.drag[2].n);
  }

  fprintf(fOut, "\n============================== Thermal fluctuation "
		"=============================\n");
  if (nAct > 0) {
	fprintf(fOut, "Thermal fluctuation of actin filaments = %s\n",
			(gTglActTherm != 0 ? "on" : "off"));
  }
  if (nAcp > 0) {
	fprintf(fOut, "Thermal fluctuation of ACPs = %s\n", 
			(gTglAcpTherm != 0 ? "on" : "off"));
  }
  if (nMot > 0) {
	fprintf(fOut, "Thermal fluctuation of motors = %s\n", 
			(gTglMotTherm != 0 ? "on" : "off"));
  }

  if (gTglMb != 0 ) {
	fprintf(fOut, "\n=========================== Conditions for a membrane "
			"==========================\n");
	fprintf(fOut, "Membrane = ");
	if (dimMbNuc == 2) {
		fprintf(fOut, "2-D cylindrical shell (axis = %c) ", 
				direc[dirNormMbNuc]);
	}
	else {
		fprintf(fOut, "3-D sphere ");
	}
	fprintf(fOut, "with the radius of %g [um] (%g)\n", 
			L_S2UM(radMb), radMb);
	fprintf(fOut, "Load membrane data from Config = %s\n",
			(gTglLoadMbNucData == 1 ? "on": "off"));
	fprintf(fOut, "Length and thickness of cylindrical membrane segments = "
			"%g [um] (%g), %g [um] (%g)\n", L_S2UM(memb.len), memb.len, 
			L_S2UM(memb.thk), memb.thk);
	fprintf(fOut, "Drag coefficient of membrane segments = %g [kg/s] (%g)\n", 
			actF.dragR * memb.dragR, memb.dragR);
	fprintf(fOut, "Bending stiffness of membrane = %g [Nm] (%g)\n",
			KB_S2NM(memb.bend.stf), memb.bend.stf);
	fprintf(fOut, "Extensional stiffness of membrane = %g [N/m] (%g)\n",
			KS_S2NPM(memb.spr.stf), memb.spr.stf);
	fprintf(fOut, "Thresholds for the membrane extension = %g, %g\n", 
			memb.spr.lo, memb.spr.hi);
	fprintf(fOut, "Extensional stiffness of a bond between membrane and bound "
			"actins = %g [N/m] (%g)\n", KS_S2NPM(memb.spr2.stf), memb.spr2.stf);
	fprintf(fOut, "Extensional stiffness of a bond between membrane and space "
			"= %g [N/m] (%g)\n", KS_S2NPM(memb.spr3.stf), memb.spr3.stf);
	fprintf(fOut, "Volume conservation of membrane = %s\n", 
			((mbVol.gTgl != 0) ? "on" : "off"));
	if (mbVol.gTgl != 0) {
		fprintf(fOut, "Stiffnesses of volume conservation = %g [Pa] (%g), %g "
				"[Pa] (%g)\n", S2PA(memb.vol.stf), memb.vol.stf, 
				S2PA(memb.vol.stf2), memb.vol.stf2);
	}
	fprintf(fOut, "Area conservation of membrane = %s\n", 
			((mbAre.gTgl != 0) ? "on" : "off"));
	if (mbAre.gTgl != 0) {
		fprintf(fOut, "Type of area conservation = %s\n", 
				(mbAre.gTgl == 1) ? "local" : "global");
		fprintf(fOut, "Stiffness of area conservation = %g [N/m] (%g)\n", 
				KS_S2NPM(memb.area.stf), memb.area.stf);
	}
	fprintf(fOut, "Thermal fluctuation of membrane = %s\n", 
			((gTglMbTherm != 0) ? "on" : "off"));
	fprintf(fOut, "Protrusion of membrane points = %s\n", 
			((mbPro.gTgl != 0) ? "on" : "off"));
	if (mbPro.gTgl != 0) {
  		fprintf(fOut, "Protrusion rate(k) [1/s] = %g [1/s]\n", K_MEMB_PROT);
		fprintf(fOut, "Factor that adjusts k = %g\n", mbPro.facK);
		fprintf(fOut, "Protrusion force = %g [pN]\n", F_S2PN(mbPro.f));
		fprintf(fOut, "Duration of protrusion = %g [s]\n", mbPro.durR);
	}
	fprintf(fOut, "Fixation of membrane points = %s\n", 
			((mbFix.gTgl != 0) ? "on" : "off"));
	if (mbFix.gTgl != 0) { 
		fprintf(fOut, "Fixation rate(k) = %g [1/s]\n", K_MEMB_FIX);
		fprintf(fOut, "Factor that adjusts k = %g\n", mbFix.facK);
	}
	fprintf(fOut, "Release of fixed membrane points = %s\n", 
			((mbDef.gTgl != 0) ? "on" : "off"));
	if (mbDef.gTgl != 0) { 
  		fprintf(fOut, "Release rate without applied force(k0) = %g [1/s]\n", 
				K0_MEMB_DEF);
		fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
				X_MEMB_DEF);
		fprintf(fOut, "Factor that adjusts k0 = %g\n", mbDef.facK0);
		fprintf(fOut, "Factor that adjusts x = %g\n", mbDef.facX);
	}
	if (nAct > 0) {
		fprintf(fOut, "Relative thickness from membrane within which actin "
				"nucleation is allowed = %g\n", thkMbActNuc);
		fprintf(fOut, "Binding of actin segments to membrane points = ");
		if (mbReb.gTgl != 0 && mbSld.gTgl == 0) {
			fprintf(fOut, "static\n");
		}
		else if (mbReb.gTgl == 0 && mbSld.gTgl != 0) {
			fprintf(fOut, "sliding\n");
		}
		else {
			fprintf(fOut, "none\n");
		}
		fprintf(fOut, "Binding of entire parts of the filaments"
				" = %s\n", ((mbReb.gTglPa != 0) ? "on" : "off"));
		if (mbReb.gTgl != 0) { 
			fprintf(fOut, "Portion of membrane points where actin segments can "
					"bind = %g\n", mbReb.por);
			fprintf(fOut, "Maximum number of actin segments which can bind to"
					" a single membrane point = %d\n", nMbAct);
			fprintf(fOut, "Binding rate(k) = %g [1/s]\n", K_MEMB_REB);
			fprintf(fOut, "Factor that adjusts k = %g\n", mbReb.facK);
			fprintf(fOut, "Unbinding of actin segments from membrane points "
					"= %s\n", ((mbUnb.gTgl != 0) ? "on" : "off"));
			if (mbUnb.gTgl != 0) { 
		  		fprintf(fOut, "Unbinding rate without applied force(k0) = %g "
						"[1/s]\n", K0_MEMB_UNB);
				fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
						X_MEMB_UNB);
				fprintf(fOut, "Factor that adjusts k0 = %g\n", mbUnb.facK0);
				fprintf(fOut, "Factor that adjusts x = %g\n", mbUnb.facX);
			}
			if (mbMat.gTgl != 0) {
		  		fprintf(fOut, "Maturation rate of a membrane point"
						"(k0) = %g [1/s]\n", K0_MEMB_MAT);
				fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
						X_MEMB_MAT);
				fprintf(fOut, "Factor that adjusts k0 = %g\n", mbMat.facK0);
				fprintf(fOut, "Factor that adjusts x = %g\n", mbMat.facX);
				fprintf(fOut, "Maximum level of maturation = %d\n", 
						mbMat.maxNFA);
			}
		}
	}
	if ((mbSld.act.gTgl != 0 && nAct> 0) || (mbSld.abp.gTgl != 0 && nAbp > 0)) {
		if (mbSld.gTgl != 0) {
			if (mbSld.act.gTgl != 0 && mbSld.act.por > 0.) {
				fprintf(fOut, "Maximum portion of actin filaments bound to "
						"membrane = %g (%d)\n", mbSld.act.por, 
						(int)(mbSld.act.por * nAct));
			}
			if (mbSld.abp.gTgl != 0 && mbSld.abp.por > 0.) {
				fprintf(fOut, "Maximum portion of ABPs bound to "
						"membrane = %g (%d)\n", mbSld.abp.por, 
						(int)(mbSld.abp.por * nAbp));
			}
			fprintf(fOut, "Distance within which elements are bound to "
			        "membrane = %g [nm] (%g)\n", L_S2NM(mbSld.critDist), 
					mbSld.critDist);
			fprintf(fOut, "Equilibrium distance between elements and "
					"membrane = %g [nm] (%g)\n", L_S2NM(mbSld.eqDist),
					mbSld.eqDist);
			fprintf(fOut, "Stiffness of elastic constraint between elements "
					"and membrane in a normal direction = %g [N/m] (%g)\n", 
					KS_S2NPM(mbSld.stf), mbSld.stf);
			fprintf(fOut, "Drag coefficient between elements and membrane "
					"in a tangential direction = %g [kg/s] (%g)\n",
					actF.dragR * mbSld.drag, mbSld.drag);
		}
	}
  }

  if (gTglNuc != 0) {
	fprintf(fOut, "\n============================ Conditions for a nucleus "
			"==========================\n");
	fprintf(fOut, "Nucleus = ");
	if (dimMbNuc == 2) {
		fprintf(fOut, "2-D cylindrical shell (axis = %c) ", 
				direc[dirNormMbNuc]);
	}
	else {
		fprintf(fOut, "3-D sphere ");
	}
	fprintf(fOut, "with the radius of %g [um] (%g)\n", 
			L_S2UM(radNuc), radNuc);
	fprintf(fOut, "Load nucleus data from Config = %s\n",
			(gTglLoadMbNucData == 1 ? "on": "off"));
	fprintf(fOut, "Length and thickness of cylindrical nucleus segments = "
			"%g [um] (%g), %g [um] (%g)\n", L_S2UM(memb.nucLen), memb.nucLen, 
			L_S2UM(memb.nucThk), memb.nucThk);
	fprintf(fOut, "Drag coefficient of nucleus segments = %g [kg/s] (%g)\n", 
			actF.dragR * memb.nucDragR, memb.nucDragR);
	fprintf(fOut, "Bending stiffness of nucleus = %g [Nm] (%g)\n",
			KB_S2NM(memb.nucBend.stf), memb.nucBend.stf);
	fprintf(fOut, "Extensional stiffness of nucleus = %g [N/m] (%g)\n",
			KS_S2NPM(memb.nucSpr.stf), memb.nucSpr.stf);
	fprintf(fOut, "Thresholds for the nucleus extension = %g, %g\n", 
			memb.nucSpr.lo, memb.nucSpr.hi);
	fprintf(fOut, "Extensional stiffness of a bond between nucleus and bound "
			"actins = %g [N/m] (%g)\n", KS_S2NPM(memb.nucSpr2.stf), 
			memb.nucSpr2.stf);
	fprintf(fOut, "Volume conservation of nucleus = %s\n", 
			((nucVol.gTgl != 0) ? "on" : "off"));
	if (nucVol.gTgl != 0) {
		fprintf(fOut, "Stiffnesses of volume conservation = %g [Pa] (%g), %g "
				"[Pa] (%g)\n", S2PA(memb.nucVol.stf), memb.nucVol.stf, 
				S2PA(memb.nucVol.stf2), memb.nucVol.stf2);
	}
	fprintf(fOut, "Area conservation of nucleus = %s\n", 
			((nucAre.gTgl != 0) ? "on" : "off"));
	if (nucAre.gTgl != 0) {
		fprintf(fOut, "Type of area conservation = %s\n", 
				(nucAre.gTgl == 1) ? "local" : "global");
		fprintf(fOut, "Stiffness of area conservation = %g [Pa] (%g)\n", 
				S2PA(memb.nucArea.stf), memb.nucArea.stf);
	}
	fprintf(fOut, "Thermal fluctuation of nucleus = %s\n", 
			((gTglNucTherm != 0) ? "on" : "off"));
	if (nAct > 0) {
		fprintf(fOut, "Binding of actin segments to nucleus points = %s\n", 
				((nucReb.gTgl != 0) ? "on" : "off"));
		if (nucReb.gTgl != 0) { 
			fprintf(fOut, "Binding of entire parts of the filaments"
					" = %s\n", ((nucReb.gTglPa != 0) ? "on" : "off"));
			fprintf(fOut, "Portion of nucleus points where actin segments can "
					"bind = %g\n", nucReb.por);
			fprintf(fOut, "Maximum number of actin segments which can bind to"
					" a single nucleus point = %d\n", nNucAct);
			fprintf(fOut, "Binding rate(k) = %g [1/s]\n", K_NUC_REB);
			fprintf(fOut, "Factor that adjusts k = %g\n", nucReb.facK);
		}
		fprintf(fOut, "Unbinding of actin segments from nucleus points = %s\n", 
				((nucUnb.gTgl != 0) ? "on" : "off"));
		if (nucUnb.gTgl != 0) { 
	  		fprintf(fOut, "Unbinding rate without applied force(k0) = %g [1/s]"
					"\n", K0_NUC_UNB);
			fprintf(fOut, "Sensitivity to applied force(x) = %g [m]\n", 
					X_NUC_UNB);
			fprintf(fOut, "Factor that adjusts k0 = %g\n", nucUnb.facK0);
			fprintf(fOut, "Factor that adjusts x = %g\n", nucUnb.facX);
		}
	}
  }

  if (gTglBead != 0) {
	fprintf(fOut, "\n============================== Conditions for a bead "
			"===========================\n");
	fprintf(fOut, "Bead = ");
	if (dir2D < 0) {
		fprintf(fOut, "3-D sphere ");
	}
	else {
		fprintf(fOut, "2-D cylindrical shell (axis = %c) ", direc[dir2D]);
	}
	fprintf(fOut, "with the radius of %g [um] (%g)\n", bead.rad2, 
			L_UM2S(bead.rad2));
	fprintf(fOut, "Number of beads = %d\n", bead.n);
	fprintf(fOut, "Strength of repulsive forces between actin filaments = %g "
			"[N/m] (%g)\n",	STF_REP_BEAD, KS_NPM2S(STF_REP_BEAD));
	fprintf(fOut, "Factor that adjusts the strength = %g\n", bead.facStfRep);
	fprintf(fOut, "Relative thickness to the bead radius that defines hardness"
			" of the repulsive force = %g\n", bead.thkRep);
	fprintf(fOut, "Type of microrheology = %s\n", (bead.rheoWay == 0) 
			? "passive" : "active");
	if (bead.rheoWay == 1) {
		if (bead.rheoSig == 0) { 
    		fprintf(fOut, "Direction of applied force = %c\n", 
					direc[bead.rheoDir]);
    		fprintf(fOut, "Amount of applied pre-force = %g [pN] (%g)\n", 
					F_S2PN(bead.rheoPreMag), bead.rheoPreMag);
    		fprintf(fOut, "Rate of applied pre-force = %g [pN/s]\n", 
					bead.rheoPreRate);
    		fprintf(fOut, "Amplitude of applied force = %g [pN] (%g)\n", 
					bead.rheoAmp, F_PN2S(bead.rheoAmp));
		}
		else {
    		fprintf(fOut, "Direction of applied displacement = %c\n", 
					direc[bead.rheoDir]);
    		fprintf(fOut, "Amount of applied pre-displacement = %g [nm] (%g)\n",
					L_S2NM(bead.rheoPreMag), bead.rheoPreMag);
    		fprintf(fOut, "Rate of applied pre-displacement = %g [nm/s]\n", 
					bead.rheoPreRate);
    		fprintf(fOut, "Amplitude of applied displacement = %g [nm] (%g)\n", 
					bead.rheoAmp, L_NM2S(bead.rheoAmp));
		}
	}
  }

  if (nAct > 0 || nAbp > 0) {
	fprintf(fOut, "\n========================= Stiffness of force models "
			"============================\n");
  }
  if (nAbp > 0) {
  fprintf(fOut, "Strength of repulsive forces between ABPs and between ABP "
		"and actin = %g [N/m] (%g)\n",	KS_S2NPM(abpF.rep.stf), 
		abpF.rep.stf);
  }
  if (nAct > 0) {
	fprintf(fOut, "------------------------------------- actin "
			"------------------------------------\n");
	fprintf(fOut, "Strength of repulsive forces between actin filaments = %g "
			"[N/m] (%g)\n",	KS_S2NPM(actF.rep.stf), actF.rep.stf);
	fprintf(fOut, "Bending stiffness of actin filaments = %g [Nm] (%g)\n",
			KB_S2NM(actF.bend.stf), actF.bend.stf);
	fprintf(fOut, "Persistence length of actin filaments = %g [um]\n", 
			KB_S2NM(actF.bend.stf) * L_SCALE_IN_M / KT_IN_J * 1e6);
	fprintf(fOut, "Extensional stiffness of actin filaments = %g [N/m] (%g)\n",
			KS_S2NPM(actF.spr.stf), actF.spr.stf);
  }
  if (nAbpDet[0] > 0) {
	fprintf(fOut, "------------------------------------- ACP^C "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^C (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[0].stf), abpF.bend[0].stf);
	fprintf(fOut, "Bending stiffness of ACP^C (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[0].stf), abpF.a90[0].stf);
	fprintf(fOut, "Extensional stiffness of ACP^C = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[0].stf), abpF.spr[0].stf);
  }
  if (nAbpDet[1] > 0) {
	fprintf(fOut, "------------------------------------- ACP^B "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of ACP^B (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[1].stf), abpF.bend[1].stf);
	fprintf(fOut, "Bending stiffness of ACP^B (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[1].stf), abpF.a90[1].stf);
	fprintf(fOut, "Extensional stiffness of ACP^B = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[1].stf), abpF.spr[1].stf);
  }
  if (nMot > 0) {
	fprintf(fOut, "------------------------------------- motor "
			"------------------------------------\n");
	fprintf(fOut, "Bending stiffness of motor (theta 1) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.bend[2].stf), abpF.bend[2].stf);
	fprintf(fOut, "Bending stiffness of motor (theta 2) = %g [Nm] (%g)\n",
			KB_S2NM(abpF.a90[2].stf), abpF.a90[2].stf);
	fprintf(fOut, "Extensional stiffness 1 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf), abpF.spr[2].stf);
	fprintf(fOut, "Extensional stiffness 2 of motor = %g [N/m] (%g)\n",
			KS_S2NPM(abpF.spr[2].stf2), abpF.spr[2].stf2);
	if (motSA.gTgl != 0) {
	  	fprintf(fOut, "Bending stiffness of motor backbone = %g [Nm] (%g)\n",
				KB_S2NM(motSA.bend.stf), motSA.bend.stf);
	  	fprintf(fOut, "Extensional stiffness of motor backbone = %g [N/m] "
				"(%g)\n", KS_S2NPM(motSA.spr.stf), motSA.spr.stf);
	}
  }
  fprintf(fOut, "\n============================ Parallel processing "
			"===============================\n");
  fprintf(fOut, "Number of cores used for the job = %d\n", nCpu);
  fprintf(fOut, "Number of subdomains: %d in x, %d in y, %d in z\n", nCell[0], 
		nCell[1], nCell[2]);
  fprintf(fOut, "Method for communication between subdomains = %s\n",
		(mpiMethod == 0) ? "normal" : "Plympton");
  fprintf(fOut, "Initial Subdomain size: x = %g [um] (%g), y = %g "
		"[um] (%g), z = %g [um] (%g)\n", L_S2UM(dimDomC[0]), 
		dimDomC[0], L_S2UM(dimDomC[1]), dimDomC[1], L_S2UM(dimDomC[2]), 
		dimDomC[2]);
  RecordInitParameterSubroutine2("Adjust the size of subdomains to make CPU "
		"loads even", &updSubdSize, fOut);
  fprintf(fOut, "Thickness of overlapping regions = %g [um] (%g)\n",
		L_S2UM(neiEdge), neiEdge);
  FOR_NDIM(k) {
	chk[k][0] = (nCell[k] > 1) ? 1 : 0; 
	chk[k][1] = (nCell[(k+1)%NDIM] > 1 && nCell[(k+2)%NDIM] > 1) ? 1 : 0;
  }
  chk2 = (nCell[0] > 1 && nCell[1] > 1 && nCell[2] > 1)  ? 1 : 0;
  fprintf(fOut, "The ratio of an overlapping region per to a subdomain = %g "
		"[percent]\n", (2. * neiEdge * (dimDomC[0] * dimDomC[1] * chk[2][0]
		+ dimDomC[1] * dimDomC[2] * chk[0][0] + dimDomC[2] * dimDomC[0] * 
		chk[1][0]) - 4. * SQR(neiEdge) * (dimDomC[0] * chk[0][1] + dimDomC[1] 
		* chk[1][1] + dimDomC[2] * chk[2][1]) + 8. * CUBE(neiEdge) * chk2) 
		/ V3PROD(dimDomC) * 100.);	

  fprintf(fOut, "\n================================= Data loading "
		"=================================\n");
  fprintf(fOut, "Load network data from Config = %s\n",
		(gTglLoadNetData == 1 ? "on": "off"));
  fprintf(fOut, "Load information of actins which are previously fixed and "
		"clamped = %s\n", (gTglLoadNetDataFix == 1 ? "on": "off"));

  fprintf(fOut, "\n=============================== Data recording "
		"=================================\n");
  fprintf(fOut, "Directory for saving data = %s\n", dataFold);
  fprintf(fOut, "Deletion of pre-existing data files = %s\n", 
		(DELETE_FILE != 0 ? "on" : "off"));
  RecordInitParameterSubroutine2("Record Output and Progress", &recProg, fOut);

  fprintf(fOut, "----------------------------- Structural information "
		"---------------------------\n");
  RecordInitParameterSubroutine2("Record Config", &recConf, fOut);
  fprintf(fOut, "Recording structural information for visualization via MATLAB "
        "= %s\n", (recConfMlb.tgl != 0) ? "on" : "off");
  if (recConfMlb.tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", recConfMlb.prdR, 
			recConfMlb.prd);
  }
  fprintf(fOut, "Recording structural information for visualization via VMD "
		"= %s\n", (recConfVmd.tgl != 0) ? "on" : "off");
  if (recConfVmd.tgl != 0) {
	fprintf(fOut, "- Period = %g [s] (%d timesteps)\n", recConfVmd.prdR, 
			recConfVmd.prd);
	if (recConfVmd.mode == 0) { fprintf(fOut, "- Multiple files\n"); }
	else { fprintf(fOut, "- Single file\n"); }
	fprintf(fOut, "- Show boundaries in the network drawn via VMD = %s\n",
			(recConfVmd.gTglBnd != 0 ? "on" : "off"));
	fprintf(fOut, "- Record information for coloring a network in VMD = %s\n",
			(recConfVmd.gTglInfo != 0 ? "on" : "off"));
  }
  RecordInitParameterSubroutine2("Record the length of actin filaments", 
		&recFilaL, fOut);
  if (motSA.gTgl != 0) {
	  RecordInitParameterSubroutine2("Record size of self-assembled motors", 
			&recMotSize, fOut);
	  RecordInitParameterSubroutine2("Record position of self-assembled motors",
			&recMotPos, fOut);
  }
  RecordInitParameterSubroutine2("Record distances between active ABPs", 
		&recCrsDist, fOut);
  RecordInitParameterSubroutine2("Record pore size of networks",
		&recPoreSize, fOut);
  RecordInitParameterSubroutine2("Record connectivity between actin filaments",
		&recConn, fOut);
  RecordInitParameterSubroutine2("Record the percolation of a network",
		&recPerc, fOut);
  RecordInitParameterSubroutine2("Find a supportive framework",
		&findSupp, fOut);
  if (findSupp.tgl != 0) {
	fprintf(fOut, "- Way to find the supportive framework = %g percent of "
			"high ABP %s forces\n", porFindSupp * 100., 
			(kindFindSupp == 0 ? "bending" : "extensional"));
  }

  fprintf(fOut, "--------------------------- Force, stress, and energy "
		"--------------------------\n");
  RecordInitParameterSubroutine2("Record longitudinal forces acting on ABPs", 
		&recLongF, fOut);
  RecordInitParameterSubroutine2("Record the mechanical energy of a network", 
		&recE, fOut);
  RecordInitParameterSubroutine2("Record internal elastic/viscous stresses",
		&recSecStre, fOut);
  if (recSecStre.tgl != 0) {
	fprintf(fOut, "- The number of measurements for the internal stresses = %d "
			"in x, %d in y, %d in z\n", nSecStreDiv[0], nSecStreDiv[1], 
			nSecStreDiv[2]);
  }

  fprintf(fOut, "--------------------- Dynamic behaviors of actin and ABPs "
		"----------------------\n");
  RecordInitParameterSubroutine2("Record actin severing", &recActSev, fOut);
  RecordInitParameterSubroutine2("Record ABP unbinding", &recAbpUnb, fOut);
  RecordInitParameterSubroutine2("Record ABP binding", &recAbpBind, fOut);
  RecordInitParameterSubroutine2("Record ABP turnover", &recAbpTurn, fOut);
  RecordInitParameterSubroutine2("Record how many times each dynamic event "
		"of ABPs occurs", &recAbpDyn, fOut);

  if (gTglMb != 0) {
	fprintf(fOut, "------------------------------- Membrane/nucleus "
			"-------------------------------\n");
	if (gTglNuc != 0) {
	  	RecordInitParameterSubroutine2("Record the center of membranes and"
			" nuclei", &recMbCen, fOut);
	  	RecordInitParameterSubroutine2("Record the radius of membranes and"
			" nuclei", &recMbDim, fOut);
	}
	else {
	  	RecordInitParameterSubroutine2("Record the center of membranes",
			&recMbCen, fOut);
	  	RecordInitParameterSubroutine2("Record the radius of membranes",
			&recMbDim, fOut);
	}
  }

  if (gTglBead != 0) {
	if (recBeadLoc.tgl != 0) {
		fprintf(fOut, "--------------------------------------- Bead "
				"-----------------------------------\n");
		RecordInitParameterSubroutine2("Record the center of beads",
				&recBeadLoc, fOut);
	}
  }

  fprintf(fOut, "---------------------------------- Miscellany "
		"----------------------------------\n");
  RecordInitParameterSubroutine2("Record the information in unit of "
		"filaments", &recInfo, fOut);
  fclose(fOut);
}

/*------------------------- Recording parameter values -----------------------*/

/*---------------------------- Recording counters ----------------------------*/

// Subroutine for RecordCounters()
void RecordCountersSubroutine(int idx, int var, FILE *fOut) {
  int *varAll;

  MALLOC(varAll,int,nCpu);
  MPI_Gather(&var, 1, MPI_INT, varAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	fprintf(fOut, "%d\t", idx);
	Fprintf1dArrayInt(fOut, varAll, nCpu, 0);
  }
  free(varAll);
}

// Record the number of particles or elements in the current subdomain.
// This is useful primarily for finding a reason causing error.
void RecordCounters(void) {
  FILE *fOut;

  if (rank == 0) { 
	fOut = fopen(GenFileName("Variable"), "a"); 
  }
  RecordCountersSubroutine(0, nActMe, fOut);
  RecordCountersSubroutine(1, nAbpMe, fOut);
  RecordCountersSubroutine(2, nActCp, fOut);
  RecordCountersSubroutine(3, nAbpCp, fOut);
  RecordCountersSubroutine(4, nMbMe, fOut);
  RecordCountersSubroutine(5, nMbCp, fOut);
  RecordCountersSubroutine(6, longCh.c, fOut);
  RecordCountersSubroutine(7, sendAbpDyn.c, fOut);
  RecordCountersSubroutine(8, noAbpDyn.c, fOut);
  RecordCountersSubroutine(9, neigh.c, fOut);
  RecordCountersSubroutine(10, neighMb.c, fOut);
  RecordCountersSubroutine(11, sizeBufMsg, fOut);
  RecordCountersSubroutine(12, neigh.siz, fOut);
  RecordCountersSubroutine(13, neighMb.siz, fOut);
  RecordCountersSubroutine(14, sendActDyn.siz, fOut);
  RecordCountersSubroutine(15, sendAbpDyn.siz, fOut);
  RecordCountersSubroutine(16, sendMbDyn.siz, fOut);
  RecordCountersSubroutine(17, longChExtMsg.siz, fOut);
  RecordCountersSubroutine(18, longChIntMsg.siz, fOut);
  int n, max = -100000;
  for(n = 0; n < cntAdjRank[0]; n++) {
	if (cpPar[n].c > max) { max = cpPar[n].c; }
  }
  RecordCountersSubroutine(19, max, fOut);
  RecordCountersSubroutine(20, memb.unit.c, fOut);
  RecordCountersSubroutine(21, memb.unitCp.c, fOut);
  if (rank == 0) { fclose(fOut); }
}

/*---------------------------- Recording counters ----------------------------*/

/*-------------------------- Recording network structure ---------------------*/

// Record network configuration. Most of the data are related to the positions
// and chain information of particles.
void RecordConfig(char *fn, int period) {
  int n, k, *nActAll, *nAbpAll, *nMbAll;
  double *rActOff, *rAbpOff, *rMbOff;
  ListInt parL;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(nMbAll,int,nCpu);
  MALLOC(rActOff,double,nActC*NDIM);
  MALLOC(rAbpOff,double,nAbpC*NDIM);
  MALLOC(rMbOff,double,nMbC*NDIM);

  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  MPI_Gather(&nMbMe, 1, MPI_INT, nMbAll, 1, MPI_INT, 0,
        MPI_COMM_WORLD);
  if (rank == 0) {
	if (SumArrInt(nActAll, nCpu) != nAct) {
		Printf("Error: the sum of collected nActMe is different from "
				"nAct!\n");
		exit(-1);
	}
	if (SumArrInt(nAbpAll, nCpu) != nAbp) {
		Printf("Error: the sum of collected nAbpMe is different from "
				"nAbp!\n");
		exit(-1);
	}
	
	fOut = fopen(fn, ((period > 0) ? "a" : "w")); 
	if (period > 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", (int)(currTimeStep / period));
	}
	fprintf(fOut, "nAct : %d\n", nAct);
	fprintf(fOut, "nAbp : %d\n", nAbp);
	fprintf(fOut, "nMb : %d\n", nMb);
	fprintf(fOut, "dimDom : %g, %g, %g\n", dimDom[0], 
			dimDom[1], dimDom[2]);
	fprintf(fOut, "nActPerSeg, nChAcX, nChAcY : %d, %d, %d\n\n", nActPerSeg, 
			nChAcX, nChAcY);
  }
  // Offset the positions of particles. This can be recovered later using the
  // information of rGrid below.
  FOR_ACTME(n) { 
	FOR_NDIM(k) {
		P2(rActOff,n,k) = P2(act.r,n,k) - rGrid[k][0];
	}
  }
  FOR_ABPME(n) { 
	FOR_NDIM(k) {
		P2(rAbpOff,n,k) = P2(abp.r,n,k) - rGrid[k][0];
	}
  }
  FOR_MBME(n) { 
	FOR_NDIM(k) {
		P2(rMbOff,n,k) = P2(memb.r,n,k) - rGrid[k][0];
	}
  }
  // Positions and chain information
  if (rank == 0) { fprintf(fOut, "## Position for actin ##\n"); }
  RecordGather2dArrayDouble(nActMe, nActAll, 3, act.id, rActOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Position for ABP ##\n"); }
  RecordGather2dArrayDouble(nAbpMe, nAbpAll, 3, abp.id, rAbpOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for actin ##\n"); }
  RecordGather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for ABP ##\n"); }
  RecordGather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Position for membrane ##\n"); }
  RecordGather2dArrayDouble(nMbMe, nMbAll, 3, memb.id, rMbOff, fOut, 0);
  if (rank == 0) { fprintf(fOut, "## Chain for membrane ##\n"); }
  RecordGather2dArrayInt(nMbMe, nMbAll, nChMb, memb.id, memb.ch, fOut, 0);
  if (rank == 0) {
	fprintf(fOut, "currTimeStep = %lld\n", currTimeStep);
	fprintf(fOut, "rGrid = %g, %g, %g, %g, %g, %g\n", rGrid[0][0], rGrid[1][0],
			rGrid[2][0], rGrid[0][nGrid[0] - 1], rGrid[1][nGrid[1] - 1], 
			rGrid[2][nGrid[2] - 1]);
	fprintf(fOut, "pbc = %d, %d, %d\n", pbc[0], pbc[1], pbc[2]);
	fprintf(fOut, "rheoWay = %d\n", rheoWay);
	fprintf(fOut, "bulkRheoWay = %d\n", bulkRheoWay);
	fprintf(fOut, "bulkRheoType = %d\n\n", bulkRheoType);
  }
  // Related to rheological measurements
  if (rheoWay > 0) {
	if (rank == 0) { 
		fprintf(fOut, "stra.accDsp = %g\n\n", stra.accDsp); 
	}
	RecordGather1dArrayIntWoIndWoCnt(meaStreParMe.c, meaStreParMe.l, 
			fOut, "meaStrePar");
	RecordGather1dArrayIntWoIndWoCnt(appStraParMe.c, appStraParMe.l, 
			fOut, "appStraPar");
  }
  else {
	if (rank == 0) { 
		fprintf(fOut, "stra.accDsp = 0\n\n");
		fprintf(fOut, "meaStrePar = 0\n\n");
		fprintf(fOut, "appStraPar = 0\n\n");
	}
  }
  if (recTraj.tgl != 0) {
	  RecordGather1dArrayIntWoIndWoCnt(recTraj.actMe.c, recTraj.actMe.l, 
			fOut, "actTraj");
  }
  else { 
	if (rank == 0) {
		fprintf(fOut, "actTraj = 0\n\n");
	}
  }
  if (recTraj.tgl2 != 0) {
	RecordGather1dArrayIntWoIndWoCnt(recTraj.abpMe.c, recTraj.abpMe.l, 
			fOut, "abpTraj");
  }
  else { 
	if (rank == 0) {
		fprintf(fOut, "abpTraj = 0\n\n");
	}
  }

  if (rank == 0) { 
	MALLOC(fixAct,int,nAct);
	if (motSA.gTgl != 0) { MALLOC(abpMotId,int,nAbp); }
  }
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  if (motSA.gTgl != 0) {
	Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.mId, abpMotId);
  }
  if (rank == 0) {
	MALLOC(parL.l,int,nAct);
	parL.c = 0;
	FOR_ACT(n) {
		CONT(fixAct[n] < 0); 
		InsertElement1dArrayWoChk(parL.l, &parL.c, n);
	} 
	fprintf(fOut, "fixAct = %d\n", parL.c);
	for(n = 0; n < parL.c; n++) {
		fprintf(fOut, "%d\t%d\n", parL.l[n], fixAct[parL.l[n]]);
	} 
	free(parL.l);
	fprintf(fOut, "\n");

	parL.c = 0;
	if (motSA.gTgl != 0) {
		MALLOC(parL.l,int,nAbp);
		FOR_ABP(n) {
			CONT(abpMotId[n] < 0); 
			InsertElement1dArrayWoChk(parL.l, &parL.c, n);
		} 
	}
	fprintf(fOut, "abp.mId = %d\n", parL.c);
	for(n = 0; n < parL.c; n++) {
		fprintf(fOut, "%d\t%d\n", parL.l[n], abpMotId[parL.l[n]]);
	} 
	free(fixAct);
	if (motSA.gTgl != 0) { 
		free(parL.l);
		free(abpMotId); 
	}
	fclose(fOut); 
  }
  free(rActOff);
  free(rAbpOff);
  free(rMbOff);
  free(nActAll);
  free(nAbpAll);
  free(nMbAll);
}

int RecordConfigVmdSubroutine(double *r1, double *r2, int *ind, int *mul) {
  int k, adjInd, fac[NDIM];
  double dr[NDIM];

  V3COPY(fac, ind);
  V3SUB(dr, r1, r2);
  FOR_NDIM(k) {
	CONT(!(fabs(dr[k]) > dimDomH[k]));
	fac[k] += (dr[k] > dimDomH[k] ? 1 : -1);
	if (fac[k] >= mul[k]) { fac[k] = 0; }
	else if (fac[k] < 0) { fac[k] = mul[k] - 1; } 
  }
  V3IND_BACK_INT(adjInd, fac, mul);

  return adjInd;
}

int RecordConfigInformation(double *recActG, double *recAbpG, double *recMbG,
		int *nActAll, int *nAbpAll, int *nMbAll, int nD) {
  int m, n, k, nRec, ch;
  double *sendRecAct, *sendRecAbp, *sendRecMb, *pntRec;

  MALLOC(sendRecAct,double,nActC * nD);
  MALLOC(sendRecAbp,double,nAbpC * nD);
  if (gTglMb != 0) {
	  MALLOC(sendRecMb,double,nMbC * nD);
  }
  FOR_ACTME(n) { 
	ch = BinaryPackActinChainArray(act.id[n]);
	V6SET(&P2A(sendRecAct,n,0,nD), recAct.len[n], 
			P2A(recAct.allF,n,NDIM,NDIM + 1), recAct.sprF[n], 
			recAct.bendF[n], ch, (double)recAct.cnt[n]);
  }
  FOR_ABPME(n) { 
	ch = BinaryPackAbpChainArray(abp.id[n]);
	V6SET(&P2A(sendRecAbp,n,0,nD), P2A(recAbp.len,n,0,recAbp.nL), 
			P2A(recAbp.allF,n,NDIM,NDIM + 1), recAbp.sprF[n], 
			recAbp.bendF[n], ch, (double)recAbp.cnt[n]);
  }
  if (gTglMb != 0) {
	FOR_MBME(n) { 
		V6SET(&P2A(sendRecMb,n,0,nD), recMb.len[n], 
				P2A(recMb.allF,n,NDIM,NDIM + 1), recMb.sprF[n], 
				recMb.bendF[n], 0, (double)recMb.cnt[n]);
	}
  }
  Gather2dArrayDouble(nActMe, nActAll, nD, act.id, sendRecAct, recActG);	
  Gather2dArrayDouble(nAbpMe, nAbpAll, nD, abp.id, sendRecAbp, recAbpG);	
  if (gTglMb != 0) {
	Gather2dArrayDouble(nMbMe, nMbAll, nD, memb.id, sendRecMb, recMbG);
  }
  if (rank == 0) {
	for(m = 0; m < (gTglMb != 0 ? 3 : 2); m++) {
		if (m == 0) {
			nRec = nAct;
			pntRec = recActG;
		}
		else if (m == 1) {
			nRec = nAbp;
			pntRec = recAbpG;
		}
		else {
			nRec = nMb;
			pntRec = recMbG;
		}	
		for(n = 0; n < nRec; n++) {
			for(k = 0; k < 4; k++) {
				if (P2A(pntRec,n,nD - 1,nD) == 0) { 
					V4SET_ALL(&P2A(pntRec,n,0,nD), 0.);
				}
				else {
					P2A(pntRec,n,k,nD) /= P2A(pntRec,n,nD - 1,nD);
	 				if (k > 0) {
						P2A(pntRec,n,k,nD) = F_S2PN(P2A(pntRec,n,k,nD));
					}
				}
			}
		}
	}
  }
  free(sendRecAct);
  free(sendRecAbp);
  if (gTglMb != 0) {
	free(sendRecMb);
  }
  return nD;
}

// Record configuration for visualization via VMD. Positions and chain 
// information are recorded in the PDB and PSF formats. Actin, ACP, and motor
// have different SEGNAMEs for distinction. 
// mode = 0: bonds between actins and ABPs are drawn to the end points of 
//           actin segments. 
//        1: actin segments are divided to sub-segments, so the bonds are
//           drawn to actual points. This will cost more memory in VMD.
void RecordConfigVmd(int mode) {
  int m, n, k, ind, ind2[NDIM], ind3, ind4, ind5, cnt, curr, amp = 1;
  int nRec, nD, ch, *pArr, *pArr2, loc[2], nCh, indMb;
  int nActVmd, nAbpVmd, mul[NDIM];
  int facInfoConfVmd, facBndConfVmd;
  int *nActAll, *nAbpAll, *chActVmd, *chAbpVmd, *iFila, *chkBond;
  int nMbVmd, *nMbAll, *chMbVmd, *rebMb, *fixMb, *chActAlt;
  int nBeadVmd;
  int bndCh[] = {0, 4, 1, 5, 2, 6, 3, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 
				2, 1, 3, 4, 6, 5, 7};
  double temp, dr[NDIM], domCen[NDIM], rPos[NDIM], sft[NDIM];
  double *rActVmd, *rAbpVmd, *occu, *mass;
  double *recActG, *recAbpG, *recMbG; //, *sendRecAct, *sendRecAbp, *sendRecMb;
  double *pntRec;
  double *rMbVmd, *rBeadVmd;
  char **segType, fn[80], fnIn[80];
  FILE *fOut, *fOut2;
  ListInt chain;

int *sldAbp;

  if (recConfVmd.mode == 0) {
	sprintf(fnIn, "ConfVmd_%d", (int)(currTimeStep / recConfVmd.prd));
  }
  else {
	sprintf(fnIn, "ConfVmd");
  }

  facInfoConfVmd = (recConfVmd.gTglInfo > 0) ? 1 : 0;
  facBndConfVmd = (recConfVmd.gTglBnd != 0) ? 8 : 0;

  MALLOC(chActAlt,int,nChAc);
  // If a direction has periodic boundary condition, the configuration is
  // automatically duplicated in the direction for a clearer view.
  FOR_NDIM(k) {
	mul[k] = (pbc[k] == 1) ? 2 : 1;
  }
  V3SET_ALL(mul, 1);

  nActVmd = nAct * V3PROD(mul) * ((mode != 0) ? nChAcX : 1);
  nAbpVmd = nAbp * V3PROD(mul);
  nMbVmd = (gTglMb != 0) ? nMb * V3PROD(mul) : 0;
  nBeadVmd = (gTglBead != 0) ? bead.n * V3PROD(mul) : 0;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (gTglMb != 0) { 
	  MALLOC(nMbAll,int,nCpu);
  }

  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(iFila,int,nAct);
	if (gTglMb != 0) {
		MALLOC(rMb,double,nMb*NDIM);
		MALLOC(chMb,int,nMb*nChMb);
		if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
			MALLOC(rebMb,int,nMb);
		}
		if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
			MALLOC(fixMb,int,nMb);
		}
		MALLOC(idxMb,int,nMb);
    	if (mbSld.abp.gTgl != 0) {
			MALLOC(sldAbp,int,nAbp * dimMbNuc);
		}
	}
  }

  // Gather the positions and chain information of particles.
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);  
  if (gTglMb != 0) {
	GatherMembChainPosition(nMbAll, chMb, rMb);
	if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
		Gather2dArrayInt(nMbMe, nMbAll, 1, memb.id, memb.reb, rebMb);  
	}
	if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
		Gather2dArrayInt(nMbMe, nMbAll, 1, memb.id, memb.fix, fixMb);  
	}
	Gather2dArrayInt(nMbMe, nMbAll, 1, memb.id, memb.idx, idxMb);  
	indMb = (dimMbNuc == 2) ? 2 : 1;
    if (mbSld.abp.gTgl != 0) {
//		Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, mbSld.abp.l, sldAbp);  
		Gather2dArrayInt(nAbpMe, nAbpAll, dimMbNuc, abp.id, mbSld.abp.l, sldAbp);  
    }
  }

if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
  if (rank == 0) {
	fOut2 = fopen(GenFileName("MembBndAbp"), "a");
	cnt = 0;
	FOR_ABP(n) {
//		CONT(sldAbp[n] == -1);
		CONT(P2A(sldAbp,n,0,dimMbNuc) == -1);
		cnt++;
	}
	fprintf(fOut2, "%d\t%lld\t0\t0\t0\n", cnt, currTimeStep);
	cnt = 0;
	FOR_ABP(n) {
		CONT(P2A(sldAbp,n,0,dimMbNuc) == -1);
		fprintf(fOut2, "%d\t%d\t%g\t%g\t%g\n", cnt, n, 
				P2(rAbp,n,0), P2(rAbp,n,1), P2(rAbp,n,2));
		for(k = 0; k < dimMbNuc; k++) {
			ind = P2A(sldAbp,n,k,dimMbNuc);
			fprintf(fOut2, "%d\t%d\t%g\t%g\t%g\n", k, ind, 
					P2(rMb,ind,0), P2(rMb,ind,1), P2(rMb,ind,2));
		}
		cnt++;
	}
	fclose(fOut2);
  }
}
  
  if (recConfVmd.gTglInfo > 0) {
	nD = 6;
	if (rank == 0) {
		MALLOC(recActG,double,nAct * nD);
		MALLOC(recAbpG,double,nAbp * nD);
		if (gTglMb != 0) {
			MALLOC(recMbG,double,nMb * nD);
		}
	}
	RecordConfigInformation(recActG, recAbpG, recMbG, 
			nActAll, nAbpAll, nMbAll, nD);
  }

  // Prepare the (duplicated) positions and chain information of particles.
  if (rank == 0) {
	MALLOC(rActVmd, double, nActVmd * NDIM);
	MALLOC(rAbpVmd, double, nAbpVmd * NDIM);
	MALLOC(chActVmd, int, nActVmd * ((mode != 0) ? nChAcY + 2 : nChAc));
	MALLOC(chAbpVmd, int, nAbpVmd * nChAb);
	MALLOC(occu, double, nAct + nAbp + nMb);
	MALLOC(mass, double, nAct + nAbp + nMb);
	MALLOC2(segType, char, nActVmd + nAbp); 
	MALLOC(chain.l, int, (nActVmd + nAbpVmd * 4 + nMbVmd 
			* ((dimMbNuc == 2) ? 2 : 6) + ((gTglMb != 0) ? nActVmd : 0) 
			+ ((recConfVmd.gTglBnd != 0) ? 12 : 0)) * 2);
	for (k = 0; k < nActVmd + nAbp; k++) {
		MALLOC(segType[k], char, 4);
	}
	memset(chActVmd, -1, sizeof(int) * nActVmd 
			* ((mode != 0) ? nChAcY + 2 : nChAc));
	memset(chAbpVmd, -1, sizeof(int) * nAbpVmd * nChAb);
	if (gTglMb != 0) {
		MALLOC(rMbVmd, double, nMbVmd*NDIM);
		MALLOC(chMbVmd, int, nMbVmd*nChMb);
		memset(chMbVmd, -1, sizeof(int)*nMbVmd*nChMb);
		if (dimMbNuc == 3) {
			MALLOC(chkBond, int, nMbVmd*6);
			memset(chkBond, 0, sizeof(int)*nMbVmd*6);
		}
	}
	if (gTglBead != 0) {
		MALLOC(rBeadVmd, double, nBeadVmd*NDIM);
	}
	// Copy position and chain information to multiplied array
	for(ind = 0; ind < V3PROD(mul); ind++) {
		V3IND_ASSIGN_INT(ind, mul, 1, ind2);
		FOR_ACT(n) {
			pArr = &P2A(chAct,n,0,nChAc);
			memset(chActAlt, -1, sizeof(int) * nChAc);
			for (m = 0; m < nChAc; m++) {
				ind5 = pArr[m];
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAct,n,0), (m < 2 ? 
						&P2(rAct,ind5,0) : &P2(rAbp,ind5,0)), ind2, mul);
				chActAlt[m] = ind5 + ind4 * ((m < 2) ? nAct : nAbp);
			}
			if (mode == 0) {
				ind3 = n + ind * nAct;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				for (m = 0; m < nChAc; m++) {
					P2A(chActVmd,ind3,m,nChAc) = chActAlt[m];
				}
			}
			else {
				ind3 = (n + ind * nAct) * nChAcX;
				VSV3ADD(&P2(rActVmd,ind3,0), &P2(rAct,n,0), dimDom, ind2);	
				pArr2 = &P2A(chActVmd,ind3,0,nChAcY + 2);
				if (pArr[1] > -1) { 
					pArr2[1] = (chActAlt[1] + ind4 * nAct + 1) * nChAcX - 1;
				}
				if (pArr[0] > -1) {
					for(m = 1; m < nChAcX; m++) {
						// Position
						CalcPosOnActSeg(&P2(rAct,n,0), &P2(rAct,pArr[0],0),
								rPos, (double)m / (double)nChAcX);
						VSV3ADD(&P2(rActVmd,ind3 + m,0), rPos, dimDom, ind2);
						// Chain
						V2SET(&P2A(pArr2,m,0,nChAcY + 2), ind3 + m + 1, 
								ind3 + m - 1);
					}
					// Chain at ends
					pArr2[0] = ind3 + 1;
					P2A(pArr2,nChAcX - 1,0,nChAcY + 2) 
							= (chActAlt[0] + ind4 * nAct) * nChAcX;

					for(m = 2; m < nChAc; m++) {
						CONT(pArr[m] < 0);
						loc[0] = (m - 2) / nChAcY;
						loc[1] = (m - 2) % nChAcY;
						P2A(pArr2,loc[0],loc[1] + 2,nChAcY + 2) 
								= chActAlt[m] + ind4 * nAbp;
					}
				}
				else {
					for(m = 1; m < nChAcX; m++) {
						V3SET_ALL(&P2(rActVmd,ind3 + m,0), 0.);
					}
				}

			}
		}
		FOR_ABP(n) {
			ind3 = n + ind * nAbp;
			VSV3ADD(&P2(rAbpVmd,ind3,0),&P2(rAbp,n,0),dimDom,ind2);	
			for (m = 0; m < 2; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAct,ind5,0), ind2, mul);
				if (mode == 0) { 
					P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAct;
				}
				else {
					loc[0] = FindAbpActinChain(ind5, n, 1);
					loc[0] = (loc[0] - 2) / nChAcY;
					P2A(chAbpVmd,ind3,m,nChAb) 
							= (ind4 * nAct + ind5) * nChAcX + loc[0];
				}
			}
			P2A(chAbpVmd,ind3,2,nChAb) = P2A(chAbp,n,2,nChAb);
			for (m = 3; m < 4; m++) {
				ind5 = P2A(chAbp,n,m,nChAb);
				CONT(ind5 < 0);
				ind4 = RecordConfigVmdSubroutine(&P2(rAbp,n,0), 
						&P2(rAbp,ind5,0), ind2, mul);
				P2A(chAbpVmd,ind3,m,nChAb) = ind5 + ind4 * nAbp;
			}
		}
		FOR_MB(n) {
			ind3 = n + ind * nMb;
			VSV3ADD(&P2(rMbVmd,ind3,0),&P2(rMb,n,0),dimDom,ind2);	
			for (m = 0; m < ((dimMbNuc == 2) ? 2 : 6); m++) {
				CONT(!(P2A(chMb,n,m,nChMb) > -1));
				ind4 = RecordConfigVmdSubroutine(&P2(rMb,n,0), 
						&P2(rMb,P2A(chMb,n,m,nChMb),0), ind2, mul);
				P2A(chMbVmd,ind3,m,nChMb) = 
						P2A(chMb,n,m,nChMb) + ind4 * nMb;
			}
		}
		for(n = 0; n < bead.n; n++) {
			ind3 = n + ind * bead.n;
			VSV3ADD(&P2(rBeadVmd,ind3,0),&P2(bead.r,n,0),dimDom,ind2);	
		}
	}
	// Fill in segment type of actin, ACP, and motor.
	for(n = 0; n < nActVmd; n++) { 
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		pArr = &P2A(chAct,ind,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) {
			sprintf(segType[n], "G");
		}
		else if (pArr[0] > -1 && pArr[1] < 0) {
			sprintf(segType[n], "FP");
			if (mode == 1) {
				if (P2A(chActVmd,n,1,nChAcY + 2) > -1) {
					sprintf(segType[n], "F");
				}
			}
		}
		else if (pArr[0] < 0 && pArr[1] > -1) {
			sprintf(segType[n], "FB");
			if (mode == 1) {
				if (n % nChAcX > 0) {
					sprintf(segType[n], "G");
				}
			}
		}
		else { sprintf(segType[n], "F"); }
	}
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		cnt = 1;
		if (pArr[0] < 0) { cnt++; }
		if (pArr[1] < 0) { cnt++; }
		if (motSA.gTgl != 0) {
			if (pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 && pArr[4] < 0) {
				cnt = 4;
			}
		}
		switch(pArr[2]) {
			case 0: sprintf(segType[n + nActVmd], "CR%d", cnt);	break;
			case 1: sprintf(segType[n + nActVmd], "BU%d", cnt);	break;
			case 2: sprintf(segType[n + nActVmd], "MO%d", cnt);	break;
		}
	}
	
	if (recConfVmd.gTglInfo > 0 && currTimeStep > 1) {
		for(m = 0; m < 2; m++) {
			nRec = (m == 0) ? nAct : nAbp;
			pntRec = (m == 0) ? recActG : recAbpG;
			for(n = 0; n < nRec; n++) {
				temp = TrimDblVal(P2A(pntRec,n,2,nD), recConfVmd.minF, 
						recConfVmd.maxF);
				mass[n + ((m == 1) ? nAct : 0)] 
						= 1. + 20. * (temp  - recConfVmd.minF) 
						/ (recConfVmd.maxF - recConfVmd.minF);
			}
		}
		for(n = 0; n < nMb; n++) {
			temp = TrimDblVal(P2A(recMbG,n,2,nD), recConfVmd.minF, 
					recConfVmd.maxF);
			mass[n + nAct + nAbp] 
					= 1. + 20. * (temp  - recConfVmd.minF) 
					/ (recConfVmd.maxF - recConfVmd.minF);
		}
	}
	else {
		for(n = 0; n < nAct; n++) {
			mass[n] = 12.;
		}
		for(n = 0; n < nAbp; n++) {
			mass[n + nAct] = 14.;
		}
		for(n = 0; n < nMb; n++) {
			mass[n + nAct + nAbp] = 13.;
		}
	}
	
	for(n = 0; n < nAct + nAbp + nMb; n++) { occu[n] = 1.; }

if (gTglMb != 0 && mbSld.abp.gTgl != 0) {
	for(n = 0; n < nAbp; n++) {
		CONT(P2A(sldAbp,n,0,dimMbNuc) == -1);
		occu[nAct + n] = 2.;
		for(k = 0; k < dimMbNuc; k++) {
			occu[nAct + nAbp + P2A(sldAbp,n,k,dimMbNuc)] = 2.;
		}
	}
}
	if (recPerc.tgl != 0) {
		for(n = 0; n < nAct + nAbp; n++) {
			if (iPerc[n] / 10 == 1) { occu[n] = 2.; }
			if (iPerc[n] % 10 == 1) { occu[n] = 3.; }
		}
	}

	// Write PDB file
	sprintf(fn, "%s.pdb", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf (fOut, "REMARK\n");

	FOR_NDIM(k) {	
  		if (bndMv.gTgl != 0 || (rheoWay > 0 && bulkRheoType == 1)) {
			domCen[k] = 0.5 * ((2 - mul[k]) * P2A(rGridInit,0,k,NDIM)
					+ mul[k] * P2A(rGridInit,1,k,NDIM));
		}
		else {
			domCen[k] = 0.5 * mul[k] * dimDom[k];
		}
	}
	// actin
	for(n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
	    VV3SUB(&P2(rActVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  C   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + 1, P2(rActVmd,n,0) * amp, 
				P2(rActVmd,n,1) * amp, P2(rActVmd,n,2) * amp, 
				occu[ind], segType[n]);
	}
	// ABP
	for(n = 0; n < nAbpVmd; n++) {
		VV3SUB(&P2(rAbpVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  N   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", n + nActVmd + 1, P2(rAbpVmd,n,0) * amp, 
				P2(rAbpVmd,n,1) * amp, P2(rAbpVmd,n,2) * amp, 
				occu[(n % nAbp) + nAct], segType[(n % nAbp) + nActVmd]);
	}
	// If a coloring method is used, a reference particle is necessary.
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "ATOM %6d  O   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", nActVmd + nAbpVmd + 1, 0., 0., 0., 1., "REF");
	}
	// Visualization of boundaries using cylinders
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			V3IND_ASSIGN_CONST_INT(n, 2, ind2);
			V3SET_ALL(sft, 0.);
			// If it is a bulk rheology with shear, the location of boundary
			// is adjusted..
			if (rheoWay > 0 && bulkRheoType == 0) {
				sft[dirStr] = (ind2[dirNoPBC] == 1) ?  stra.acc 
						* dimDom[dirNoPBC] : 0.;
			}
			FOR_NDIM(k) {
				rPos[k] = rGrid[k][(ind2[k] == 0) ? 0 : nGrid[k] - 1] 
						- domCen[k] + sft[k] 
						+ ((ind2[k] == 1) ? dimDom[k] * (mul[k] - 1) : 0.); 
			}
			fprintf(fOut, "ATOM %6d  H   LYS     1    %8.3f%8.3f%8.3f  %4.2f"
					"  0.00      %s\n", nAbpVmd + nActVmd + n + facInfoConfVmd 
					+ 1, rPos[0], rPos[1], rPos[2], 1., "BND");
		}
	}
	for(n = 0; n < nMbVmd; n++) {
	    VV3SUB(&P2(rMbVmd,n,0), domCen);
		for(k = 0; k < indMb; k++) {
			fprintf(fOut, "ATOM %6d  P   LYS     1    ", nAbpVmd + nActVmd 
					+ facInfoConfVmd + 1 + facBndConfVmd + n * indMb + k);
			for(m = 0; m < NDIM; m++) {
				if (dimMbNuc == 2) { 
					fprintf(fOut, "%8.3f", ((m == dirNormMbNuc) ? domCen[m] 
							* amp * (k * 2 - 1) :  P2(rMbVmd,n,m) * amp));
				}
				else {
					fprintf(fOut, "%8.3f", P2(rMbVmd,n,m) * amp);
				}
			}
			fprintf(fOut, "  %4.2f  0.00      %s\n", occu[n % nMb + nAct + nAbp], 
					(ISNUC(idxMb[n % nMb]) ? "NUC" : "MB")); 
		}
	}
	for(n = 0; n < nBeadVmd; n++) {
	    VV3SUB(&P2(rBeadVmd,n,0), domCen);
		fprintf(fOut, "ATOM %6d  S   LYS     1    %8.3f%8.3f%8.3f  %4.2f  "
				"0.00      %s\n", nActVmd + nAbpVmd + facInfoConfVmd + 1 
				+ facBndConfVmd + nMbVmd * indMb + n, P2(rBeadVmd,n,0) * amp, 
				P2(rBeadVmd,n,1) * amp, P2(rBeadVmd,n,2) * amp, 1., "BD");
	}

	fprintf(fOut, "ENDMDL\n");
	fclose(fOut);

	// Write PSF file
	sprintf(fn, "%s.psf", fnIn);
	fOut = fopen(GenFileName(fn), ((recConfVmd.mode == 0) ? "w" : "a"));
	if (recConfVmd.mode != 0) {
		fprintf(fOut, "------------------------------------ %6d ------------"
				"------------------------\n", 
				(int)(currTimeStep / recConfVmd.prd));
	}
	fprintf(fOut, "PSF CMAP\n\n");
	fprintf(fOut, "%8d !NTITLE\n\n", 7);
	fprintf(fOut, "%8d !NATOM\n", nActVmd + nAbpVmd + nMbVmd 
			* indMb + facInfoConfVmd + facBndConfVmd + nBeadVmd);
	// actin
	for (n = 0; n < nActVmd; n++) {
		ind = (int)(n / ((mode != 0) ? nChAcX : 1)) % nAct;
		fprintf(fOut, "%8d %-4s 1    LYS  C      11%11.5f%14.4f%11d\n", 
				n + 1, segType[n], ((iSupp[ind] > -1) ? 1. : 0.), 
				mass[ind], 0);
	}
	// ABP
	for (n = 0; n < nAbpVmd; n++) {
		fprintf(fOut, "%8d %-4s 1    LYS  N      38%11.5f%14.4f%11d\n", 
				n + nActVmd + 1, segType[(n % nAbp) + nActVmd], 
				((iSupp[(n % nAbp) + nAct] > -1) ? 1. : 0.),
				mass[(n % nAbp) + nAct], 0);
	}
	// A reference particle for coloring method
	if (recConfVmd.gTglInfo > 0) {
		fprintf(fOut, "%8d %-4s 1    LYS  O      38%11.5f%14.4f%11d\n", 
				nAbpVmd + nActVmd + 1, "REF", 0., 16., 0);
	}
	// Visualization of boundaries
	if (recConfVmd.gTglBnd != 0) {
		for(n = 0; n < 8; n++) {
			fprintf(fOut, "%8d %-4s 1    LYS  H      38%11.5f%14.4f%11d\n", 
					nAbpVmd + nActVmd + n + 1 + facInfoConfVmd, 
					"BND", 0., 1., 0);
		}
	}
	// Membrane
	for (n = 0; n < nMbVmd * indMb; n++) {
		if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) {
			cnt = (rebMb[(n / indMb) % nMb] > -1) ? 2 : 1;
			if (cnt == 2) {
				for(k = nChMb - nMbAct; k < nChMb; k++) {
					if (P2A(chMb,(n / indMb) % nMb,k,nChMb) > -1) { cnt++; }
				}
			}
		}
		else { cnt = 1; }
		fprintf(fOut, "%8d %-4s 1    LYS  P      38%11.5f%14.4f%11d\n", 
				nAbpVmd + nActVmd + 1 + facInfoConfVmd + facBndConfVmd + n,
				(ISNUC(idxMb[(n / indMb) % nMb]) ? "NUC" : "MB"), 
				0., mass[(n / indMb) % nMb + nAct + nAbp], 0); 
	}
	for(n = 0; n < nBeadVmd; n++) {
	    VV3SUB(&P2(rBeadVmd,n,0), domCen);
		fprintf(fOut, "%8d %-4s 1    LYS  S      38%11.5f%14.4f%11d\n", 
				nActVmd + nAbpVmd + facInfoConfVmd + 1 + facBndConfVmd 
				+ nMbVmd * indMb + n, "BD", 0., 18., 0);
	}
	chain.c = 0;
	// Chain information from actin (chActVmd)
	nCh = (mode != 0) ? nChAcY + 2 : nChAc;
	for (n = 0; n < nActVmd; n++) {
		CONT(!(P2A(chActVmd,n,0,nCh) > -1 && P2A(chActVmd,n,1,nCh) < 0));
		curr = n;
		while(P2A(chActVmd,curr,0,nCh) > -1) {
			V3SUB(dr, &P2(rActVmd,curr,0), 
					&P2(rActVmd,P2A(chActVmd,curr,0,nCh),0));
			FOR_NDIM(k) {
				BREAK(pbc[k] != 0 && fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = curr + 1;
				P2A(chain.l,chain.c,1,2) = 
						P2A(chActVmd,curr,0,nCh) + 1;
				chain.c++; 
			}
			curr = P2A(chActVmd,curr,0,nCh);	
		}
	}
	// Chain information from ABP (chAbpVmd)
	for (n = 0; n < nAbpVmd; n++) {
		pArr = &P2A(chAbpVmd,n,0,nChAb);
		for(m = 0; m < nChAb - 1; m++) {
			CONT(m == 2);
			CONT(pArr[m] < 0);
			CONT(!(ISMTF(pArr[2])) && m > 2);
			if (m < 2) { 
				V3SUB(dr, &P2(rAbpVmd,n,0), &P2(rActVmd,pArr[m],0));
			}
			else {
				V3SUB(dr, &P2(rAbpVmd,n,0),	&P2(rAbpVmd,pArr[m],0));
			}
			FOR_NDIM(k) {
				BREAK(fabs(dr[k]) > dimDomH[k])
			}
			if (k == NDIM) { 
				P2A(chain.l,chain.c,0,2) = n + nActVmd + 1;
				P2A(chain.l,chain.c,1,2) = pArr[m] + 1 
						+ ((m < 2) ? 0 : nActVmd);
				(chain.c)++; 
			}
		}
    }
	// Chain information is also required for cylinders showing boundaries.
	if (recConfVmd.gTglBnd != 0) {
		ind = nActVmd + nAbpVmd + 1 + facInfoConfVmd;
		for(n = 0; n < 12; n++) {
			P2A(chain.l,chain.c,0,2) = P2A(bndCh,n,0,2) + ind;
			P2A(chain.l,chain.c,1,2) = P2A(bndCh,n,1,2) + ind;
			(chain.c)++; 
		}	
	}
	// Chains for membrane
	ind = nActVmd + nAbpVmd + facInfoConfVmd + 1 + facBndConfVmd;
	for (n = 0; n < nMbVmd; n++) {
		for(m = 0; m < ((dimMbNuc == 2) ? 1 : 6); m++) {
			ind3 = P2A(chMbVmd,n,m,nChMb);
			CONT(!(ind3 > -1));
			if (dimMbNuc == 3) {
				CONT(P2A(chkBond,n,m,6) != 0);
			}
			V3SUB(dr, &P2(rMbVmd,n,0), &P2(rMbVmd,ind3,0));
			FOR_NDIM(k) {
				BREAK(fabs(dr[k]) > dimDomH[k])
			}
			CONT(k != NDIM);
			for(k = 0; k < indMb; k++) {
				P2A(chain.l,chain.c,0,2) = n * indMb + k + ind;
				P2A(chain.l,chain.c,1,2) = ind3 * indMb + k + ind;
				(chain.c)++; 
			}
			if (dimMbNuc == 3) {
		        loc[0] = FindElementArray(&P2A(chMbVmd,ind3,0,
						nChMb), 6, n, 0, 1);
				P2A(chkBond,ind3,loc[0],6) = 1;
			}
		}
	}
	fprintf(fOut, "\n%8d !NBOND: bonds\n", chain.c);
	for(n = 0; n < chain.c; n++) {
		fprintf(fOut, "%8d%8d", P2A(chain.l,n,0,2), P2A(chain.l,n,1,2));
		CONT((n + 1) % 4 != 0);
		fprintf(fOut, "\n"); 
	}
	if (chain.c % 4 != 0) {
		fprintf(fOut, "\n"); 
	}	
   	fclose(fOut);

	free(rActVmd);
	free(rAbpVmd);
	free(chActVmd);
	free(chAbpVmd);
	free(chain.l);
	free(occu);
	free(mass);
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
	free(iFila);
	if (gTglMb != 0) {
		free(rMbVmd);
		free(chMbVmd);
		free(rMb);
		free(chMb);
		if (mbReb.gTgl != 0 || mbUnb.gTgl != 0) { free(rebMb); }
		if (mbFix.gTgl != 0 || mbDef.gTgl != 0) { free(fixMb); }
		free(idxMb);
		if (dimMbNuc == 3) {
			free(chkBond);
		}
    	if (mbSld.abp.gTgl != 0) {
			free(sldAbp);
		}
	}
	if (gTglBead != 0) {
		free(rBeadVmd);
	}
	for (n = 0; n < nActVmd + nAbp; n++) { free(segType[n]); }
	free(segType);
	if (recConfVmd.gTglInfo > 0) {
		free(recActG);
		free(recAbpG);
		if (gTglMb != 0)  {
			free(recMbG);
		}
	}
  }

  free(nActAll);
  free(nAbpAll);
  if (gTglMb != 0) {
	free(nMbAll);
  }
  free(chActAlt);
}

void RecordConfigMatlab(void) {
  int m, n, k, ind, cnt, *pArr, *nActAll, *nAbpAll;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
  }
  // Gather information about the positions and chains of actins and ABPs
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);

  if (rank == 0) {
	fOut = fopen(GenFileName("MlbActPos"), "a");
	if (currTimeStep == recConfMlb.prd) {
		fprintf(fOut, "%d\t0\t0\t0\n", nAct);
	}
	FOR_ACT(n) {
		fprintf(fOut, "%d\t%g\t%g\t%g\n", n, P2(rAct,n,0), P2(rAct,n,1), 
				P2(rAct,n,2));
	}
	fclose(fOut);

	fOut = fopen(GenFileName("MlbActCh"), "a");
	for(m = 0; m < 2; m++) {
		cnt = 0;
		FOR_ACT(n) {
			pArr = &P2A(chAct,n,0,nChAc);
			CONT(pArr[0] < 0 && pArr[1] < 0);
			ind = (pArr[0] > -1) ? pArr[0] : n;
			FOR_NDIM(k) {
				if (fabs(P2(rAct,ind,k) - P2(rAct,n,k)) > dimDomH[k] 
						&& pbc[k] == 1) {
					break;
				}
			}
			CONT(k < NDIM);
			if (m == 1) {
				fprintf(fOut, "%d\t%g\t%g\t%g\t", cnt, P2(rAct,n,0), 
						P2(rAct,n,1), P2(rAct,n,2));
				fprintf(fOut, "%g\t%g\t%g\n", P2(rAct,ind,0), P2(rAct,ind,1), 
						P2(rAct,ind,2));
			}
			cnt++;
		}
		if (m == 0) {
			fprintf(fOut, "%d\t0\t0\t0\t0\t0\t0\n", cnt);
		}
	}
	fclose(fOut);

	if (nAbp - nMot > 0) {
		fOut = fopen(GenFileName("MlbAcpPos"), "a");
		for(m = 0; m < 2; m++) {
			cnt = 0;
			FOR_ABP(n) {
				CONT(P2A(chAbp,n,2,nChAb) == 2);
				if (m == 1) {
					fprintf(fOut, "%d\t%g\t%g\t%g\n", cnt, P2(rAbp,n,0), 
							P2(rAbp,n,1), P2(rAbp,n,2));
				}
				cnt++;
			}
			if (m == 0) {
				fprintf(fOut, "%d\t0\t0\t0\n", cnt);
			}
		}
		fclose(fOut);
	}

	if (currTimeStep == recConfMlb.prd) {
		fOut = fopen(GenFileName("MlbMotPos"), "a");
		for(m = 0; m < 2; m++) {
			cnt = 0;
			FOR_ABP(n) {
				CONT(P2A(chAbp,n,2,nChAb) != 2);
				if (m == 1) {
					fprintf(fOut, "%d\t%g\t%g\t%g\n", cnt, P2(rAbp,n,0), 
							P2(rAbp,n,1), P2(rAbp,n,2));
				}
				cnt++;
			}
			if (m == 0) {
				fprintf(fOut, "%d\t0\t0\t0\n", cnt);
			}
		}
		fclose(fOut);
	}
  }
  free(nActAll);
  free(nAbpAll);
  if (rank == 0) {
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
  }
}

/*-------------------------- Recording network structure ---------------------*/

/*------------------- Dynamic behaviors of actin and ABPs --------------------*/

void RecordAbpDynamics(int period) {
  int n, *nAbpAll, *abpDynCntU, *abpDynCntB, *abpDynCntW;
  FILE *fOut;
  
  MALLOC(nAbpAll, int, nCpu);
  if (rank == 0) {
	MALLOC(abpMotId,int,nAbp);
	MALLOC(abpDynCntU,int,nAbp);
	MALLOC(abpDynCntB,int,nAbp);
	MALLOC(abpDynCntW,int,nAbp);
	MALLOC(chAbp,int,nAbp * nChAb);
  }
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.mId, abpMotId);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abpDyn.cntU, abpDynCntU);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abpDyn.cntB, abpDynCntB);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abpDyn.cntW, abpDynCntW);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);

  if (rank == 0) {
	fOut = fopen(GenFileName("AbpDyn"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t0\t0\t0\t0\n", nAbp);
	}
	FOR_ABP(n) {
		if (P2A(chAbp,n,2,nChAb) == 2) {
			fprintf(fOut, "%d\t%d\t%d\t%d\t%d\n", n, abpDynCntU[n],
					abpDynCntB[n], abpDynCntW[n], abpMotId[n]); 
		}
		else {
			fprintf(fOut, "%d\t%d\t%d\t-1\t-1\n", n, abpDynCntU[n],
					abpDynCntB[n]);
		}
	}  
	fclose(fOut);

	free(abpMotId);
	free(abpDynCntU);
	free(abpDynCntB);
	free(abpDynCntW);
	free(chAbp);
  }
  free(nAbpAll);
}

void RecordActinSeverEvent(void) {
  int sumCnt, *cntAll, chk[12];
  double *recvArr;
  FILE *fOut;

  MALLOC(cntAll,int,nCpu); 
  MPI_Gather(&sevLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvArr,double,sumCnt * 12);
  }
  Gather2dArrayDoubleWoIndWoSort(sevLog.c, cntAll, sevLog.l, recvArr, 12);
  if (rank == 0) {
	qsort(recvArr, sumCnt, 12 * sizeof(double), CompDbl);
	fOut = fopen(GenFileName("ActSev"), "a");
	memset(chk, 0, sizeof(int) * 12);
	chk[0] = chk[1] = chk[2] = chk[3] = 1;
	Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 12, chk);
	fclose(fOut);
	free(recvArr);
  }
  free(cntAll);
  sevLog.c = 0;
}

void RecordAbpBindEvent(int abpInd, int side, int mode) {
  int n, locAbpInd, actInd, locActInd, sumCnt, nextActInd, actSide;
  int chk[26], *cntAll;
  double *recvArr, rPos[NDIM], dr[NDIM];
  FILE *fOut;  

  if (mode < 2) {
	locAbpInd = iAbp[abpInd];
	V5SET(&P2A(bindLog.l,bindLog.c,0,26), mode, (double)currTimeStep, 
			(double)abpInd, K_ABP(locAbpInd), (double)side);
	V3COPY(&P2A(bindLog.l,bindLog.c,5,26), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			CalcUnitVec(dr, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(bindLog.l,bindLog.c,8 + 9 * n,26), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(bindLog.l,bindLog.c,10 + 9 * n,26), rPos);
			V3COPY(&P2A(bindLog.l,bindLog.c,13 + 9 * n,26), dr);
			P2A(bindLog.l,bindLog.c,16 + 9 * n,26) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(bindLog.l,bindLog.c,8 + 9 * n,26), 
					9, 0.);
		}
	}
	(bindLog.c)++;
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&bindLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt * 26);
	}
	Gather2dArrayDoubleWoIndWoSort(bindLog.c, cntAll, bindLog.l, recvArr, 26);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 26 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpBind"), "a");
		memset(chk, 0, sizeof(int) * 26);
		chk[0] = chk[2] = chk[3] = chk[8] = chk[16] = chk[17] = chk[25] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 26, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	bindLog.c = 0;
  }
}

void RecordAbpUnbindEvent(int abpInd, int side, int mode) {
  int n, locAbpInd, actInd, locActInd, sumCnt, nextActInd, actSide;
  int chk[29], *cntAll;
  double *recvArr, rPos[NDIM], dr[NDIM];
  FILE *fOut;  

  double len, len2, dr2[NDIM], drAxis[NDIM], fac, rPos2[NDIM], fInst[2];
  double fi[NDIM], f, f2;
  int kind, k;

  if (mode == 0) {
	locAbpInd = iAbp[abpInd];
	V4SET(&P2A(unbLog.l,unbLog.c,0,29), 0, (double)currTimeStep, (double)abpInd,
			(double)side);
	V3COPY(&P2A(unbLog.l,unbLog.c,4,29), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			CalcUnitVec(dr, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(unbLog.l,unbLog.c,9 + 11 * n,29), rPos);
			V3COPY(&P2A(unbLog.l,unbLog.c,12 + 11 * n,29), dr);
			V2COPY(&P2A(unbLog.l,unbLog.c,15 + 11 * n,29), 
					&P2A(recInstSprFabp,locAbpInd,side * 2,4));
			P2A(unbLog.l,unbLog.c,17 + 11 * n,29) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), 
					11, 0.);
		}
	}
	(unbLog.c)++;
  }
  else if (mode == 1) {
	locAbpInd = iAbp[abpInd];
	V4SET(&P2A(unbLog.l,unbLog.c,0,29), 1, (double)currTimeStep, (double)abpInd,
			(double)side);
	V3COPY(&P2A(unbLog.l,unbLog.c,4,29), &P2(abp.r,locAbpInd,0));
	for(n = 0; n < 2; n++) {
		actInd = P2A(abp.ch,locAbpInd,(n == 0 ? side : 1 - side),nChAb);
		if (actInd > -1) {
			locActInd = iAct[actInd];
			nextActInd = P2A(act.ch,locActInd,0,nChAc);
			kind = K_ABP(locAbpInd);

			// Find the location of ABP on the actin segment.
			actSide = FindAbpActinChain(locActInd, abpInd, 0);
			CalcPosOnActSegSide(&P2(act.r,locActInd,0), 
					&P2(act.r,iAct[nextActInd],0), rPos, actSide);
			len = CalcVecDist(dr, rPos, &P2(abp.r,locAbpInd,0), 0);
			if (ISMTF(kind)) {
				CalcVec(drAxis, &P2(act.r,locActInd,0), 
						&P2(act.r,iAct[nextActInd],0));
				fac = V3DOT(dr, drAxis) / V3LEN_SQ(drAxis);
				VS3SUB(rPos2, rPos, drAxis, fac);
				ApplyBoundCondVector(rPos2, -1, 0);
				len = CalcVecDist(dr, rPos2, &P2(abp.r,locAbpInd,0), 0);
				len2 = CalcVecDist(dr2, rPos, rPos2, 0);
			}	
			f = SPRING(abpF.spr[kind].stf, len, abpF.spr[kind].eq);
			if (ISMTF(kind)) {
				f2 = SPRING(abpF.spr[kind].stf2, len2, 0.);
			}

			f /= len;
			if (ISMTF(kind) && len2 > 0) { f2 /= len2; }
			FOR_NDIM(k) {
				fi[k] = f * dr[k];
				if (ISMTF(kind)) { fi[k] += f2 * dr2[k]; }
			}
			fInst[0] = REVSIGN(f * len);
			if (ISMTF(kind)) {
				fInst[1] = f2 * len2 * ((fac < 0) ? -1. : 1.);
			}
			else {
				if (f < 0) {
					CalcUnitVec(drAxis, &P2(act.r,locActInd,0),
							&P2(act.r,iAct[nextActInd],0));
					fInst[1] = V3DOT(fi, drAxis);
				}
				else {
					fInst[1] = 0;
				}
			}
			CalcUnitVec(drAxis, &P2(act.r,iAct[nextActInd],0), 
					&P2(act.r,locActInd,0));
			V2SET(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), (double)actInd, 
					(double)((int)((actSide - 2) / nChAcY)) / (double)nChAcX);
			V3COPY(&P2A(unbLog.l,unbLog.c,9 + 11 * n,29), rPos);
			V3COPY(&P2A(unbLog.l,unbLog.c,12 + 11 * n,29), drAxis);
			V2COPY(&P2A(unbLog.l,unbLog.c,15 + 11 * n,29), fInst);
			P2A(unbLog.l,unbLog.c,17 + 11 * n,29) = (double)act.iF[locActInd];
		}
		else {
			SetAllValue1dArrayDouble(&P2A(unbLog.l,unbLog.c,7 + 11 * n,29), 
					11, 0.);
		}
	}
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&unbLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt * 29);
	}
	Gather2dArrayDoubleWoIndWoSort(unbLog.c, cntAll, unbLog.l, recvArr, 29);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 29 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpUnb"), "a");
		memset(chk, 0, sizeof(int) * 29);
		chk[0] = chk[2] = chk[3] = chk[7] = chk[17] = chk[18] = chk[28] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 29, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	unbLog.c = 0;
  }
}

// Record the turnover of ABPs: to find the displacement from the location 
// of unbinding to the location of following binding and whether unbound ABPs
// bind to the same filament or not.
// [Columns of abpTurn]
// 0: the moment of unbinding 
// 1, 2, 3: previous position, 4: force applied at the moment of unbinding
// 5: previous shear strain, 6: previous filament number 
// mode = 0: called at the moment of unbinding, 1: called at that of binding
void RecordAbpTurnover(int abpInd, int actInd, int side, int mode) {
  int locAbpInd, sameFila, *cntAll, sumCnt, chk[14];
  double dr[NDIM], len, *recvArr;
  FILE *fOut;

  locAbpInd = iAbp[abpInd];
  // At the moment of unbinding
  if (mode == 0) {
	P2A(abpTurn,locAbpInd,0,7) = (double)currTimeStep;
	V3COPY(&P2A(abpTurn,locAbpInd,1,7), &P2(abp.r,locAbpInd,0));
	V3SET(&P2A(abpTurn,locAbpInd,4,7), stra.acc, 
			P2A(recInstSprFabp,locAbpInd,side * 2,4),
			(double)act.iF[iAct[actInd]]);
  }
  // At the moment of binding
  else if (mode == 1) {
	V3SUB(dr, &P2(abp.r,locAbpInd,0), &P2A(abpTurn,locAbpInd,1,7));
	if (rheoWay > 0 && bulkRheoType == 0) {
		dr[dirStr] -= (stra.acc - P2A(abpTurn,locAbpInd,4,7)) 
				* P2A(abpTurn,locAbpInd,dirNoPBC + 1,7);
	}
    ApplyBoundCondVecDiff(dr);
	len = V3LEN(dr);
	sameFila = (act.iF[iAct[actInd]] 
			== (int)P2A(abpTurn,locAbpInd,6,7)) ? 1 : 0;
	V6SET(&P2A(toLog.l,toLog.c,0,14), P2A(abpTurn,locAbpInd,0,7), 
			(double)currTimeStep, (double)abpInd, 
			F_S2PN(P2A(abpTurn,locAbpInd,5,7)), L_S2NM(len), (double)sameFila);
	V3COPY(&P2A(toLog.l,toLog.c,6,14), &P2(abp.r,locAbpInd,0));
	V4COPY(&P2A(toLog.l,toLog.c,9,14), &P2A(abpTurn,locAbpInd,1,7));
	P2A(toLog.l,toLog.c,13,14) = stra.acc;
	(toLog.c)++;
	SetAllValue1dArrayDouble(&P2A(abpTurn,locAbpInd,0,7), 7, -1.);
  }
  else {
	MALLOC(cntAll,int,nCpu); 
	MPI_Gather(&toLog.c, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		sumCnt = SumArrInt(cntAll, nCpu);
		MALLOC(recvArr,double,sumCnt*14);
	}
	Gather2dArrayDoubleWoIndWoSort(toLog.c, cntAll, toLog.l, recvArr, 14);
	if (rank == 0) {
		qsort(recvArr, sumCnt, 14 * sizeof(double), CompDbl);
		fOut = fopen(GenFileName("AbpTurnover"), "a");
		memset(chk, 0, sizeof(int) * 14);
		chk[0] = chk[1] = chk[2] = chk[5] = 1;
		Fprintf2dArrayIntDouble(fOut, recvArr, sumCnt, 14, chk);
		fclose(fOut);
		free(recvArr);
	}
	free(cntAll);
	toLog.c = 0;
  }
}

/*------------------- Dynamic behaviors of actin and ABPs --------------------*/

/*------------------------- Motor thick filaments ----------------------------*/

void RecordMotorPosition(void) {
  int n, k, len, curr, next, nMot, *nAbpAll, *chAbp, *iMot, *ageAbp, *ageMot;
  double *rMot, *rAbp;
  FILE *fOut;

  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(iMot,int,nAbp*motSA.nMotPerTF);
	MALLOC(rMot,double,nAbp*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAbp,int,nAbp*nChAb);
	if (abpAge.gTgl != 0) {
		MALLOC(ageMot,int,nAbp);
		MALLOC(ageAbp,int,nAbp);
	}
	memset(iMot, -1, sizeof(int) * nAbp * motSA.nMotPerTF);
  }
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM, abp.id, abp.r, rAbp);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);
  if (abpAge.gTgl != 0) {
	Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.age, ageAbp);
  }

  if (rank == 0) {
	fOut = fopen(GenFileName("MotPos"), "a");
	nMot = 0;
	FOR_ABP(n) {
		CONT(P2A(chAbp,n,2,nChAb) != 2);
		CONT(!(P2A(chAbp,n,3,nChAb) > -1 && P2A(chAbp,n,4,nChAb) < 0));
		len = 0;
		curr = n;
		while(curr > -1) {
			len++;
			curr = P2A(chAbp,curr,3,nChAb);
		}
		curr = n;
		for(k = 0; k < len / 2 - 1; k++) {
			curr = P2A(chAbp,curr,3,nChAb);
		}
		next = P2A(chAbp,curr,3,nChAb);
		V3AVG(&P2(rMot,nMot,0), &P2(rAbp,curr,0), &P2(rAbp,next,0));
		FOR_NDIM(k) {
			if (fabs(P2(rAbp,next,k) - P2(rAbp,curr,k)) > dimDomH[k]
					&& pbc[k] == 1) {
				P2(rMot,nMot,k) += dimDomH[k];
				if (P2(rMot,nMot,k) >= dimDom[k]) {
					P2(rMot,nMot,k) -= dimDom[k];
				}
			}
		}
		if (abpAge.gTgl != 0) {
			curr = n;
			ageMot[nMot] = ageAbp[n];
			for(k = 0; k < len; k++) {
				if (ageAbp[curr] < ageMot[nMot]) {
					ageMot[nMot] = ageAbp[curr];
				}
				P2A(iMot,nMot,k,motSA.nMotPerTF) = curr;
				curr = P2A(chAbp,curr,3,nChAb);
			}
		}
		nMot++;
	}
	fprintf(fOut, "%d\t", nMot);
	fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dFillerInt(fOut, 0, motSA.nMotPerTF + NDIM, 0);
	for(n = 0; n < nMot; n++) {
		fprintf(fOut, "%d\t", n);
		Fprintf1dArrayDoubleWoRet(fOut, &P2(rMot,n,0), NDIM);
		Fprintf1dArrayIntWoRet(fOut, &P2A(iMot,n,0,motSA.nMotPerTF), 
				motSA.nMotPerTF);
		if (abpAge.gTgl != 0) {
			fprintf(fOut, "%d", ageMot[n]);
		}
		fprintf(fOut, "\n");
	}
	free(iMot);
	free(rMot);
	free(rAbp);
	free(chAbp);
	if (abpAge.gTgl != 0) {
		free(ageMot);
		free(ageAbp);
	}
	fclose(fOut);
  }
  free(nAbpAll);
}

// Record the distribution of actin filaments
void RecordMotorFilamentSize(void) {
  int n, k, curr, maxLen, CS;
  int len[2], cntFila[3], **filaLen, *nAbpAll, *pArr;
  double avgFilaLen[3], stdFilaLen[3];
  FILE *fOut;

  maxLen = NEG_LARGE_VALUE;
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) { 
	MALLOC(chAbp,int,nAbp*nChAb); 
	MALLOC(abpMotId,int,nAbp);
  }
  // Gather information about chains.
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp); 
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, abp.mId, abpMotId); 

  if (rank == 0) {
	MALLOC2(filaLen,int,3);
	for(n = 0; n < 3; n++) {
		MALLOC(filaLen[n],int,nAbp);
		memset(filaLen[n], 0, sizeof(int) * nAbp);
	}
	V3SET_ALL(cntFila, 0); 
	V3SET_ALL(avgFilaLen, 0.);
	V3SET_ALL(stdFilaLen, 0.);
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		CONT(abpMotId[n] < 0);
		CONT(pArr[3] < 0 && pArr[4] < 0);
		CONT(!(pArr[3] > -1 && pArr[4] < 0));
		CS = 1;
		V2SET_ALL(len, 0);
		curr = n;
		while (curr > -1) {
			len[0]++;
			len[1]++;
			if (abpMotId[curr] == 0 && CS == 1) {
				filaLen[1][len[1]]++;
				len[1] = 0;
				CS = 0;
			}
			curr = P2A(chAbp,curr,3,nChAb);
		} 
		filaLen[0][len[0]]++;
		filaLen[2][len[1]]++;
		if (maxLen < len[0]) { maxLen = len[0]; }
	}
	if (maxLen < 0) { maxLen = 0; }
	fOut = fopen(GenFileName("MotFilaLen"), "a");
	for(n = 1; n < maxLen + 1; n++) {
		for(k = 0; k < 3; k++) {
		    cntFila[k] += filaLen[k][n];
			avgFilaLen[k] += (double)(n * filaLen[k][n]);
			stdFilaLen[k] += (double)(SQR(n) * filaLen[k][n]);
		}
	}
	// Record the number of lines recorded.
	fprintf(fOut, "%d\t0\t%d\t%d\t%d\n", maxLen + 2, cntFila[0], cntFila[1], 
			cntFila[2]);
	// Calculate average length and standard deviation.
	for(n = 1; n < maxLen + 1; n++) {
	    fprintf(fOut, "%d\t%d\t%d\t%d\t%d\n", n - 1, n, filaLen[0][n], 
				filaLen[1][n], filaLen[2][n]);
	}
	for(n = 0; n < 3; n++) {
		if (cntFila[n] == 0) { 
			avgFilaLen[n] = 0.; 
			stdFilaLen[n] = 0.;
		}
		else { 
			avgFilaLen[n] /= (double)cntFila[n]; 
			stdFilaLen[n] = sqrt(stdFilaLen[n] / (double)cntFila[n] 
					- SQR(avgFilaLen[n]));
		}
	}
	fprintf(fOut, "%d\t0\t%g\t%g\t%g\n", maxLen, avgFilaLen[0], avgFilaLen[1], 
			avgFilaLen[2]);
	fprintf(fOut, "%d\t0\t%g\t%g\t%g\n", maxLen + 1, stdFilaLen[0], 
			stdFilaLen[1], stdFilaLen[2]);
	fclose(fOut);
	for(n = 0; n < 3; n++) {
		free(filaLen[n]);
	}
	free(filaLen);
	free(chAbp);
	free(abpMotId);
  }
  free(nAbpAll);
}

/*------------------------- Motor thick filaments ----------------------------*/

/*--------------------- Structural properties of network ---------------------*/

// Record the distribution of actin filaments
void RecordFilamentLength(int period) {
  int n, len, curr, cntFila, maxLen, *filaLen, *nActAll, *pArr;
  double avgFilaLen, stdFilaLen, len2;
  FILE *fOut;

  maxLen = -1;
  MALLOC(nActAll,int,nCpu);
  if (rank == 0) { MALLOC(chAct,int,nAct*nChAc); }
  // Gather information about chains.
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct); 

  if (rank == 0) {
	MALLOC(filaLen,int,nAct+1);
	memset(filaLen, 0, sizeof(int) * (nAct + 1) );
	cntFila = 0; 
	avgFilaLen = 0.; 
	stdFilaLen = 0.;

	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if (pArr[0] < 0 && pArr[1] < 0) { 
			filaLen[0]++; 
		} 
		else if (pArr[0] > -1 && pArr[1] < 0) { 
			len = 1; 
			curr = pArr[0];
			while (P2A(chAct,curr,0,nChAc) > -1) {
				curr = P2A(chAct,curr,0,nChAc);
				len++;
			}
			filaLen[len]++;
			if (maxLen < len) { maxLen = len; }
		}
	}
	fOut = fopen(GenFileName("ActFilaLen"), "a");
	// Record the number of lines recorded.
    Fprintf1dArrayIntWFil(fOut, &maxLen, 1, 0, 3);
	// Calculate average length and standard deviation.
	for(n = 1; n < maxLen + 1; n++) {
		len2 = L_S2UM(n);
	    fprintf(fOut, "%d\t%g\t%d\n", n - 1, len2, filaLen[n]);
	    cntFila += filaLen[n];
		avgFilaLen += (double)(len2 * filaLen[n]);
		stdFilaLen += (double)(SQR(len2) * filaLen[n]);
	}
	if (cntFila == 0) { 
		avgFilaLen = 0.; 
		stdFilaLen = 0.;
	}
	else { 
		avgFilaLen /= (double)cntFila; 
		stdFilaLen = sqrt(stdFilaLen / (double)cntFila - SQR(avgFilaLen));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxLen - 1, avgFilaLen, stdFilaLen);
	fclose(fOut);
	free(filaLen);
	free(chAct);
  }
  free(nActAll);
}

// Record the distribution of angles formed by two crosslinked filaments
void RecordCrosslinkAngle(void) {
  double dr1[NDIM], dr2[NDIM], mag1, mag2; 
  double c, angIntv, avgAng, stdAng;
  int n, *crsAng, cntCrossAng, min, max, cntAng;
  int *nActAll, *nAbpAll, *pArr;
  FILE *fOut;

  angIntv = 5;
  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
  }
  // Gather information about the positions and chains of actins and ABPs
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);

  if (rank == 0) {
	cntCrossAng = (int)(180. / (double)angIntv)+ 1;
	MALLOC(crsAng, int, cntCrossAng);
	memset(crsAng, 0, sizeof(int) * cntCrossAng);

	fOut = fopen(GenFileName("CrossAng"), "a");
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		CONT(ISMTF(pArr[2]));
		CONT(!(pArr[0] > -1 && pArr[1] > -1));
		mag1 = CalcVecDist(dr1, &P2(rAct,P2A(chAct,pArr[0],1,nChAc),0),
			&P2(rAct,P2A(chAct,pArr[0],0,nChAc),0), 0);
		mag2 = CalcVecDist(dr2, &P2(rAct,P2A(chAct,pArr[1],1,nChAc),0),
			&P2(rAct,P2A(chAct,pArr[1],0,nChAc),0), 0);
		c = V3DOT(dr1, dr2) / (mag1 * mag2);
		crsAng[(int)(RAD2DEG(Acos(c)) / angIntv)]++;
	}
	for(n = 0; n < cntCrossAng; n++) {
		CONT(crsAng[n] == 0);
		min = n;
		break;
	}
	for(n = cntCrossAng - 1; n >= 0; n--) {
		CONT(crsAng[n] == 0);
		max = n;
		break;
	}
	// Calculate average and standard deviation.
	cntAng = 0;
	avgAng = 0.;
	stdAng = 0.;
	for(n = min; n < max + 1; n++) {
		cntAng += crsAng[n];	
	}
	fprintf(fOut, "%d\t%d\t0\n", max - min + 2, cntAng);	
	for(n = min; n < max + 1; n++) {
		fprintf(fOut, "%d\t%g\t%d\n", n - min, (n + 0.5) * angIntv, 
				crsAng[n]);
		avgAng += ((n + 0.5) * angIntv) * (double)crsAng[n];
		stdAng += SQR(((n + 0.5) * angIntv)) * (double)crsAng[n];
	}
	if (cntAng == 0) { 
		avgAng = 0.; 
		stdAng = 0.; 
	}
	else {
		avgAng /= (double)cntAng;
		stdAng = sqrt(stdAng / (double)cntAng - SQR(avgAng));
	}
	fprintf(fOut, "%d\t%g\t%g\n", max - min + 1, avgAng, stdAng);
	fclose(fOut);
  }
  free(nActAll);
  free(nAbpAll);
  if (rank == 0) {
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
	free(crsAng);
  }
}

// Record the distribution of crosslinking distance which corresponds to 
// distance between active ABPs on a single filament.
void RecordCrosslinkDistance(void) {
  int n, k, l, curr, len, maxLen, actInd, abpInd, cntDist, CS, cnt, prevIMot;
  int *nActAll, *nAbpAll, *crossDist, *pArr, *pArr2, *iMot, *iFila;
  double avgCrossDist, stdCrossDist, len2;
  FILE *fOut;
  ListInt motIFila;

  maxLen = NEG_LARGE_VALUE;
  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	if (motSA.gTgl != 0) {
		MALLOC(motIFila.l,int,nAbp*2);
		MALLOC(iFila, int, nAct);
		MALLOC(iMot, int, nAbp);
		memset(iMot, -1, sizeof(int) * nAbp);
	}
  }

  // Gather information about chains of actin
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, nChAc, act.id, act.ch, chAct);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp);
  if (motSA.gTgl != 0) {
	Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  }

  if (rank == 0) {
	MALLOC(crossDist, int, nAct * nChAcX);
	fOut = fopen(GenFileName("CrossDist"),"a");
	memset(crossDist, 0, sizeof(int) * nAct * nChAcX);
	avgCrossDist = 0.;
	cntDist = 0;

	cnt = 0;
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		CONT(!(ISMTF(pArr[2])));
		CONT(!(pArr[3] > -1 && pArr[4] < 0));
		motIFila.c = 0;
		curr = n;
		while(curr > -1) { 
			for(k = 0; k < 2; k++) {
				actInd = P2A(chAbp,curr,k,nChAb);
				CONT(!(actInd > -1));
				InsertElement1dArrayWChk(motIFila.l, &motIFila.c, 
						iFila[actInd]);
			}
			curr = P2A(chAbp,curr,3,nChAb);
		}
		if (motIFila.c > 1) {
			curr = n;
			while(curr > -1) { 
				iMot[curr] = cnt;
				curr = P2A(chAbp,curr,3,nChAb);
			}	
			cnt++;
		}
	}
	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		CONT(!(pArr[0] > -1 && pArr[1] < 0));
		len = 0;
		CS = 0;
		prevIMot = -1;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = P2A(chAct,curr,2 + k * nChAcY + l,nChAc);
					CONT(abpInd < 0);
					pArr2 = &P2A(chAbp,abpInd,0,nChAb);
					if (ISMTF(pArr2[2])) {
						CONT(iMot[abpInd] < 0);
						if (iMot[abpInd] == prevIMot) {
							len = 0;
							continue;
						}
						prevIMot = iMot[abpInd];
					}
					else {
						CONT(!(pArr2[0] > -1 && pArr2[1] > -1));
						prevIMot = -1;
					}
					// if it is not the end of filament
					if (CS != 0) {
						crossDist[len]++;
						if (maxLen < len) { maxLen = len; }
					}
					len = 0;
					CS = 1;
				}
				len++;
			}
			curr = P2A(chAct,curr,0,nChAc);
		}
	}
	if (maxLen < 0) { maxLen = -1; }
	for (n = 0; n < maxLen + 1; n++) {
		cntDist += crossDist[n];
	}
	fprintf(fOut, "%d\t%d\t0\n", maxLen + 2, cntDist);
	// Calculate average and standard deviation.
	for (n = 0; n < maxLen + 1; n++) {
		len2 = L_S2NM((double)n / (double)nChAcX);
		fprintf(fOut, "%d\t%g\t%d\n", n, len2, crossDist[n]);
		avgCrossDist += (double)(crossDist[n] * len2);
		stdCrossDist += (double)(crossDist[n] * SQR(len2));
	}
	if (cntDist == 0) {
		avgCrossDist = 0.;
		stdCrossDist = 0.;
	}
	else { 
		avgCrossDist /= (double)cntDist; 
		stdCrossDist = sqrt(stdCrossDist / (double)cntDist - SQR(avgCrossDist));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxLen + 1, avgCrossDist, stdCrossDist);
	fclose(fOut);
	free(crossDist);
	free(chAct);
	free(chAbp);
	if (motSA.gTgl != 0) { 
		free(motIFila.l);
		free(iFila);
		free(iMot); 
	}
  }
  free(nActAll);
  free(nAbpAll);
}

void RecordConnectivitySubroutine1(int *chkIFila, int *iFila, int *nFila) {
  int n, ind;

  if (rank == 0) {
	// Mark the number of each filament
	memset(chkIFila, -1, nAct * sizeof(int));
	*nFila = 0;
	FOR_ACT(n) {
		CONT(!(iFila[n] > -1));
		InsertElement1dArrayWChk(chkIFila, nFila, iFila[n]);
	}
	FOR_ACT(n) {
		if (iFila[n] > -1) {
			ind = FindElementArray(chkIFila, *nFila, iFila[n], 0, 1);
			iFila[n] = ind;
		}
		else { iFila[n] = -1; }
	}
  }
  MPI_Bcast(nFila, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void RecordConnectivitySubroutine2(int **perc, int *cntPerc, 
		int *iFila, int nFila, int thres) {
  int m, n, k, mm, *pArr, div = 1, actInd, curr;
  ListInt motIFila;

  if (motSA.gTgl != 0) {
	MALLOC(motIFila.l,int,nAbp);
  }
  // Fill in the perc[0] array using "chAbp"
  if (rank == 0) {
	memset(perc[0], 0, sizeof(int) * nFila * nFila);
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		if (ISMTF(pArr[2])) { 
			CONT(!(pArr[3] > -1 && pArr[4] < 0));
			motIFila.c = 0;
			curr = n;
			while(curr > -1) { 
				for(k = 0; k < 2; k++) {
					actInd = P2A(chAbp,curr,k,nChAb);
					CONT(!(actInd > -1));
					InsertElement1dArrayWChk(motIFila.l, &motIFila.c, 
							iFila[actInd]);
				}
				curr = P2A(chAbp,curr,3,nChAb);
			}
			for(k = 0; k < motIFila.c; k++) {
				for(m = k + 1; m < motIFila.c; m++) {
					P2A(perc[0],iFila[motIFila.l[k]],iFila[motIFila.l[m]],nFila)
							+= 1;
				}
			}			
		}
		else {
			CONT(!(pArr[0] > -1 && pArr[1] > -1));
			P2A(perc[0],iFila[pArr[0]],iFila[pArr[1]],nFila) += 1;
			P2A(perc[0],iFila[pArr[1]],iFila[pArr[0]],nFila) += 1;
		}
	}
	for(n = 0; n < nFila; n++) {
		for(k = 0; k < nFila; k++) {
			P2A(perc[0],n,k,nFila) = (P2A(perc[0],n,k,nFila) >= thres) ? 1 : 0;
		}
	}
  }
  MPI_Bcast(perc[0], nFila * nFila, MPI_INT, 0, MPI_COMM_WORLD);
  cntPerc[0] = SumArrInt(perc[0], nFila * nFila);
  if (cntPerc[0] == 0) {
	return; 
  }
  for(m = 0; m < nFila / div - 2; m++) {
	// gradual allocation of arrays (not at once due to memory limit)
	MALLOC(perc[m + 1], int, nFila * nFila);
	MultiplySqMatrixIntParallel(perc[m], perc[0], perc[m + 1], nFila);
	// Calculate the sum of all elements in perc using the property of
	// a symmetric array
	for(n = 0; n < nFila; n++) {
		P2A(perc[m + 1],n,n,nFila) = 0;
		for(k = 0; k < nFila; k++) {
			for(mm = 0; mm < m + 1; mm++) {
				CONT(!(P2A(perc[mm],n,k,nFila) > 0));
				P2A(perc[m + 1],n,k,nFila) = 0;
				P2A(perc[m + 1],k,n,nFila) = 0;
				break;
			}
			cntPerc[m + 1] += P2A(perc[m + 1],n,k,nFila);
		}
	}
	cntPerc[m + 1] *= 2;
	BREAK(cntPerc[m + 1] == 0);
  }
  if (motSA.gTgl != 0) { 
	free(motIFila.l);
  }
}

// Record the connectivity level of networks using the adjacent matrix concept.
// thres: the number of crosslinks required for a percolation between 
//        two filaments. 0: one crosslink, 1: two crosslinks, etc..
void RecordConnectivity(int thres) {
  int n, nFila, maxPerc, div = 1;
  int **perc, *cntPerc, *nActAll, *nAbpAll, *iFila, *chkIFila;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(iFila,int,nAct);
	MALLOC(chkIFila,int,nAct);
  }
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nAbpMe, nAbpAll, nChAb, abp.id, abp.ch, chAbp); 
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);

  RecordConnectivitySubroutine1(chkIFila, iFila, &nFila);
  // perc: connectivity information
  // Allocation and initialization of arrays
  MALLOC2(perc, int, nFila / div);
  MALLOC(perc[0], int, nFila * nFila);
  MALLOC(cntPerc, int, nFila / div);
  memset(cntPerc, 0, sizeof(int) * nFila / div);
  RecordConnectivitySubroutine2(perc, cntPerc, iFila, nFila, thres);
  maxPerc = 0;
  while(cntPerc[maxPerc] > 0) {
    maxPerc++;
  }
  if (rank == 0) {
	fOut = fopen(GenFileName("Connect"), "a");
	fprintf(fOut, "%d\t%d\t0\n", maxPerc, nFila);
	for (n = 0; n < maxPerc; n++) {
		fprintf(fOut, "%d\t%d\t%d\n", n, n + 2, cntPerc[n] / 2);
	}
	fclose(fOut);
	free(chAbp);
	free(iFila);
	free(chkIFila);
  }
  for (n = 0; n < maxPerc + 1; n++) {
	free(perc[n]);
  }
  free(perc);
  free(cntPerc);
  free(nActAll);
  free(nAbpAll);
}

// Subroutine for RecordPoreSize()
double RecordPoreSizeSubroutine1(double *rPnt, double *rPnt1, double *rPnt2,
		double *dr) {
  int k;
  double rCen[NDIM], rDiff[NDIM], rOff[3][NDIM], len, ratio, r1[NDIM], r2[NDIM];

  V3COPY(r1, rPnt1);
  V3COPY(r2, rPnt2);
  if (rheoWay > 0) {
	ApplyBoundCondVector(r1, -1, 0);
	ApplyBoundCondVector(r2, -1, 0);
  }

  V3AVG(rCen, r1, r2);
  V3SUB(rDiff, r1, r2);
  FOR_NDIM(k) {
	CONT(!(fabs(rDiff[k]) > dimDomH[k]));
	if (rCen[k] > dimDomH[k] + rGrid[k][0]) { rCen[k] -= dimDomH[k]; }
	else { rCen[k] += dimDomH[k]; }
  }
  OffsetSegEndPnt(rOff[0], rCen, rPnt);
  OffsetSegEndPnt(rOff[1], r1, rOff[0]);
  OffsetSegEndPnt(rOff[2], r2, rOff[0]);
  len = CalcSegPntDist(&rOff[1], rPnt, dr, &ratio);
  return len;
}
void RecordPoreSizeSubroutine2(ListDbl *poreAll, int *k) {
  int n;
  for(n = *k; n < poreAll->c - 1; n++) {
	V7COPY(&P2A(poreAll->l,n,0,7), &P2A(poreAll->l,n + 1,0,7));
  }
  (poreAll->c)--;
  (*k)--;
}

void RecordPoreSizeSubroutine3(double *currPore, double *oldPore, 
		int *nInitPnt, double *thres, int *CS) {
  int k, cnt;
  double dr[NDIM];

  // compare the final center of the pore with other previous pores
  CalcVec(dr, currPore, oldPore);
  if (fabs(dr[0]) < thres[1] && fabs(dr[1]) < thres[1]
		&& fabs(dr[2]) < thres[1]) {
	*CS = 0;
	return;
  }
  // if the final pore radius is similar to that of 
  // previous pores, compare their initial positions.
  if (fabs(currPore[3] - oldPore[3]) < thres[2]) {
	V3SUB(dr, &currPore[4], &oldPore[4]);

	cnt = 0;
	FOR_NDIM(k) {
		CONT(!(fabs(dr[k]) < thres[0] || fabs(dr[k]) > nInitPnt[k] - thres[0]));
		cnt++;
	}
	if (cnt == NDIM) { *CS = 0; }
  }
}

// Record the distribution of pore sizes in networks. 
// In this measurement method, many spheres are uniformly placed within 
// the network. Then, the size of the spheres is increased until they cannot
// become larger due to steric effects caused by actin filaments. Then, the 
// final diameter of the sphere is the size of the pore.
void RecordPoreSize(void) {
  int m, n, k, CS, histIntv, actInd, cnt, totNInitPnt;
  int poreSize, maxPoreSize, minPoreSize, cntOverl, maxCntOverl, nPore;
  int ind[NDIM], nInitPnt[NDIM], *nActAll, *poreCall, *actCh;
  double initIntv, totPoreVol, len, radPore, avgPoreSize, stdPoreSize;
  double initRadPore, minDimDomH;
  double dr[NDIM], rPore[NDIM], drPore[NDIM], thres[3], currPore[7];
  ListDbl pore, poreAll;
  ListInt actCyl, histPore;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(rAct,double,nAct*NDIM);
  MALLOC(chAct,int,nAct*2);
  MALLOC(actCh,int,nActMe*2);
  for(n = 0; n < nActMe; n++) {
	V2COPY(&P2A(actCh,n,0,2), &P2A(act.ch,n,0,nChAc));
  }

  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nActMe, nActAll, NDIM, act.id, act.r, rAct); 
  Gather2dArrayInt(nActMe, nActAll, 2, act.id, actCh, chAct); 

  MPI_Bcast(rAct, nAct * NDIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chAct, nAct * 2, MPI_INT, 0, MPI_COMM_WORLD);

  actCyl.c = 0;
  MALLOC(actCyl.l,int,nAct*2);
  FOR_ACT(n) {
	actInd = P2A(chAct,n,0,2);
	CONT(!(actInd > -1));
	V2SET(&P2A(actCyl.l,actCyl.c,0,2), n, actInd);
	(actCyl.c)++;
  }

  // initIntv is a kind of resolution in a dimensionless unit.
  initRadPore = 0.5;
  initIntv = 1.2;
  histIntv = 20; // in nm
  minPoreSize = POS_LARGE_VALUE;
  maxPoreSize = NEG_LARGE_VALUE;
  maxCntOverl = 10;

  // thres[0]: threshold for initial pore positions
  // thres[1]: threshold for final pore positions
  // thres[2]: threshold for pore radius
  V3SET(thres, 1.2, 1.2, 0.5);

  FOR_NDIM(k) {
	nInitPnt[k] = (int)(dimDom[k] / initIntv);
  }
  totNInitPnt = V3PROD(nInitPnt);
  MALLOC(pore.l,double,totNInitPnt*7);
  pore.c = 0;
	
  minDimDomH = POS_LARGE_VALUE;
  FOR_NDIM(k) {
 	if (dimDomH[k] < minDimDomH) { minDimDomH = dimDomH[k]; }
  }
  histPore.siz = (int)(minDimDomH * L_SCALE_IN_NM / (double)histIntv) + 1;
  MALLOC(histPore.l,int,histPore.siz);
  memset(histPore.l, 0, sizeof(int) * histPore.siz); 
	
  for(m = 0; m < totNInitPnt; m++) {
 	V3IND_ASSIGN_INT(m, nInitPnt, 1, ind);
	// Determine the initial locations of spheres
	VS3COPY(rPore, ind, initIntv);	
	FOR_NDIM(k) { rPore[k] += rGrid[k][0]; }

	cnt = 0;
	FOR_NDIM(k) {
		CONT(!(rPore[k] >= P2A(bnd.r,0,k,NDIM) 
				&& rPore[k] < P2A(bnd.r,1,k,NDIM)));
		cnt++;
	}
	CONT(cnt != NDIM);

	if (rheoWay > 0 && bulkRheoType == 0) {
		rPore[dirStr] += rPore[dirNoPBC] * (stra.accDsp / dimDom[dirNoPBC]);
	}
	CS = 1;
	// Check whether the initial sphere overlaps with actins
	CONT(CS != 1);
	for (n = 0; n < actCyl.c; n++) {
		len = RecordPoreSizeSubroutine1(rPore, &P2(rAct,P2A(actCyl.l,n,0,2),0), 
				&P2(rAct,P2A(actCyl.l,n,1,2),0), dr);
		if (len < AVG2(initRadPore * 2., actF.dia)) { 
			CS = 0; 
			break;
		}
	}
	CONT(CS != 1);
	radPore = initRadPore;
	CS = 0;
	// Increase the size of spheres until it cannot become
	// larger anymore.
	do {
		V3ZERO(drPore);	
		cntOverl = 0;
		for (n = 0; n < actCyl.c; n++) {
			len = RecordPoreSizeSubroutine1(rPore, 
					&P2(rAct,P2A(actCyl.l,n,0,2),0), 
					&P2(rAct,P2A(actCyl.l,n,1,2),0), dr);
			CONT(!(len < radPore + 0.5 * actF.dia));
			VVS3SUB(drPore, dr, 0.1 / len);		
			cntOverl++;
		}
		VV3ADD(rPore, drPore);
		ApplyBoundCondVector(rPore, 0, 0);
		if (cntOverl > maxCntOverl) { CS = 1; }
		BREAK(radPore > dimDomH[0] || radPore > dimDomH[1] 
				|| radPore > dimDomH[2])
		radPore += 0.1;
	} while (CS == 0);
	// Eliminate possible double-counting of the same pores
	for (n = 0; n < pore.c; n++) {
		V3COPY(currPore, rPore);
		currPore[3] = radPore;
		V3COPY(&currPore[4], ind);
		RecordPoreSizeSubroutine3(currPore, &P2A(pore.l,n,0,7), nInitPnt,
				thres, &CS);
		BREAK(CS == 0);
	}
	if (CS == 1) {
		V3COPY(&P2A(pore.l,pore.c,0,7), rPore);
		P2A(pore.l,pore.c,3,7) = radPore;
		V3COPY(&P2A(pore.l,pore.c,4,7), ind);
		(pore.c)++;
	}
  }
  MALLOC(poreCall,int,nCpu);
  if (rank == 0) { MALLOC(poreAll.l,double,totNInitPnt*7); }
  MPI_Gather(&pore.c, 1, MPI_INT, poreCall, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDoubleWoIndWoSort(pore.c, poreCall, pore.l, poreAll.l, 7);

  if (rank == 0) {
	poreAll.c = SumArrInt(poreCall, nCpu);
	for(n = 0; n < poreAll.c; n++) {
		for(k = n + 1; k < poreAll.c; k++) {
			CS = 1;
			RecordPoreSizeSubroutine3(&P2A(poreAll.l,n,0,7), 
					&P2A(poreAll.l,k,0,7), nInitPnt, thres, &CS);
			if (CS == 0) {
				RecordPoreSizeSubroutine2(&poreAll, &k);
			}
		}
	}	
	for(n = 0; n < poreAll.c; n++) {
		poreSize = (int)(P2A(poreAll.l,n,3,7) * L_SCALE_IN_NM 
				/ (double)histIntv);
		histPore.l[poreSize]++;
		if (poreSize > maxPoreSize) { maxPoreSize = poreSize; }
		if (poreSize < minPoreSize) { minPoreSize = poreSize; }
	}
	// Calculate average and standard deviation.
	totPoreVol = 0.;
	avgPoreSize = 0.;
	nPore = 0;
	for(n = minPoreSize; n < maxPoreSize + 1; n++) {
		nPore += histPore.l[n];
	}
	fOut = fopen(GenFileName("PoreSize"), "a");
	fprintf(fOut, "%d\t%d\t0\n", maxPoreSize - minPoreSize + 3, nPore);
	for(n = minPoreSize; n < maxPoreSize + 1; n++) {
	    radPore = ((double)n + 0.5) * histIntv;
	    fprintf(fOut, "%d\t%g\t%d\n", n - minPoreSize, 
				radPore * 2., histPore.l[n]);
	    totPoreVol += 4. / 3. * PI * CUBE(radPore) * histPore.l[n];
		avgPoreSize += ((double)n + 0.5) * histIntv * 2. * histPore.l[n];
		stdPoreSize += SQR(((double)n + 0.5) * histIntv * 2.) * histPore.l[n];
	}
	if (nPore == 0) { avgPoreSize = 0.; stdPoreSize = 0.; }
	else { 
		avgPoreSize /= (double)nPore; 
		stdPoreSize = sqrt(stdPoreSize / (double)nPore - SQR(avgPoreSize));
	}
	fprintf(fOut, "%d\t%g\t%g\n", maxPoreSize - minPoreSize + 1, 
			avgPoreSize, stdPoreSize);
	fprintf(fOut, "%d\t%g\t0\n", maxPoreSize - minPoreSize + 2,
			totPoreVol / (V3PROD(dimDom) * CUBE(L_SCALE_IN_NM)) * 100.); 
	fclose(fOut); 

	free(poreAll.l); 
  }
  free(poreCall);
  free(pore.l);
  free(nActAll);
  free(rAct);
  free(chAct);
  free(actCh);
  free(actCyl.l);
  free(histPore.l);
}

void FilamentEndPositionForceInBundle(void) {
  int n, *nActAll, ind;
  double f, rOri, stf, pos[2], force[2];
  FILE *fOut;

  stf = KS_NPM2S(1e-3);

  FOR_ACTME(n) {
	CONT(!(P2A(act.ch,n,0,nChAc) < 0 && P2A(act.ch,n,1,nChAc) > -1));
	rOri = (act.id[n] < nAct / 2) ? 0 : dimDom[2] - 0.1;
	f = SPRING(stf, P2(act.r,n,2), rOri);
	P2(act.f,n,2) += f;
  }

  if (currTimeStep % (int)(1e-4 / dtReal) == 0) {
	MALLOC(nActAll,int,nCpu);
	if (rank == 0) {
	    MALLOC(rAct,double,nAct*NDIM);
	    MALLOC(chAct,int,nAct*nChAc);
	}
	MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	GatherActChainPosition(nActAll, chAct, rAct);
	if (rank == 0) {
		fOut = fopen(GenFileName("BunFilaEnd"), "a");
		FOR_ACT(n) {
			if (P2A(chAct,n,0,nChAc) < 0 && P2A(chAct,n,1,nChAc) > -1) {
				ind = (n < nAct / 2) ? 0 : 1;
				rOri = (n < nAct / 2) ? 0 : dimDom[2] - 0.1;
				force[ind] = F_S2PN(stf * fabs(P2(rAct,n,2) - rOri));
			}
			if (P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0) {
				ind = (n < nAct / 2) ? 0 : 1;
				pos[ind] = P2(rAct,n,2);
			}
		}
		for(n = 0; n < 2; n++) {
			fprintf(fOut, "%g\t%g\t%g\n", currTimeStep * dtReal * 100,
					pos[n], force[n]);
		}
		fclose(fOut);
		free(rAct);
		free(chAct);
	}
	free(nActAll);
  }
}

/*--------------------- Structural properties of network ---------------------*/

/*----------------------------- Energy calculation ---------------------------*/

// mode = 0: measure energy of entire network
// mode = 1: measure energy of selected portion of the network
//			 to the other 
void RecordMechEnergy(int mode) {
  int m, n, k, locActInd, abpInd, *pArr, side, kind;
  double len, lenEq, ang, fac;
  double E[4], Esum[4], *Eall, dr[2][NDIM + 1];
  double rPos[NDIM], rPos2[NDIM], *rPnt[2];
  FILE *fOut;
  ListInt *pAbpL;
 
  if (rank == 0) { MALLOC(Eall,double,nCpu*4); }
  V4SET_ALL(E, 0.);

  // Bending of actin filament
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
    CONT(!(pArr[0] > -1 && pArr[1] > -1));
	CONT((mode == 1 && iPerc[act.id[n]] == 0) 
			|| (mode == 2 && iSupp[act.id[n]] < 0));
	CalcVec(dr[0], &P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0));
	CalcVec(dr[1], &P2(act.r,iAct[pArr[1]],0), &P2(act.r,n,0));
	ang = V3ANG(dr[0], dr[1]);
	E[0] += 0.5 * actF.bend.stf * SQR(ang);
  }

  // Extension of actin filament (double-counting)
  FOR_ACTME(n) {
	pArr = &P2A(act.ch,n,0,nChAc);
	CONT(pArr[0] < 0);
	CONT((mode == 1 && iPerc[act.id[n]] == 0) 
			|| (mode == 2 && iSupp[act.id[n]] < 0));
	len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[pArr[0]],0), 0);
	E[1] += 0.5 * actF.spr.stf * SQR(len - actF.spr.eq);
  }

  // Bending of ABP (actin-ABP-actin)
  FOR_ABPME(n) {
	pArr = &P2A(abp.ch,n,0,nChAb);
	CONT(!(pArr[0] > -1 && pArr[1] > -1));
	kind = pArr[2];
	CONT(ISMTF(kind));
	CONT((mode == 1 && iPerc[nAct + abp.id[n]] == 0) 
			|| (mode == 2 && iSupp[nAct + abp.id[n]] < 0));
	for(k = 0; k < 2; k++) {
		dr[k][NDIM] = CalcVecDistActinAbp(dr[k], pArr[k], abp.id[n], 0);
	}
	V3REVSIGN(dr[1]);
	ang = V3ANG(dr[0], dr[1]);
	E[2] += 0.5 * abpF.bend[kind].stf * SQR(ang - abpF.bend[kind].eq);
  }
  // Bending of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		abpInd = abp.id[n];
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
		CONT((mode == 1 && iPerc[nAct + abpInd] == 0) 
				|| (mode == 2 && iSupp[nAct + abpInd] < 0));
	    CONT(!(pArr[3] > -1 && pArr[4] > -1));
		CalcVec(dr[0], &P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0));
		CalcVec(dr[1], &P2(abp.r,iAbp[pArr[4]],0), &P2(abp.r,n,0));
		ang = V3ANG(dr[0], dr[1]);
		E[2] += 0.5 * motSA.bend.stf * SQR(ang);
	}
  }

  // Extension of ABP
  FOR_ABPME(n){
	CONT(ISABPM(n));
	pArr = &P2A(abp.ch,n,0,nChAb);
	kind = pArr[2];
	CONT(mode == 1 && iPerc[nAct + abp.id[n]] == 0);
	for(k = 0; k < 2; k++) {
		CONT(pArr[k] < 0);
		locActInd = iAct[pArr[k]];
		side = FindAbpActinChain(locActInd, abp.id[n], 0);
		rPnt[0] = &P2(act.r,locActInd,0);
		rPnt[1] = &P2(act.r,iAct[P2A(act.ch,locActInd,0,nChAc)],0);
		CalcPosOnActSegSide(rPnt[0], rPnt[1], rPos, side);
  		CalcVec(dr[0], rPos, &P2(abp.r,n,0));
		if (ISMTF(kind)) {
			CalcVec(dr[1], rPnt[0], rPnt[1]);
			fac = V3DOT(dr[0], dr[1]) / V3LEN_SQ(dr[1]);
			VS3SUB(rPos2, rPos, dr[1], fac);
			ApplyBoundCondVector(rPos2, -1, 0);
			CalcVec(dr[0], rPos2, &P2(abp.r,n,0));
			CalcVec(dr[1], rPos, rPos2);			
		}
		len = V3LEN(dr[0]);
		E[3] += 0.5 * abpF.spr[kind].stf 
				* SQR(len - abpF.spr[kind].eq);
		if (ISMTF(kind)) {	
			len = V3LEN(dr[1]);
			E[3] += 0.5 * abpF.spr[kind].stf2 * SQR(len);
		}
	}
  }
  
  // Extension of motor-assembled structures
  if (motSA.gTgl != 0) {
	for(n = 0; n < nAbpMe; n++) {
		abpInd = abp.id[n];
		pArr = &P2A(abp.ch,n,0,nChAb);
		CONT(pArr[2] != 2);
		CONT((mode == 1 && iPerc[nAct + abpInd] == 0) 
				|| (mode == 2 && iSupp[nAct + abpInd] < 0));
		CONT(pArr[3] < 0);
		len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[pArr[3]],0), 0);
		lenEq = (abp.mId[n] == 0 && abp.mId[iAbp[pArr[3]]] == 0) 
				? motSA.cenDist : motSA.spr.eq;
		E[3] += 0.5 * motSA.spr.stf * SQR(len - lenEq);
	}
  }

  MPI_Gather(E, 4, MPI_DOUBLE, Eall, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    V4SET_ALL(Esum, 0.); 
    for(n = 0; n < nCpu; n++) {
        for(k = 0; k < 4; k++) {
            Esum[k] += Eall[4 * n + k];
        }
    }
	if (mode == 0) {
	    fOut = fopen(GenFileName("MechEall"), "a");
	}
	else if (mode == 1) {
	    fOut = fopen(GenFileName("MechEperc"), "a");
	}
	else {
	    fOut = fopen(GenFileName("MechEsupp"), "a");
	}
    fprintf(fOut, "%lld\t", currTimeStep);
	Fprintf1dArrayDouble(fOut, Esum, 4, 0);
    fclose(fOut);
    free(Eall);
  }
}

/*----------------------------- Energy calculation ---------------------------*/

/*------------------------ Internal stress measurement -----------------------*/

void RecordElasViscStressSubroutine1(double *rPnt, double *widDiv, int *dir, 
		double minDimDomH, int *indDomDiv, int *nSec, int mode) {
  int k;
  double dr[NDIM], ang, len;
  
  V3SET_ALL(indDomDiv, -1);
  if (mode == 0) {  
	FOR_NDIM(k) { 
		indDomDiv[k] = (int)((rPnt[k] - rGrid[k][0]) /  widDiv[k]);
		if (indDomDiv[k] < 0) { 
			indDomDiv[k] = (pbc[k] == 1) ? nSec[k] - 1 : 0;
		}
		else if (indDomDiv[k] >= nSec[k]) { 
			indDomDiv[k] = (pbc[k] == 1) ? 0 : nSec[k] - 1;
		}
	}
  }
  else {
	FOR_NDIM(k) {
		dr[k] = rPnt[k] - rGrid[k][0];
		CONT(k == dir[2]); 
		dr[k] -= dimDomH[k];
	}
	len = sqrt(SQR(dr[dir[0]]) + SQR(dr[dir[1]]));
	if (len >= minDimDomH) {
		return;
	}
	indDomDiv[0] = (int)(len / widDiv[0]);
	TrimIntVal(indDomDiv[0], 0, nSec[0] - 1);

	ang = atan(dr[dir[1]] / dr[dir[0]]);
	if (dr[dir[0]] < 0) { ang += PI; }
	if (ang < 0) { ang += 2 * PI; }
	indDomDiv[1] = (int)(ang / widDiv[1]);
	TrimIntVal(indDomDiv[1], 0, nSec[1] - 1);

	indDomDiv[2] = (int)(dr[dir[2]] / widDiv[2]);
	if (indDomDiv[2] < 0) { 
		indDomDiv[2] = (pbc[dir[2]] == 1) ? nSec[2] - 1 : 0;
	}
	else if (indDomDiv[2] >= nSec[2]) { 
		indDomDiv[2] = (pbc[dir[2]] == 1) ? 0 : nSec[2] - 1;
	} 
  }
}

void RecordElasViscStressSubroutine2(double stiff, double lenEq, 
		int *indDomDiv, double *rPnt, double *rPnt2, double *widDiv, int *dir, 
		double minDimDomH, int *nSec, int mode) {
  int n, k, dir2, CS, mm, l, idx, sft[NDIM];
  int ind[NDIM], indAdj[NDIM], indDomDivAdj[NDIM], begin[NDIM], end[NDIM];
  double len, len2, slope, ang, aa, bb, cc, dd, diff, low[NDIM], hi[NDIM];
  double dr[NDIM], dR[2][NDIM], rSec[NDIM], rAdj[2][NDIM], rCrs[NDIM];

  RecordElasViscStressSubroutine1(rPnt2, widDiv, dir, minDimDomH, 
		indDomDivAdj, nSec, mode);
  FOR_NDIM(k) {
	BREAK(indDomDiv[k] != indDomDivAdj[k]);
  }
  if (k == NDIM) { return; }
  FOR_NDIM(k) {
	rAdj[1][k] = rPnt2[k] - rGrid[k][0];
	rAdj[0][k] = rPnt[k] - rGrid[k][0];
	diff = rPnt2[k] - rPnt[k];
	if (pbc[k] == 1 && fabs(diff) > dimDomH[k]) {
		sft[k] = 1;
		if (diff >= dimDomH[k]) { rAdj[1][k] -= dimDom[k]; }
		else if (diff < -1. * dimDomH[k]) { rAdj[0][k] -= dimDom[k]; }
	}
	else { sft[k] = 0; }
	dr[k] = rAdj[1][k] - rAdj[0][k];
  }	

  FOR_NDIM(k) {
	if (((mode == 0 || (mode == 1 && k == dir[2])) && sft[k] == 1) 
			|| (mode == 1 && k == dir[1] && (double)abs(indDomDivAdj[k] 
			- indDomDiv[k]) >= (double)nSec[k] / 2.)) {
		begin[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDivAdj[k] 
				: indDomDiv[k];
		end[k] = ((indDomDiv[k] <= indDomDivAdj[k]) ? indDomDiv[k] 
				: indDomDivAdj[k]) + nSec[k];
	}
	else {
		begin[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDiv[k] 
				: indDomDivAdj[k];
		end[k] = (indDomDiv[k] <= indDomDivAdj[k]) ? indDomDivAdj[k] 
				: indDomDiv[k];
	}
  }
  for(ind[0] = begin[0]; ind[0] <= end[0]; ind[0]++) {
	for(ind[1] = begin[1]; ind[1] <= end[1]; ind[1]++) {
		for(ind[2] = begin[2]; ind[2] <= end[2]; ind[2]++) {
			FOR_NDIM(k) {
				indAdj[k] = ind[k] % nSec[k];
			}
			V3IND_BACK_INT(idx, indAdj, nSec);
			V3MUL(rSec, indAdj, widDiv);
			if (mode == 0) {
				FOR_NDIM(k) {
					CONT(pbc[k] == 0 && indAdj[k] == 0);
					CONT((rAdj[0][k] - rSec[k]) * (rAdj[1][k] - rSec[k]) >= 0);
					CS = 1;
					for(n = 1; n < 3; n++) {
						dir2 = (k + n) % NDIM;
						rCrs[dir2] = (rSec[k] - rAdj[0][k]) * dr[dir2] / dr[k] 
								+ rAdj[0][dir2];
						if (rCrs[dir2] < rSec[dir2] || rCrs[dir2] >= rSec[dir2] 
								+ widDiv[dir2]) {
							CS = 0;
							break;
						}
					}
					CONT(CS == 0);
					rCrs[k] = rSec[k];
					len = V3LEN(dr);
					P2A(elasStre,idx,k,NDIM) += stiff * (len - lenEq) 
							* fabs(dr[k]) / len;

				}  
			}	
			else {
				FOR_NDIM(k) {
					for(n = 0; n < 2; n++) {
						dR[n][dir[k]] = rAdj[n][dir[k]] 
								- (k == dir[2] ? 0 : dimDomH[dir[k]]);
					}
					if (dR[0][dir[k]] < dR[1][dir[k]]) {
						low[k] = dR[0][dir[k]];
						hi[k] = dR[1][dir[k]];
					}
					else {
						low[k] = dR[1][dir[k]];
						hi[k] = dR[0][dir[k]];
					}
				}
				dd = dR[1][dir[0]] * dR[0][dir[1]] 
						- dR[0][dir[0]] * dR[1][dir[1]];
				FOR_NDIM(k) {
					CS = 1;
					if (k == 0) {
						len = sqrt(SQR(dR[0][dir[0]]) + SQR(dR[0][dir[1]]));
						len2 = sqrt(SQR(dR[1][dir[0]]) + SQR(dR[1][dir[1]]));
						CONT((len - rSec[k]) * (len2 - rSec[k]) >= 0);
						aa = 1 + SQR(dr[dir[1]] / dr[dir[0]]);
						bb = dr[dir[1]] * dd / SQR(dr[dir[0]]);
						cc = SQR(dd / dr[dir[0]]) - SQR(rSec[k]);
						for(n = 0; n < 2; n++) {
							rCrs[dir[0]] = (-1. * bb + (n == 0 ? 1. : -1.) 
									* sqrt(bb * bb - aa * cc)) / aa;
							BREAK(rCrs[dir[0]] >= low[dir[0]] 
									&& rCrs[dir[0]] < hi[dir[0]]);
						}
						for(n = 1; n < 3; n++) {
							rCrs[dir[n]] = dr[dir[n]] / dr[dir[0]] 
									* (rCrs[dir[0]] - dR[0][dir[0]]) 
									+ dR[0][dir[n]];
						}
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						if (ang < rSec[1] || ang >= rSec[1] + widDiv[1]
								|| rCrs[dir[2]] < rSec[2] 
								|| rCrs[dir[2]] >= rSec[2] + widDiv[2]) {
							CS = 0;
						}
					}
					else if (k == 1) {
						if (fabs(rSec[1] - 0.5 * PI) < 0.01 
								|| fabs(rSec[1] - 1.5 * PI) < 0.01) {
							CONT(fabs(dr[dir[0]]) < 0.01);
							CONT(dR[0][dir[0]] * dR[1][dir[0]] >= 0.);
							rCrs[dir[0]] = 0.;
						}
						else {
							slope = tan(rSec[1]);
							CONT(fabs(slope - dr[dir[1]] / dr[dir[0]]) 
									< 0.01);
							rCrs[dir[0]] = dd / (slope * dr[dir[0]] 
									- dr[dir[1]]);
							CONT(rCrs[dir[0]] < low[dir[0]] 
									|| rCrs[dir[0]] >= hi[dir[0]]);
						}
						for(n = 1; n < 3; n++) {
							rCrs[dir[n]] = dr[dir[n]] / dr[dir[0]] 
									* (rCrs[dir[0]] - dR[0][dir[0]]) 
									+ dR[0][dir[n]];
						}
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						CONT(ang < rSec[1] || ang >= rSec[1] + widDiv[1]);

						len = sqrt(SQR(rCrs[dir[0]]) + SQR(rCrs[dir[1]]));
						if (len < rSec[0] || len >= rSec[0] + widDiv[0]
								|| rCrs[dir[2]] < rSec[2] 
								|| rCrs[dir[2]] >= rSec[2] + widDiv[2]) {
							CS = 0;
						}
					}
					else {
						CONT(pbc[k] == 0 && indAdj[k] == 0);
						CONT((rAdj[0][k] - rSec[k]) 
								* (rAdj[1][k] - rSec[k]) >= 0);
						CS = 1;
						rCrs[dir[0]] = (rSec[k] - rAdj[0][k]) * dr[dir[0]] 
								/ dr[k] + dR[0][dir[0]];
						rCrs[dir[1]] = (rSec[k] - rAdj[0][k]) * dr[dir[1]] 
								/ dr[k] + dR[0][dir[1]];
						len = sqrt(SQR(rCrs[dir[0]]) + SQR(rCrs[dir[1]]));
						ang = atan(rCrs[dir[1]] / rCrs[dir[0]]);
						if (rCrs[dir[0]] < 0) { ang += PI; }
						if (ang < 0) { ang += 2 * PI; }
						if (len < rSec[0] || len >= rSec[0] + widDiv[0] 
								|| ang < rSec[1] || ang >= rSec[1] 
								+ widDiv[1]) {
							CS = 0;
						}
					}
					CONT(CS == 0);
					len = V3LEN(dr);
					P2A(elasStre,idx,k,NDIM) += stiff * (len - lenEq) 
							* fabs(dr[k]) / len;
				}  
			}
		}
	}
  }
}

// mode = 0: Cartesian, 1: cylindrical
void RecordElasViscStress(int period, int *nSec, int mode) {
  int mm, m, n, k, end, nSecP, nSecPN, CS, side;
  int abpInd, *pArr, ind, indDomDiv[NDIM], dir[NDIM];
  double fac, minDimDomH, widDiv[NDIM];
  double dr[NDIM], dr2[NDIM], drAxis[NDIM], rPos[NDIM], rPos2[NDIM];
  double *rPnt, *fPnt, *rPntAdj, *rPntAdjAdj, *streAll, *streSum, *sendArr;
  FILE *fOut;

  minDimDomH = POS_LARGE_VALUE;
  nSecP = V3PROD(nSec);
  if (mode == 0) {
	V3DIV(widDiv, dimDom, nSec);
  }
  else {
	V3SET(dir, 0, 1, 2);
	FOR_NDIM(k) {
		CONT(!(k != dir[2] && dimDomH[k] < minDimDomH));
		minDimDomH = dimDomH[k];
	}
	widDiv[0] = minDimDomH / nSec[0];
	widDiv[1] = 2 * PI / nSec[1];
	widDiv[2] = dimDom[dir[2]] / nSec[2];
  }

  for(n = 0; n < nActMe + nAbpMe; n++) {
	CS = 1;
	if (n < nActMe) {
		CONT(ISACTM(n));
		fPnt = &P2(act.f,n,0);
		rPnt = &P2(act.r,n,0);
		if (act.fix[n] > -1 && bndReb.gTglSpr == 0) { CS = -1; }
	}
	else {
		abpInd = n - nActMe;
		CONT(ISABPM(abpInd));
		pArr = &P2A(abp.ch,abpInd,0,nChAb);
		fPnt = &P2(abp.f,abpInd,0);
		rPnt = &P2(abp.r,abpInd,0);
	}
	RecordElasViscStressSubroutine1(rPnt, widDiv, dir, minDimDomH, 
			indDomDiv, nSec, mode);
	V3IND_BACK_INT(ind, indDomDiv, nSec);
	CONT(indDomDiv[0] < 0);
	// Viscous
	if (CS == 1 && isnan(fPnt[0]) == 0) {
		VV3ADD(&P2A(viscStre,ind,0,NDIM), fPnt);
	}
	// Elastic
	if (n < nActMe) {
		end = 1;
	}
	else {
		end = (motSA.gTgl != 0) ? nChAb - 1 : 2;
	}
	for(m = 0; m < end; m++) {
		if (n < nActMe) {
			CONT(P2A(act.ch,n,0,nChAc) < 0);
			rPntAdj = &P2(act.r,iAct[P2A(act.ch,n,0,nChAc)],0);
		}
		else {
			CONT(m == 2);
			CONT(pArr[m] < 0);
			CONT(pArr[2] != 2 && m > 2);
			rPntAdj = (m < 2) ? &P2(act.r,iAct[pArr[m]],0) 
					: &P2(abp.r,iAbp[pArr[m]],0);
		}
		if (n >= nActMe && m < 2) {
			side = FindAbpActinChain(iAct[pArr[m]], abp.id[abpInd], 0);
			rPntAdjAdj = &P2(act.r,
					iAct[P2A(act.ch,iAct[pArr[m]],0,nChAc)],0);
			CalcPosOnActSegSide(rPntAdj, rPntAdjAdj, rPos, side);
  			CalcVec(dr, rPos, rPnt);
			if (ISMTF(pArr[2])) {
				CalcVec(drAxis, rPntAdj, rPntAdjAdj);
				fac = V3DOT(dr, drAxis) / V3LEN_SQ(drAxis);
				VS3SUB(rPos2, rPos, drAxis, fac);
				ApplyBoundCondVector(rPos2, -1, 0);
			}
			RecordElasViscStressSubroutine2(abpF.spr[pArr[2]].stf,  
					abpF.spr[pArr[2]].eq, indDomDiv, rPnt, ((ISMTF(pArr[2]))
					? rPos2 : rPos), widDiv, dir, minDimDomH, nSec, mode);
			if (ISMTF(pArr[2])) {
				RecordElasViscStressSubroutine2(abpF.spr[pArr[2]].stf2,  
						0., indDomDiv, rPnt, rPos, widDiv, dir, minDimDomH,
						nSec, mode);
			}
		}
		else {
			if (n < nActMe) { 
				RecordElasViscStressSubroutine2(actF.spr.stf, actF.spr.eq, 
						indDomDiv, rPnt, rPntAdj, widDiv, dir, minDimDomH, 
						nSec, mode);
			}
			else {
				RecordElasViscStressSubroutine2(motSA.spr.stf, motSA.spr.eq,
						indDomDiv, rPnt, rPntAdj, widDiv, dir, minDimDomH, 
						nSec, mode);
			}
		}
	}
  }
  nSecPN = nSecP * NDIM;
  if (currTimeStep % period == 0) {
    MALLOC(sendArr,double,nSecPN * 2);
    if (rank == 0) {
        MALLOC(streAll,double,nSecPN * 2 * nCpu);
        MALLOC(streSum,double,nSecPN * 2);
    }
    for(n = 0; n < nSecPN; n++) {
        sendArr[n] = elasStre[n];
        sendArr[n + nSecPN] = viscStre[n];
    }
    MPI_Gather(sendArr, nSecPN * 2, MPI_DOUBLE, streAll,
            nSecPN * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for(n = 0; n < nSecPN * 2; n++) {
            streSum[n] = 0.;
            for(k = 0; k < nCpu; k++) {
                streSum[n] += streAll[k * nSecPN * 2 + n];
            }
            streSum[n] /= (double)period;
            streSum[n] = F_S2PN(streSum[n]);
        }
        fOut = fopen(GenFileName("ElasViscStre"), "a");
        if (currTimeStep / period == 1) {
            fprintf(fOut, "%d\t", nSecP * 2);
            Fprintf1dArrayIntWoRet(fOut, nSec, NDIM);
            FOR_NDIM(k) {
                fprintf(fOut, "%g\t", ((mode == 1 && k == 1)
                        ? RAD2DEG(widDiv[k]) : L_S2UM(widDiv[k])));
            }
            fprintf(fOut, "\n");
        }

        for(n = 0; n < nSecP * 2; n++) {
            V3IND_ASSIGN_INT((n % nSecP), nSec, 1, indDomDiv);
            fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\n", n, indDomDiv[0],
                    indDomDiv[1], indDomDiv[2], P2(streSum,n,0),
                    P2(streSum,n,1), P2(streSum,n,2));
        }
        fclose(fOut);
        free(streSum);
        free(streAll);
    }
    free(sendArr);
    for(n = 0; n < nSecPN; n++) {
        viscStre[n] = 0.;
        elasStre[n] = 0.;
    }
  }
}
/*------------------------ Internal stress measurement -----------------------*/


/*-------------------------- Instantaneous information -----------------------*/

// Record the instantaneous length of all chains
void RecordInstantChainLength(void) {
  int m, n, k, ind, begin, end, deg;
  int *cntChLenAct, *cntChLenAbp, *pCnt, *nActAll, *nAbpAll, *nMbAll, *pArr;
  double len, min, max, *chLen, *pArr2; 
  FILE *fOut;

  deg = 50;
  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(nMbAll,int,nCpu);
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(rMb,double,nMb*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(chMb,int,nMb*nChMb);
  }
  // Gather information about the positions and chains of actins and ABPs
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  GatherMembChainPosition(nMbAll, chMb, rMb);

  if (rank == 0) {
	MALLOC(cntChLenAct, int, deg + 1);
	MALLOC(cntChLenAbp, int, deg + 1);
	MALLOC(chLen, double, nAct * nChAc);
	memset(cntChLenAct, 0, sizeof(int) * (deg + 1));
	memset(cntChLenAbp, 0, sizeof(int) * (deg + 1));

	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
	    for (k = 0; k < nChAc; k++) {
	        if (pArr[k] > -1) {
				if (k < 2) {
					len = CalcDist(&P2(rAct,n,0), &P2(rAct,pArr[k],0), 0);
				}
				else {
					len = CalcDistActinAbp(n, pArr[k], 1);
				}
	            P2A(chLen,n,k,nChAc) = len;
	        }
			else {
	            P2A(chLen,n,k,nChAc) = 0;
			}
	    }
	} 
	// Record the distribution of chain length of actin and ABPs
	for(m = 0; m < 2; m++) {
		if (m == 0) { 
			begin = 0;
			end = 2;
			pCnt = cntChLenAct;
			fOut = fopen(GenFileName("ActInstLenS"), "a");
		}
		else {
			begin = 2;
			end = nChAc;
			pCnt = cntChLenAbp;
			fOut = fopen(GenFileName("AbpInstLenS"), "a");
		}
		fprintf(fOut, "%d\t0\t0\n", deg + 1);
	    min = POS_LARGE_VALUE;	
		max = NEG_LARGE_VALUE; 
		FOR_ACT(n) {
			pArr2 = &P2A(chLen,n,0,nChAc);
			for(k = begin; k < end; k++) {
				if (pArr2[k] > 0.) {
					if (max < pArr2[k]) { max = pArr2[k]; }
					if (min > pArr2[k]) { min = pArr2[k]; }
				}
			}
		}
		FOR_ACT(n) {
			for(k = begin; k < end; k++) {
				if (P2A(chLen,n,k,nChAc) > 0.) {
					ind = (int)((P2A(chLen,n,k,nChAc) - min) 
							/ (max - min) * deg);
					pCnt[ind]++;
				}
			}
		}
		for (n = 0; n < deg + 1; n++) {
			if (m == 0) { pCnt[n] /= 2; } 
			fprintf(fOut, "%d\t%g\t%d\n", n, ((double)n / (double)deg 
					* (max - min) + min) * L_SCALE_IN_NM, pCnt[n]);
		}
		fclose(fOut);
	}

	// Record individual chain length

	// actin
	fOut = fopen(GenFileName("ActInstLenD"), "a");
	fprintf(fOut, "%d\t", nAct);
	Fprintf1dFillerInt(fOut, 0, nChAc, 0);
	Fprintf2dArrayDouble(fOut, chLen, nAct, nChAc, 0, 0);
	fclose(fOut);

	// ABP
	fOut = fopen(GenFileName("AbpInstLenD"), "a");
	fprintf(fOut, "%d\t", nAbp);
	Fprintf1dFillerInt(fOut, 0, 4, 0);
	FOR_ABP(n) {
		fprintf(fOut, "%d\t", n);
		for(k = 0; k < 2; k++) {
			ind = P2A(chAbp,n,k,nChAb);
			if (ind > -1) {
				len = CalcDistActinAbp(ind, n, 1);
				fprintf(fOut, "%g\t", len);
			}
			else { fprintf(fOut, "%g\t", 0.); }
		}
		for(k = 0; k < 2; k++) {
			ind = P2A(chAbp,n,k + 3,nChAb);
			if (ind > -1 && P2A(chAbp,n,2,nChAb) == 2) {
				len = CalcDist(&P2(rAbp,n,0), &P2(rAbp,ind,0), 0);
				fprintf(fOut, "%g\t", len);
			}
			else { fprintf(fOut, "%g\t", 0.); }
		}	
		fprintf(fOut, "\n");
	}
	fclose(fOut);

	if (gTglMb != 0) {
		// Membrane
		fOut = fopen(GenFileName("MembInstLenD"), "a");
		fprintf(fOut, "%d\t", nMb);
		Fprintf1dFillerInt(fOut, 0, nChMb, 0);
		FOR_MB(n) {
			fprintf(fOut, "%d\t", n);
			for(k = 0; k < nChMb - nMbAct; k++) {
				ind = P2A(chMb,n,k,nChMb);
				if (ind > -1) {
					len = CalcDist(&P2(rMb,n,0), &P2(rMb,ind,0), 0);
					fprintf(fOut, "%g\t", len);
				}
				else { fprintf(fOut, "%g\t", 0.); }
			}
			for(k = 0; k < nMbAct; k++) {
				ind = P2A(chMb,n,k + nChMb - nMbAct,nChMb);
				if (ind > -1) {
					len = CalcDist(&P2(rMb,n,0), &P2(rAct,ind,0), 0);
					fprintf(fOut, "%g\t", len);
				}
				else { fprintf(fOut, "%g\t", 0.); }
			}	
			fprintf(fOut, "\n");
		}
		fclose(fOut);
	}

	free(cntChLenAct);
	free(cntChLenAbp);
	free(chLen);
	free(rAct);
	free(rAbp);
	free(rMb);
	free(chAct);
	free(chAbp);
	free(chMb);
  }
  free(nActAll);
  free(nAbpAll);
  free(nMbAll);
}

// Record the instantaneous forces actin on actin and ABP
void RecordInstantForces(void) {
  FILE *fOut;
  int *nActAll, *nAbpAll;
  double *fActin, *fAbp;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(fActin,double,nAct*NDIM);
	MALLOC(fAbp,double,nAbp*NDIM);
  }
  // Gather information about forces
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nActMe, nActAll, NDIM, act.id, act.f, fActin);
  Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM, abp.id, abp.f, fAbp);

  if (rank == 0) {
	fOut = fopen (GenFileName("AllInstFor"), "a");
	fprintf(fOut, "%d\t%d\t", nAct, nAbp);
	Fprintf1dFillerInt(fOut, 0, 2, 0);
	Fprintf2dArrayDouble(fOut, fActin, nAct, NDIM, 0, 0);
	Fprintf2dArrayDouble(fOut, fAbp, nAbp, NDIM, nAct, 0);
	fclose (fOut);
	free(fActin);
	free(fAbp);
  }
  free(nActAll);
  free(nAbpAll);
}

// Record the orientation of cylindrical segments of actin and ABP
void RecordActinAbpInstantOrient(int period) {
  FILE *fOut;
  int n, k, *nActAll, *nAbpAll, cntAct, cntAbp, *pArr, nCol;
  double dr[2][NDIM + 1], dr2[NDIM], angBend, ang90;

  nCol = (NDIM + 2) * 2 + 3;
  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
  }
  // Gather information about positions and chains.
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);

  if (rank == 0) {
	fOut = fopen(GenFileName("AllInstOri"), "a");

	cntAct = 0;
	FOR_ACT(n) {
		CONT(P2A(chAct,n,0,nChAc) < 0);
		cntAct++;
	}
	cntAbp = 0;
	FOR_ABP(n) {
		CONT(!(P2A(chAbp,n,0,nChAb) > -1 || P2A(chAbp,n,1,nChAb) > -1));
		cntAbp++;
	}
	fprintf(fOut, "%d\t%d\t", cntAct, cntAbp);
	Fprintf1dFillerInt(fOut, 0, nCol - 2, 0);
	// actin
	cntAct = 0;
	FOR_ACT(n) {
		CONT(P2A(chAct,n,0,nChAc) < 0);
		dr[0][NDIM] = CalcVecDist(dr[0], &P2(rAct,n,0), 
				&P2(rAct,P2A(chAct,n,0,nChAc),0), 0);
		fprintf(fOut, "%d\t", cntAct);
		Fprintf1dArrayDoubleWFil(fOut, dr[0], 4, 0, nCol - 1);
		cntAct++;
	}
	// ABP
	cntAbp = 0;
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		CONT(!(pArr[0] -1 || pArr[1] > -1));
		fprintf(fOut, "%d\t", cntAct + cntAbp);
		for(k = 0; k < 2; k++) {
			if (pArr[k] -1) {
				dr[k][NDIM] = CalcVecDistActinAbp(dr[k], pArr[k], n, 1);
				CalcVec(dr2, &P2(rAct,P2A(chAct,pArr[k],0,nChAc),0), 
						&P2(rAct,P2A(chAct,pArr[k],1,nChAc),0));
				ang90 = V3ANG(dr[k],dr2);
			}
			else { 
				V4SET_ALL(dr[k], 0.); 
				ang90 = 0.;
			}
			Fprintf1dArrayDoubleWoRet(fOut, dr[k], NDIM + 1);
			fprintf(fOut, "%g\t", RAD2DEG(ang90));
		}
		if (V3LEN_SQ(dr[0]) != 0 && V3LEN_SQ(dr[1]) != 0) {
			angBend = V3ANG(dr[0], dr[1]);
		}
		else { angBend = 0.; }
		fprintf(fOut, "%g\t%d\n", RAD2DEG(angBend), P2A(chAbp,n,2,nChAb));
		cntAbp++;
	}
	fclose (fOut);
	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
  }
  free(nActAll);
  free(nAbpAll);
}

/*-------------------------- Instantaneous information -----------------------*/

/*--------------------------- Accumulated information ------------------------*/

// Accumulate forces acting on actin and ABP. This information is not recorded 
// to file, but used for coloring methods via VMD.
void RecordAccuLengthForces(void) {
  int n, k, cntCh;
  int actInd, nextInd;
  double lenAccu, len;
  // actin
  FOR_ACTME(n) {
	CONT(ISACTM(n));
	VV3ADD(&P2A(recAct.allF,n,0,NDIM + 1), &P2(act.f,n,0));
	P2A(recAct.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(act.f,n,0));
	nextInd = P2A(act.ch,n,0,nChAc);
	if (nextInd > -1) {
		len = CalcDist(&P2(act.r,n,0), &P2(act.r,iAct[nextInd],0), 0);
		recAct.len[n] += len;
		recAct.cnt[n] += 1;
	}
	else {
		recAct.len[n] = 0.;
		recAct.cnt[n] = 0;
	}
  } 
  
  // ABP
  FOR_ABPME(n) {
	CONT(ISABPM(n));
	VV3ADD(&P2A(recAbp.allF,n,0,NDIM + 1), &P2(abp.f,n,0));
	P2A(recAbp.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(abp.f,n,0));
	cntCh = 0;
	lenAccu = 0.;
	for(k = 0; k < 2; k++) {
		actInd = P2A(abp.ch,n,k,nChAb);
		CONT(actInd < 0);
		len = CalcDistActinAbp(actInd, abp.id[n], 0);
		lenAccu += len / abpF.spr[K_ABP(n)].eq;
		cntCh++;
	}
	if (cntCh > 0) {
		P2A(recAbp.len,n,0,recAbp.nL) += lenAccu / (double)cntCh;
		recAbp.cnt[n] += 1;
	}
	else {
		P2A(recAbp.len,n,0,recAbp.nL) = 0.;
		if (!(ISMTF(K_ABP(n)))) {
			recAbp.cnt[n] = 0;
		}
	}
	CONT(!(ISMTF(K_ABP(n))));
	for(k = 0; k < 2; k++) {
		nextInd = P2A(abp.ch,n,k + 3,nChAb);
		if (nextInd > -1) {
			len = CalcDist(&P2(abp.r,n,0), &P2(abp.r,iAbp[nextInd],0), 0);
			P2A(recAbp.len,n,1,recAbp.nL) += len;
			if (cntCh == 0) { 
				recAbp.cnt[n] += 1;
			}
		}
		else {
			P2A(recAbp.len,n,1,recAbp.nL) = 0.;
			if (cntCh == 0) {
				recAbp.cnt[n] = 0;
			}
		}
	}
  }

  // Membrane
  FOR_MBME(n) {
	VV3ADD(&P2A(recMb.allF,n,0,NDIM + 1), &P2(memb.f,n,0));
	P2A(recMb.allF,n,NDIM,NDIM + 1) += V3LEN(&P2(memb.f,n,0));
	cntCh = 0;
	lenAccu = 0.;
	for(k = 0; k < nChMb - nMbAct; k++) {	
		nextInd = P2A(memb.ch,n,k,nChMb);
		CONT(nextInd < 0);
		len = CalcDist(&P2(memb.r,n,0), &P2(memb.r,iMb[nextInd],0), 0);
		lenAccu += len;
		cntCh++;
	}
	recMb.len[n] += lenAccu / (double)cntCh;
	recMb.cnt[n] += 1;
  } 
}

void ResetAccuLengthForces(void) {
  int n;

  FOR_ACTC(n) {
	recAct.sprF[n] = 0.;
	recAct.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAct.allF,n,0,NDIM + 1), 0.);
	recAct.len[n] = 0.;
	recAct.cnt[n] = 0.;
  }
  FOR_ABPC(n) {
	recAbp.sprF[n] = 0.;
	recAbp.bendF[n] = 0.;
	V4SET_ALL(&P2A(recAbp.allF,n,0,NDIM + 1), 0.);
	SetAllValue1dArrayDouble(&P2A(recAbp.len,n,0,recAbp.nL), recAbp.nL, 0.);
	recAbp.cnt[n] = 0.;
  }
}

// Calculate and record the longitudinal forces acting on ACP or motor. This
// information is very useful to analyze the activity of motors since 
// the longitudinal forces are known to affect the rate of walking events.
void RecordAccuLongSpringForces(int period) {
  int n, k, *nAbpAll;
  double *f, *fAbp;
  FILE *fOut;

  // Accmulate the calculated longitudinal forces.
  FOR_ABPME(n) {
	CONT(ISABPM(n));
    for (k = 0; k < 2; k++) {
		CONT(P2A(abp.ch,n,k,nChAb) < 0);
        // find force
        f = &P2A(recInstSprFabp,n,k * 2,4);
		P2A(recLongSprFabp,n * 2 + k,0,2) += f[0];
		P2A(recLongSprFabp,n * 2 + k,1,2) += f[1];
	}
  }
  if (currTimeStep % period == 0) {
	MALLOC(nAbpAll,int,nCpu);
	if (rank == 0) { MALLOC(fAbp,double,nAbp * 2 * 2); }
	MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
	Gather2dArrayDouble(nAbpMe, nAbpAll, 4, abp.id, recLongSprFabp, fAbp);	
    for(n = 0; n < nAbpC * 4; n++) {
        recLongSprFabp[n] = 0.;
    }

	if (rank == 0) {
		for(n = 0; n < nAbp * 4; n++) {
			fAbp[n] = F_S2PN(fAbp[n]) / (double)period;
		}
		fOut = fopen(GenFileName("AbpLongFor"), "a");
		if (currTimeStep / period == 1) {
    		Fprintf1dArrayIntWFil(fOut, &nAbp, 1, 0, 5);
		}
		Fprintf2dArrayDouble(fOut, fAbp, nAbp, NDIM + 1, 0, 0);
		fclose(fOut);
		free(fAbp);
	}
	free(nAbpAll);
  }
}

/*--------------------------- Accumulated information ------------------------*/

/*------------------ Information for filaments and segments  -----------------*/

void RecordIndvSegFilaInformation(int period) {
  int n, k, l, curr, ind, iPercSum, iSuppSum, cnt, cntFor, CS, side;
  int abpInd, actInd, prevAbpInd, prevActInd, sft[NDIM], cntXlink[3];
  int *iFila, *nActAll, *nAbpAll, *cntAct, *cntAbp, *iMot;
  int *binAct, *binAbp, *binActAll, *binAbpAll, *pArr, *pArr2;
  double len, lenCtr, lenEE, lenSeg, ang, angSum[NDIM], crsAng[2];
  double tensFactSum, tensFabpSum, bendFactSum, bendFabpSum;
  double dr[NDIM], dr2[NDIM], drCrs[NDIM], drSum[NDIM];
  double viscFsum[NDIM], rPos[2][NDIM], ratio[2];
  double *fAct, *fAbp, *fBendAct, *fBendAbp;
  double *lAct, *lAbp, *tensF, *bendF, *viscF;
  FILE *fOut;

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  MALLOC(binAct,int,nActC);
  MALLOC(binAbp,int,nAbpC);
  MPI_Gather(&nActMe, 1, MPI_INT, nActAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nAbpMe, 1, MPI_INT, nAbpAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0) {
	MALLOC(fAct,double,nAct * (NDIM + 1));
	MALLOC(fAbp,double,nAbp * (NDIM + 1));
	MALLOC(fBendAct,double,nAct);
	MALLOC(fBendAbp,double,nAbp);
	MALLOC(rAct,double,nAct * NDIM);
	MALLOC(rAbp,double,nAbp * NDIM);
	MALLOC(lAct,double,nAct);
	MALLOC(lAbp,double,nAbp * recAbp.nL);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(cntAct,int,nAct);
	MALLOC(cntAbp,int,nAbp);
	MALLOC(tensF,double,nAct + nAbp);
	MALLOC(bendF,double,nAct + nAbp);
	MALLOC(viscF,double,(nAct + nAbp) * NDIM);
	MALLOC(fixAct,int,nAct);
	MALLOC(iFila,int,nAct);
	MALLOC(binActAll,int,nAct);
	MALLOC(binAbpAll,int,nAbp);
	if (motSA.gTgl != 0) {
		MALLOC(iMot,int,nAbp);
	}
  }
  FOR_ACTME(n) { 
	binAct[n] = BinaryPackActinChainArray(act.id[n]);
  }
  FOR_ABPME(n) { 
	binAbp[n] = BinaryPackAbpChainArray(abp.id[n]);
  }
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  Gather2dArrayDouble(nActMe, nActAll, NDIM + 1, act.id, recAct.allF, fAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, NDIM + 1, abp.id, recAbp.allF, fAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.bendF, fBendAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, 1, abp.id, recAbp.bendF, fBendAbp);
  Gather2dArrayDouble(nActMe, nActAll, 1, act.id, recAct.len, lAct);
  Gather2dArrayDouble(nAbpMe, nAbpAll, recAbp.nL, abp.id, recAbp.len, lAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, recAct.cnt, cntAct);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, recAbp.cnt, cntAbp);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.iF, iFila);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, act.fix, fixAct);
  Gather2dArrayInt(nActMe, nActAll, 1, act.id, binAct, binActAll);
  Gather2dArrayInt(nAbpMe, nAbpAll, 1, abp.id, binAbp, binAbpAll);

  if (rank == 0) {
	if (motSA.gTgl != 0) {
		memset(iMot, -1, sizeof(int) * nAbp);
		cnt = 0;
		FOR_ABP(n) {
			pArr = &P2A(chAbp,n,0,nChAb);
			CONT(!(pArr[3] > -1 && pArr[4] < 0));
			curr = n;
			while(P2A(chAbp,curr,3,nChAb) > -1) {
				iMot[curr] = cnt;
				curr = P2A(chAbp,curr,3,nChAb);
			}
			iMot[curr] = cnt;
			cnt++;
		}
	}
	fOut = fopen(GenFileName("InfoIndv"), "a");
	if (currTimeStep / period == 1) {
		fprintf(fOut, "%d\t%d\t%d\t", nAct + nAbp, nAct, nAbp);
		Fprintf1dFillerInt(fOut, 0, 16, 0);
	}
	FOR_ACT(n) {
		pArr = &P2A(chAct,n,0,nChAc);
		if ((pArr[0] < 0 && pArr[1] < 0) || cntAct[n] == 0) {
			fprintf(fOut, "%d\t", n);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 7, 0);
			tensF[n] =	0.;
			bendF[n] =	0.;
			V3SET_ALL(&P2(viscF,n,0), 0.);
			continue;
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		for(k = 0; k < 2; k++) {
			CONT(pArr[k] < 0);
			CalcUnitVec(dr, &P2(rAct,n,0), &P2(rAct,pArr[k],0));
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);
		tensF[n] = F_S2PN(actF.spr.stf * (lAct[n] / (double)cntAct[n] 
				- actF.spr.eq));
		bendF[n] = F_S2PN(fBendAct[n] / (double)cntAct[n]);
		VS3COPY(&P2(viscF,n,0), &P2A(fAct,n,0,NDIM + 1), INV(cntAct[n]));
		FOR_NDIM(k) { P2(viscF,n,k) = F_S2PN(P2(viscF,n,k)); }
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n, P2(rAct,n,0), 
				P2(rAct,n,1), P2(rAct,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[n], bendF[n], P2(viscF,n,0), P2(viscF,n,1), P2(viscF,n,2),
				pArr[0], pArr[1], iPerc[n], iSupp[n], fixAct[n], iFila[n], 
				binActAll[n]);
	}
	FOR_ABP(n) {
		pArr = &P2A(chAbp,n,0,nChAb);
		if ((ISMTF(pArr[2]) && pArr[0] < 0 && pArr[1] < 0 && pArr[3] < 0 
				&& pArr[4] < 0) || (!(ISMTF(pArr[2])) && pArr[0] < 0 
				&& pArr[1] < 0) || cntAbp[n] == 0) {
			fprintf(fOut, "%d\t", n + nAct);
			Fprintf1dFillerInt(fOut, 0, 11, 1);
			Fprintf1dFillerInt(fOut, -1, 7, 0);
			tensF[n + nAct] = 0.;
			bendF[n + nAct] = 0.;
			V3SET_ALL(&P2(viscF,n + nAct,0), 0.);
			continue;
		}
		cnt = 0;
		V3SET_ALL(drSum, 0.);
		ind = (ISMTF(pArr[2])) ? 3 : 0;
		for(k = 0; k < 2; k++) {
			CONT(pArr[k + ind] < 0);
			if (ind == 0) {
				CalcUnitVecActinAbp(dr, pArr[k], n, 1);
			}
			else {
				CalcUnitVec(dr, &P2(rAbp,pArr[k + ind],0), &P2(rAbp,n,0)); 
			}
			if (k == 1) { V3REVSIGN(dr); }
			VV3ADD(drSum, dr);
			cnt++;
		}
		NormVec(drSum);
		ind = n + nAct;
		tensF[ind] = F_S2PN(abpF.spr[pArr[2]].stf * abpF.spr[pArr[2]].eq 
				* (P2A(lAbp,n,0,recAbp.nL) / (double)cntAbp[n] - 1.));
		bendF[ind] = F_S2PN(fBendAbp[n] / (double)cntAbp[n]);
		VS3COPY(&P2(viscF,ind,0), &P2A(fAbp,n,0,NDIM + 1), INV(cntAbp[n]));
		FOR_NDIM(k) { P2(viscF,ind,k) = F_S2PN(P2(viscF,ind,k)); }
		for(k = 0; k < 2; k++) {
            if (pArr[k] > -1) {
                side = FindAbpActinChain(pArr[k], n, 1);
                ratio[k] = (side > -1) ? (double)((side - 2) / nChAcY)
                        / (double)nChAcX : -1.;
            }
            else {
                ratio[k] = -1.;
            }
		}
		fprintf(fOut, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t"
				"%g\t%d\t%d\t%d\t%d\t%g\t%g\t%d\n", ind, P2(rAbp,n,0), 
				P2(rAbp,n,1), P2(rAbp,n,2), drSum[0], drSum[1], drSum[2], 
				tensF[ind], bendF[ind], P2(viscF,ind,0), P2(viscF,ind,1), 
				P2(viscF,ind,2), pArr[0], pArr[1], iPerc[ind], iSupp[ind], 
				ratio[0], ratio[1],	binAbpAll[n]);
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFila"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		cnt++;
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 26);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		tensFactSum = 0.;
		tensFabpSum = 0.;
		bendFactSum = 0.;
		bendFabpSum = 0.;
		lenSeg = 0.;
		lenCtr = 0.;
		iPercSum = 0;
		iSuppSum = 0;
		V3SET_ALL(angSum, 0.);
		V3SET_ALL(cntXlink, 0);
		V3SET_ALL(drSum, 0.);
		V3SET_ALL(viscFsum, 0.);
		V3SET_ALL(sft, 0);
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			lenCtr += len;
			CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
			if (pArr[1] > -1) {
				// Calculate bending angles
				CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
				V3CROSS(drCrs, dr, dr2);
				NormVec(drCrs);
				ang = V3ANG(dr, dr2);
				VVS3ADD(angSum, drCrs, ang);
			}
			VVS3ADD(drSum, dr, INV(len));
			tensFactSum += tensF[curr];
			bendFactSum += bendF[curr];
			for(k = 2; k < nChAc; k++) {
				abpInd = pArr[k];
				CONT(abpInd < 0);
				pArr2 = &P2A(chAbp,abpInd,0,nChAb);
				CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
				tensFabpSum += tensF[nAct + abpInd];
				bendFabpSum += bendF[nAct + abpInd];
				iSuppSum += (iSupp[nAct + abpInd] > -1) ? 1 : 0;
				cntXlink[(int)(pArr2[2] / 2)]++; 
				CONT(ISMTF(pArr2[2]));
				actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
				CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
						&P2(rAct,actInd,0));
				if (V3DOT(dr, dr2) > 0.) { cntXlink[2]++; }
			}
			VV3ADD(viscFsum, &P2(viscF,curr,0));
			iPercSum += iPerc[curr] % 10;
			iSuppSum += (iSupp[curr] > -1) ? 1 : 0;
			lenSeg += 1.;
			curr = pArr[0];
		}
		// Normalize the orientation vector.
		NormVec(drSum);
		// Measure the end-to-end distance of the filament.
		V3SUB(dr, &P2(rAct,curr,0), &P2(rAct,n,0));
		VSV3ADD(dr, dr, sft, dimDom);
		lenEE = V3LEN(dr);
		if (lenSeg > 0.) {
			tensFactSum /= lenSeg;
			bendFactSum /= lenSeg;
			V3SCALE(viscFsum, INV(lenSeg));
		}
		if (cntXlink[0] + cntXlink[1] > 0) {
			tensFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
			bendFabpSum /= (double)(cntXlink[0] + cntXlink[1]);
		}
		fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
				"\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\t%d\n", cnt, 
				iFila[n], n, curr, drSum[0], drSum[1], drSum[2], tensFactSum, 
				tensFabpSum, bendFactSum, bendFabpSum, L_S2UM(lenSeg), 
				L_S2UM(lenCtr), L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
				viscFsum[0], viscFsum[1], viscFsum[2], cntXlink[0], cntXlink[1],
				cntXlink[2], iPercSum, iPerc[n] / 10, iSuppSum);
		cnt++;
	}
	fclose(fOut);

	// Record information in unit of filaments
	fOut = fopen(GenFileName("InfoFilaSeg"), "a");
	// Count and print the number of actin filaments
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = P2A(chAct,curr,2 + k * nChAcY + l,nChAc);
					CONT(abpInd < 0);
					CONT(!(ISMTF(P2A(chAbp,abpInd,2,nChAb)))
							&& P2A(chAbp,abpInd,1,nChAb) < 0);
					if (CS == 1) { 
						cnt++;
					}
					CS = 1;
					prevAbpInd = abpInd;
				}
			}
			curr = P2A(chAct,curr,0,nChAc);
		}
	}
    Fprintf1dArrayIntWFil(fOut, &cnt, 1, 0, 22);
	// Record the information
	cnt = 0;
	FOR_ACT(n) {
		CONT(!(P2A(chAct,n,0,nChAc) > -1 && P2A(chAct,n,1,nChAc) < 0));
		CS = 0;
		curr = n;
		while(P2A(chAct,curr,0,nChAc) > -1) {
			pArr = &P2A(chAct,curr,0,nChAc);
			len = CalcVecDist(dr, &P2(rAct,pArr[0],0), &P2(rAct,curr,0), 0);
			if (CS == 1) {
				VVS3ADD(drSum, dr, INV(len));
				tensFactSum += tensF[curr];
				bendFactSum += bendF[curr];
				VV3ADD(viscFsum, &P2(viscF,curr,0));
				cntFor++;
				iPercSum += iPerc[curr] % 10;
				iSuppSum += (iSupp[curr] > -1) ? 1 : 0;
			}
			for(k = 0; k < nChAcX; k++) {
				for(l = 0; l < nChAcY; l++) {
					abpInd = pArr[2 + k * nChAcY + l];
					CONT(abpInd < 0);
					pArr2 = &P2A(chAbp,abpInd,0,nChAb);
					CONT(!(ISMTF(pArr2[2])) && pArr2[1] < 0);
					CalcPosOnActSeg(&P2(rAct,curr,0), &P2(rAct,pArr[0],0), 
							rPos[1], (double)k / (double)nChAcX);
					crsAng[1] = -1.;
					if (!(ISMTF(pArr2[2]))) {
						actInd = pArr2[(pArr2[0] == curr) ? 1 : 0];
						CalcVec(dr2, &P2(rAct,P2A(chAct,actInd,0,nChAc),0), 
								&P2(rAct,actInd,0));
						crsAng[1] = V3ANG(dr, dr2);
					}
					if (CS == 1) { 
						if (curr != prevActInd) {
							CheckCrossBound(sft, &P2(rAct,curr,0), rPos[1]);
						}
						V3SUB(dr2, rPos[1], rPos[0]);
						VSV3ADD(dr2, dr2, sft, dimDom);
						lenEE = V3LEN(dr2);
						if (cntFor > 0) {
							tensFactSum /= (double)cntFor;
							bendFactSum /= (double)cntFor;
							V3SCALE(viscFsum, INV((double)cntFor));
						}
						NormVec(drSum);
						fprintf(fOut, "%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t"
								"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\t"
								"%g\n", cnt, iFila[n], prevActInd, curr, 
								drSum[0], drSum[1], drSum[2], tensFactSum, 
								bendFactSum, L_S2UM(lenSeg), L_S2UM(lenCtr), 
								L_S2UM(lenEE), angSum[0], angSum[1], angSum[2],
								viscFsum[0], viscFsum[1], viscFsum[2],
								prevAbpInd * 10 + P2A(chAbp,prevAbpInd,2,nChAb),
								abpInd * 10 + pArr2[2], crsAng[0], crsAng[1]);
						cnt++;
					}
					lenCtr = 0.;
					lenSeg = 0.;
					tensFactSum = tensF[curr];
					bendFactSum = bendF[curr];
					V3COPY(viscFsum, &P2(viscF,curr,0));
					cntFor = 1;
					V3SET_ALL(angSum, 0.);
					VS3COPY(drSum, dr, INV(len));
					iPercSum = iPerc[curr] % 10;
					iSuppSum = (iSupp[curr] > -1) ? 1 : 0;
					V3SET_ALL(sft, 0);
					prevAbpInd = abpInd;
					prevActInd = curr;
					V3COPY(rPos[0], rPos[1]);
					crsAng[0] = crsAng[1];
					CS = 1;
				}
				if (CS == 1) {
					lenCtr += len / (double)nChAcX;
					lenSeg += 1. / (double)nChAcX;
				}
			}
			if (CS == 1) { 
				CheckCrossBound(sft, &P2(rAct,curr,0), &P2(rAct,pArr[0],0));
				if (pArr[1] > -1) {
					// Calculate bending angles
					CalcVec(dr2, &P2(rAct,curr,0), &P2(rAct,pArr[1],0));
					V3CROSS(drCrs, dr, dr2);
					NormVec(drCrs);
					ang = V3ANG(dr, dr2);
					VVS3ADD(angSum, drCrs, ang);
				}
			}
			curr = pArr[0];
		}
	}
	fclose(fOut);

  }

  if (rank == 0) {
	free(rAct);
	free(rAbp);
	free(fAct);
	free(fAbp);
	free(fBendAct);
	free(fBendAbp);
	free(lAct);
	free(lAbp);
	free(chAct);
	free(chAbp);
	free(cntAct);
	free(cntAbp);
	free(tensF);
	free(bendF);
	free(viscF);
	free(fixAct);
	free(iFila);
	free(binActAll);
	free(binAbpAll);
	if (motSA.gTgl != 0) {
		free(iMot);
	}
  }
  free(nActAll);
  free(nAbpAll);
  free(binAct);
  free(binAbp);
}

/*------------------ Information for filaments and segments  -----------------*/

void RecordConfigParaview(void) {
  int m, n, k, l, cnt, ind, nD, seq, ind2[NDIM], loc, gTglCenAdj;
  int cntKind[3];
  int *nActAll, *nAbpAll, *nMbAll, *pI, *iAct, *iAbp, *iMb, *iMbSld, *sldAbp;
  double *recActG, *recAbpG, *recMbG;
  double sft[NDIM], rBnd[8][NDIM], avgPos[NDIM];
  FILE *fOut;

  gTglCenAdj = 0;
  nD = 6;
  seq = (int)(currTimeStep / recConfPV.prd);

  MALLOC(nActAll,int,nCpu);
  MALLOC(nAbpAll,int,nCpu);
  if (gTglMb != 0) { 
	  MALLOC(nMbAll,int,nCpu);
  }
  if (rank == 0) {
	MALLOC(rAct,double,nAct*NDIM);
	MALLOC(rAbp,double,nAbp*NDIM);
	MALLOC(chAct,int,nAct*nChAc);
	MALLOC(chAbp,int,nAbp*nChAb);
	MALLOC(iAct,int,nAct);
	MALLOC(iAbp,int,nAbp);
	MALLOC(recActG,double,nAct * nD);
	MALLOC(recAbpG,double,nAbp * nD);
	if (gTglMb != 0) {
		MALLOC(rMb,double,nMb*NDIM);
		MALLOC(chMb,int,nMb*nChMb);
		MALLOC(iMb,int,nMb*nChMb);
		MALLOC(recMbG,double,nMb * nD);	
    	if (mbSld.abp.gTgl != 0) {
			MALLOC(sldAbp,int,nAbp * dimMbNuc);
		}
		if (mbSld.act.gTgl != 0 || mbSld.abp.gTgl != 0) {
			MALLOC(iMbSld,int,nMb * nChMb);
			memset(iMbSld, -1, sizeof(int) * nMb * nChMb);
		}
	}
  }
  // Gather the positions and chain information of particles.
  GatherActChainPosition(nActAll, chAct, rAct);
  GatherAbpChainPosition(nAbpAll, chAbp, rAbp);
  if (gTglMb != 0) {
	GatherMembChainPosition(nMbAll, chMb, rMb);
    if (mbSld.abp.gTgl != 0) {
		Gather2dArrayInt(nAbpMe, nAbpAll, dimMbNuc, abp.id, mbSld.abp.l, 
				sldAbp);  
	}
  }
  RecordConfigInformation(recActG, recAbpG, recMbG, 
		nActAll, nAbpAll, nMbAll, nD);

  if (rank == 0) {
	V3SET_ALL(cntKind, 0);
	FOR_ABP(n) {
		cntKind[P2A(chAbp,n,2,nChAb)]++;
	}
	fOut = fopen(GenFileName("ConfPV"), "a");
	if (seq == 1) {
		fprintf(fOut, "%d %d %d %d %d %d\n", (nAct > 0 ? 1 : 0), 
				(cntKind[0] > 0 ? 1 : 0), (cntKind[1] > 0 ? 1 : 0), 
				(cntKind[2] > 0 ? 1 : 0), (gTglMb > 0 ? 1 : 0),
				(gTglCenAdj != 1 ? 1 : 0));
	}

	if (gTglCenAdj == 1 && gTglMb != 0) {
		V3SET_ALL(avgPos, 0.);
		cnt = 0;
		FOR_MB(n) {
			VV3ADD(avgPos, &P2(rMb,n,0));
			cnt++;
		}
		V3SCALE(avgPos, 1. / (double)cnt);
		FOR_NDIM(k) {
			sft[k] = rGrid[k][0] + dimDomH[k] - avgPos[k];
		}

		FOR_ACT(n) {
			VV3ADD(&P2(rAct,n,0), sft);
		}
		FOR_ABP(n) {
			VV3ADD(&P2(rAbp,n,0), sft);
		}
		FOR_MB(n) {
			VV3ADD(&P2(rMb,n,0), sft);
		}
	}

	fprintf(fOut, "------------------------------------ %6d ------------"
			"------------------------\n", seq);
	if (nAct > 0) {
		// Actin
	  	cnt = 0;
		memset(iAct, -1, sizeof(int) * nAct);
		fprintf(fOut, "# vtk DataFile Version 2.0\n");
		fprintf(fOut, "ACTIN\n");
		fprintf(fOut, "ASCII\n");
		fprintf(fOut, "DATASET UNSTRUCTURED_GRID\n");
		FOR_ACT(n) {
			ind = P2A(chAct,n,0,nChAc);
			CONT(ind < 0);
			FOR_NDIM(k) {
				CONT(pbc[k] != 1);
				BREAK(fabs(P2(rAct,ind,k) - P2(rAct,n,k)) 
						> 0.5 * dimDom[k]);
			}
			CONT(k != NDIM);
			iAct[n] = 1;
			cnt++;
		}
		fprintf(fOut, "POINTS %d double\n", cnt * 2);
		FOR_ACT(n) {
			CONT(iAct[n] != 1);
			ind = P2A(chAct,n,0,nChAc);
			fprintf(fOut, "%lf %lf %lf\n", P2(rAct,n,0), P2(rAct,n,1), 
					P2(rAct,n,2)); 
			fprintf(fOut, "%lf %lf %lf\n", P2(rAct,ind,0), 
					P2(rAct,ind,1),	P2(rAct,ind,2));
		}
	
		fprintf(fOut, "CELLS %d %d\n", cnt, 3 * cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "2 %d %d\n", 2 * n, 2 * n + 1);
		}
		fprintf(fOut, "CELL_TYPES %d\n", cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "3\n");
		}
		fprintf(fOut, "CELL_DATA %d\n", cnt);
	 	fprintf(fOut, "SCALARS force double 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ACT(n) {
			CONT(iAct[n] != 1);	
			ind = P2A(chAct,n,0,nChAc);
			fprintf(fOut, "%f\n", 0.5 * (P2A(recActG,n,2,nD) 
				+ P2A(recActG,ind,2,nD)));
		}

		fprintf(fOut, "SCALARS filaEnd int 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ACT(n) {
			CONT(iAct[n] != 1);	
			pI = &P2A(chAct,n,0,nChAc);
			if (pI[1] < 0) {
				fprintf(fOut, "2\n");
			} 
			else {
				if (P2A(chAct,pI[0],0,nChAc) < 0) {
					fprintf(fOut, "0\n");
				}
				else {
					fprintf(fOut, "1\n");
				}
			}
		}
		fprintf(fOut, "#\n");
	}

	for(m = 0; m < (motSA.gTgl != 0 ? 2 : 3); m++) {
		CONT(cntKind[m] == 0);
		// ACP
		cnt = 0;
		memset(iAbp, -1, sizeof(int) * nAbp);
		fprintf(fOut, "# vtk DataFile Version 2.0\n");
		if (m == 0) {
			fprintf(fOut, "ACPC\n");
		}
		else if (m == 1) {
			fprintf(fOut, "ACPB\n");
		}
		else {
			fprintf(fOut, "MOTOR\n");
		}
		fprintf(fOut, "ASCII\n");
		fprintf(fOut, "DATASET UNSTRUCTURED_GRID\n");
		FOR_ABP(n) {
			pI = &P2A(chAbp,n,0,nChAb);
			CONT(pI[2] != m || (pI[0] < 0 && pI[1] < 0));
			iAbp[n] = 1;
			cnt++;
		}
		fprintf(fOut, "POINTS %d double\n", cnt);
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);
			fprintf(fOut, "%lf %lf %lf\n", P2(rAbp,n,0), P2(rAbp,n,1), 
					P2(rAbp,n,2)); 
		}

		fprintf(fOut, "CELLS %d %d\n", cnt, 2 * cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "1 %d\n", n);
		}
		fprintf(fOut, "CELL_TYPES %d\n", cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "1\n");
		}
		fprintf(fOut, "CELL_DATA %d\n", cnt);
	 	fprintf(fOut, "SCALARS force double 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);	
			fprintf(fOut, "%f\n", P2A(recAbpG,n,2,nD));
		}

		fprintf(fOut, "SCALARS state int 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);	
			pI = &P2A(chAbp,n,0,nChAb);
			if (pI[0] < 0 && pI[1] < 0) {
				fprintf(fOut, "0\n");
			} 
			else if (pI[0] > -1 && pI[1] < 0) {
				fprintf(fOut, "1\n");
			}
			else {
				fprintf(fOut, "2\n");
			}
		}
		if (gTglMb != 0 && mbSld.abp.gTgl != 0 && m < 2) {
			fprintf(fOut, "SCALARS membBound int 1\n");
			fprintf(fOut, "LOOKUP_TABLE default\n");
			FOR_ABP(n) {
				CONT(iAbp[n] != 1);	
				if (P2A(sldAbp,n,0,dimMbNuc) == -1) {
					fprintf(fOut, "0\n");
				}
				else {
					fprintf(fOut, "1\n");
					for(k = 0; k < dimMbNuc; k++) {
						ind = P2A(sldAbp,n,k,dimMbNuc);
						for(l = 1; l < dimMbNuc; l++) {
		        			loc = FindElementArray(&P2A(chMb,ind,0,nChMb), 
									((dimMbNuc == 2) ? 2 : 6), 
									P2A(sldAbp,n,(k + l) % dimMbNuc,dimMbNuc),
									0, 1);
							P2A(iMbSld,ind,loc,nChMb) = 1;
						}
					}
				}
			}
		}
		fprintf(fOut, "#\n");
	}

	// Motor
	if (motSA.gTgl != 0 && cntKind[2] > 0) {
	  	cnt = 0;
		memset(iAbp, -1, sizeof(int) * nAbp);
		fprintf(fOut, "# vtk DataFile Version 2.0\n");
		fprintf(fOut, "MOTOR\n");
		fprintf(fOut, "ASCII\n");
		fprintf(fOut, "DATASET UNSTRUCTURED_GRID\n");
		FOR_ABP(n) {
			pI = &P2A(chAbp,n,0,nChAb);
			ind = pI[3];
			CONT(pI[2] != 2 || ind < 0);
			FOR_NDIM(k) {
				CONT(pbc[k] != 1);
				BREAK(fabs(P2(rAbp,ind,k) - P2(rAbp,n,k)) 
					> 0.5 * dimDom[k]);
			}
			CONT(k != NDIM);
			iAbp[n] = 1;
			cnt++;
		}
		fprintf(fOut, "POINTS %d double\n", cnt * 2);
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);
			ind = P2A(chAbp,n,3,nChAb);
			fprintf(fOut, "%lf %lf %lf\n", P2(rAbp,n,0), P2(rAbp,n,1), 
					P2(rAbp,n,2)); 
			fprintf(fOut, "%lf %lf %lf\n", P2(rAbp,ind,0), 
					P2(rAbp,ind,1),	P2(rAbp,ind,2));
		}
		fprintf(fOut, "CELLS %d %d\n", cnt, 3 * cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "2 %d %d\n", 2 * n, 2 * n + 1);
		}
		fprintf(fOut, "CELL_TYPES %d\n", cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "3\n");
		}
		fprintf(fOut, "CELL_DATA %d\n", cnt);
	 	fprintf(fOut, "SCALARS force double 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);	
			fprintf(fOut, "%f\n", P2A(recAbpG,n,2,nD));
		}
		fprintf(fOut, "SCALARS state int 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_ABP(n) {
			CONT(iAbp[n] != 1);	
			pI = &P2A(chAbp,n,0,nChAb);
			if (pI[0] < 0 && pI[1] < 0) {
				fprintf(fOut, "0\n");
			} 
			else if (pI[0] > -1 && pI[1] < 0) {
				fprintf(fOut, "1\n");
			}
			else {
				fprintf(fOut, "2\n");
			}
		}
		fprintf(fOut, "#\n");
	}
	// Membrane
	if (gTglMb != 0) {
	  	cnt = 0;
		memset(iMb, -1, sizeof(int) * nMb * nChMb);
		fprintf(fOut, "# vtk DataFile Version 2.0\n");
		fprintf(fOut, "MEMBRANE\n");
		fprintf(fOut, "ASCII\n");
		fprintf(fOut, "DATASET UNSTRUCTURED_GRID\n");
		FOR_MB(n) {
			pI = &P2A(chMb,n,0,nChMb);
			for(k = 0; k < nChMb; k++) {
				ind = pI[k];
				CONT(ind < 0);
				CONT(n > ind);
				FOR_NDIM(l) {
					CONT(pbc[l] != 1);
					BREAK(fabs(P2(rMb,ind,l) - P2(rMb,n,l)) 
							> 0.5 * dimDom[l]);
				}
				CONT(l != NDIM);
				P2A(iMb,n,k,nChMb) = 1;
				cnt++;
			}
		}
		fprintf(fOut, "POINTS %d double\n", cnt * 2);
		FOR_MB(n) {
			for(k = 0; k < nChMb; k++) {
				CONT(P2A(iMb,n,k,nChMb) != 1);
				ind = P2A(chMb,n,k,nChMb);
				fprintf(fOut, "%lf %lf %lf\n", P2(rMb,n,0), 
						P2(rMb,n,1), P2(rMb,n,2)); 
				fprintf(fOut, "%lf %lf %lf\n", P2(rMb,ind,0), 
						P2(rMb,ind,1),	P2(rMb,ind,2));
			}
		}

		fprintf(fOut, "CELLS %d %d\n", cnt, 3 * cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "2 %d %d\n", 2 * n, 2 * n + 1);
		}
		fprintf(fOut, "CELL_TYPES %d\n", cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "3\n");
		}
		fprintf(fOut, "CELL_DATA %d\n", cnt);
	 	fprintf(fOut, "SCALARS force double 1\n");
		fprintf(fOut, "LOOKUP_TABLE default\n");
		FOR_MB(n) {
			pI = &P2A(chMb,n,0,nChMb);
			for(k = 0; k < nChMb; k++) {
				ind = pI[k];
				CONT(P2A(iMb,n,k,nChMb) != 1);	
				fprintf(fOut, "%f\n", 0.5 * (P2A(recMbG,n,2,nD) 
					+ P2A(recMbG,ind,2,nD)));
			}
		}
		if (mbSld.abp.gTgl != 0) {
			fprintf(fOut, "SCALARS acpBound int 1\n");
			fprintf(fOut, "LOOKUP_TABLE default\n");
			FOR_MB(n) {
				for(k = 0; k < nChMb; k++) {
					CONT(P2A(iMb,n,k,nChMb) != 1);	
					if (P2A(iMbSld,n,k,nChMb) == 1) {
						fprintf(fOut, "1\n");
					}
					else {
						fprintf(fOut, "0\n");
					}
				}
			}
		}
		fprintf(fOut, "#\n");
	}

	if (gTglCenAdj != 1) {
		// Boundaries
	  	cnt = 12;
		fprintf(fOut, "# vtk DataFile Version 2.0\n");
		fprintf(fOut, "BOUNDARY\n");
		fprintf(fOut, "ASCII\n");
		fprintf(fOut, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(fOut, "POINTS %d double\n", cnt * 2);
	
		for(n = 0; n < 8; n++) {
			V3IND_ASSIGN_CONST_INT(n, 2, ind2);
			V3SET_ALL(sft, 0.);
			// If it is a bulk rheology with shear, the location of boundary
			// is adjusted..
	
			if (rheoWay > 0 && bulkRheoType == 0) {
				sft[dirStr] = (ind2[dirNoPBC] == 1) ? stra.acc 
						* dimDom[dirNoPBC] : 0.;
			}
			FOR_NDIM(k) {
				rBnd[n][k] = rGrid[k][(ind2[k] == 0) ? 0 : nGrid[k] - 1] 
						+ sft[k];
			}
		}
		// In x direction
		for(n = 0; n < 4; n++) {
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n][0], rBnd[n][1], rBnd[n][2]);
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n + 4][0], rBnd[n + 4][1], 
					rBnd[n + 4][2]);
		}
		// In z direction
		for(n = 0; n < 4; n++) {
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n * 2][0], rBnd[n * 2][1], 
					rBnd[n * 2][2]);
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n * 2 + 1][0], 
					rBnd[n * 2 + 1][1], rBnd[n * 2 + 1][2]);
		}
		// In y direction
		for(n = 0; n < 2; n++) {
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n][0], rBnd[n][1], rBnd[n][2]);
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n + 2][0], rBnd[n + 2][1], 
					rBnd[n + 2][2]);
		}
		for(n = 0; n < 2; n++) {
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n + 4][0], rBnd[n + 4][1], 
					rBnd[n + 4][2]);
			fprintf(fOut, "%lf %lf %lf\n", rBnd[n + 6][0], rBnd[n + 6][1], 
					rBnd[n + 6][2]);
		}
		fprintf(fOut, "CELLS %d %d\n", cnt, 3 * cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "2 %d %d\n", 2 * n, 2 * n + 1);
		}
		fprintf(fOut, "CELL_TYPES %d\n", cnt);
		for(n = 0; n < cnt; n++){
			fprintf(fOut, "3\n");
		}
		fprintf(fOut, "#\n");
	}

	fclose(fOut);

	free(rAct);
	free(rAbp);
	free(chAct);
	free(chAbp);
	free(iAct);
	free(iAbp);
	free(recActG);
	free(recAbpG);
	if (gTglMb != 0) {
		free(rMb);
		free(chMb);
		free(iMb);
		free(recMbG);
    	if (mbSld.abp.gTgl != 0) {
			free(sldAbp);
		}
		if (mbSld.act.gTgl != 0 || mbSld.abp.gTgl != 0) {
			free(iMbSld);
		}
	}
  }
  free(nActAll);
  free(nAbpAll);
  if (gTglMb != 0) {
	free(nMbAll);
  }
}
