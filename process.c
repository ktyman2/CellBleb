// This file contains detailed processes.

void InitSingleProcessSubroutine(double *arr, int cnt) {
  double *ip;

  ip = arr;
  while (ip < arr + cnt) {
    *ip = 0.;
    ip++;
  }
}

// Initiate each single step 
void InitSingleProcess(void) {
  int n;

  currTimeStep++;
  // Forces are set to zero before SingleProcess().
  // forces of actins.
  InitSingleProcessSubroutine(act.f, (nActMe + nActCp) * NDIM);
  // Brownian forces of actins.
  InitSingleProcessSubroutine(act.fBr, (nActMe + nActCp) * NDIM);
  // forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.f, (nAbpMe + nAbpCp) * NDIM);
  // Brownian forces of ACPs and motors.
  InitSingleProcessSubroutine(abp.fBr, (nAbpMe + nAbpCp) * NDIM);
  // forces of membranes.
  if (gTglMb != 0) {
	InitSingleProcessSubroutine(memb.f, (nMbMe + nMbCp) * NDIM);
  }
  if (gTglBead != 0) {
	InitSingleProcessSubroutine(bead.f, bead.n * NDIM);
  }
  // forces actin on boundaries.
  if (bndMv.gTgl != 0 || rheoWay > 0) {
	for(n = 0; n < NDIM * 2 * NDIM; n++) { bnd.f[n] = 0.; }
  }
  if (actDgd.gTgl != 0) {
	actDgd.c = 0;
  }
  if (tglActFormDyn != 0) {
	cntNucAssAnn = 0;
  }

  if (gTglMb != 0) {
	FOR_ACTME(n) { 
		act.distMb[n] = radMb;
	}
	FOR_ABPME(n) { 
		abp.distMb[n] = radMb;
	}
  }

  if (currTimeStep == netForm.dur && motActiv.dur > 0) {
	Printf0("============================ Motor activation phase "
			"============================\n");
  }
  else if (currTimeStep == netForm.dur + motActiv.dur && rheo.dur > 0) {
	Printf0("=============================== Main measurement "
			"===============================\n");
  }
}

void SingleProcess(void) {
  // Initialize each step
  InitSingleProcess ();
  // Adjust size of subdomains depending on the number of particles in each 
  // Measure displacements of all particles to determine whether
  // neighboring list should be updated or not at a current time step.
  MeasureDisplacementForNL();
  // Calculate and apply stress or strain to a network for bulk 
  // rheology measurement
  if (rheoWay > 0 && currTimeStep >= netForm.dur + motActiv.dur) { 
	ApplyStressStrain(); 
  }
  // Repulsive forces between actin segments
  if (actF.rep.facStf > 0. || abpF.rep.facStf > 0.) {
	CalcRepulsiveForces(); 
  }
  // Spring forces for actin, ACP, and motor
  CalcSpringForces();
  // Bending forces for actin filaments
  if (actF.bend.facStf > 0.) {
	CalcFilaBendForces();
  }
  // Brownian forces for actin, ACP, and motor
  CalcBrownForces();
  // Two kinds of bending forces for ACPs
  if (!(motSA.gTgl != 0 && nAbp - nMot == 0)) {
	CalcAbpBendForces();
  }
  // Spring and bending forces acting on the backbone of multimerized motors
  if (motSA.gTgl != 0 && nMot > 0) { 
	CalcMotorBackboneForces();
  }

  /*-------------------------------- Membrane --------------------------------*/
  // Adjust size of subdomains depending on the number of particles in each 
  // subdomain. The aim is to maintain the number at similar level.
  if (gTglMb != 0) {
	UpdateMembraneUnitArea();
	CalcMembraneRepulsiveForces();
	if (mbSld.gTgl != 0) { 
		if (currTimeStep % recProg.prd == 0 || currTimeStep == 1) {
			UpdateMembraneSlideList();
		}
		CalcMembraneSlideForces();
	}
	CalcMembraneSpringForces();
	if (memb.bend.facStf > 0.) {
		if (dimMbNuc == 2) {
			Calc2dMembraneBendForces();
		}
		else {
			Calc3dMembraneBendForces();
		}
	}
	if ((mbVol.gTgl != 0 || nucVol.gTgl != 0) && currTimeStep > 1) {
		CalcMembraneVolumeForces();
	}

	if ((mbAre.gTgl != 0 || nucAre.gTgl != 0) && currTimeStep > 1 
			&& dimMbNuc == 3) {
		CalcMembraneAreaForces();
	}
	if (mbPro.gTgl != 0 && currTimeStep >= netForm.dur) {
		CalcMembraneProtrusiveForces();
	} 
  }
  /*----------------------------------- Bead ---------------------------------*/
  if (gTglBead != 0) {
	CalcBeadRepulsiveForces();
	if (bead.rheoWay == 1 
			&& currTimeStep >= netForm.dur + motActiv.dur) {
		ApplyBeadSinuSignal();
	}
	if (beadBind.gTgl != 0) {
		CalcBeadActinBindingForces();
  		if (!(actDis.tgl == 0 && actBst.tgl == 0 
				&& currTimeStep >= netForm.dur + motActiv.dur)) {
			UpdateBeadActinBinding();
		}
	}
  }
  /*--------------------------------- Boundary -------------------------------*/
  if (bndReb.gTgl != 0) {
	CalcBoundSpringForces();
  }
  if (bnd.gTglRnd != 0) { 
	CalcBoundRepulsiveForces();
  }

  /*--------------------------- Dynamics of actins ---------------------------*/
  // Check whether severing occurs at a previous time step
  if (nCpu > 1 && actSev.tgl != 0) {
	UpdateActinSeverEvents(); 
  }
  // Capping and uncapping of the ends of actin filaments
  if (actCap.tgl != 0 || actUnc.tgl != 0 
		|| (actSev.tgl != 0 && gTglActSevCap != 0)) {
	UpdateActinCapUncap();
  }
  // Polymerization of actin monomers on filaments
  if (actAss.tgl != 0 && actM.c > 0) {
	UpdateActinAssembly();
  }
  // Depolymerization of actin monomers on filaments
  if (actDis.tgl != 0 || actBst.tgl != 0) {
	UpdateActinDisassembly();
  }
  // Depolymerization of actin monomers on filaments
  if (actBst.tgl != 0) {
	UpdateActinBurstDisassembly();
  }
  // Nucleation of actin
  if (actNuc.tgl != 0 && actM.c > 1) {
	if (!(gTglMb != 0 && ((currTimeStep < netForm.dur && gTglMbActNuc == 0) 
			|| currTimeStep >= netForm.dur))) {
		UpdateActinNucleation();
	}
  }
  // Annealing of actin filaments
  if (actAnn.tgl != 0) { 
	UpdateActinAnnealing();
  }
  // Severing of actin filaments
  if (actSev.tgl != 0) { 
	UpdateActinSevering();
  }
  // Eliminate expired elements in noActDyn.l
  if (tglActDyn != 0) {
	UpdateNoActinDynamicsList();
  }
  // If actin disassembly/severing is allowed with bound ABPs, this function 
  // should be executed here to let other adjacent subdomains know about ABP 
  // unbinding events caused by disassembly/severing.
  if (nCpu > 1 && tglActDisBstSev != 0) {
	UpdateActinAbpDynamicsEvents(0); 
  }
  if (nCpu > 1 && gTglMb != 0 && actDis.gTgl != 0 && actDgd.gTgl != 0) {
	UpdateMembraneDynamicsEvents();   
  }
  /*----------------------- Dynamics of ACPs or motors -----------------------*/
  // Binding event of monomeric ACP or motor
  if ((acpMoBind.tgl != 0 && acpM.c > 0 && gTglImpAcpM != 0) 
		|| (motMoBind.tgl != 0 && motM.c > 0 && gTglImpMotM != 0)) {
	UpdateAbpMonomerBinding();
  }
  if (motSA.gTgl != 0) {
	if (motSA.to.gTgl != 0 && currTimeStep >= netForm.dur) {
		UpdateMotorTurnover();
	}
	if ((motM.c > 0 && gTglImpMotM != 0) || gTglImpMotM == 0) {
		UpdateMotorAssembly();
	}
  }
  // Unbinding event of inactive ACP or motor to monomer
  if (acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) {
	UpdateInactAbpUnbinding();
  }
  // Binding event of inactive ACP or motor
  if (acpReb.tgl != 0 || motReb.tgl != 0) {
	UpdateAbpBinding();
  }
  // Unbinding event of active ACP or motor
  if (acpUnb.tgl != 0 || motUnb.tgl != 0) { 
	UpdateActiveAbpUnbinding();
  }
  // Walking event of active motor
  if (motWalk.tgl != 0) {
	UpdateMotorWalking();
  }
  // Eliminate expired elements in noAbpDyn.l 
  if (tglAbpAcInaDyn != 0 || acpMoBind.tgl != 0 || motMoBind.tgl != 0 
		|| tglActDisBstSevAbp != 0) { 
	UpdateNoAbpUnbindList();
  }
  // transfer actin severing/disassembly and ABP unbinding/binding events 
  // stored in sendAbpDyn.l to adjacent CPUs
  if (nCpu > 1 && tglAbpAcInaDyn != 0) { 
	UpdateActinAbpDynamicsEvents(1);   
  }
  if (updMono.tgl != 0) {
	if (nCpu > 1 && currTimeStep % updMono.prd  == 0  
			&& (tglAbpInaMoDyn != 0 || tglActMoDyn != 0) 
			&& !(gTglMb != 0 && gTglLoadNetData == 0 && (currTimeStep < netForm.dur / 2 
			|| currTimeStep >= netForm.dur))) {
		UpdateActinAbpMonomerList();
	}
  }
  if (motSA.gTgl != 0) {
	if (nCpu > 1 && currTimeStep % updMotN.prd == 0 && updMotN.tgl != 0) {
		UpdateMotorNucleCounter();
	}
  }
  /*-------------------------------- Boundary --------------------------------*/
  // Unbinding and binding of filaments from boundaries
  if (bnd.gotF != 0) {
	if ((gTglLoadNetData != 0 && currTimeStep >= netForm.dur) 
			|| (gTglLoadNetData == 0 && currTimeStep < netForm.dur)) {
		if (bndUnb.gTgl != 0) {
			UpdateBoundaryActinUnbindMature();
		}
		if (bndReb.gTgl != 0) {
			UpdateBoundaryActinBinding();
		}
	}
  }
  /*--------------------------------- Membrane -------------------------------*/
  // Dynamics related to membranes
  if (gTglMb != 0) {
	if (mbUnb.gTgl != 0 || mbReb.gTgl != 0) {
		UpdateMembraneUnbindMature();
	}
	if (mbReb.gTgl != 0) {
		UpdateMembraneBinding();
	}
	if (nCpu > 1 && (mbUnb.gTgl != 0 || mbReb.gTgl != 0)) {
		UpdateMembraneDynamicsEvents();   
	}
	if (mbDef.gTgl != 0) {
		UpdateMembraneRelease();
	}
	if (currTimeStep >= netForm.dur + pres.dur + motActiv.dur) {
		if (mbFix.gTgl != 0) {
			UpdateMembraneFixation();
		}
	}
	if (mbPro.gTgl != 0 && currTimeStep >= netForm.dur) {
		UpdateMembraneProtrusion();
		if (gTglMbCont != 0) {
			UpdateMembraneProteinRecover();
		}
	} 
	if ((recMbDim.tgl != 0 && currTimeStep == 1) 
			|| ((mbVol.gTgl != 0 || nucVol.gTgl != 0 
			|| (mbPro.gTgl != 0 && mbPro.gTglCT != 0)))) {
		UpdateMembraneCenterPosition();
		UpdateMembraneVolume();
	}
	if ((recMbDim.tgl != 0 && currTimeStep == 1) 
			|| ((mbAre.gTgl == 2 || nucAre.gTgl == 2))) {
		UpdateMembraneTotalArea();
	}
	if (tglMbNucDyn != 0) {
		UpdateNoMembraneDynamicsList();
	}
  }
  /*--------------------------------------------------------------------------*/
  // Calculate new positions of all particles based on forces 
  UpdateNewLocation();
  if (gTglBead != 0 && ((currTimeStep >= netForm.dur && bead.rheoWay == 0)
		|| (currTimeStep >= netForm.dur && currTimeStep < netForm.dur 
		+ motActiv.dur && bead.rheoWay == 1))) {
	UpdateBeadNewLocation();
  }
  // Apply the periodic boundary condition to all particles
  ApplyBoundCondAll();
  // Only for a case with more than 1 CPU
  if (nCpu > 1) { 
	// Adjust the locations of boundaries following stress-strain relationship
	if (bndMv.gTgl != 0 && currTimeStep >= netForm.dur 
			+ pres.dur + motActiv.dur) {
		UpdateBoundaryLocation();		
	}
	// Adjust size of subdomains depending on the number of particles in each 
	// subdomain. The aim is to maintain the number at similar level.
	// MoveParticles() and CopyParticles() should be excuted right after 
	// UpdateSubdomSectionLocation().
	if (updSubdSize.tgl != 0) {
		if ((currTimeStep == netForm.dur && gTglLoadNetData == 0) 
				|| (currTimeStep % updSubdSize.prd == 0
				&& !(currTimeStep < netForm.dur && gTglLoadNetData == 0))) {
			UpdateSubdomSectionLocation(0);
			if (rank == 0) { 
				RecordSubdomSectionLocation(); 
			}
			MoveParticles();
			CopyParticles();
			UpdateNeighborList();
		}
	}
	// Process longCh.l in either normal or Plympton way
	if (mpiMethod == 0) { 
		UpdateLongChainNormal(); 
	}
	else { 
		UpdateLongChainPlympton(); 
	}
	// Move particles escaping from a current subdomain to adjacent subdomains
	MoveParticles();
	// Copy particles located near boundaries to adjacent subdomains 
	// for synchronization
	CopyParticles();
  }
}


void InitToggleParameters(int mode) {
  // Actin assembly or disassembly is on or off 
  // depending on parameters in 'condition'.
  if (nAct > 0 && (mode == 2 || (mode == 1 && gTglActDynPres != 0)
		|| (mode == 0 && gTglActDynNF != 0))) {
	actAss.tgl = actAss.gTgl;
	actDis.tgl = actDis.gTgl;
	actBst.tgl = actBst.gTgl;
	actSev.tgl = actSev.gTgl;
	actAnn.tgl = actAnn.gTgl;
	actNuc.tgl = actNuc.gTgl;
	actCap.tgl = actCap.gTgl;
	actUnc.tgl = actUnc.gTgl;
  }
  else {
	actAss.tgl = 0;
	actDis.tgl = 0;
	actBst.tgl = 0;
	actSev.tgl = 0;
	actAnn.tgl = 0;
	actNuc.tgl = 0;
	actCap.tgl = 0;
	actUnc.tgl = 0;
  }
  // Binding of monomeric ACPs or motors is on or off.
  if ((mode == 2 && nAbp > 0) || (mode == 1 && ((gTglAcpDynPres != 0 
		&& nAbp - nMot > 0) || (gTglMotUnbRebPres != 0 && nMot > 0)))
		|| (mode == 0 && ((gTglAcpDynNF != 0 && nAbp - nMot > 0) 
		|| (gTglMotUnbRebNF != 0 && nMot > 0)))) { 
	acpMoBind.tgl = acpMoBind.gTgl;
	motMoBind.tgl = motMoBind.gTgl;
  }
  else { 
	acpMoBind.tgl = 0; 
	motMoBind.tgl = 0; 
  }
  // Unbinding of active ACPs or motors and 
  // unbinding and binding of inactive ACPs of motors are on or off.
  if (nAbp - nMot > 0 && (mode == 2 || (mode == 1 && gTglAcpDynPres != 0)
		|| (mode == 0 && gTglAcpDynNF != 0))) {
	acpUnb.tgl = acpUnb.gTgl;
	acpReb.tgl = acpReb.gTgl;
	acpInaUnb.tgl = acpInaUnb.gTgl;
  }
  else {
	acpUnb.tgl = 0;
	acpReb.tgl = 0;
	acpInaUnb.tgl = 0;
  }
  // Unbinding of active motors and the unbinding or binding of inactive
  // motors are on or off. 
  if (nMot > 0 && (mode == 2 || (mode == 1 && gTglMotUnbRebPres != 0)
		|| (mode == 0 && gTglMotUnbRebNF != 0))) {
	motUnb.tgl = motUnb.gTgl;
	motReb.tgl = motReb.gTgl;
	motInaUnb.tgl = motInaUnb.gTgl;
  }
  else {
	motUnb.tgl = 0;
	motReb.tgl = 0;
	motInaUnb.tgl = 0;
  }
  // Walking of active motors is on or off.
  if (nMot > 0 && (mode == 2 || (mode == 1 && gTglMotWalkPres != 0)
		|| (mode == 0 && gTglMotWalkNF != 0))) {
	motWalk.tgl = motWalk.gTgl;
  }
  else { 
	motWalk.tgl = 0; 
  }
  // Toggle parameters representing multiple toggles are on or off.
  tglAcpAcInaDyn = (acpUnb.tgl != 0 || acpReb.tgl != 0 
		|| acpInaUnb.tgl != 0) ? 1 : 0;
  tglMotAcInaDyn = (motUnb.tgl != 0 || motReb.tgl != 0 || motWalk.tgl != 0 
		|| motInaUnb.tgl != 0) ? 1 : 0;
  tglAbpAcInaDyn = (tglAcpAcInaDyn != 0 || tglMotAcInaDyn != 0) ? 1 : 0;
  tglAbpInaMoDyn = (acpMoBind.tgl != 0 || motMoBind.tgl != 0
		|| acpInaUnb.tgl != 0 || motInaUnb.tgl != 0) ? 1 : 0;
  tglActDisBstSev = (actBst.tgl != 0 || actDis.tgl != 0 
		|| actSev.tgl != 0) ? 1 : 0;
  tglActDisBstSevAbp = ((actDis.tgl != 0 && actDis.facKWA > 0.) 
		|| (actBst.tgl != 0 && actBst.facKWA > 0.) 
		|| (actSev.tgl != 0 && actSev.facKWA > 0.)) ? 1 : 0;
  tglActFormDyn = (actAss.tgl != 0 || actNuc.tgl != 0 || actAnn.tgl != 0) 
		? 1 : 0;
  tglActMoDyn = (actAss.tgl != 0 || actDis.tgl != 0 || actBst.tgl != 0 
		|| actNuc.tgl != 0) ? 1 : 0;
  tglActDyn = (tglActMoDyn != 0 || actSev.tgl != 0 || actAnn.tgl != 0) ? 1 : 0;
  tglMbDyn = (gTglMb != 0 && (mbFix.gTgl != 0 || mbDef.gTgl != 0 
		|| mbReb.gTgl != 0 || mbUnb.gTgl != 0 )) ? 1 : 0;
  tglMbNucDyn = (gTglMb != 0 && (tglMbDyn != 0 
		|| nucReb.gTgl != 0 || nucUnb.gTgl != 0)) ? 1 : 0;
  tglMbNucLocAre = (gTglMb != 0 && dimMbNuc == 3 
		&& (mbAre.gTgl == 1 || nucAre.gTgl == 1)) ? 1 : 0;
  tglNeiAbpSC = (gTglMb != 0 || acpReb.gTgl != 0 || motReb.gTgl != 0 
		|| abpF.rep.facStf > 0.) ? 1 : 0;
  tglNeiAbpDC = (gTglMb != 0 || abpF.rep.facStf > 0.) ? 1 : 0;
}

// Main process
// There are two separate parts for each rheology measurement.
// rheoWay = 0 -> segment-rheology, 1 -> bulk rheology
void MainProcess(void) {
  // Synchronization
  MPI_Barrier(MPI_COMM_WORLD);
  // Initial output.
  if (rank == 0) { 
	InitOutput(); 
  }
  InitToggleParameters(0);
  // This loop is iterated until currTimeStep reaches a last time step.
  // "netForm.dur" is the time step required for network formation.
  // "pres.dur" is the time step required to apply prestrain/prestress.
  // "rheo.dur" is the time of rheological measurement in "condition".
  while (currTimeStep <= netForm.dur + rheo.dur + pres.dur + motActiv.dur) {
	if (currTimeStep == netForm.dur) {
		InitToggleParameters(1);
	}
	if (currTimeStep == netForm.dur + pres.dur + motActiv.dur) {
		InitToggleParameters(2);
	}
	// Excute the main part of the process run every time step.
	SingleProcess();
	// Check and record things after SingleProcess().
	TaskAfterSingleProcess();
  }
  Printf0("All process is done!\n"); 
  // Finalize MPI processes
  MPI_Finalize();
}

void TaskAfterSingleProcess(void) {
  char fn[80];

  /*------------------------------ Record progress ---------------------------*/
  // Record progress information and the mechanical energy of networks.
  if (currTimeStep % recProg.prd == 0 || currTimeStep == 1) {
    RecordProgress(); 
  }
  /*------------------------- Rheological measurement ------------------------*/
  // Choose actin segments which will be traced for microrheology
  if (currTimeStep > netForm.dur) {
	if (rheoWay == 0 || rheoWay == 2) {
		// Trajectory of actin segments is recorded at the interval of
		// "recTraj.prd"
	  	if (currTimeStep == netForm.dur + 1 && recTraj.gTglCho != 0) {
			ChooseTrajectoryList(); 
		}
	    if ((currTimeStep - netForm.dur) % recTraj.prd == 0 
				|| currTimeStep == netForm.dur + 1) {
			if (recTraj.tgl != 0 || recTraj.tgl2 != 0) {
		        RecordTrajectory(recTraj.prd);
			}
		}
	}
	// In bulk rheology, trajectory of actin segments is recorded at the 
	// interval of "recStre.prd".
	if (rheoWay > 0) {
		RecordStress(recStre.prd);
	}
  }
  if (recSecStre.tgl != 0) {
	// Measure internal elastic and viscous stress
	RecordElasViscStress(recSecStre.prd, nSecStreDiv, recSecStre.mode);
  }
  /*----------------------- Boundaries and subdomains ------------------------*/
  if ((bndUnb.gTgl != 0 || bndReb.gTgl != 0) && recBndUnbReb.tgl != 0) {
	if (currTimeStep % recBndUnbReb.prd == 0) {
		RecordBoundaryActinUnbindBind();
	}
  }
  if (bndMv.gTgl != 0 && recBndLoc.tgl != 0 && rank == 0) {
	if (currTimeStep % recBndLoc.prd == 0) {
		RecordBoundaryLocation();
	}
  }
  if (bndUnb.gTgl != 0 && bndMat.gTgl != 0 && recBndActMat.tgl != 0) {
	if (currTimeStep % recBndActMat.prd == 0) {
		RecordBoundaryActinMature();
	}
  }
  if (bndReb.gTgl != 0 && bndReb.gTglSpr == 1 && recBndTracF.tgl != 0) {
	if (currTimeStep % recBndTracF.prd == 0 && currTimeStep >= netForm.dur) {
		RecordBoundaryTractionForce();
	}
  }
  /*--------------------------- Network information --------------------------*/
  if (recConf.tgl != 0) {
	// Record network information including the chains and positions of 
	// particles.
	if (currTimeStep % recConf.prd == 0 || currTimeStep == 1) {
		RecordConfig(GenFileName("AccuConf"), recConf.prd);
	}
  }
  // Update iSupp[]
  if (findSupp.tgl != 0) { 
	if (currTimeStep % findSupp.prd == 0) {
		FindSupportiveFramework(); 
	}
  }
  if (recConfVmd.tgl != 0) {
	// Record networks for visualization via VMD. 
	if (currTimeStep % recConfVmd.prd == 0 || currTimeStep == 1) {
		RecordConfigVmd(0);
	}
  }
  if (recConfPV.tgl != 0) {
	if (currTimeStep % recConfPV.prd == 0) {
		RecordConfigParaview();
	}
  }
  if (recConfMlb.tgl != 0) {
	if (currTimeStep % recConfMlb.prd == 0) {
		RecordConfigMatlab();
	}
  }
  /*--------------------------- Energy calculation ---------------------------*/
  if (recE.tgl != 0) {
	if (currTimeStep % recE.prd == 0 || currTimeStep == 1) {
	    RecordMechEnergy(0);
		if (recPerc.tgl != 0) { 
			RecordMechEnergy(1); 
		}	
		if (findSupp.tgl != 0) {
			RecordMechEnergy(2);
		}
	}
  }
  /*-------- Information in unit of elements, segments, and filaments --------*/
  if (recInfo.tgl != 0) { 
	if (currTimeStep % recInfo.prd == 0) {
		RecordIndvSegFilaInformation(recInfo.prd);
	}
  }
  if (confVmdInfo.tgl != 0) {
	// If a specific coloring method is used, the accumulated values of chain  
	// lengths and forces have to recorded.
    RecordAccuLengthForces();
	if (currTimeStep % confVmdInfo.prd == 0) {
		ResetAccuLengthForces();
	}
  }
  if (recLongF.tgl != 0) {
	RecordAccuLongSpringForces(recLongF.prd);
  }
  /*------------------ Dynamics behaviors of actin and ABPs ------------------*/

  if (currTimeStep % recProg.prd == 0) {
	RecordCounters();
  }
  if (recAbpDyn.tgl != 0) {
	if (currTimeStep % recAbpDyn.prd == 0) {
		RecordAbpDynamics(recAbpDyn.prd);
	}
  }
  if (recAbpUnb.tgl != 0) {
	if (currTimeStep % recAbpUnb.prd == 0) {
  		RecordAbpUnbindEvent(-1, -1, 2);
	}
  }
  if (recAbpBind.tgl != 0) {
	if (currTimeStep % recAbpBind.prd == 0) {
  		RecordAbpBindEvent(-1, -1, 2);
	}
  }
  if (recAbpTurn.tgl != 0) {
	if (currTimeStep % recAbpTurn.prd == 0) {
  		RecordAbpTurnover(-1, -1, -1, 2);
	}
  }
  if (recActSev.tgl != 0) {
	if (currTimeStep % recActSev.prd == 0) {
  		RecordActinSeverEvent();
	}
  }   
  if (motSA.gTgl != 0) { 
	if (recMotSize.tgl != 0) {
		if (currTimeStep % recMotSize.prd == 0 || currTimeStep == 1) {
			RecordMotorFilamentSize();
		}
	}
	if (recMotPos.tgl != 0) {
		if (currTimeStep >= netForm.dur
				&& (currTimeStep - netForm.dur) % recMotPos.prd == 0) {
			RecordMotorPosition();
		}
	} 
  }
  /*-------------------- Structural properties of network --------------------*/
  if (recCrsDist.tgl != 0) {
	if (currTimeStep % recCrsDist.prd == 0 || currTimeStep == 1) {
		RecordCrosslinkDistance();
	}
  }
  if (recFilaL.tgl != 0) {
	if (currTimeStep % recFilaL.prd == 0) {
		RecordFilamentLength(recFilaL.prd);
	}
  }
  if (recConn.tgl != 0) {
	if (currTimeStep % recConn.prd == 0 ) {
		RecordConnectivity(1);
	}
  }
  if (recPerc.tgl != 0) { 
	if (currTimeStep % recPerc.prd == 0) {
		FindPercolatedFilaments(dirRecPerc, recPerc.prd); 
	}
  }

  if (recPoreSize.tgl != 0) {
	if (currTimeStep % recPoreSize.prd == 0) {
		RecordPoreSize();
	}
  }   
  /*-------------------------------- Membrane --------------------------------*/
  if (gTglMb != 0) {
	if (mbFix.gTgl != 0 && (currTimeStep % (recProg.prd / 10) == 0 
			|| currTimeStep == 1)) {
		RecordMembraneAdhesion();
	}
	if (recMbDim.tgl != 0) { 
		if (currTimeStep % recMbDim.prd == 0) {
			RecordMembraneDimension(recMbDim.prd);
		}
	}
	if (recMbCen.tgl != 0) { 
		if (currTimeStep % recMbCen.prd == 0) {
			RecordMembraneCenter(recMbCen.prd);
		}
	}
  }
  /*---------------------------------- Bead ----------------------------------*/
  if (gTglBead != 0) {
	if (currTimeStep >= netForm.dur + motActiv.dur) {
		RecordBeadForceLocation();
	}
	if (beadBind.gTgl != 0 && currTimeStep % recProg.prd == 0) {
		RecordBeadActinBinding();
	}
  }

  // If an error signal appears, the run shoulde be terminated.
  if (stopSig != 0) {
	printf("Due to the large force, the code is stopped at currTimeStep = "
		"%lld!\n", currTimeStep);
	while(1);
		exit(-1);
	}
}

void PrepareStateWoNetworkData(void) {
  int n, nIFilaP;

  nAct = 0;
  nActMe = 0;
  nActCp = 0;
  actM.c = 0; 
  nActFilaMe = 0;

  nAbp = 0;
  nAbpMe = 0;
  nAbpCp = 0; 
  if (gTglImpAcpM != 0) { acpM.c = 0; }
  if (gTglImpMotM != 0) { motM.c = 0; }

  recTraj.actMe.c = 0; 
  recTraj.abpMe.c = 0; 
  if (recTraj.gTglCho == 0 && (rheoWay == 0 || rheoWay == 2)) {
	Printf0("Warning: new components won't be chosen for tracking trajectory "
			"without loaded network data. Nothing will be traced!\n\n");
  }

  if (rheoWay > 0) { 
	meaStreParMe.c = 0;
	appStraParMe.c = 0;
  }
  neigh.c = 0;
  if (gTglMb != 0 ) { 
	  neighMb.c = 0;
  }
  if (actSev.gTgl != 0 || actNuc.gTgl != 0 || actAnn.gTgl != 0) {
	nIFilaP = (int)(nActGoal / nCpu);
	for(n = 0; n < nIFilaP; n++) {
		iFilaP.l[n] = rank * nIFilaP + n;
	}
	iFilaP.c = nIFilaP;
  }
  UpdateChainList();
  
  if (rheoWay > 0) {
	if (bndReb.gTgl == 0) {
		Printf0("Error: if Config is not loaded, the binding of actins "
				"on boundaries should be allowed for bulk rheology!\n\n");
		exit(-1);
	}
  }
}
