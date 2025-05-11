// This file contains the main process which calls many functions to perform
// a simulation.

#include "common.c"

int main (int argc, char *argv[]) { 
  // The general initialization for MP
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nCpu);

  // Check whether 'condition' and 'parallel' exist.
  InitFileCheck();
  // Initialize the number of cells in x,y,z directions based on 'parallel'
  InitCellNumber();
  // Load initial parameters from two files, "condition" and "Config"
  LoadInitParameter(); 
  // Initialize a run - initialize the random-number generator, 
  // assign values to variables, define arrays,and assign values to the arrays
  InitRun();    
  if (gTglMb != 0) { InitMembraneInformation(); }
  if (gTglBead != 0) {
	InitBeadInformation();
  }
  // If network data are loaded from given files
  if (gTglLoadNetData != 0) {
	// Load configuration file from "Config". (Primarily chain and position)
	LoadConfig(1);	
	if (gTglMb != 0 && gTglLoadMbNucData == 0) { TrimNetworkForMembrane(); }
	// Prepare procedures for bulk rheology measurement.
	// (Sever filaments crossing boundaries, eliminate free actins)
	if (rheoWay > 0) { PrepareBulkRheology(); }
	else { PrepareNotBulkRheology(); }
	// Delete a portion of ABPs depending on the given ABP concentrations
	DeleteAbp();
	// Eliminate free actin segments
	if (actAss.gTgl == 0 && actDis.gTgl == 0 && actBst.gTgl == 0) { 
		DeleteFreeActin(); 
	}
	// Extract information for a current subdomain from a whole network
	ExtractConfig();
	// Inspect the integrity of data loaded from Config
	InspectChainListAndLength();
  }
  else {
	PrepareStateWoNetworkData();
  }
  // If necessary, add free actins and ABPs
  AddFreeActin();
  AddFreeAbp();
  AssignMembraneSlideList();

  // Update the neighboring list
  UpdateNeighborList();
  // Record the loaded parameter values
  if (rank == 0) { RecordInitParameter(); }
  // Main process
  MainProcess();
  return 0;
}
