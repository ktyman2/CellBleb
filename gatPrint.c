// This file contains functions gathering and printing data.

// Gather all information about the positions and chains of actin to main CPU
void GatherActChainPosition(int *cntAll, int *chain, double *posi) {
  MPI_Gather(&nActMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nActMe, cntAll, NDIM, act.id, act.r, posi); 
  Gather2dArrayInt(nActMe, cntAll, nChAc, act.id, act.ch, chain); 
}

// Gather all information about the positions and chains of ABP to main CPU
void GatherAbpChainPosition(int *cntAll, int *chain, double *posi) {
  MPI_Gather(&nAbpMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nAbpMe, cntAll, NDIM, abp.id, abp.r, posi); 
  Gather2dArrayInt(nAbpMe, cntAll, nChAb, abp.id, abp.ch, chain); 
}

// Gather all information about the positions and chains of ABP to main CPU
void GatherMembChainPosition(int *cntAll, int *chain, double *posi) {
  MPI_Gather(&nMbMe, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  Gather2dArrayDouble(nMbMe, cntAll, NDIM, memb.id, memb.r, posi); 
  Gather2dArrayInt(nMbMe, cntAll, nChMb, memb.id, memb.ch, chain); 
}


// print string on screen and also
// save it into a file
void Printf0(const char *str, ...) {
  FILE *fOut;
  va_list args;

  if (rank == 0) {
	va_start(args, str);
	vprintf(str, args);
	va_end(args);
	fOut = fopen(GenFileName(OUTPUT_FILE), "a");
	va_start(args, str);
	vfprintf(fOut, str, args);
	va_end(args);
	fclose(fOut);
  }
}

void Printf(const char *str, ...) {
  FILE *fOut;
  fOut = fopen(GenFileName(OUTPUT_FILE), "a");
  va_list args;
  va_start(args, str);
  vprintf(str, args);
  va_end(args);
  va_start(args, str);
  vfprintf(fOut, str, args);
  va_end(args);
  fclose(fOut);
}

// Record 1-D array
// mode = 0: print elements horizontally, 1: print elements vertically
// For integer array
void Fprintf1dArrayInt(FILE *fOut, int *arr, int cnt, int mode) {
  int k;

  for (k = 0; k < cnt; k++) { 
	fprintf(fOut, "%d", arr[k]); 
	fprintf(fOut, ((mode == 0) ? "\t" : "\n"));
  }
  if (mode == 0) { fprintf(fOut, "\n"); }
}

// For double array
void Fprintf1dArrayDouble(FILE *fOut, double *arr, int nCol, int mode) {
  int k;

  for (k = 0; k < nCol; k++) { 
	fprintf(fOut, "%g", arr[k]); 
	fprintf(fOut, ((mode == 0) ? "\t" : "\n"));
  }
  if (mode == 0) { fprintf(fOut, "\n"); }
}

// Record 1-D array horizontally without returning at end.
// For integer array
void Fprintf1dArrayIntWoRet(FILE *fOut, int *arr, int cnt) {
  int k;

  for (k = 0; k < cnt; k++) { 
	fprintf(fOut, "%d\t", arr[k]); 
  }
}

void Fprintf1dArrayDoubleWoRet(FILE *fOut, double *arr, int cnt) {
  int k;

  for (k = 0; k < cnt; k++) { 
	fprintf(fOut, "%g\t", arr[k]); 
  }
}

// Record 1-D array horizontally and fill the rest of columns with a certain
// number.
// For integer array
void Fprintf1dArrayIntWFil(FILE *fOut, int *arr, int nCol, 
		int fil, int nTotCol) {
  Fprintf1dArrayIntWoRet(fOut, arr, nCol);
  Fprintf1dFillerInt(fOut, fil, nTotCol - nCol, 0);
}

// For double array
void Fprintf1dArrayDoubleWFil(FILE *fOut, double *arr, int nCol, 
		int fil, int nTotCol) {
  Fprintf1dArrayDoubleWoRet(fOut, arr, nCol);
  Fprintf1dFillerInt(fOut, fil, nTotCol - nCol, 0);
}


// Record 2-D array with mixed contents in the format of double and integer. 
// "chk" is a marker for indicating the format.
void Fprintf2dArrayIntDouble(FILE *fOut, double *arr, int nRow, 
		int nCol, int *chk) {
  int n, k;

  for(n = 0; n < nRow; n++) {
	for(k = 0; k < nCol; k++) {
		if (chk[k] == 1) {
			fprintf(fOut, "%d\t", (int)P2A(arr,n,k,nCol));
		}
		else {
			fprintf(fOut, "%g\t", P2A(arr,n,k,nCol));
		}
	}
	fprintf(fOut, "\n");
  }
}

// Record 2-D array.
// mode = -1: no index, >= 0: index beginning from "mode"
// mode2 = 0: no blank line at end, 1: a blank line at end
// For interger array
void Fprintf2dArrayInt(FILE *fOut, int *arr, int nRow, int nCol, 
		int mode, int mode2) {
  int n;

  for(n = 0; n < nRow; n++) {
    if (mode > -1) { fprintf(fOut, "%d\t", mode + n); }
	Fprintf1dArrayInt(fOut, &P2A(arr,n,0,nCol), nCol, 0);
  }
  if (mode2 != 0) { fprintf(fOut, "\n"); }
}

// For double array
void Fprintf2dArrayDouble(FILE *fOut, double *arr, int nRow, int nCol, 
		int mode, int mode2) {
  int n;

  for(n = 0; n < nRow; n++) {
    if (mode > -1) { fprintf(fOut, "%d\t", mode + n); }
	Fprintf1dArrayDouble(fOut, &P2A(arr,n,0,nCol), nCol, 0); 
  }
  if (mode2 != 0) { fprintf(fOut, "\n"); }
}

// Record 2-D array with the index beginning from "startInd" and fills the 
// rest of columns with a certain number.
// mode = -1: no index, >= 0: index beginning from "mode"
// For interger array
void Fprintf2dArrayIntWFil(FILE *fOut, int *arr, int nRow,
        int nCol, int mode, int fil, int nTotCol) {
  int n;

  for(n = 0; n < nRow; n++) {
    fprintf(fOut, "%d\t", mode + n);
	Fprintf1dArrayIntWFil(fOut, arr, nCol, fil, nTotCol);
  }
}

// For double array
void Fprintf2dArrayDoubleWFil(FILE *fOut, double *arr, int nRow,
        int nCol, int mode, int fil, int nTotCol) {
  int n;
  for(n = 0; n < nRow; n++) {
    fprintf(fOut, "%d\t", mode + n);
	Fprintf1dArrayDoubleWFil(fOut, arr, nCol, fil, nTotCol);
  }
}

// Record 2-D array having chain information. At the beginning of each line,
// the index of each row is recorded. 
// mode = 1: print a label, 0: don't print it.
// For integer array
void RecordChainArrayInt(FILE *fOut, const char *name, int *arr, 
		int nRow, int nCol, int mode) {
  if (mode == 1) { fprintf(fOut, "%s\n\n", name); }
  Fprintf2dArrayInt(fOut, arr, nRow, nCol, 0, 1);
}

// For double array
void RecordChainArrayDouble(FILE *fOut, const char *name, double *arr, 
		int nRow, int nCol, int mode) {
  if (mode == 1) { fprintf(fOut, "%s\n\n", name); }
  Fprintf2dArrayDouble(fOut, arr, nRow, nCol, 0, 1);
}

// Gather local 1-D array (sendArr) to "recvData" at the main node.
// "cntLocal" is the number of elements in "sendArr" in each CPU, which is
// gathered to "cntAll". 
// For interger array.
void Gather1dArrayIntWoInd(int cntLocal, int *cntAll, int *sendArr, 
		int *recvData) {
  int n, cnt, sumCnt, tag = 0;

  if (rank > 0) {
    MPI_Send(sendArr, cntLocal, MPI_INT, 0, tag, MPI_COMM_WORLD);
  }
  else {
	sumCnt = SumArrInt(cntAll, nCpu);
	for(n = 0; n < cntLocal; n++) {
		recvData[n] = sendArr[n];
	}
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
        MPI_Recv(recvData + cnt, cntAll[n], MPI_INT, n, tag, 
			MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
    qsort(recvData, sumCnt, sizeof(int), CompInt);
  }
}

// For double array.
void Gather1dArrayDoubleWoInd(int cntLocal, int *cntAll, double *sendArr, 
		double *recvData) {
  int n, cnt, sumCnt, tag = 0;

  if (rank > 0) {
    MPI_Send(sendArr, cntLocal, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
  else {
	sumCnt = SumArrInt(cntAll, nCpu);
	for(n = 0; n < cntLocal; n++) {
		recvData[n] = sendArr[n];
	}
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
        MPI_Recv(recvData + cnt, cntAll[n], MPI_DOUBLE, n, tag, 
			MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
    qsort(recvData, sumCnt, sizeof(double), CompDbl);
  }
}

// Gather local 1-D array (sendArr) at the main node. Then, record the data
// to a file (fOut).
// "cntLocal" is the number of elements in "sendArr" in each CPU.
// For interger array
void RecordGather1dArrayIntWoInd(int cntLocal, int *cntAll, int *sendArr, 
		FILE *fOut) {
  int *recvData, sumCnt;

  if (rank == 0) { 
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvData,int,sumCnt); 
  }
  Gather1dArrayIntWoInd(cntLocal, cntAll, sendArr, recvData);
  if (rank == 0) {
	Fprintf1dArrayInt(fOut, recvData, sumCnt, 1);
	fprintf(fOut, "\n");
    free(recvData);
  }
}

// For double array
void RecordGather1dArrayDoubleWoInd(int cntLocal, int *cntAll, double *sendArr, 
		FILE *fOut) {
  int sumCnt;
  double *recvData;

  if (rank == 0) { 
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvData,double,sumCnt); 
  }
  Gather1dArrayDoubleWoInd(cntLocal, cntAll, sendArr, recvData);
  if (rank == 0) {
	Fprintf1dArrayDouble(fOut, recvData, sumCnt, 1);
	fprintf(fOut, "\n");
    free(recvData);
  }
}

// Gather local 1-D array (sendArr) at the main node. Then, record the data
// to a file (fOut).
// "cntLocal" is the number of elements in "sendArr" in each CPU, which is
// gathered to "cntAll". 
// For interger array.
void RecordGather1dArrayIntWoIndWoCnt(int cntLocal, int *sendArr, 
		FILE *fOut, const char *str) {
  int *cntAll;

  MALLOC(cntAll,int,nCpu); 
  MPI_Gather(&cntLocal, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0 && strcmp(str, "none") != 0) {
	fprintf(fOut, "%s = %d\n", str, SumArrInt(cntAll, nCpu)); 
  }
  RecordGather1dArrayIntWoInd(cntLocal, cntAll, sendArr, fOut);
  free(cntAll);
}

// For double array.
void RecordGather1dArrayDoubleWoIndWoCnt(int cntLocal, double *sendArr, 
		FILE *fOut, const char *str) {
  int *cntAll;

  MALLOC(cntAll,int,nCpu); 
  MPI_Gather(&cntLocal, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0 && strcmp(str, "none") != 0) {
	fprintf(fOut, "%s = %d\n", str, SumArrInt(cntAll, nCpu)); 
  }
  RecordGather1dArrayDoubleWoInd(cntLocal, cntAll, sendArr, fOut);
  free(cntAll);
}

// Gather local 2-D array (sendArr) in "recvData" at the main node. 
// Also, the data in "recvData" are sorted by "index".
// "cntLocal" is the number of elements in "sendArr" in each CPU, which is
// gathered to "cntAll". 
// For interger array.
void Gather2dArrayInt(int cntLocal, int *cntAll, int nCol, int *index, 
		int *sendArr, int *recvData) {
  int n, k, cnt, sumCnt, tag = 0;
  int *sendBuf, *recvBuf;

  if (rank > 0) {
	MALLOC(sendBuf,int,cntLocal*(nCol+1));  
    for(n = 0; n < cntLocal; n++) {
        P2A(sendBuf,n,0,nCol + 1) = index[n];
		for(k = 0; k < nCol; k++) {
	       P2A(sendBuf,n,k + 1,nCol + 1) = P2A(sendArr,n,k,nCol);
		}
    }
    MPI_Send(sendBuf, (nCol + 1) * cntLocal, MPI_INT, 0, tag, 
		MPI_COMM_WORLD);
	free(sendBuf);
  }
  else {
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvBuf,int,sumCnt*(nCol+1));
    for(n = 0; n < cntLocal; n++) {
        P2A(recvBuf,n,0,nCol + 1) = index[n];
		for(k = 0; k < nCol; k++) {
	       P2A(recvBuf,n,k + 1,nCol + 1) = P2A(sendArr,n,k,nCol);
		}
    }
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
        MPI_Recv(recvBuf + cnt * (nCol + 1), cntAll[n] * 
				(nCol + 1), MPI_INT, n, tag, MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
    qsort(recvBuf, sumCnt, (nCol + 1) * sizeof(int), 
			CompInt);
	for(n = 0; n < sumCnt; n++) {
		for(k = 0; k < nCol; k++) {
			P2A(recvData,n,k,nCol) = P2A(recvBuf,n,k + 1,nCol + 1);
		}
	}
	free(recvBuf);
  }
}

// For double array.
void Gather2dArrayDouble(int cntLocal, int *cntAll, int nCol, int *index, 
		double *sendArr, double *recvData) {
  int n, k, cnt, sumCnt, tag = 0;
  double *sendBuf, *recvBuf;

  if (rank > 0) {
	MALLOC(sendBuf,double,cntLocal*(nCol+1));  
    for(n = 0; n < cntLocal; n++) {
        P2A(sendBuf,n,0,nCol + 1) = (double)index[n];
		for(k = 0; k < nCol; k++) {
	       P2A(sendBuf,n,k + 1,nCol + 1) = P2A(sendArr,n,k,nCol);
		}
    }
    MPI_Send(sendBuf, (nCol + 1) * cntLocal, MPI_DOUBLE, 0, tag, 
			MPI_COMM_WORLD);
	free(sendBuf);
  }
  else {
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvBuf,double,sumCnt*(nCol+1));
    for(n = 0; n < cntLocal; n++) {
        P2A(recvBuf,n,0,nCol + 1) = (double)index[n];
		for(k = 0; k < nCol; k++) {
	       P2A(recvBuf,n,k + 1,nCol + 1) = P2A(sendArr,n,k,nCol);
		}
    }
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
        MPI_Recv(recvBuf + cnt * (nCol + 1), cntAll[n] * 
				(nCol + 1), MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
    qsort(recvBuf, sumCnt, (nCol + 1) * sizeof(double), CompDbl);
	for(n = 0; n < sumCnt; n++) {
		for(k = 0; k < nCol; k++) {
			P2A(recvData,n,k,nCol) = P2A(recvBuf,n,k + 1,nCol + 1);
		}
	}
	free(recvBuf);
  }
}

// Gather and record local 2-D array (sendArr) at the main node. 
// Also, the gathered data are sorted by "index".
// "cntLocal" is the number of elements in "sendArr" in each CPU, which is
// gathered to "cntAll". 
// mode = -1: no index, >= 0: index beginning from "mode"
// For interger array.
void RecordGather2dArrayInt(int cntLocal, int *cntAll, int nCol, 
		int *index, int *sendArr, FILE *fOut, int mode) {
  int sumCnt, *recvData;

  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu); 
	MALLOC(recvData,int,sumCnt*nCol);
  }
  Gather2dArrayInt(cntLocal, cntAll, nCol, index, sendArr, recvData);
  if (rank == 0) { 
	Fprintf2dArrayInt(fOut, recvData, sumCnt, nCol, mode, 1); 
	free(recvData);
  }
}

// For double array.
void RecordGather2dArrayDouble(int cntLocal, int *cntAll, int nCol, 
		int *index, double *sendArr, FILE *fOut, int mode) {
  int sumCnt;
  double *recvData;

  if (rank == 0) {
	sumCnt = SumArrInt(cntAll, nCpu);
	MALLOC(recvData,double,sumCnt*nCol);
  }
  Gather2dArrayDouble(cntLocal, cntAll, nCol, index, sendArr, recvData);
  if (rank == 0) {
	Fprintf2dArrayDouble(fOut, recvData, sumCnt, nCol, mode, 1); 
	free(recvData);
  }
}

// Gather and record local 2-D array (sendArr) at the main node. 
// Also, the gathered data are sorted by "index".
// "cntLocal" is the number of elements in "sendArr" in each CPU.
// mode = -1: no index, >= 0: index beginning from "mode"
// For interger array.
void RecordGather2dArrayIntWoCnt(int cntLocal, int nCol,
        int *index, int *sendArr, FILE *fOut, int mode) {
  int *cntAll;

  MALLOC(cntAll,int,nCpu); 
  MPI_Gather(&cntLocal, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  RecordGather2dArrayInt(cntLocal, cntAll, nCol, index, sendArr, fOut, mode);
  free(cntAll);
}

// For double array.
void RecordGather2dArrayDoubleWoCnt(int cntLocal, int nCol,
        int *index, double *sendArr, FILE *fOut, int mode) {
  int *cntAll;

  MALLOC(cntAll,int,nCpu); 
  MPI_Gather(&cntLocal, 1, MPI_INT, cntAll, 1, MPI_INT, 0, MPI_COMM_WORLD);
  RecordGather2dArrayDouble(cntLocal, cntAll, nCol, index, sendArr, fOut, mode);
  free(cntAll);
}

void Gather2dArrayIntWoIndWoSort(int cntLocal, int *cntAll, int *sendArr, 
		int *recvData, int col) {
  int n, cnt, tag = 0;

  if (rank > 0) {
	if (cntLocal > 0) { 
	    MPI_Send(sendArr, cntLocal * col, MPI_INT, 0, tag, MPI_COMM_WORLD);
	}
  }
  else {
	for(n = 0; n < cntLocal * col; n++) {
		recvData[n] = sendArr[n];
	}
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
		CONT(cntAll[n] == 0);
        MPI_Recv(&P2A(recvData,cnt,0,col), cntAll[n] * col, MPI_INT, n, tag, 
			MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
  }
}

void Gather2dArrayDoubleWoIndWoSort(int cntLocal, int *cntAll, double *sendArr, 
		double *recvData, int col) {
  int n, cnt, tag = 0;

  if (rank > 0) {
	if (cntLocal > 0) { 
	    MPI_Send(sendArr, cntLocal * col, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
  }
  else {
	for(n = 0; n < cntLocal * col; n++) {
		recvData[n] = sendArr[n];
	}
    cnt = cntAll[0];
    for(n = 1; n < nCpu; n++) {
		CONT(cntAll[n] == 0);
        MPI_Recv(&P2A(recvData,cnt,0,col), cntAll[n] * col, MPI_DOUBLE, n, tag, 
			MPI_COMM_WORLD, &status);
        cnt += cntAll[n];
    }
  }
}

void Fprintf1dFillerInt(FILE *fOut, int fil, int cnt, int mode) {
  int k;

  for (k = 0; k < cnt; k++) { 
	fprintf(fOut, "%d\t", fil); 
  }
  if (mode == 0) { 
	fprintf(fOut, "\n");
  }
}

void Fprintf1dFillerDouble(FILE *fOut, double fil, int cnt, int mode) {
  int k;

  for (k = 0; k < cnt; k++) { 
	fprintf(fOut, "%g\t", fil); 
  }
  if (mode == 0) { 
	fprintf(fOut, "\n");
  }
}
