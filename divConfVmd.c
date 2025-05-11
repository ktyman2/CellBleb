#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  int n, k, mode, cnt, tgl;
  char str[100], fn[80];
  FILE *fIn, *fOut;
  
  if (argc != 2) {
	printf("Error: there should be only one argument!\n");
	exit(-1);
  }
  if (strcmp(argv[1], "all") == 0) {
	mode = -1;
  }
  else {
	mode = atoi(argv[1]);
  }
  fIn = fopen("ConfVmd.pdb", "r");
  cnt = 0;
  tgl = 0;
  while(fgets(str, 80, fIn) != NULL) {
	if (str[0] == 'R' && str[1] == 'E' && str[2] == 'M') {
		if (mode == -1 || mode == cnt) {
			sprintf(fn, "ConfVmd_%d.pdb", cnt);
			fOut = fopen(fn, "w");
		}
		tgl = 1;
	}
	if (tgl == 1 && (mode == -1 || mode == cnt)) {
		fputs(str, fOut);
	}
	if (str[0] == 'E' && str[1] == 'N' && str[2] == 'D') {
		if (mode == -1 || mode == cnt) {
			fclose(fOut);
		}
		tgl = 0;
		cnt++;
	}
  }
  fclose(fIn);

  fIn = fopen("ConfVmd.psf", "r");
  cnt = 0;
  tgl = 0;
  while(fgets(str, 80, fIn) != NULL) {
	if (str[0] == 'P' && str[1] == 'S' && str[2] == 'F') {
		if (mode == -1 || mode == cnt) {
			sprintf(fn, "ConfVmd_%d.psf", cnt);
			fOut = fopen(fn, "w");
		}
		tgl = 1;
	}
	if (str[0] == '-' && str[1] == '-' && tgl == 1) {
		if (mode == -1 || mode == cnt) {
			fclose(fOut);
		}
		tgl = 0;
		cnt++;
	}
	if (tgl == 1 && (mode == -1 || mode == cnt)) {
		fputs(str, fOut);
	}
  }
  fclose(fIn);

  if (mode > cnt) {
	printf("Error: the frame which you want is not included in these files!\n");
	exit(-1);
  }
}


