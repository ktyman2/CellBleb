#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  int n, k, mode, cnt, tag[6];
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

  fIn = fopen("ConfPV", "r");
  fscanf(fIn, "%d %d %d %d %d %d\n", &tag[0], &tag[1], &tag[2], &tag[3], 
		&tag[4], &tag[5]);

  cnt = 0;
  while(fgets(str, 100, fIn) != NULL) {
	for(n = 0; n < 6; n++) {
		if (tag[n] == 0) {
			continue;
		}
		if (n == 0) {
			sprintf(fn, "pv_actin_%d.vtk", cnt);
		}
		else if (n == 1) {
			sprintf(fn, "pv_ACPC_%d.vtk", cnt);
		}
		else if (n == 2) {
			sprintf(fn, "pv_ACPB_%d.vtk", cnt);
		}
		else if (n == 3) {
			sprintf(fn, "pv_motor_%d.vtk", cnt);
		}
		else if (n == 4) {
			sprintf(fn, "pv_memb_%d.vtk", cnt);
		}
		else {
			sprintf(fn, "pv_bound_%d.vtk", cnt);
		}
  		fgets(str, 100, fIn);
		if (mode == cnt || mode == -1) {
			fOut = fopen(fn, "w");
			fputs(str, fOut);
		}
		while(1) {
  			fgets(str, 100, fIn);
			if (str[0] == '#') {
				 break;
			}
			if (mode == cnt || mode == -1) {
				fputs(str, fOut);
			}
		}
		if (mode == cnt || mode == -1) {
			fclose(fOut);
		}
	}
	cnt++;
  }
  fclose(fIn);

  if (mode > cnt) {
	printf("Error: the frame which you want is not included in these files!\n");
	exit(-1);
  }
}


