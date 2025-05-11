// This file contains functions used as tools.

// Find maximum and minimum values in an integer array.
// mode = 0: minimum, 1: maximum
int FindMaxMin1dArrayInt(int *arr, int cnt, int mode) {
  int n, max = NEG_LARGE_VALUE, min = POS_LARGE_VALUE;
  if (mode == 0) {
	for(n = 0; n < cnt; n++) {
		if (min > arr[n]) { min = arr[n]; }
	}
    return min;
  }
  else if (mode == 1) {
	for(n = 0; n < cnt; n++) {
		if (max < arr[n]) { max = arr[n]; }
	}
    return max;
  }
  return 0;
}

// Find maximum and minimum values in a double array.
// mode = 0: minimum, 1: maximum
double FindMaxMin1dArrayDbl(double *arr, int cnt, int mode) {
  int n;
  double max = NEG_LARGE_VALUE, min = POS_LARGE_VALUE;
  if (mode == 0) {
	for(n = 0; n < cnt; n++) {
		if (min > arr[n]) { min = arr[n]; }
	}	
    return min;
  }
  else if (mode == 1) {
	for(n = 0; n < cnt; n++) {
		if (max < arr[n]) { max = arr[n]; }
	}
    return max;
  }
  return 0;
}

void SetAllValue1dArrayInt(int *arr, int sizeArr, int val) {
  int n;

  for(n = 0; n < sizeArr; n++) { 
	arr[n] = val;
  }
}

void SetAllValue1dArrayDouble(double *arr, int sizeArr, double val) {
  int n;

  for(n = 0; n < sizeArr; n++) { 
	arr[n] = val;
  }
}

void Copy1dArrayInt(int *arrRecv, int *arrSrc, int sizeArr) {
  int n;

  for(n = 0; n < sizeArr; n++) { 
	arrRecv[n] = arrSrc[n];
  }
}

void Copy1dArrayDouble(double *arrRecv, double *arrSrc, int sizeArr) {
  int n;

  for(n = 0; n < sizeArr; n++) { 
	arrRecv[n] = arrSrc[n];
  }
}

// Search for (ele) in (pos)th column of array having (cntArr) rows and 
// (col) columns. If found, the number of a row having (ele) is returned.
// If not, -1 is returned.
int FindElementArray(int *arr, int cntArr, int ele, int pos, int col) {
  int n, CS = -1;
  for(n = 0; n < cntArr; n++) {
	if (ele == P2A(arr,n,pos,col)) {
		CS = n;
		break;
	}
  }
  return CS;
}

// Search for two elements (ele1 and ele2) in (start)~(start+1) columns 
// of array having (cntArr) rows and (col) columns. If found, the number 
// of a row having (ele) is returned. If not, -1 is returned.
int Find2ElementArray(int *arr, int cntArr, int ele1, int ele2, 
		int start, int col) {
  int n, CS = -1;
  for(n = 0; n < cntArr; n++) {
	if (ele1 == P2A(arr,n,start,col) && ele2 == P2A(arr,n,start + 1,col)) {
		CS = n;
		break;
	}
  }
  return CS;
}

// Search for two elements (ele1, ele2, and ele3) in (start)~(start+2) 
// columns of array having (cntArr) rows and (col) columns. If found, the 
// number of a row having (ele) is returned. If not, -1 is returned.
int Find3ElementArray(int *arr, int cntArr, int ele1, int ele2, int ele3,
		int start, int col) {
  int n, CS = -1;
  for(n = 0; n < cntArr; n++) {
	if (ele1 == P2A(arr,n,start,col) && ele2 == P2A(arr,n,start + 1,col)
			&& ele3 == P2A(arr,n,start + 2,col)) {
		CS = n;
		break;
	}
  }
  return CS;
}

// First, search for an element in array with one column. If not found,
// (ele) is added behind the last element of the array.
int InsertElement1dArrayWChk(int *arr, int *cntArr, int ele) {
  int CS;
  CS = FindElementArray(arr, *cntArr, ele, 0, 1);
  if (CS == -1) {
	arr[*cntArr] = ele;
	(*cntArr)++;
  }
  return CS;
}

// (ele) is added behind the last element of the array without lookup.
void InsertElement1dArrayWoChk(int *arr, int *cntArr, int ele) {
  arr[*cntArr] = ele;
  (*cntArr)++;
}

// Insert a row in an array by a row number
void InsertElementArrayByIndex(int *arr, int *cntArr, int *ele, 
		int idx, int col)
{
  int n;
  switch(col) {
	case 1:
	  	for(n = *cntArr; n >= idx + 1; n--) { 
			arr[n] = arr[n - 1]; 
		}
		arr[idx] = *ele;
		break;
	case 2:
	  	for(n = *cntArr; n >= idx + 1; n--) { 
			V2COPY(&P2A(arr,n,0,2),&P2A(arr,n - 1,0,2));
		}
		V2COPY(&P2A(arr,idx,0,2),ele);
		break;
	case 3:
	  	for(n = *cntArr; n >= idx + 1; n--) { 
			V3COPY(&P2A(arr,n,0,3),&P2A(arr,n - 1,0,3));	
		}
		V3COPY(&P2A(arr,idx,0,3),ele);
		break;
	case 4:
	  	for(n = *cntArr; n >= idx + 1; n--) { 
			V4COPY(&P2A(arr,n,0,4),&P2A(arr,n - 1,0,4));	
		}
		V4COPY(&P2A(arr,idx,0,4),ele);
		break;
	case 5:
	  	for(n = *cntArr; n >= idx + 1; n--) { 
			V5COPY(&P2A(arr,n,0,5),&P2A(arr,n - 1,0,5));	
		}
		V5COPY(&P2A(arr,idx,0,5),ele);
		break;
  }
  (*cntArr)++;
}

// First, search for an element in array with one column. If found,
// delete the row, and the number of a row is returned. If not, -1 is returned.
void DeleteElement1dArray(int *arr, int *cntArr, int ele) {
  int CS;
  CS = FindElementArray(arr, *cntArr, ele, 0, 1);
  if (CS > -1) { DeleteElementArrayByIndex(arr, cntArr, CS, 1); }
}

// First, search for two elements (ele1 and ele2) in (start)~(start+1) 
// columns of array having (cntArr) rows and (col) columns. If found, 
// delete the row, and the number of a row is returned. If not, -1 is returned.
int DeleteTwoElementArray(int *arr, int *cntArr, int ele1, int ele2, 
		int start, int col) {  
  int CS;
  CS = Find2ElementArray(arr, *cntArr, ele1, ele2, start, col);
  if (CS > -1) { DeleteElementArrayByIndex(arr, cntArr, CS, col); }
  return CS;
}

// Delete a row in an array by a row number
void DeleteElementArrayByIndex(int *arr, int *cntArr, int idx, int col) {
  int n;
  (*cntArr)--;
  switch(col) {
	case 1:
	  	for(n = idx; n < *cntArr; n++) { arr[n] = arr[n + 1]; }
		break;
	case 2:
		for(n = idx; n < *cntArr; n++) {
			V2COPY(&P2A(arr,n,0,2),&P2A(arr,n + 1,0,2));	
		}
		break;
	case 3:
		for(n = idx; n < *cntArr; n++) {
			V3COPY(&P2A(arr,n,0,3),&P2A(arr,n + 1,0,3));	
		}
		break;
	case 4:
		for(n = idx; n < *cntArr; n++) {
			V4COPY(&P2A(arr,n,0,4),&P2A(arr,n + 1,0,4));	
		}
		break;
	case 5:
		for(n = idx; n < *cntArr; n++) {
			V5COPY(&P2A(arr,n,0,5),&P2A(arr,n + 1,0,5));	
		}
		break;
	case 6:
		for(n = idx; n < *cntArr; n++) {
			V6COPY(&P2A(arr,n,0,6),&P2A(arr,n + 1,0,6));	
		}
		break;
  }
  if (col > 6) {
	for(n = idx; n < *cntArr; n++) {
		Copy1dArrayInt(&P2A(arr,n,0,col), &P2A(arr,n + 1,0,col), col);	
	}
  }
}

// Find the sum of values in array having integer numbers.
int SumArrInt(int *arr, int cntArr) {
  int n, sum = 0;
  for(n = 0; n < cntArr; n++) {
	sum += arr[n];
  }
  return sum;
}

// Find the sum of values in array having double numbers.
double SumArrDbl(double *arr, int cntArr) {
  int n;
  double sum = 0.;
  for(n = 0; n < cntArr; n++) {
	sum += arr[n];
  }
  return sum;
}

// Find the mean of values in array having integer numbers.
double AvgArrInt(int *arr, int cntArr) {
  int sum;
  sum = SumArrInt(arr, cntArr);
  return ((double)sum / cntArr);
}

// Find the mean of values in array having double numbers.
double AvgArrDbl(double *arr, int cntArr) {
  double sum;
  sum = SumArrDbl(arr, cntArr);
  return ((double)sum / cntArr);
}

void CheckArraySize(ListInt *arr, int *sizeArr, int col, int mode) {
  int n, *arrBack;

  if (arr->c * col >= (*sizeArr) / 2) {
	if (mode == 1) {
		MALLOC(arrBack,int,arr->c * col);
		for(n = 0; n < arr->c * col; n++) {
			arrBack[n] = arr->l[n];
		}
	}
	free(arr->l);
	(*sizeArr) *= 2;
	MALLOC(arr->l,int,*sizeArr);
	if (mode == 1) {
		for(n = 0; n < arr->c * col; n++) {
			arr->l[n] = arrBack[n];
		}
		free(arrBack);
	}
  } 
}

void Check2dArraySize(ListInt2 *arr, int *sizeArr, int col, int mode) {
  int n, k, *arrBack;
 
  if (arr->c >= (*sizeArr) / 2) {
	if (mode == 1) {
		MALLOC(arrBack,int,arr->c * col);
		for(n = 0; n < arr->c; n++) {
			for(k = 0; k < col; k++) {
				P2A(arrBack,n,k,col) = arr->l[n][k];
			}
		}
	}
	for(n = 0; n < arr->c; n++) {
		free(arr->l[n]);
	}
	free(arr->l);
	(*sizeArr) *= 2;
	MALLOC2(arr->l,int,*sizeArr);
	for(n = 0; n < arr->c; n++) {
		MALLOC(arr->l[n],int,col);
		if (mode == 1) {
			for(k = 0; k < col; k++) {
				arr->l[n][k] = P2A(arrBack,n,k,col);
			}
		}
	}
	if (mode == 1) { free(arrBack); }
  }
}

int CheckTriaPntDist(double rT[][NDIM], double *rP, double thres) {
  int CS;
  double len, cen[NDIM], dr[NDIM];
  len = CalcTriaCentroidRadius(rT, cen);

  V3SUB(dr, cen, rP);
  CS = (V3LEN(dr) > thres + len) ? 0 : 1;
  return CS;
}

int CheckTriaSegDist(double rT[][NDIM], double rS[][NDIM], double thres) {
  int CS;
  double lenT, lenS, dr[NDIM], cenT[NDIM], cenS[NDIM];

  lenT = CalcTriaCentroidRadius(rT, cenT);

  V3AVG(cenS,rS[1],rS[0]);
  V3SUB(dr,rS[1],rS[0]);
  lenS = V3LEN(dr);

  V3SUB(dr, cenT, cenS);
  CS = (V3LEN(dr) > thres + lenT + 0.5 * lenS) ? 0 : 1;
  return CS;
}

int CheckTriaTriaDist(double rT1[][NDIM], double rT2[][NDIM], double thres) {
  int CS;
  double lenT1, lenT2, dr[NDIM], cenT1[NDIM], cenT2[NDIM];
  lenT1 = CalcTriaCentroidRadius(rT1, cenT1);
  lenT2 = CalcTriaCentroidRadius(rT2, cenT2);

  V3SUB(dr, cenT1, cenT2);
  CS = (V3LEN(dr) > thres + lenT1 + lenT2) ? 0 : 1;
  return CS;
}

// dr: a vector of [point on the triangle] - [point]
double CalcTriaPntDist(double rT[][NDIM], double *rP, double *dr,
        double *ratio) {
  double s, t, det, distSq, dr2[3][NDIM], lenSq[3], dot[3];
  int n;

  CalcVec(dr2[0], rT[0], rP);
  CalcVec(dr2[1], rT[1], rT[0]);
  CalcVec(dr2[2], rT[2], rT[0]);
  for(n = 0; n < 3; n++) {
    lenSq[n] = V3LEN_SQ(dr2[n]);
    dot[n] = V3DOT(dr2[(n + 1) % 3], dr2[(n + 2) % 3]);
  }
  det = fabs(lenSq[1] * lenSq[2] - SQR(dot[0]));
  s = dot[0] * dot[1] - lenSq[2] * dot[2];
  t = dot[0] * dot[2] - lenSq[1] * dot[1];
  if (s + t <= det) {
    if (s <= 0.) {
        s = 0.;
        if (t <= 0.)    {  // region 4
            t = 0.;
        }
        else {  // region 3
            t = -1. * dot[1] / lenSq[2];
            t = TrimDblVal(t, 0., 1.);
        }
    }
    // s > 0
    else {
        if (t <= 0.) {  // region 5
            s = -1 * dot[2] / lenSq[1];
            t = 0.;
            s = TrimDblVal(s, 0., 1.);
        }
        else {  // region 0
            // minimum at interior point
            s /= det;
            t /= det;
        }
    }
  }
  // s+t > det
  else {
    if (s <= 0.) {  // region 2
        s = 0.;
        t = 1.;
    }
    else {
        if (t <= 0.) {  // region 6
            s = 1.;
            t = 0.;
        }
        else {  // region 1
            s = (lenSq[2] + dot[1] - dot[0] - dot[2]) 
					/ (lenSq[1] - 2 * dot[0] + lenSq[2]);
            t = 1. - s;
        }
    }
  }
  ratio[0] = s;
  ratio[1] = t;
  distSq = s * (lenSq[1] * s + dot[0] * t + 2 * dot[2])
        + t * (dot[0] * s + lenSq[2] * t + 2 * dot[1]) + lenSq[0];
  VSS3ADD(dr, dr2[1], dr2[2], s, t);
  VV3ADD(dr, dr2[0]);
  V3REVSIGN(dr);
  // Account for numerical round-off error.
  if (distSq < 0) { distSq = 0; }
  return sqrt(distSq);
}

// dr: a vector of [point on the triangle] - [point on the segment]
double CalcTriaSegDist(double rT[][NDIM], double rS[][NDIM], double *dr,
		double *ratio) {
  int n, ind;
  double e1[NDIM], e2[NDIM], crsProd[NDIM], dr2[NDIM], u[NDIM], v[NDIM];
  double cen[NDIM], dir[NDIM], rS2[2][NDIM];
  double crsddir, ude1, ude2, uddr, vde1, vde2, vddr, dde1, dde2, dddr;
  double b0, b1, b2, invLen, invDet, minDist, dist, len, ratio2[2];

  V3AVG(cen, rS[0], rS[1]);
  len = CalcUnitVec(dir, rS[0], rS[1]);

  CalcVec(e1,rT[1],rT[0]);
  CalcVec(e2,rT[2],rT[0]);
  V3CROSS(crsProd,e1,e2);
  NormVec(crsProd);
  // Check whether triangle is parallel to line
  crsddir = V3DOT(crsProd,dir);
  minDist = POS_LARGE_VALUE;
  if (fabs(crsddir) > 0) {
	CalcVec(dr2,cen,rT[0]);
    if (fabs(dir[0]) >= fabs(dir[1])) {
        // W.x or W.z is the largest magnitude component, swap them
        invLen = 1 / sqrt(SQR(dir[0]) + SQR(dir[2]));
        V3SET(u, -1. * dir[2] * invLen, 0., dir[0] * invLen);
    }
	else {
        // W.y or W.z is the largest magnitude component, swap them
        invLen = 1 / sqrt(SQR(dir[1]) + SQR(dir[2]));
        V3SET(u, 0., dir[2] * invLen, -1. * dir[1] * invLen);
    }
	V3CROSS(v, dir, u);
	ude1 = V3DOT(u,e1);
	ude2 = V3DOT(u,e2);
	uddr = V3DOT(u,dr2);
	vde1 = V3DOT(v,e1);
	vde2 = V3DOT(v,e2);
	vddr = V3DOT(v,dr2);
	invDet = 1 / (ude1 * vde2 - ude2 * vde1);
	// Barycentric coordinates for the point of intersection.
	b1 = (vde2 * uddr - ude2 * vddr) * invDet;
	b2 = (ude1 * vddr - vde1 * uddr) * invDet;
	b0 = 1 - b1 - b2;
	if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
		dde1 = V3DOT(dir,e1);
		dde2 = V3DOT(dir,e2);
		dddr = V3DOT(dir,dr2);
		V2SET(ratio, b1, b2);
		ratio[2] = -(b1 * dde1 + b2 * dde2 - dddr) / len + 0.5;
		if (ratio[2] < 0.) {
			dist = -1. * ratio[2]  * len;
			VS3COPY(dr, dir, dist);
			ratio[2] = 0.;
		}
		else if (ratio[2] > 1.) {
			dist = (ratio[2] - 1.) * len;
			VS3COPY(dr, dir, -1. * dist);
			ratio[2] = 1.;
		}
		else { 
			dist = 0.; 
			V3COPY(dr, dir);
		}
		minDist = SQR(dist);
	}
  }
  if (!(minDist < POS_SMALL_VALUE)) { 
	for (n = 0; n < 3; n++){
		ind = (n + 1) % 3;
		V3COPY(rS2[0], rT[n]);
		V3COPY(rS2[1], rT[ind]);
		dist = CalcSegSegDist(rS, rS2, dr2, ratio2, 1);
		if (dist < minDist) {
		   ratio[2] = ratio2[0];
		   if (n == 0) {
			   V2SET(ratio, ratio2[1], 0.);
		   }
		   else if (n == 1) {
			   V2SET(ratio, 1. - ratio2[1], ratio2[1]);
		   }
		   else {
			   V2SET(ratio, 0., 1. - ratio2[1]);
   			}
			minDist = dist;
			V3COPY(dr, dr2);
			if (minDist < POS_SMALL_VALUE) { break; }
		}
	}
  }
  V3REVSIGN(dr);
  return sqrt(minDist);
}

// dr: a vector of [point on the triangle 1] - [point on the triangle 2]
double CalcTriaTriaDist(double rT1[][NDIM], double rT2[][NDIM], double *dr,
		double *ratio) {
  int m, n, ind;
  double (*pT)[NDIM], (*pT2)[NDIM], rS[2][NDIM], dr2[NDIM], minDist, dist;
  double ratio2[3];

  minDist = POS_LARGE_VALUE;
  for(m = 0; m < 2; m++) {
	pT = (m == 0) ? rT1 : rT2;
	pT2 = (m == 0) ? rT2 : rT1;
	for(n = 0; n < 3; n++) {
		ind = (n + 1) % 3;
		V3COPY(rS[0], pT2[n]);
		V3COPY(rS[1], pT2[ind]);
		dist = CalcTriaSegDist(pT, rS, dr2, ratio2);
		dist = V3LEN_SQ(dr2);
		if (dist < minDist) {
			V3COPY(ratio, ratio2);
			if (n == 0) {
				ratio[3] = 0.;
			}
			else if (n == 1) {
				ratio[3] = ratio[2];
				ratio[2] = 1. - ratio[3];
			}
			else {
				ratio[3] = 1. - ratio[2];
				ratio[2] = 0.;
			}
			if (m == 1) {
				SwapDbl(&ratio[0], &ratio[2]);
				SwapDbl(&ratio[1], &ratio[3]);
				V3REVSIGN(dr2);
			}
			minDist = dist;
			V3COPY(dr, dr2);
		}
		if (minDist < POS_SMALL_VALUE) { break; }
	}
	if (minDist < POS_SMALL_VALUE) { break; }
  }  
  return sqrt(minDist);
}

// Calculate a distance between segment and point.
// x0: a point, x1, x2: end points of a segment
// ratio: (distance between x1 and a point where the distance is measured on the
// segment) / (distance between x1 and x2)
// dr = x0 - x1
double CalcSegPntDist(double x1[][NDIM], double *x0, double *dr, 
		double *ratio) {
  double a[NDIM], b[NDIM], t;

  V3SUB(a, x1[1], x1[0]);
  V3SUB(b, x1[0], x0);
  ApplyBoundCondVecDiff(a);  
  ApplyBoundCondVecDiff(b);
 
  t = - V3DOT(b, a) / V3LEN_SQ(a);
 
  if ( t >= 0 && t <= 1 ) {
	VS3ADD(dr, b, a, t);
	*ratio = t;
  }
  else if ( t > 1 ) {
	V3ADD(dr, a, b);
	*ratio = 1.;
  }
  else {
	V3COPY(dr, b);
	*ratio = 0.;
  }
  V3REVSIGN(dr);

  return V3LEN(dr);
}

int CheckSegSegDist(double x1[][NDIM], double x2[][NDIM], double thres) {
  int CS;
  double u[NDIM], v[NDIM], dr[3][NDIM];

  V3AVG(u,x1[1],x1[0]);
  V3AVG(v,x2[1],x2[0]);
  V3SUB(dr[0],x1[1],x1[0]);
  V3SUB(dr[1],x2[1],x2[0]);
  V3SUB(dr[2], u, v);
  CS = (V3LEN(dr[2]) > thres + 0.5 * (V3LEN(dr[0]) + V3LEN(dr[1]))) ? 0 : 1;
  return CS;
}


// Calculate a distance between two segments.
// x1: end points of a segment, x3: end points of the other segment
// ratio: (distance between x1 and x2) / (distance between x3 and x4) 
// mode = 0: length, 1: the square of length
// dr: x2 - x1
double CalcSegSegDist(double x1[][NDIM], double x2[][NDIM], double *dr, 
		double *ratio, int mode) {
  double u[NDIM], v[NDIM], w[NDIM], a, b, c, d ,e , D;
  double sc, sN, sD, tc, tN, tD;
  int k;
 
  V3SUB(u,x1[1],x1[0]);
  V3SUB(v,x2[1],x2[0]);
  V3SUB(w,x1[0],x2[0]);

  a = V3DOT(u,u);
  b = V3DOT(u,v);
  c = V3DOT(v,v);
  d = V3DOT(u,w);
  e = V3DOT(v,w);
  D = a * c - b * b;
  sD = D;
  tD = D;
 
  // Compute the line parameters of the two closest points
  if (D < POS_SMALL_VALUE) {  
    sN = 0.;          
    sD = 1.;          
    tN = e;
    tD = c;
  }
  else {                
    sN = b * e - c * d;
    tN = a * e - b * d;
    if (sN < 0.) {      
        sN = 0.;
        tN = e;
        tD = c;
    }
    else if (sN > sD) {  
        sN = sD;
        tN = e + b;
        tD = c;
    }
  }
 
  if (tN < 0.) {         
    tN = 0.;
    if (-d < 0.) { sN = 0.; }
    else if (-d > a) { sN = sD; }
    else { sN = -d; sD = a; }
  }
  else if (tN > tD) {     
    tN = tD;
    if ((-d + b) < 0.) { sN = 0; }
    else if ((-d + b) > a) { sN = sD; }
    else { sN = (-d + b); sD = a; }
  }
  sc = sN / sD;
  tc = tN / tD;
 
  FOR_NDIM(k) {
    dr[k] = -(w[k] + sc * u[k] - tc * v[k]);  // = S1(sc) - S2(tc)
  }
  ratio[0] = sc;
  ratio[1] = tc; 

  if (mode == 0) { return V3LEN(dr); }
  else { return V3LEN_SQ(dr); }
}

double CalcTriaCentroidRadius(double rT[][NDIM], double *cen) {
  int k;
  double maxLenSq, lenSq;

  V3SET_ALL(cen, 0.);
  for(k = 0; k < 3; k++) {
 	VV3ADD(cen, rT[k]);
  }
  V3SCALE(cen, INV(3.));
  maxLenSq = NEG_LARGE_VALUE;
  for(k = 0; k < 3; k++) {
	lenSq = CalcDist(cen, rT[k], 1);
	if (lenSq > maxLenSq) { maxLenSq = lenSq; }
  }
  return sqrt(maxLenSq);
}

void CalcTriaCircumCenter(double rT[][3], double *cen) {
  double denom, lenSq[2], dr[3][NDIM], crsProd[NDIM];

  lenSq[0] = CalcVecDist(dr[0], rT[1], rT[0], 1);
  lenSq[1] = CalcVecDist(dr[1], rT[2], rT[0], 1);
  V3CROSS(crsProd, dr[0], dr[1]);
  denom = 2. * V3LEN_SQ(crsProd);
  VSS3ADD(dr[2], dr[1], dr[0], lenSq[0], REVSIGN(lenSq[1]));
  V3CROSS(cen, dr[2], crsProd);
  V3SCALE(cen, INV(denom));
  VV3ADD(cen, rT[0]);
}

double CalcTriaArea(double *r1, double *r2, double *r3) {
  double len, area,  dr[3][NDIM];

  CalcVec(dr[0], r1, r2);
  CalcVec(dr[1], r1, r3);
  V3CROSS(dr[2], dr[0], dr[1]);
  area = 0.5 * V3LEN(dr[2]);
  return area;
}

double CalcTetraVol(double *r1, double *r2, double *r3, double *r4) {
  double dr[3][NDIM], crosProd[NDIM], dotProd, vol;

  CalcVec(dr[0], r2, r1);
  CalcVec(dr[1], r3, r1);
  CalcVec(dr[2], r4, r1);
  V3CROSS(crosProd, dr[1], dr[2]);
  dotProd = V3DOT(dr[0], crosProd);
  vol = fabs(dotProd) / 6.;
  return vol;
}

void DistributeForceOnTria(double (*rPnt)[NDIM], double **fPnt, 
		double *ratio, double *fi) {
  int k, ind[2];
  double denom, frac, dr[2][NDIM], crsProd[NDIM], vec[3][NDIM], len[3];
  
  CalcVec(dr[0], rPnt[1], rPnt[0]);
  CalcVec(dr[1], rPnt[2], rPnt[0]);
  FOR_NDIM(k) {
	vec[0][k] = -1. * ratio[0] * rPnt[1][k] - ratio[1] * rPnt[2][k];
	vec[1][k] = -1. * rPnt[0][k] + (1. - ratio[0]) * rPnt[1][k] 
		- ratio[1] * rPnt[2][k];
	vec[2][k] = -1. * rPnt[0][k] - ratio[0] * rPnt[1][k] 
		+ (1. - ratio[1]) * rPnt[2][k];
  }
  for(k = 0; k < 3; k++) {
	ind[0] = (k + 1) % 3;
	ind[1] = (k + 2) % 3;
	V3CROSS(crsProd, vec[ind[0]], vec[ind[1]]);
	len[k] = V3LEN(crsProd);
  }
  denom = V3SUM(len);
  for(k = 0; k < 3; k++) {
	frac = len[k] / denom;
	VS3COPY(fPnt[k], fi, frac);
  } 
}

// Calculate the difference of two vectors and apply boundary condition.
void CalcVec(double *dr, double *r1, double *r2) {
  V3SUB(dr, r1, r2);
  ApplyBoundCondVecDiff(dr);
}

// Calculate the difference of two vectors and apply boundary condition.
// Also, calculate the norm of the difference vector. 
// If mode is 0, return the norm. If mode is 1, return the square of the norm.
double CalcVecDist(double *dr, double *r1, double *r2, int mode) {
  double dist;
  CalcVec(dr, r1, r2);
  dist = (mode == 0) ? V3LEN(dr) : V3LEN_SQ(dr); 
  return dist;
}

// Calculate the norm of difference of two vectors with periodic boundary 
// condition.
// If mode is 0, return the norm. If mode is 1, return the square of the norm.
double CalcDist(double *r1, double *r2, int mode) {
  double dr[NDIM], dist;
  dist = CalcVecDist(dr, r1, r2, mode);
  return dist;
}

// Calculate the norm of difference of two vectors with periodic boundary 
// condition. Make it a unit vector too.
double CalcUnitVec(double *dr, double *r1, double *r2) {
  double len;
  CalcVec(dr, r1, r2);
  len = V3LEN(dr);
  V3SCALE(dr, INV(len));
  return len;
}

// Normalize a vector
void NormVec(double *dr) {
  double len;
  len = V3LEN(dr);
  // Check whether it is a zero vector
  if (len > 0) { 
	V3SCALE(dr, INV(len));
  }
}

double CalcDistActinAbp(int actInd, int abpInd, int mode) {
  double dr[NDIM], len;

  CalcVecActinAbp(dr, actInd, abpInd, mode);
  len = V3LEN(dr);
  return len;
}

double CalcVecDistActinAbp(double *dr, int actInd, int abpInd, 
		int mode) {
  double len;

  CalcVecActinAbp(dr, actInd, abpInd, mode);
  len = V3LEN(dr);
  return len;
}

void CalcUnitVecActinAbp(double *dr, int actInd, int abpInd,
		int mode) {
  double invLen;
 
  CalcVecActinAbp(dr, actInd, abpInd, mode);
  invLen = INV(V3LEN(dr));
  V3SCALE(dr, invLen);
}

double CalcVecActinAbp(double *dr, int actInd, int abpInd, int mode) {
  int side;
  double rPos[NDIM], *rActPnt, *rAdjActPnt, *rAbpPnt, ratio;

  if (mode == 0) { 
	side = FindAbpActinChain(iAct[actInd], abpInd, 0);
	rActPnt = &P2(act.r,iAct[actInd],0);
	rAdjActPnt = &P2(act.r,iAct[P2A(act.ch,iAct[actInd],0,nChAc)],0);
	rAbpPnt = &P2(abp.r,iAbp[abpInd],0);
  }
  else {
	side = FindAbpActinChain(actInd, abpInd, 1);
	rActPnt = &P2(rAct,actInd,0);
	rAdjActPnt = &P2(rAct,P2A(chAct,actInd,0,nChAc),0);
	rAbpPnt = &P2(rAbp,abpInd,0);
  }
  CalcPosOnActSegSide(rActPnt, rAdjActPnt, rPos, side);
  CalcVec(dr, rPos, rAbpPnt);
  ratio = (double)((int)((side - 2) / nChAcY)) / (double)nChAcX;

  return ratio;
}

void CalcPosOnActSegSide(double *r1, double *r2, double *r3, int side) {
  double ratio;

  ratio = (double)((int)((side - 2) / nChAcY)) / (double)nChAcX;
  CalcPosOnActSeg(r1, r2, r3, ratio);
}

void CalcPosOnActSeg(double *r1, double *r2, double *r3, double ratio) {
  int k, ind;
  double rPnt[2][NDIM];

  if (ratio == 0.) {
	V3COPY(r3, r1);
  }
  else {
	V3COPY(rPnt[0], r1);
	V3COPY(rPnt[1], r2);
	FOR_NDIM(k) {
		if (pbc[k] == 1	&& fabs(rPnt[0][k] - rPnt[1][k]) > 1.5 * dimDomH[k]) {
			ind = (rPnt[0][k] < rPnt[1][k]) ? 0 : 1;
			rPnt[ind][k] += dimDom[k];
		}
		r3[k] = rPnt[0][k] * (1 - ratio) + rPnt[1][k] * ratio;
	}
	ApplyBoundCondVector(r3, -1, 0);
  }
}

void OffsetSegEndPnt(double *rPnt, double *rOri, double *rCen) {
  int k;
  FOR_NDIM(k) {
	if (pbc[k] == 1 && fabs(rCen[k] - rOri[k]) > dimDomH[k]) {
		if (rCen[k] < rOri[k]) {
			rPnt[k] = rOri[k] - dimDom[k];
		}
		else { rPnt[k] = rOri[k] + dimDom[k]; }
	}
	else { rPnt[k] = rOri[k]; }
  }
}

void CalcInactAbpPosition(double *r, double *rOri, double *rC, double *rU, 
		int kind, int side) {
  double ang, len, drDirec[NDIM], drAxis[NDIM], drPos[NDIM];

  CalcPosOnActSegSide(rC, rU, rOri, side);
  ang = genrand_real3() * 2 * PI;
  V3SET(drDirec,cos(ang),sin(ang),0.);
  V3SUB(drAxis,rU,rC);
  ApplyBoundCondVecDiff(drAxis);
  V3CROSS(drPos,drAxis,drDirec);
  len = V3LEN(drPos);
  V3SCALE(drPos,abpF.spr[kind].eq / len);
  V3ADD(r,rOri,drPos);
  ApplyBoundCondVector(r, -1, 0);
}

void CalcAbpArmEndPos(double *rPos, int actInd, int abpInd) {
  int locAbpInd;
  double mag, fac, dr[NDIM];

  locAbpInd = iAbp[abpInd];
  mag = CalcVecDistActinAbp(dr, actInd, abpInd, 0);
  V3REVSIGN(dr);
  fac = abpF.len[K_ABP(locAbpInd)].n / mag;
  VS3ADD(rPos, &P2(abp.r,locAbpInd,0), dr, fac);
}

void CalcMembUnitNormalDirec(int *cylInd, double *normDir) {
  double dr[2][NDIM], v[NDIM];

  if (dimMbNuc == 2) {
	V3SET_ALL(v, 0.);
	v[dirNormMbNuc] = 1.;
	CalcUnitVec(dr[0], &P2(memb.r,iMb[cylInd[0]],0), 
			&P2(memb.r,iMb[cylInd[1]],0));
	V3CROSS(normDir, dr[0], v);
  }
  else {
	CalcVec(dr[0],&P2(memb.r,iMb[cylInd[1]],0), &P2(memb.r,iMb[cylInd[0]],0));
	CalcVec(dr[1],&P2(memb.r,iMb[cylInd[2]],0), &P2(memb.r,iMb[cylInd[0]],0));
	V3CROSS(normDir, dr[0], dr[1]);
	NormVec(normDir);
  }
}

int FindAbpActinChain(int ind, int ele, int mode) {
  int loc, *pArr;
 
  if (mode == 0) { pArr = act.ch; }
  else { pArr = chAct; }
  loc = FindElementArray(&P2A(pArr,ind,2,nChAc), nChAc - 2, ele, 0, 1);
  if (loc > -1) { loc += 2; }
  return loc;
}

int HowManyAbpActinChain(int ind, int mode) {
  int n, *pArr, cnt;
 
  if (mode == 0) { pArr = act.ch; }
  else { pArr = chAct; }
  cnt = 0;
  for(n = 2; n < nChAc; n++) { 
	CONT(!(P2A(pArr,ind,n,nChAc) > -1));
	cnt++;
  }
  return cnt;
}

void GenRandPosSubdom(double *r) {
  int k;

  FOR_NDIM(k) {
	r[k] = P2A(bnd.r,0,k,NDIM) + (P2A(bnd.r,1,k,NDIM)
			- P2A(bnd.r,0,k,NDIM)) * genrand_real3();
  }
}

void GenRandDirecVec(double *dr) {
  double ang, z, mag;

  ang = 2. * PI * genrand_real3();
  if (dir2D < 0) {
	z = 2. * genrand_real3() - 1.;
	mag = Sqrt(1. - SQR(z));
	V3SET(dr, mag * cos(ang), mag * sin(ang), z);
  }
  else {
	dr[dir2D] = 0.;
	dr[(dir2D + 1) % NDIM] = cos(ang);
	dr[(dir2D + 2) % NDIM] = sin(ang);
  }
}

int GenRandIntIndex(int nInd) {
  int ind;

  ind = (int)(genrand_real3() * (double)nInd - 0.1); 
  if (ind < 0) { ind = 0; }
  return ind;
}

// Calculate iCell[] from position
void CalcIndMolecule(double *pos, int *ind) {
  int n, k;
  double r[NDIM];

  V3COPY(r, pos);
  if (rheoWay > 0 && bulkRheoType == 0) { ConvertRectDomainVector(r, 0); }
  FOR_NDIM(k) {
	for(n = 0; n < nCell[k]; n++) {
		if (pbc[k] == 1 || (pbc[k] != 1 && n > 0 && n < nCell[k] - 1)) {
			if (r[k] >= rGrid[k][n] && r[k] < rGrid[k][n + 1]) {
				ind[k] = n;
				break;
			}
		}
    	else if ((n == 0 && r[k] < rGrid[k][1]) 
				|| (n == nCell[k] - 1 && r[k] >= rGrid[k][nCell[k] - 1])) {
			ind[k] = n;
			break;
		}
	}
  }
}

// Calculate rank from position
int CalcRankMolecule(double *pos) {
  int ind[NDIM], calRank;
  CalcIndMolecule(pos, ind);
  V3IND_BACK_INT(calRank, ind, nCell);
  return calRank;
}

// Calculate iCell[] and rank from position
int CalcRankIndMolecule(double *pos, int *ind) {
  int calRank;
  CalcIndMolecule(pos, ind);
  V3IND_BACK_INT(calRank, ind, nCell);
  return calRank;
}

// This function accelerates or decelerates the rates of dynamics of actins 
// and ABPs depending on the phase of simulations (e.g. network formation, 
// application of prestress/prestrain) and given factors in 'condition'.
double AdjustDynamicsRate(double p) {
  double pAdj;

  pAdj = p;
  if (currTimeStep < netForm.dur) { 
	pAdj *= netForm.facK;
  }
  else if (currTimeStep >= netForm.dur 
		&& currTimeStep < netForm.dur + motActiv.dur) {
	pAdj *= motActiv.facK;
  }
  else if (currTimeStep >= netForm.dur + motActiv.dur
		&& currTimeStep < netForm.dur + pres.dur + motActiv.dur) {
	pAdj *= pres.facK;
  }
  return pAdj;
}

// Simple matrix multiplication
// matC = matA * matB
void MultiplyMatrixIntSerial(int *matA, int *matB, int *matC, 
		int nDimA, int nDimAB, int nDimB) {
  int m, n, k;
  
  for(n = 0; n < nDimA * nDimB; n++) {
	matC[n] = 0.;
  }
  for(m = 0; m < nDimA; m++) {
	for(n = 0; n < nDimB; n++) {
		for(k = 0; k < nDimAB; k++) {
			P2A(matC,m,n,nDimB) += P2A(matA,m,k,nDimAB) * P2A(matB,k,n,nDimB);
		}
	}
  }
}

void MultiplySqMatrixIntParallel(int *matA, int *matB, int *matC, int nDim) {
  int n, nDimLoc, *matCloc;

  nDimLoc  = (int)(nDim / nCpu);
  MALLOC(matCloc,int,nDim * (nDimLoc + ((rank == nCpu - 1) 
		? nDim - nDimLoc * nCpu : 0)));

  MultiplyMatrixIntSerial(matA + rank * nDimLoc * nDim, matB, matCloc, 
		nDimLoc + ((rank == nCpu - 1) ? nDim - nDimLoc * nCpu : 0), nDim, nDim);
  MPI_Allgather(matCloc, nDimLoc * nDim, MPI_INT, matC, nDimLoc * nDim, 
		MPI_INT, MPI_COMM_WORLD);
  if (nDim % nCpu != 0) {
	if (rank == nCpu - 1) {
		for(n = 0; n < (nDim % nCpu) * nDim; n++) {
			matC[nDimLoc * nCpu * nDim + n] = matCloc[nDimLoc * nDim + n];
		}
	}
	MPI_Bcast(matC + nDimLoc * nCpu * nDim, (nDim % nCpu) * nDim, 
			MPI_INT, nCpu - 1, MPI_COMM_WORLD);
  }
  free(matCloc);
}

// Convert a file name to store in a data folder.
char *GenFileName(const char *fnIn) {
  sprintf(fnOut, "%s/%s", dataFold, fnIn);
  return(fnOut);
}


// Swap two integer numbers.
void SwapInt(int *a, int *b) {
  int temp;
  temp = *a;
  *a = *b;
  *b = temp;
}

// Swap two double numbers.
void SwapDbl(double *a, double *b) {
  double temp;
  temp = *a;
  *a = *b;
  *b = temp;
}

// Return upper limit if a value is larger than the upper limit.
// Return lower limit if a value is smaller than the lower limit.
double TrimDblVal(double val, double loBnd, double hiBnd) {
  if (val < loBnd) { val = loBnd; }
  else if (val > hiBnd) { val = hiBnd; }
  return val;
}

int TrimIntVal(int val, int loBnd, int hiBnd) {
  if (val < loBnd) { val = loBnd; }
  else if (val > hiBnd) { val = hiBnd; }
  return val;
}

int SetKind(int ind, int n1, int n2) {
  int kind;

  if (ind < n1) { kind = 0; }
  else if (ind >= n1 && ind < n1 + n2) { kind = 1; }
  else { kind = 2; }
  return kind;
}

// Calculate acos(input) using precalculated values to save computational cost.
double Acos(double val) {
  val = TrimDblVal(val, -1., 1.);
  return arrAcos[(int)( (val + 1.) * (double)DEG_ARR_ACOS)];
}

double Sqrt(double val) {
  if (val < 0.) { val = 0.; }
  return sqrt(val);
}

// Compare two integer numbers.
int CompInt(const void *arg1, const void *arg2) {
  return *(int *)arg1 - *(int *)arg2;
}

// Compare two double numbers.
int CompDbl(const void *arg1, const void *arg2) {
  return *(double *)arg1 - *(double *)arg2;
}

// Return 1 if input is positive, -1 if negative, 0 is zero.
int SignDbl(double val) {
  int sign;
  if (val > 0.) { sign = 1; }
  else if (val < 0.) { sign = -1; }
  else { sign = 0; }
  return sign;
}

// Convert an integer value to string.
char* IntToStr(int val, int base) {
  static char buf[32] = {0};
  int n = 30;
  for(; val && n ; --n, val /= base)
	buf[n] = "0123456789abcdef"[val % base];
  return &buf[n+1];
} 
