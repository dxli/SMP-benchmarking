/*
p 200MHz, 1.09 
1.6 P-M, 15.1 
p4 2.8 pc 533      18.1
p4 3.0 pc 400 ddr; 19.2
p4 3.06-m 333 ddr 16.1
14.8892 P/S on the Barton 2500 on ddr333

duron 1394.066 MHz, pc2100, 8.6 P/s
athlon64 3000+(2.0GHz) ddr 400, 17.6 P/S
Speed= 16.016 P/S with the Barton 2500 now
<marienz> my xp 2600+ is doing reasonably well with a speed of 14.3 approximately then
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#define dsteps 25.
#define v_steplength 0.0014
#define rc 30.
#define rv 31.
#define gtemp 175.
#define psteps1 2
#define psteps2 6
#define intersteps 4.
#define gp2d 330.
#define rrange 2.5
typedef struct neighbour neighbour, *nblink;
typedef struct neighbours neighbours;
typedef struct
{
  double re;
  double im;
} complex;
struct neighbour
{
  int ni;
  nblink np;
  nblink nadd;
  short flag;
  double arg;
  double rsquared;
};
struct neighbours
{
  nblink np;
  double arg;
};
#define upotential(ppointer) pti=(int) dx;dx-=pti;dy=(ppointer[pti]+(ppointer[pti+1]-ppointer[pti])*dx);
void newvlist (double *), df_pattern (void), df_number (int, int *);
void newproperties (double *), writeproperties (double *);
double new1_energy (double *, int), new1_exy (double *, int), gtpth (double *,
								     double,
								     int);
double *MonteCarlo (double *);
neighbour *buflist, *pnbtv2, *pnbtv3, **pnbv1, **pnbv2;
neighbours *ncp;
complex *g6list, *g6buf;
unsigned defect5, defect7, *g2list, g2steps = 0, *zdist, g2length;
int j2, l, k, kp1, lp1, n, number1, number2, *vlist1, df_save[7][3], *dflag,
  nc[3], MCsteps =
  50000, **pv1, **pv2, *nlist, **bufcells, **bufcells2, **nbufcells,
  **nbufcells2, **pcells2, **lncells, **lncells2, ***cells, ***ncells,
  ****plncells, ***buflncells, ***buflncells2, *cellsbuf, *cellsbuf2,
  lncellsend = -1, gs0;
double opsq = 0., op = 0., tpsq = 0., tp =
  0., lx, ly, halflx, halfly, c, lx_new, ly_new, halflx_new, halfly_new, d =
  0.0276452, dl, r22s, deltar, v_origin, *pt_c1, *pt_c2, g_pressure,
  temperature =
  200., *oldc, *pne, **pesave, **pesave2, *uij, **peaddb, ***peadd, ***peadd2,
  *pvsave, *pt1, *pt2, *pt3, *pt9, *vz, rvs, tkpb, *pcoord1, *pcoord2, xi, yi,
  zi, q0, vznew, nrx, nry, upot, usquared, xis0, xis1, xis2, *sxy, gth0 =
  0., gth =
  0., oprun, tprun, gtpth (double *p, double gth, int n1),
trans_op (double *p, int n1);
char *filename1, random_buf[256 + 1], fn5[35];
FILE *f1, *f2;
struct timeval time0, time1;
int
main (int argc, char *argv[])
{
  int i, j;
  int ptis;
  double r, x, y, dx0, dy0, dx1, dy1;
  static char string1[35] = "c2.dat";
/* Initilize the random generator*/

  filename1 = random_buf;
  initstate (1, random_buf, 256);
      for (i = 0; i < 256; i++)
	*filename1++ = (char) i;
/* */
  temperature = gtemp;
  number1 = 4012;
  lx = 380.596;
  ly = 381.309;
  g_pressure = gp2d / temperature;
  pt2 = (double *) malloc (10002 * sizeof (double));
  pt9 = pt2;
  i = 0;
  while (i < 10000)
    {
      y = x = 37.3 / (0.25 * i);
      y *= x * x;
      *(pt9++) = 4450. * y * (y - 1.) / temperature;
      i++;
    }
  ptis = (int) (1.5 * 3.14159 * rv * rv * number1 / lx / ly) + 10;
  halflx = lx / 2.;
  halfly = ly / 2.;
  n = number1 + number2;
  uij = ((double *) malloc (ptis * (number1 + 1) * sizeof (double)));
  peaddb = (double **) malloc (ptis * (number1 + 1) * sizeof (double *));
  vlist1 = (int *) malloc (ptis * (number1 + 1) * sizeof (int));
  pcoord1 = (double *) malloc (2 * (number1 + 1) * sizeof (double));
  pcoord2 = (double *) malloc (2 * (number1 + 1) * sizeof (double));
  pt9 = pcoord1;
  v_origin = 380.596 * 381.309;
  r = 1. / sqrt (sqrt (3.) / 2 * number1 / (lx * ly));
  dx0 = r * sqrt (3.);
  dy0 = r;
  dx1 = dx0 / 2.;
  dy1 = dy0 / 2.;
  x = 0.01;
  for (i = 0; i++ < 34; x += dx0)
    {
      y = 0.01;
      for (j = 0; j++ < 59; y += dy0)
	{
	  *pt9++ = x;
	  *pt9++ = y;
	  *pt9++ = x + dx1;
	  *pt9++ = y + dy1;
	}
    }
  vz = (double *) malloc ((number1 + 1) * sizeof (double));
  pv1 = (int **) malloc ((number1 + 1) * sizeof (int *));
  pv2 = (int **) malloc ((number1 + 1) * sizeof (int *));
  dflag = (int *) malloc ((number1 + 1) * sizeof (int));
  pesave2 = (double **) malloc ((number1 + 1) * sizeof (double *));
  pesave = (double **) malloc ((number1 + 1) * sizeof (double *));
  peadd = (double ***) malloc ((number1 + 1) * sizeof (double **));
  peadd2 = (double ***) malloc ((number1 + 1) * sizeof (double **));
  pne = (double *) malloc ((ptis + 1) * sizeof (double));
  buflist = (neighbour *) malloc (ptis * (number1 + 1) * sizeof (neighbour));
  pnbv1 = (neighbour * *)malloc ((number1 + 1) * sizeof (neighbour *));
  pnbv2 = (neighbour * *)malloc ((number1 + 1) * sizeof (neighbour *));
  g6buf = (complex *) malloc ((number1 + 1) * sizeof (complex));
  ncp = (neighbours *) malloc (ptis * sizeof (neighbours));
  nlist = (int *) malloc ((number1 + 1) * sizeof (int));
  g2steps = k = i = 0;
  gs0 = i;
  j = 0;
  for (i = 0; i < number1; i++)
    {
      *(pv1 + i) = vlist1 + j;
      *(pesave + i) = uij + j;
      *(peadd + i) = peaddb + j;
      pnbv1[i] = buflist + j;
      j += ptis;
    }
  k = (int) (lx / rv);
  l = (int) (ly / rv);
  kp1 = k;
  lp1 = l;
  nrx = k / lx;
  nry = l / ly;
  j2 = (int) (number1 / (0.6 * k * l)) + 10;
  k++;
  l++;
  bufcells = (int **) malloc (k * l * sizeof (int *));
  cells = (int ***) malloc (l * sizeof (int **));
  nbufcells = (int **) malloc (k * l * sizeof (int *));
  ncells = (int ***) malloc (l * sizeof (int **));
  lncells = (int **) malloc (k * l * 9 * sizeof (int *));
  plncells = (int ****) malloc (l * sizeof (int ***));
  buflncells = (int ***) malloc (k * l * sizeof (int **));
  cellsbuf = (int *) malloc (j2 * k * l * sizeof (int));
  cellsbuf2 = cellsbuf;
  bufcells2 = bufcells;
  nbufcells2 = nbufcells;
  lncells2 = lncells;
  buflncells2 = buflncells;
  k--;
  l--;
  k--;
  l--;
  for (j = 0; j <= l; j++)
    {
      cells[j] = bufcells2;
      ncells[j] = nbufcells2;
      bufcells2 += (k + 1);
      nbufcells2 += (k + 1);
      for (i = 0; i <= k; i++)
	{
	  cells[j][i] = cellsbuf2;
	  cellsbuf2 += j2;
	}
    }
  for (j = 0; j <= l; j++)
    {
      plncells[j] = buflncells2;
      buflncells2 += (k + 1);
      for (i = 0; i <= k; i++)
	{
	  plncells[j][i] = lncells2;
	  if (j == 0)
	    {
	      *(lncells2++) = cells[l][i];
	      *(lncells2++) = cells[j + 1][i];
	      if (i < k)
		{
		  *(lncells2++) = cells[l][i + 1];
		  *(lncells2++) = cells[j][i + 1];
		  *(lncells2++) = cells[j + 1][i + 1];
		  if (i > 0)
		    {
		      *(lncells2++) = cells[j + 1][i - 1];
		      *(lncells2++) = cells[l][i - 1];
		    }
		  else
		    {
		      *(lncells2++) = cells[l][k];
		      *(lncells2++) = cells[j][k];
		      *(lncells2++) = cells[1][k];
		    }
		}
	      else
		{
		  *(lncells2++) = cells[j + 1][i - 1];
		  *(lncells2++) = cells[l][i - 1];
		  *(lncells2++) = cells[1][0];
		  *(lncells2++) = cells[l][0];
		}

	    }
	  else
	    {
	      if (j < l)
		{
		  *(lncells2++) = cells[j + 1][i];
		  if (i < k)
		    {
		      *(lncells2++) = cells[j][i + 1];
		      *(lncells2++) = cells[j + 1][i + 1];
		      if (i > 0)
			{
			  *(lncells2++) = cells[j + 1][i - 1];
			}
		      else
			{
			  *(lncells2++) = cells[j][k];
			  *(lncells2++) = cells[j + 1][k];
			}
		    }
		  else
		    {
		      *(lncells2++) = cells[j + 1][i - 1];
		      *(lncells2++) = cells[j + 1][0];
		    }
		}
	      else
		{
		  if (i < k)
		    {
		      *(lncells2++) = cells[j][i + 1];
		      if (i == 0)
			*(lncells2++) = cells[j][k];
		    }
		}
	    }
	  *(lncells2++) = NULL;
	}
    }
  c = pow (2., sizeof (int) * 8 - 1);
  deltar = rv - rc;
  deltar = deltar * deltar;
  r22s = rv * rv;
  rvs = rrange * sqrt (lx * ly / number1 * 2. / sqrt (3.));
  rvs *= rvs;
  MonteCarlo (pcoord1);
  return (0);
}


void
newvlist (double *p)
{
  double yi0;
  register double *pt7, dx, dy, xi0;
  register int ii, pti, *nptc;
  int kk, ll, mm, ll3, mm3, jj, *nptc2, **lncells3;
  lncells2 = nbufcells;
  lncells3 = bufcells;
  jj = kp1 * lp1;
  ii = 0;
  while (ii++ < jj)
    *(lncells2++) = *lncells3++;
  kk = 0;
  pt7 = p;
  for (ii = 0; ii < number1; ii++)
    {
      dx = *pt7++;
      dy = *pt7++;
      *(ncells[((int) (dy * nry)) % lp1][((int) (dx * nrx)) % kp1]++) = ii;
      peadd2[ii] = peadd[ii];
      pesave2[ii] = pesave[ii];
      pv2[ii] = pv1[ii];
    }
  lncells2 = nbufcells;
  ii = 0;
  kk = ii;
  while (ii++ < jj)
    *(*lncells2++) = -1;
  while (kk <= l)
    {
      ii = 0;
      while (ii <= k)
	{
	  nptc2 = cells[kk][ii];
	  while (*(nptc = nptc2++) >= 0)
	    {
	      ll = *nptc;
	      ll3 = ll << 1;
	      pt7 = p + ll3;
	      xi0 = *pt7++;
	      yi0 = *pt7;
	      while ((mm = *(++nptc)) >= 0)
		{
		  mm3 = mm << 1;
		  pt7 = p + mm3;
		  dx = halflx - fabs (halflx - fabs (xi0 - *pt7++));
		  dx *= dx;
		  dy = halfly - fabs (halfly - fabs (yi0 - *pt7));
		  dx += dy * dy;
		  if (dx <= r22s)
		    {
		      *(pv2[ll]++) = mm3;
		      *(pv2[mm]++) = ll3;
		      dx *= intersteps;
		      upotential (pt2);
		      *(*(peadd2[ll]++) = (pesave2[mm]++)) = dy;
		      *(*(peadd2[mm]++) = (pesave2[ll]++)) = dy;
		    }
		}
	      lncells2 = plncells[kk][ii];
	      while (*lncells2 != NULL)
		{
		  nptc = *lncells2;
		  while ((mm = *nptc++) >= 0)
		    {
		      mm3 = mm << 1;
		      pt7 = p + mm3;
		      dx = halflx - fabs (halflx - fabs (xi0 - *pt7++));
		      dx *= dx;
		      dy = halfly - fabs (halfly - fabs (yi0 - *pt7));
		      dx += dy * dy;
		      if (dx <= r22s)
			{
			  *(pv2[ll]++) = mm3;
			  *(pv2[mm]++) = ll3;
			  dx *= intersteps;
			  upotential (pt2);
			  *(*(peadd2[ll]++) = (pesave2[mm]++)) = dy;
			  *(*(peadd2[mm]++) = (pesave2[ll]++)) = dy;
			}
		    }
		  lncells2++;
		}
	    }
	  ii++;
	}
      kk++;
    }
  ii = 0;
  while (ii < number1)
    *pv2[ii++] = -1;
}


void
newproperties (double *p)
{
  int kk, ll, mm, ll3, mm3, jj, *nptc2, **lncells3, ir, level, nlength;
  neighbours *pnbtv1, rra;
  double *pt8, yi0, xi1, yi1, xi2, yi2, index;
  double opre0, opim0;
  register double *pt7, dx, dy, xi0;
  register int ii, pti, *nptc;
  defect5 = defect7 = 0;
  lncells2 = nbufcells;
  lncells3 = bufcells;

  jj = lp1 * kp1;
  ii = 0;
  while (ii++ < jj)
    *(lncells2++) = *lncells3++;
  kk = 0;
  pt7 = p;
  for (ii = 0; ii < number1; ii++)
    {
      dx = *pt7++;
      dy = *pt7++;
      *(ncells[((int) (dy * nry)) % lp1][((int) (dx * nrx)) % kp1]++) = ii;
      peadd2[ii] = peadd[ii];
      pesave2[ii] = pesave[ii];
      pv2[ii] = pv1[ii];
      pnbv2[ii] = pnbv1[ii];
    }
  lncells2 = nbufcells;
  ii = 0;
  kk = ii;
  while (ii++ < jj)
    *(*lncells2++) = -1;
  while (kk <= l)
    {
      ii = 0;
      while (ii <= k)
	{
	  nptc2 = cells[kk][ii];
	  while (*(nptc = nptc2++) >= 0)
	    {
	      ll = *nptc;
	      ll3 = ll << 1;
	      pt7 = p + ll3;
	      xi0 = *pt7++;
	      yi0 = *pt7;
	      while ((mm = *(++nptc)) >= 0)
		{
		  mm3 = mm << 1;
		  pt7 = p + mm3;
		  dx = halflx - fabs (halflx - fabs (xi0 - *pt7++));
		  dy = halfly - fabs (halfly - fabs (yi0 - *pt7));
		  dx = dy * dy + dx * dx;
		  if (dx <= r22s)
		    {
		      *(pv2[ll]++) = mm3;
		      *(pv2[mm]++) = ll3;
		      if (dx <= rvs)
			{
			  (((pnbv2[ll]->nadd) = pnbv2[mm])->ni) = ll3;
			  (((pnbv2[mm]->nadd) = pnbv2[ll])->ni) = mm3;
			  pnbv2[ll]->rsquared = (pnbv2[mm]->rsquared = dx);
			  (pnbv2[mm]++)->flag = ((pnbv2[ll]++)->flag = 0);
			}
		      dx *= intersteps;
		      upotential (pt2);
		      *(*(peadd2[ll]++) = (pesave2[mm]++)) = dy;
		      *(*(peadd2[mm]++) = (pesave2[ll]++)) = dy;
		    }
		}
	      lncells2 = plncells[kk][ii];
	      while (*lncells2 != NULL)
		{
		  nptc = *lncells2;
		  while ((mm = *nptc++) >= 0)
		    {
		      mm3 = mm << 1;
		      pt7 = p + mm3;
		      dx = halflx - fabs (halflx - fabs (xi0 - *pt7++));
		      dy = halfly - fabs (halfly - fabs (yi0 - *pt7));
		      dx = dy * dy + dx * dx;
		      if (dx <= r22s)
			{
			  *(pv2[ll]++) = mm3;
			  *(pv2[mm]++) = ll3;
			  if (dy <= rvs)
			    {
			      (((pnbv2[ll]->nadd) = pnbv2[mm])->ni) = ll3;
			      (((pnbv2[mm]->nadd) = pnbv2[ll])->ni) = mm3;
			      pnbv2[ll]->rsquared = (pnbv2[mm]->rsquared =
						     dx);
			      (pnbv2[mm]++)->flag = ((pnbv2[ll]++)->flag = 0);
			    }
			  dx *= intersteps;
			  upotential (pt2);
			  *(*(peadd2[ll]++) = (pesave2[mm]++)) = dy;
			  *(*(peadd2[mm]++) = (pesave2[ll]++)) = dy;
			}
		    }
		  lncells2++;
		}
	    }
	  ii++;
	}
      kk++;
    }
  /* start the newproperties construction */
/* sorting for each atom*/
  for (ii = 0; ii < number1; ii++)
    {
      pnbv2[ii]->ni = -1;
      pnbtv3 = pnbv1[ii];
      pnbtv1 = ncp;
      nlength = -1;
      jj = ii << 1;
      xi0 = p[jj++];
      yi0 = p[jj];
      while ((jj = pnbtv3->ni) >= 0)
	{
	  if ((pnbtv3->flag) >= 0)
	    {
	      nlength++;
	      pnbtv1->np = pnbtv3;
	      dx = xi0 - p[jj++];
	      dx -= ((int) (dx / halflx)) * lx;
	      dy = yi0 - p[jj];
	      dy -= ((int) (dy / halfly)) * ly;
	      (pnbtv3->arg) = (pnbtv1++)->arg = atan2 (dy, dx);
	    }
	  pnbtv3++;
	}
      if (nlength <= 0)
	continue;
      level = (nlength >> 1) + 1;
      ir = nlength;
      for (;;)
	{
	  if (level > 0)
	    {
	      rra = ncp[--level];
	    }
	  else
	    {
	      rra = ncp[ir];
	      ncp[ir] = *ncp;
	      if (--ir == 0)
		{
		  *ncp = rra;
		  break;
		}
	    }
	  pti = level;
	  jj = pti << 1;
	  jj++;
	  while (jj <= ir)
	    {
	      if (jj < ir && ncp[jj].arg < ncp[jj + 1].arg)
		jj++;
	      if (rra.arg < ncp[jj].arg)
		{
		  ncp[pti] = ncp[jj];
		  pti = jj;
		  jj <<= 1;
		  jj++;
		}
	      else
		jj = ir + 1;
	    }
	  ncp[pti] = rra;
	}
      pnbtv3 = ncp->np;
      jj = 0;
      pti = nlength;
      while (jj++ < pti)
	{
	  pnbtv3 = (pnbtv3->np = ncp[jj].np);
	}
      pnbtv3->np = (pnbtv2 = ncp->np);


/* newproperties Construction */
      dx = p[pnbtv3->ni] - xi0;
      if (dx > halflx)
	{
	  dx -= lx;
	}
      else
	{
	  if (dx + halflx < 0)
	    dx += lx;
	}
      dy = p[pnbtv3->ni + 1] - yi0;
      if (dy > halfly)
	{
	  dy -= ly;
	}
      else
	{
	  if (dy + halfly < 0)
	    dy += ly;
	}
      xi1 = p[pnbtv2->ni] - xi0;
      if (xi1 > halflx)
	{
	  xi1 -= lx;
	}
      else
	{
	  if (xi1 + halflx < 0)
	    xi1 += lx;
	}
      yi1 = p[pnbtv2->ni + 1] - yi0;
      if (yi1 > halfly)
	{
	  yi1 -= ly;
	}
      else
	{
	  if (yi1 + halfly < 0)
	    yi1 += ly;
	}
      pnbtv2 = pnbtv2->np;
      jj = 0;
      do
	{
	  xi2 = p[pnbtv2->ni] - xi0;
	  if (xi2 > halflx)
	    {
	      xi2 -= lx;
	    }
	  else
	    {
	      if (xi2 + halflx < 0)
		xi2 += lx;
	    }
	  yi2 = p[pnbtv2->ni + 1] - yi0;
	  if (yi2 > halflx)
	    {
	      yi2 -= ly;
	    }
	  else
	    {
	      if (yi2 + halflx < 0)
		yi2 += ly;
	    }
	  if (pnbtv3->np->flag < 1)
	    {
	      index =
		(pnbtv3->rsquared * (yi2 * xi1 - yi1 * xi2) +
		 pnbtv2->rsquared * (dx * yi1 - dy * xi1)) / (dx * yi2 -
							      dy * xi2);
	      if (index > 0 && index <= pnbtv3->np->rsquared)
		{
		  pnbtv3->np->nadd->flag = -1;
		  pnbtv3->np = pnbtv2;
		  jj = 0;
		  pti--;
		}
	      else
		{
		  pnbtv3 = pnbtv3->np;
		  jj++;
		  dx = xi1;
		  dy = yi1;
		}
	    }
	  else
	    {
	      pnbtv3 = pnbtv3->np;
	      jj++;
	      dx = xi1;
	      dy = yi1;
	    }
	  pnbtv2 = pnbtv2->np;
	  xi1 = xi2;
	  yi1 = yi2;
	}
      while (jj <= pti);
      for (jj = 0; jj <= pti; jj++)
	{
	  pnbtv3->nadd->flag = 1;
	  pnbtv3 = pnbtv3->np;
	}
      pnbv2[ii] = pnbtv3;
      pti++;

/* identify defects */
      switch (pti)
	{
	case 5:
	  defect5++;
	  break;
	case 7:
	  defect7++;
	  break;
	  /*      default: */
	}

      nlist[ii] = pti;
    }

  opre0 = opim0 = 0.;
  for (ii = 0; ii < number1; *(pv2[ii++]) = -2)
    {
      g6buf[ii].re = cos (dx = 6 * (pnbtv2 = pnbv2[ii])->arg);
      g6buf[ii].im = sin (dx);
      pti = nlist[ii];
      for (jj = 1; jj < pti; jj++)
	{
	  g6buf[ii].re += cos (dx = 6 * (pnbtv2 = pnbtv2->np)->arg);
	  g6buf[ii].im += sin (dx);
	}
      opre0 += (g6buf[ii].re /= pti);
      opim0 += (g6buf[ii].im /= pti);
    }				/* 
				   putchar('\n');  */
  oprun = sqrt (opre0 * opre0 + opim0 * opim0) / number1;
  pt7 = p;
  for (ii = 0; ii < number1; ii++)
    {
      xi0 = *pt7++;
      yi0 = *pt7++;
      pt8 = pt7;
      dx = g6buf[ii].re * g6buf[ii].re;
      dx += (g6buf[ii].im * g6buf[ii].im);
      for (jj = ii + 1; jj < number1; jj++)
	{
	  dx = halflx - fabs (halflx - fabs (xi0 - *pt8++));
	  dy = halfly - fabs (halfly - fabs (yi0 - *pt8++));
	  pti = (int) (sqrt (dx * dx + dy * dy) / 0.02);
	}
    }
}


double
new1_energy (double *p, int nptmp)
{
  register int *pv3, pti;
  register double *pt8, *pt8z, dx, dy, a = 0.;
  pv3 = pv1[nptmp];
  pt8z = pne;
  while (*pv3 >= 0)
    {
      pt8 = p + *pv3++;
      dx = halflx - fabs (halflx - fabs (xi - *pt8++));
      dx *= dx;
      dy = halfly - fabs (halfly - fabs (yi - *pt8));
      dx = (dx + dy * dy) * intersteps;
      upotential (pt2);
      a += (*pt8z++ = dy);
    }
  return (a);
}




double *
MonteCarlo (double *p)
{
  int j, *pv3;
  double fx, fy, *pt5, *p1c;
  int i, kj, n1_1, count_v1 = 0, count_v2 = 0, i_count = 0;
  double acc = 0., dl100, dl_lv, *ps1, **pt6, pe_tmp, dlnv, dlf;
  struct rusage t_start, t_end;
  /* initial parameters */
  c = pow (2., 8 * sizeof (int) - 1);
  dl100 = sqrt (lx * ly / number1) / dsteps / c;
  dl_lv = v_steplength / c;
  n1_1 = number1 + 1;
  i = 0;
  /* Monte Carlo */
  pt_c1 = pcoord1;
  pt_c2 = pcoord2;
  newvlist (pt_c1);
  j = 0;
  getrusage (RUSAGE_SELF, &t_start);
  while (1)
    {
      kj = 0;
      while (kj++ < number1)
	{
	  j = random () % n1_1;
	  if (j == number1)
	    {
	      fx = 0.;
	      j = 0;
	      while (j < number1)
		{
		  pt5 = pesave[j];
		  pv3 = pv1[j++];
		  while (*pv3++ >= 0)
		    {
		      fx += *(pt5++);
		    }
		}
	      pe_tmp = fx;
/*lnV change*/
	      dlf = exp (dlnv = ((random () << 1) * dl_lv));
	      ly_new = ly;
	      lx_new = lx;
	      p1c = pt_c1;
	      ps1 = pt_c2;
	      if ((random () & 0x1) || (lx > 1.08 * ly && dlnv > 0.)
		  || (lx < 0.92 * ly && dlnv < 0.))
		{
		  fy = dlf;
		  ly *= fy;
		  nry = lp1 / ly;
		  for (j = 0; j < number1; j++)
		    {
		      *ps1++ = *p1c++;
		      *ps1++ = fy * (*p1c++);
		    }
		}
	      else
		{
		  fx = dlf;
		  lx *= fx;
		  nrx = kp1 / lx;
		  for (j = 0; j < number1; j++)
		    {
		      *ps1++ = fx * (*p1c++);
		      *ps1++ = *p1c++;
		    }
		}
	      halfly = 0.5 * ly;
	      halflx = 0.5 * lx;
/*Potential after the lnV change*/
	      newvlist (pt_c2);
	      j = 0;
	      fx = j;
	      while (j < number1)
		{
		  pt5 = pesave[j];
		  pv3 = pv1[j++];
		  while (*pv3++ >= 0)
		    {
		      fx += *(pt5++);
		    }
		}
	      fy =
		n1_1 * dlnv - (fx - pe_tmp) * 0.5 - g_pressure * (dlf -
								  1) *
		(lx_new * ly_new);
	      count_v2++;
	      if (fy > 0 || (random () < c * exp (fy)))
		{
		  count_v1++;
		  p1c = pt_c2;
		  pt_c2 = pt_c1;
		  pt_c1 = p1c;
		}
	      else
		{
		  lx = lx_new;
		  ly = ly_new;
		  halflx = 0.5 * lx;
		  halfly = 0.5 * ly;
		  nrx = kp1 / lx;
		  nry = lp1 / ly;
		  newvlist (pt_c1);
		}
	      continue;
	    }
	  p1c = pt_c1 + (j << 1);
	  fx = *(p1c++) + ((double) (random () << 1)) * dl100;
	  xi=fx=fmod(fx+lx,lx);
	  fy = *p1c + ((double) (random () << 1)) * dl100;
	  yi = fmod(fy+ly,ly);
	  fy = 0.;
	  pt5 = pesave[j];
	  pv3 = pv1[j];
	  while (*pv3++ >= 0)
	    {
	      fy += *(pt5++);
	    }
	  if (random () < c * exp (fy - new1_energy (pt_c1, j)))
	    {
	      acc++;
	      *(p1c--) = yi;
	      *p1c = xi;
	      pt5 = pesave[j];
	      pt6 = peadd[j];
	      p1c = pne;
	      pv3 = pv1[j];
	      while (*pv3++ >= 0)
		{
		  *(*pt6++) = (*pt5++ = *p1c++);
		}
	    }
	}
      i++;
      if ((i >> (psteps2)) > i_count)
	{
	  i_count++;
	  if ((i_count >> (psteps1)) << (psteps1) == i_count)
	    {
	      newproperties (pt_c1);
	      getrusage (RUSAGE_SELF, &t_end);
	      fx =
		t_end.ru_utime.tv_sec - t_start.ru_utime.tv_sec +
		1.e-6 * (t_end.ru_utime.tv_usec - t_start.ru_utime.tv_usec);
              exit(0);
	      printf ("%d\t%g s\tSpeed= %g P/S\n", i, fx, i / fx);
	    }
	  else
	    {
	      getrusage (RUSAGE_SELF, &t_end);
	      fx =
		t_end.ru_utime.tv_sec - t_start.ru_utime.tv_sec +
		1.e-6 * (t_end.ru_utime.tv_usec - t_start.ru_utime.tv_usec);

	    }
	}
    }
  return (pt_c1);
}
