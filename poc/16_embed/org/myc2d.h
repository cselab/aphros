#define NOT_ZERO 1.e-30

/*-----------------------------------------------------*
 *MYC - Mixed Youngs and Central Scheme (2D)           *
 *-----------------------------------------------------*/
coord mycs (Point point, scalar c)
{
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;
  
  /* top, bottom, right and left sums of c values */
  c_t = c[-1,1] + c[0,1] + c[1,1];
  c_b = c[-1,-1] + c[0,-1] + c[1,-1];
  c_r = c[1,-1] + c[1,0] + c[1,1];
  c_l = c[-1,-1] + c[-1,0] + c[-1,1];

  /* consider two lines: sgn(my) Y =  mx0 X + alpha,
     and: sgn(mx) X =  my0 Y + alpha */ 
  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);

  /* minimum coefficient between mx0 and my0 wins */
  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }

  /* Youngs' normal to the interface */
  mm1 = c[-1,-1] + 2.0*c[-1,0] + c[-1,1];
  mm2 = c[1,-1] + 2.0*c[1,0] + c[1,1];
  mx1 = mm1 - mm2 + NOT_ZERO;
  mm1 = c[-1,-1] + 2.0*c[0,-1] + c[1,-1];
  mm2 = c[-1,1] + 2.0*c[0,1] + c[1,1];
  my1 = mm1 - mm2 + NOT_ZERO;

  /* choose between the best central and Youngs' scheme */ 
  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
	
  /* normalize the set (mx0,my0): |mx0|+|my0|=1 and
     write the two components of the normal vector  */
  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};
  
  return n;
}
