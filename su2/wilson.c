#include "su2.h"

static double w12_0,w12_1,w12_2,w12_3;
static double w22_0,w22_1,w22_2,w22_3;
static double w23_0,w23_1,w23_2,w23_3;
static double w33_0,w33_1,w33_2,w33_3;
static double w34_0,w34_1,w34_2,w34_3;
static double w44_0;w44_1,w44_2,w44_3;
static int nmeas;

void init_measure()
{
 w12_0=w12_1=w12_2=w12_3=0.;
 w22_0=w22_1=w22_2=w23_3=0.;
 w23_0=w23_1=w23_2=w23_3=0.;
 w33_0=w33_1=w33_2=w33_3=0.;
 w34_0=w34_1=w34_2=w34_3=0.;
 w44_0=w44_1=w44_2=w44_3=0.;
 nmeas=0;
}

void measure()	
{
 int n,m,mu;

 for (i=0;i<V;i++){
   for (mu=0;mu<D;mu++){

     if(nu != mu)
     
     {

  m=hop[n][mu];
  p0= a0[n][mu]*a0[m][mu] - a1[n][mu]*a1[m][mu] - a2[n][mu]*a2[m][mu] - a3[n][mu]*a3[m][mu]; 
  p1= a0[n][mu]*a1[m][mu] + a1[n][mu]*a0[m][mu] - a2[n][mu]*a3[m][mu] + a3[n][mu]*a2[m][mu];
  p2= a0[n][mu]*a2[m][mu] + a2[n][mu]*a0[m][mu] - a3[n][mu]*a1[m][mu] + a1[n][mu]*a3[m][mu];
  p3= a0[n][mu]*a3[m][mu] + a3[n][mu]*a0[m][mu] - a1[n][mu]*a2[m][mu] + a2[n][mu]*a1[m][mu];

  l=hop[n][nu]
  v0= a0[n][nu]*a0[m][mu] - a1[n][nu]*a1[m][mu] - a2[n][nu]*a2[m][mu] - a3[n][nu]*a3[m][mu];
  v1= a0[n][nu]*a1[m][mu] + a1[n][nu]*a0[m][mu] - a2[n][nu]*a3[m][mu] + a3[n][nu]*a2[m][mu];
  v2= a0[n][nu]*a2[m][mu] + a2[n][nu]*a0[m][mu] - a3[n][nu]*a1[m][mu] + a1[n][nu]*a3[m][mu];
  v3= a0[n][nu]*a3[m][mu] + a3[n][nu]*a0[m][mu] - a1[n][nu]*a2[m][mu] + a2[n][nu]*a1[m][mu];

  k=hop[m][mu];
  j=hop[l][mu];
  q0= a0[k][nu]*a0[j][mu] + a1[k][nu]*a1[j][mu] + a2[k][nu]*a2[j][mu] + a3[k][nu]*a3[j][mu];
  q1=-a0[k][nu]*a1[j][mu] + a1[k][nu]*a0[j][mu] + a2[k][nu]*a3[j][mu] - a3[k][nu]*a2[j][mu];
  q2=-a0[k][nu]*a2[j][mu] + a2[k][nu]*a0[j][mu] + a3[k][nu]*a1[j][mu] - a1[k][nu]*a3[j][mu];
  q3=-a0[k][nu]*a3[j][mu] + a3[k][nu]*a0[j][mu] + a1[k][nu]*a2[j][mu] - a2[k][nu]*a1[j][mu];
        
  r0= p0*q0 - p1*q1 - p2*q2 - p3*q3;
  r1= p0*q1 + p1*q0 - p2*q3 + p3*q2;
  r2= p0*q2 + p2*q0 - p3*q1 + p1*q3;
  r3= p0*q3 + p3*q0 - p1*q2 + p2*q1;

  s0= r0*v0 + r1*v1 + r2*v2 + r3*v3;
  s1=-r0*v1 - r1*v0 + r2*v3 - r3*v2;
  s2=-r0*v2 - r2*v0 + r3*v1 - r1*v3;
  s3=-r0*v3 - r3*v0 + r1*v2 - r2*v1;
	
  w12_0+=s0;
  w12_0+=s1;
  w12_0+=s2;
  w12_0+=s3;    
        
  l=hop[m][mu];
  k=hop[l][nu];
  q0= a0[l][nu]*a0[k][nu] - a1[l][nu]*a1[k][nu] - a2[l][nu]*a2[k][nu] - a3[l][nu]*a3[k][nu];
  q1= a0[l][nu]*a1[k][nu] + a1[l][nu]*a0[k][nu] - a2[l][nu]*a3[k][nu] + a3[l][nu]*a2[k][nu];
  q2= a0[l][nu]*a2[k][nu] + a2[l][nu]*a0[k][nu] - a3[l][nu]*a1[k][nu] + a1[l][nu]*a3[k][nu];
  q3= a0[l][nu]*a3[k][nu] + a3[l][nu]*a0[k][nu] - a1[l][nu]*a2[k][nu] + a2[l][nu]*a1[k][nu];

  j=hop[hop[k][nu]][D+mu];
  i=hop[j][D+mu];
  r0= a0[i][mu]*a0[j][mu] - a1[i][mu]*a1[j][mu] - a2[i][mu]*a2[j][mu] - a3[i][mu]*a3[j][mu];
  r1= a0[i][mu]*a1[j][mu] + a1[i][mu]*a0[j][mu] - a2[i][mu]*a3[j][mu] + a3[i][mu]*a2[j][mu];
  r2= a0[i][mu]*a2[j][mu] + a2[i][mu]*a0[j][mu] - a3[i][mu]*a1[j][mu] + a1[i][mu]*a3[j][mu];
  r3= a0[i][mu]*a3[j][mu] + a3[i][mu]*a0[j][mu] - a1[i][mu]*a2[j][mu] + a2[i][mu]*a1[j][mu];
  
  h=hop[n][nu];
  s0= a0[n][nu]*a0[h][nu] - a1[n][nu]*a1[h][nu] - a2[n][nu]*a2[h][nu] - a3[n][nu]*a3[h][nu];
  s1= a0[n][nu]*a1[h][nu] + a1[n][nu]*a0[h][nu] - a2[n][nu]*a3[h][nu] + a3[n][nu]*a2[h][nu];
  s2= a0[n][nu]*a2[h][nu] + a2[n][nu]*a0[h][nu] - a3[n][nu]*a1[h][nu] + a1[n][nu]*a3[h][nu];
  s3= a0[n][nu]*a3[h][nu] + a3[n][nu]*a0[h][nu] - a1[n][nu]*a2[h][nu] + a2[n][nu]*a1[h][nu];
   
  u0= t0*r0 + t1*r1 + t2*r2 + t3*r3;
  u1=-t0*r1 - t1*r0 + t2*r3 - t3*r2;
  u2=-t0*r2 - t2*r0 + t3*r1 - t1*r3;
  u3=-t0*r3 - t3*r0 + t1*r2 - t2*r1;

  w0= u0*s0 + u1*s1 + u2*s2 + u3*s3;
  w1=-u0*s1 - u1*s0 + u2*s3 - u3*s2;
  w2=-u0*s2 - u2*s0 + u3*s1 - u1*s3;
  w3=-u0*s3 - u3*s0 + u1*s2 - u2*s1;

  w22_0+=w0;
  w22_1+=w1;
  w22_2+=w2;
  w22_3+=w3;
	
  i=hop[k][nu];
  h=hop[hop[i][nu]][D+mu];
  r0= a0[i][nu]*a0[h][mu] + a1[i][nu]*a1[h][mu] + a2[i][nu]*a2[h][mu] + a3[i][nu]*a3[h][mu];
  r1=-a0[i][nu]*a1[h][mu] + a1[i][nu]*a0[h][mu] + a2[i][nu]*a3[h][mu] - a3[i][nu]*a2[h][mu];
  r2=-a0[i][nu]*a2[h][mu] + a2[i][nu]*a0[h][mu] + a3[i][nu]*a1[h][mu] - a1[i][nu]*a3[h][mu];
  r3=-a0[i][nu]*a3[h][mu] + a3[i][nu]*a0[h][mu] + a1[i][nu]*a2[h][mu] - a2[i][nu]*a1[h][mu];

  g=hop[h][D+mu];
  f=hop[g][D+nu];
  t0= a0[f][nu]*a0[g][mu] - a1[f][nu]*a1[g][mu] - a2[f][nu]*a2[g][mu] - a3[n][nu]*a3[h][nu];
  t1= a0[f][nu]*a1[g][mu] + a1[f][nu]*a0[g][mu] - a2[f][nu]*a3[g][mu] + a3[n][nu]*a2[h][nu];
  t2= a0[f][nu]*a2[g][mu] + a2[f][nu]*a0[g][mu] - a3[f][nu]*a1[g][mu] + a1[n][nu]*a3[h][nu];
  t3= a0[f][nu]*a3[g][mu] + a3[f][nu]*a0[g][mu] - a1[f][nu]*a2[g][mu] + a2[n][nu]*a1[h][nu];
   
  e=hop[f][D+nu];
  s0= a0[n][nu]*a0[h][nu] - a1[n][nu]*a1[h][nu] - a2[n][nu]*a2[h][nu] - a3[n][nu]*a3[h][nu];
  s1= a0[n][nu]*a1[h][nu] + a1[n][nu]*a0[h][nu] - a2[n][nu]*a3[h][nu] + a3[n][nu]*a2[h][nu];
  s2= a0[n][nu]*a2[h][nu] + a2[n][nu]*a0[h][nu] - a3[n][nu]*a1[h][nu] + a1[n][nu]*a3[h][nu];
  s3= a0[n][nu]*a3[h][nu] + a3[n][nu]*a0[h][nu] - a1[n][nu]*a2[h][nu] + a2[n][nu]*a1[h][nu];

}

void 3_3(void)
{
   m=hop[n][mu];
   p0= a0[n][mu]*a0[m][mu] - a1[n][mu]*a1[m][mu] - a2[n][mu]*a2[m][mu] - a3[n][mu]*a3[m][mu];
   p1= a0[n][mu]*a1[m][mu] + a1[n][mu]*a0[m][mu] - a2[n][mu]*a3[m][mu] + a3[n][mu]*a2[m][mu];
   p2= a0[n][mu]*a2[m][mu] + a2[n][mu]*a0[m][mu] - a3[n][mu]*a1[m][mu] + a1[n][mu]*a3[m][mu];
   p3= a0[n][mu]*a3[m][mu] + a3[n][mu]*a0[m][mu] - a1[n][mu]*a2[m][mu] + a2[n][mu]*a1[m][mu];

   l=hop[m][mu];
   k=hop[l][nu];
   q0= a0[l][nu]*a0[k][nu] - a1[l][nu]*a1[k][nu] - a2[l][nu]*a2[k][nu] - a3[l][nu]*a3[k][nu];
   q1= a0[l][nu]*a1[k][nu] + a1[l][nu]*a0[k][nu] - a2[l][nu]*a3[k][nu] + a3[l][nu]*a2[k][nu];
   q2= a0[l][nu]*a2[k][nu] + a2[l][nu]*a0[k][nu] - a3[l][nu]*a1[k][nu] + a1[l][nu]*a3[k][nu];
   q3= a0[l][nu]*a3[k][nu] + a3[l][nu]*a0[k][nu] - a1[l][nu]*a2[k][nu] + a2[l][nu]*a1[k][nu];

   j=hop[k][nu];
   i=hop[j][nu];
   q0= a0[l][nu]*a0[k][nu] - a1[l][nu]*a1[k][nu] - a2[l][nu]*a2[k][nu] - a3[l][nu]*a3[k][nu];
   q1= a0[l][nu]*a1[k][nu] + a1[l][nu]*a0[k][nu] - a2[l][nu]*a3[k][nu] + a3[l][nu]*a2[k][nu];
   q2= a0[l][nu]*a2[k][nu] + a2[l][nu]*a0[k][nu] - a3[l][nu]*a1[k][nu] + a1[l][nu]*a3[k][nu];
   q3= a0[l][nu]*a3[k][nu] + a3[l][nu]*a0[k][nu] - a1[l][nu]*a2[k][nu] + a2[l][nu]*a1[k][nu];

   h=hop[hop[i][nu]][D+mu];
   g=hop[h][D+mu];
   s0= a0[g][nu]*a0[h][mu] - a1[g][nu]*a1[h][nu] - a2[n][nu]*a2[h][nu] - a3[n][nu]*a3[h][nu];
   s1= a0[g][nu]*a1[h][mu] + a1[g][nu]*a0[h][nu] - a2[n][nu]*a3[h][nu] + a3[n][nu]*a2[h][nu];
   s2= a0[g][nu]*a2[h][mu] + a2[g][nu]*a0[h][nu] - a3[n][nu]*a1[h][nu] + a1[n][nu]*a3[h][nu];
   s3= a0[g][nu]*a3[h][mu] + a3[g][nu]*a0[h][nu] - a1[n][nu]*a2[h][nu] + a2[n][nu]*a1[h][nu];

   f=hop[g][D+mu];
   e=hop[f][D+nu];
   
   d=hop[e][D+nu];
   s0= a0[n][nu]*a0[h][nu] - a1[n][nu]*a1[h][nu] - a2[n][nu]*a2[h][nu] - a3[n][nu]*a3[h][nu];
   s1= a0[n][nu]*a1[h][nu] + a1[n][nu]*a0[h][nu] - a2[n][nu]*a3[h][nu] + a3[n][nu]*a2[h][nu];
   s2= a0[n][nu]*a2[h][nu] + a2[n][nu]*a0[h][nu] - a3[n][nu]*a1[h][nu] + a1[n][nu]*a3[h][nu];
   s3= a0[n][nu]*a3[h][nu] + a3[n][nu]*a0[h][nu] - a1[n][nu]*a2[h][nu] + a2[n][nu]*a1[h][nu];
}

void 4_3(void)
{
 int n;

   m=hop[n][mu];
   p0= a0[n][mu]*a0[m][mu] - a1[n][mu]*a1[m][mu] - a2[n][mu]*a2[m][mu] - a3[n][mu]*a3[m][mu];
   p1= a0[n][mu]*a1[m][mu] + a1[n][mu]*a0[m][mu] - a2[n][mu]*a3[m][mu] + a3[n][mu]*a2[m][mu];
   p2= a0[n][mu]*a2[m][mu] + a2[n][mu]*a0[m][mu] - a3[n][mu]*a1[m][mu] + a1[n][mu]*a3[m][mu];
   p3= a0[n][mu]*a3[m][mu] + a3[n][mu]*a0[m][mu] - a1[n][mu]*a2[m][mu] + a2[n][mu]*a1[m][mu];

   l=hop[m][mu];
   k=hop[l][mu];
   q0= a0[l][nu]*a0[k][nu] - a1[l][nu]*a1[k][nu] - a2[l][nu]*a2[k][nu] - a3[l][nu]*a3[k][nu];
   q1= a0[l][nu]*a1[k][nu] + a1[l][nu]*a0[k][nu] - a2[l][nu]*a3[k][nu] + a3[l][nu]*a2[k][nu];
   q2= a0[l][nu]*a2[k][nu] + a2[l][nu]*a0[k][nu] - a3[l][nu]*a1[k][nu] + a1[l][nu]*a3[k][nu];
   q3= a0[l][nu]*a3[k][nu] + a3[l][nu]*a0[k][nu] - a1[l][nu]*a2[k][nu] + a2[l][nu]*a1[k][nu];

   j=hop[k][mu];
   i=hop[j][nu];

   h=hop[i][nu];
   g=hop[hop[h][nu]][D+mu];

   f=hop[g][D+mu];
   e=hop[f][D+mu];

   d=hop[e][D+mu];
   c=hop[d][D+nu];

   b=hop[c][D+nu];
}

void 4_4(void)
{
   m=hop[n][mu];

   l=hop[m][mu];
   k=hop[l][mu];

   j=hop[k][mu];
   i=hop[j][nu];

   h=hop[i][nu];
   g=hop[h][nu];

   f=hop[hop[g][nu]][D+mu];
   e=hop[e][D+mu];

   d=hop[e][D+mu];
   c=hop[d][D+mu];

   b=hop[c][D+nu];
   a=hop[b][D+nu];

   z=hop[a][D+nu];

}
