#include "su2.h"

static double w11_0;
/*static double w12_0;*/
static double w22_0;
/*static double w23_0;*/
static double w33_0;
/*static double w34_0;*/
static double w44_0;
static int nmeas;

void init_measure()
{
 w11_0=0.;
/* w12_0=0.;*/
 w22_0=0.;
/* w23_0=0.;*/
 w33_0=0.;
/* w34_0=0.;*/
 w44_0=0.;
 nmeas=0;
}

void measure()
{
 int a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x;
 int mu,nu;
 int number_loop;

 double l12_0,l22_0,l23_0,l33_0,l34_0,l44_0;
 double loop12,loop22,loop23,loop33,loop34,loop44;
 double ac0,ac1,ac2,ac3;
 double ce0,ce1,ce2,ce3;
 double eg0,eg1,eg2,eg3;
 double gy0,gy1,gy2,gy3;
 double yj0,yj1,yj2,yj3;
 double jl0,jl1,jl2,jl3;
 double ln0,ln1,ln2,ln3;
 double na0,na1,na2,na3;
 double cp0,cp1,cp2,cp3;
 double pa0,pa1,pa2,pa3;
 double cx0,cx1,cx2,cx3;
 double xn0,xn1,xn2,xn3;
 double xv0,xv1,xv2,xv3;
 double vn0,vn1,vn2,vn3;
 double cr0,cr1,cr2,cr3;
 double rt0,rt1,rt2,rt3;
 double tv0,tv1,tv2,tv3;
 double gt0,gt1,gt2,gt3;
 double ap0,ap1,ap2,ap3;
 double ax0,ax1,ax2,ax3;
 double xa0,xa1,xa2,xa3;
 double av0,av1,av2,av3;
 double va0,va1,va2,va3;
 double ar0,ar1,ar2,ar3;
 double ta0,ta1,ta2,ta3;
 double ae0,ae1,ae2,ae3;
 double	ag0,ag1,ag2,ag3;	
 double ga0,ga1,ga2,ga3;
 double	ay0,ay1,ay2,ay3;
 double	yl0,yl1,yl2,yl3;
 double la0,la1,la2,la3; 
 double ya0,ya1,ya2,ya3;
 double at0,at1,at2,at3;
 double P;
 double kai11,kai22,kai33,kai44;

 number_loop=D*(D-1)*V;
 
 l12_0=0.0;
 l22_0=0.0;
 l23_0=0.0;
 l33_0=0.0;
 l34_0=0.0;
 l44_0=0.0;

 for (a=0;a<V;a++){
     for (mu=0;mu<D;mu++){
	 for (nu=0;nu<D;nu++){
	     if (nu != mu ){

  b=hop[a][mu];
  c=hop[b][mu];
  d=hop[c][mu];
  e=hop[d][mu];
  f=hop[e][nu];
  g=hop[f][nu];
  h=hop[g][nu];
  i=hop[hop[h][nu]][D+mu];
  j=hop[i][D+mu];
  k=hop[j][D+mu];
  l=hop[k][D+mu];
  m=hop[l][D+nu];
  n=hop[m][D+nu];
  o=hop[n][D+nu];
  p=hop[o][mu];
  q=hop[p][mu];
  r=hop[q][mu];
  s=hop[r][nu];
  t=hop[s][nu];
  u=hop[t][D+mu];
  v=hop[u][D+mu];
  w=hop[v][D+nu];
  x=hop[w][mu];

  ac0= a0[a][mu]*a0[b][mu] - a1[a][mu]*a1[b][mu] - a2[a][mu]*a2[b][mu] - a3[a][mu]*a3[b][mu];
  ac1= a0[a][mu]*a1[b][mu] + a1[a][mu]*a0[b][mu] - a2[a][mu]*a3[b][mu] + a3[a][mu]*a2[b][mu];
  ac2= a0[a][mu]*a2[b][mu] + a2[a][mu]*a0[b][mu] - a3[a][mu]*a1[b][mu] + a1[a][mu]*a3[b][mu];
  ac3= a0[a][mu]*a3[b][mu] + a3[a][mu]*a0[b][mu] - a1[a][mu]*a2[b][mu] + a2[a][mu]*a1[b][mu];
 
  ce0= a0[c][mu]*a0[d][mu] - a1[c][mu]*a1[d][mu] - a2[c][mu]*a2[d][mu] - a3[c][mu]*a3[d][mu];
  ce1= a0[c][mu]*a1[d][mu] + a1[c][mu]*a0[d][mu] - a2[c][mu]*a3[d][mu] + a3[c][mu]*a2[d][mu];
  ce2= a0[c][mu]*a2[d][mu] + a2[c][mu]*a0[d][mu] - a3[c][mu]*a1[d][mu] + a1[c][mu]*a3[d][mu];
  ce3= a0[c][mu]*a3[d][mu] + a3[c][mu]*a0[d][mu] - a1[c][mu]*a2[d][mu] + a2[c][mu]*a1[d][mu];

  eg0= a0[e][nu]*a0[f][nu] - a1[e][nu]*a1[f][nu] - a2[e][nu]*a2[f][nu] - a3[e][nu]*a3[f][nu];
  eg1= a0[e][nu]*a1[f][nu] + a1[e][nu]*a0[f][nu] - a2[e][nu]*a3[f][nu] + a3[e][nu]*a2[f][nu];
  eg2= a0[e][nu]*a2[f][nu] + a2[e][nu]*a0[f][nu] - a3[e][nu]*a1[f][nu] + a1[e][nu]*a3[f][nu];
  eg3= a0[e][nu]*a3[f][nu] + a3[e][nu]*a0[f][nu] - a1[e][nu]*a2[f][nu] + a2[e][nu]*a1[f][nu];

  gy0= a0[g][nu]*a0[h][nu] - a1[g][nu]*a1[h][nu] - a2[g][nu]*a2[h][nu] - a3[g][nu]*a3[h][nu];
  gy1= a0[g][nu]*a1[h][nu] + a1[g][nu]*a0[h][nu] - a2[g][nu]*a3[h][nu] + a3[g][nu]*a2[h][nu];
  gy2= a0[g][nu]*a2[h][nu] + a2[g][nu]*a0[h][nu] - a3[g][nu]*a1[h][nu] + a1[g][nu]*a3[h][nu];
  gy3= a0[g][nu]*a3[h][nu] + a3[g][nu]*a0[h][nu] - a1[g][nu]*a2[h][nu] + a2[g][nu]*a1[h][nu];

  yj0= a0[j][mu]*a0[i][mu] - a1[j][mu]*a1[i][mu] - a2[j][mu]*a2[i][mu] - a3[j][mu]*a3[i][mu];
  yj1= a0[j][mu]*a1[i][mu] + a1[j][mu]*a0[i][mu] - a2[j][mu]*a3[i][mu] + a3[j][mu]*a2[i][mu];
  yj2= a0[j][mu]*a2[i][mu] + a2[j][mu]*a0[i][mu] - a3[j][mu]*a1[i][mu] + a1[j][mu]*a3[i][mu];
  yj3= a0[j][mu]*a3[i][mu] + a3[j][mu]*a0[i][mu] - a1[j][mu]*a2[i][mu] + a2[j][mu]*a1[i][mu];

  jl0= a0[l][mu]*a0[k][mu] - a1[l][mu]*a1[k][mu] - a2[l][mu]*a2[k][mu] - a3[l][mu]*a3[k][mu];
  jl1= a0[l][mu]*a1[k][mu] + a1[l][mu]*a0[k][mu] - a2[l][mu]*a3[k][mu] + a3[l][mu]*a2[k][mu];
  jl2= a0[l][mu]*a2[k][mu] + a2[l][mu]*a0[k][mu] - a3[l][mu]*a1[k][mu] + a1[l][mu]*a3[k][mu];
  jl3= a0[l][mu]*a3[k][mu] + a3[l][mu]*a0[k][mu] - a1[l][mu]*a2[k][mu] + a2[l][mu]*a1[k][mu];

  ln0= a0[n][nu]*a0[m][nu] - a1[n][nu]*a1[m][nu] - a2[n][nu]*a2[m][nu] - a3[n][nu]*a3[m][nu];
  ln1= a0[n][nu]*a1[m][nu] + a1[n][nu]*a0[m][nu] - a2[n][nu]*a3[m][nu] + a3[n][nu]*a2[m][nu];
  ln2= a0[n][nu]*a2[m][nu] + a2[n][nu]*a0[m][nu] - a3[n][nu]*a1[m][nu] + a1[n][nu]*a3[m][nu];
  ln3= a0[n][nu]*a3[m][nu] + a3[n][nu]*a0[m][nu] - a1[n][nu]*a2[m][nu] + a2[n][nu]*a1[m][nu];

  na0= a0[a][nu]*a0[o][nu] - a1[a][nu]*a1[o][nu] - a2[a][nu]*a2[o][nu] - a3[a][nu]*a3[o][nu];
  na1= a0[a][nu]*a1[o][nu] + a1[a][nu]*a0[o][nu] - a2[a][nu]*a3[o][nu] + a3[a][nu]*a2[o][nu];
  na2= a0[a][nu]*a2[o][nu] + a2[a][nu]*a0[o][nu] - a3[a][nu]*a1[o][nu] + a1[a][nu]*a3[o][nu];
  na3= a0[a][nu]*a3[o][nu] + a3[a][nu]*a0[o][nu] - a1[a][nu]*a2[o][nu] + a2[a][nu]*a1[o][nu];

  cp0= a0[c][nu]*a0[p][mu] + a1[c][nu]*a1[p][mu] + a2[c][nu]*a2[p][mu] + a3[c][nu]*a3[p][mu];
  cp1=-a0[c][nu]*a1[p][mu] + a1[c][nu]*a0[p][mu] + a2[c][nu]*a3[p][mu] - a3[c][nu]*a2[p][mu];
  cp2=-a0[c][nu]*a2[p][mu] + a2[c][nu]*a0[p][mu] + a3[c][nu]*a1[p][mu] - a1[c][nu]*a3[p][mu];
  cp3=-a0[c][nu]*a3[p][mu] + a3[c][nu]*a0[p][mu] + a1[c][nu]*a2[p][mu] - a2[c][nu]*a1[p][mu];

  pa0= a0[a][nu]*a0[o][mu] - a1[a][nu]*a1[o][mu] - a2[a][nu]*a2[o][mu] - a3[a][nu]*a3[o][mu];
  pa1= a0[a][nu]*a1[o][mu] + a1[a][nu]*a0[o][mu] - a2[a][nu]*a3[o][mu] + a3[a][nu]*a2[o][mu];
  pa2= a0[a][nu]*a2[o][mu] + a2[a][nu]*a0[o][mu] - a3[a][nu]*a1[o][mu] + a1[a][nu]*a3[o][mu];
  pa3= a0[a][nu]*a3[o][mu] + a3[a][nu]*a0[o][mu] - a1[a][nu]*a2[o][mu] + a2[a][nu]*a1[o][mu];

  cx0= a0[c][nu]*a0[q][nu] - a1[c][nu]*a1[q][nu] - a2[c][nu]*a2[q][nu] - a3[c][nu]*a3[q][nu];
  cx1= a0[c][nu]*a1[q][nu] + a1[c][nu]*a0[q][nu] - a2[c][nu]*a3[q][nu] + a3[c][nu]*a2[q][nu];
  cx2= a0[c][nu]*a2[q][nu] + a2[c][nu]*a0[q][nu] - a3[c][nu]*a1[q][nu] + a1[c][nu]*a3[q][nu];
  cx3= a0[c][nu]*a3[q][nu] + a3[c][nu]*a0[q][nu] - a1[c][nu]*a2[q][nu] + a2[c][nu]*a1[q][nu];

  xn0= a0[n][mu]*a0[w][mu] - a1[n][mu]*a1[w][mu] - a2[n][mu]*a2[w][mu] - a3[n][mu]*a3[w][mu];
  xn1= a0[n][mu]*a1[w][mu] + a1[n][mu]*a0[w][mu] - a2[n][mu]*a3[w][mu] + a3[n][mu]*a2[w][mu];
  xn2= a0[n][mu]*a2[w][mu] + a2[n][mu]*a0[w][mu] - a3[n][mu]*a1[w][mu] + a1[n][mu]*a3[w][mu];
  xn3= a0[n][mu]*a3[w][mu] + a3[n][mu]*a0[w][mu] - a1[n][mu]*a2[w][mu] + a2[n][mu]*a1[w][mu];

  xv0= a0[x][nu]*a0[v][mu] + a1[x][nu]*a1[v][mu] + a2[x][nu]*a2[v][mu] + a3[x][nu]*a3[v][mu];
  xv1=-a0[x][nu]*a1[v][mu] + a1[x][nu]*a0[v][mu] + a2[x][nu]*a3[v][mu] - a3[x][nu]*a2[v][mu];
  xv2=-a0[x][nu]*a2[v][mu] + a2[x][nu]*a0[v][mu] + a3[x][nu]*a1[v][mu] - a1[x][nu]*a3[v][mu];
  xv3=-a0[x][nu]*a3[v][mu] + a3[x][nu]*a0[v][mu] + a1[x][nu]*a2[v][mu] - a2[x][nu]*a1[v][mu];

  vn0= a0[n][nu]*a0[m][mu] - a1[n][nu]*a1[m][mu] - a2[n][nu]*a2[m][mu] - a3[n][nu]*a3[m][mu];
  vn1= a0[n][nu]*a1[m][mu] + a1[n][nu]*a0[m][mu] - a2[n][nu]*a3[m][mu] + a3[n][nu]*a2[m][mu];
  vn2= a0[n][nu]*a2[m][mu] + a2[n][nu]*a0[m][mu] - a3[n][nu]*a1[m][mu] + a1[n][nu]*a3[m][mu];
  vn3= a0[n][nu]*a3[m][mu] + a3[n][nu]*a0[m][mu] - a1[n][nu]*a2[m][mu] + a2[n][nu]*a1[m][mu];

  cr0= a0[c][mu]*a0[d][nu] - a1[c][mu]*a1[d][nu] - a2[c][mu]*a2[d][nu] - a3[c][mu]*a3[d][nu];
  cr1= a0[c][mu]*a1[d][nu] + a1[c][mu]*a0[d][nu] - a2[c][mu]*a3[d][nu] + a3[c][mu]*a2[d][nu];
  cr2= a0[c][mu]*a2[d][nu] + a2[c][mu]*a0[d][nu] - a3[c][mu]*a1[d][nu] + a1[c][mu]*a3[d][nu];
  cr3= a0[c][mu]*a3[d][nu] + a3[c][mu]*a0[d][nu] - a1[c][mu]*a2[d][nu] + a2[c][mu]*a1[d][nu];

  rt0= a0[r][nu]*a0[s][nu] - a1[r][nu]*a1[s][nu] - a2[r][nu]*a2[s][nu] - a3[r][nu]*a3[s][nu];
  rt1= a0[r][nu]*a1[s][nu] + a1[r][nu]*a0[s][nu] - a2[r][nu]*a3[s][nu] + a3[r][nu]*a2[s][nu];
  rt2= a0[r][nu]*a2[s][nu] + a2[r][nu]*a0[s][nu] - a3[r][nu]*a1[s][nu] + a1[r][nu]*a3[s][nu];
  rt3= a0[r][nu]*a3[s][nu] + a3[r][nu]*a0[s][nu] - a1[r][nu]*a2[s][nu] + a2[r][nu]*a1[s][nu];

  tv0= a0[v][mu]*a0[u][mu] - a1[v][mu]*a1[u][mu] - a2[v][mu]*a2[u][mu] - a3[v][mu]*a3[u][mu];
  tv1= a0[v][mu]*a1[u][mu] + a1[v][mu]*a0[u][mu] - a2[v][mu]*a3[u][mu] + a3[v][mu]*a2[u][mu];
  tv2= a0[v][mu]*a2[u][mu] + a2[v][mu]*a0[u][mu] - a3[v][mu]*a1[u][mu] + a1[v][mu]*a3[u][mu];
  tv3= a0[v][mu]*a3[u][mu] + a3[v][mu]*a0[u][mu] - a1[v][mu]*a2[u][mu] + a2[v][mu]*a1[u][mu];

  gt0= a0[g][nu]*a0[t][mu] + a1[g][nu]*a1[t][mu] + a2[g][nu]*a2[t][mu] + a3[g][nu]*a3[t][mu];
  gt1=-a0[g][nu]*a1[t][mu] + a1[g][nu]*a0[t][mu] + a2[g][nu]*a3[t][mu] - a3[g][nu]*a2[t][mu];
  gt2=-a0[g][nu]*a2[t][mu] + a2[g][nu]*a0[t][mu] + a3[g][nu]*a1[t][mu] - a1[g][nu]*a3[t][mu];
  gt3=-a0[g][nu]*a3[t][mu] + a3[g][nu]*a0[t][mu] + a1[g][nu]*a2[t][mu] - a2[g][nu]*a1[t][mu];

  ap0= ac0*cp0 - ac1*cp1 - ac2*cp2 - ac3*cp3;
  ap1= ac0*cp1 + ac1*cp0 - ac2*cp3 + ac3*cp2;
  ap2= ac0*cp2 + ac2*cp0 - ac3*cp1 + ac1*cp3;
  ap3= ac0*cp3 + ac3*cp0 - ac1*cp2 + ac2*cp1;

  l12_0+= ap0*pa0 + ap1*pa1 + ap2*pa2 + ap3*pa3;
/*w12*/

  ax0= ac0*cx0 - ac1*cx1 - ac2*cx2 - ac3*cx3;
  ax1= ac0*cx1 + ac1*cx0 - ac2*cx3 + ac3*cx2;
  ax2= ac0*cx2 + ac2*cx0 - ac3*cx1 + ac1*cx3;
  ax3= ac0*cx3 + ac3*cx0 - ac1*cx2 + ac2*cx1;

  xa0= na0*xn0 - na1*xn1 - na2*xn2 - na3*xn3;
  xa1= na0*xn1 + na1*xn0 - na2*xn3 + na3*xn2;
  xa2= na0*xn2 + na2*xn0 - na3*xn1 + na1*xn3;
  xa3= na0*xn3 + na3*xn0 - na1*xn2 + na2*xn1;

  l22_0+= ax0*xa0 + ax1*xa1 + ax2*xa2 + ax3*xa3;
/*w22*/

  av0= ax0*xv0 - ax1*xv1 - ax2*xv2 - ax3*xv3;
  av1= ax0*xv1 + ax1*xv0 - ax2*xv3 + ax3*xv2;
  av2= ax0*xv2 + ax2*xv0 - ax3*xv1 + ax1*xv3;
  av3= ax0*xv3 + ax3*xv0 - ax1*xv2 + ax2*xv1;
  
  va0= na0*vn0 - na1*vn1 - na2*vn2 - na3*vn3;
  va1= na0*vn1 + na1*vn0 - na2*vn3 + na3*vn2;
  va2= na0*vn2 + na2*vn0 - na3*vn1 + na1*vn3;
  va3= na0*vn3 + na3*vn0 - na1*vn2 + na2*vn1;

  l23_0+= av0*va0 + av1*va1 + av2*va2 + av3*va3;
/*w23*/ 

  ar0= ac0*cr0 - ac1*cr1 - ac2*cr2 - ac3*cr3;
  ar1= ac0*cr1 + ac1*cr0 - ac2*cr3 + ac3*cr2;
  ar2= ac0*cr2 + ac2*cr0 - ac3*cr1 + ac1*cr3;
  ar3= ac0*cr3 + ac3*cr0 - ac1*cr2 + ac2*cr1;

  at0= ar0*rt0 - ar1*rt1 - ar2*rt2 - ar3*rt3;
  at1= ar0*rt1 + ar1*rt0 - ar2*rt3 + ar3*rt2;
  at2= ar0*rt2 + ar2*rt0 - ar3*rt1 + ar1*rt3;
  at3= ar0*rt3 + ar3*rt0 - ar1*rt2 + ar2*rt1;

  ta0= va0*tv0 - va1*tv1 - va2*tv2 - va3*tv3;
  ta1= va0*tv1 + va1*tv0 - va2*tv3 + va3*tv2;
  ta2= va0*tv2 + va2*tv0 - va3*tv1 + va1*tv3;
  ta3= va0*tv3 + va3*tv0 - va1*tv2 + va2*tv1;

  l33_0+= at0*ta0 + at1*ta1 + at2*ta2 + at3*ta3;
/*w33*/
 
  ae0= ac0*ce0 - ac1*ce1 - ac2*ce2 - ac3*ce3;
  ae1= ac0*ce1 + ac1*ce0 - ac2*ce3 + ac3*ce2;
  ae2= ac0*ce2 + ac2*ce0 - ac3*ce1 + ac1*ce3;
  ae3= ac0*ce3 + ac3*ce0 - ac1*ce2 + ac2*ce1;

  ag0= ae0*eg0 - ae1*eg1 - ae2*eg2 - ae3*eg3;
  ag1= ae0*eg1 + ae1*eg0 - ae2*eg3 + ae3*eg2;
  ag2= ae0*eg2 + ae2*eg0 - ae3*eg1 + ae1*eg3;
  ag3= ae0*eg3 + ae3*eg0 - ae1*eg2 + ae2*eg1;

  ga0= gt0*ta0 + gt1*ta1 + gt2*ta2 + gt3*ta3;
  ga1=-gt0*ta1 + gt1*ta0 + gt2*ta3 - gt3*ta2;
  ga2=-gt0*ta2 + gt2*ta0 + gt3*ta1 - gt1*ta3;
  ga3=-gt0*ta3 + gt3*ta0 + gt1*ta2 - gt2*ta1;

  l34_0+= ag0*ga0 - ag1*ga1 - ag2*ga2 - ag3*ga3;
/*w34*/ 

  ay0= ag0*gy0 - ag1*gy1 - ag2*gy2 - ag3*gy3;
  ay1= ag0*gy1 + ag1*gy0 - ag2*gy3 + ag3*gy2;
  ay2= ag0*gy2 + ag2*gy0 - ag3*gy1 + ag1*gy3;
  ay3= ag0*gy3 + ag3*gy0 - ag1*gy2 + ag2*gy1;

  yl0= jl0*yj0 - jl1*yj1 - jl2*yj2 - jl3*yj3;
  yl1= jl0*yj1 + jl1*yj0 - jl2*yj3 + jl3*yj2;
  yl2= jl0*yj2 + jl2*yj0 - jl3*yj1 + jl1*yj3;
  yl3= jl0*yj3 + jl3*yj0 - jl1*yj2 + jl2*yj1;

  la0= na0*ln0 - na1*ln1 - na2*ln2 - na3*ln3;
  la1= na0*ln1 + na1*ln0 - na2*ln3 + na3*ln2;
  la2= na0*ln2 + na2*ln0 - na3*ln1 + na1*ln3;
  la3= na0*ln3 + na3*ln0 - na1*ln2 + na2*ln1;

  ya0= la0*yl0 - la1*yl1 - la2*yl2 - la3*yl3;
  ya1= la0*yl1 + la1*yl0 - la2*yl3 + la3*yl2;
  ya2= la0*yl2 + la2*yl0 - la3*yl1 + la1*yl3;
  ya3= la0*yl3 + la3*yl0 - la1*yl2 + la2*yl1;

  l44_0+= ay0*ya0 + ay1*ya1 + ay2*ya2 + ay3*ya3;
/*w44*/
             }
	 }
     }
 }
  
  loop12=l12_0/number_loop;
  loop22=l22_0/number_loop;
  loop23=l23_0/number_loop;
  loop33=l33_0/number_loop;
  loop34=l34_0/number_loop;
  loop44=l44_0/number_loop;

  P=plaquette();

  kai11=-log(P);
  kai22=-log(loop22*P/(loop12*loop12));
  kai33=-log(loop33*loop22/(loop23*loop23));
  kai44=-log(loop44*loop33/(loop34*loop34));

  printf("%e %e %e %e",kai11,kai22,kai33,kai44);
  printf("\n");
  
/*  w11_0+=P;
  w12_0+=loop12;
  w22_0+=loop22;
  w23_0+=loop23;
  w33_0+=loop33;
  w34_0+=loop34;
  w44_0+=loop44;*/
 
  w11_0+=kai11;
  w22_0+=kai22;
  w33_0+=kai33;
  w44_0+=kai44;
  nmeas++;
}

void print_meas()
{
/*    printf("MEAS %e %e %e %e %e %e %e",w11_0/nmeas,w12_0/nmeas,w22_0/nmeas,w23_0/nmeas,w33_0/nmeas,w34_0/nmeas,w44_0/nmeas);
*/
      printf("MEAS %e %e %e %e",w11_0/nmeas,w22_0/nmeas,w33_0/nmeas,w44_0/nmeas);	
      printf("\n");
    /*    init_measure(); */
}  







