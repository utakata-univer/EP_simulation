


void move_phi(double eps)
{
  int i;
  for (i=0;i<V,i++) phi[i]+=mom[i]*eps
}


void move_mom(double eps)
{
  int i,mu;
  double J, force;

  for (i=0;i<V;i++)
  {
     J=0;
     for (mu=0;mu<2*D;mu++) J+=phi[hop[i][mu]];

     force=2*kappa*J-2*phi[i]-lambda*4*(phi[i]*phi[i]-1)*phi[i];
     mom[i]+=force*eps;
  }
}


double hamiltonian(void)
{
  int i,mu;
  double phin,H,phi2,mom2;
  double kappa =act_params.kappa;
  double lambda=act_params.lambda;

  H=0;

  for (i=0;i<V;i++)
  {
      phin=0;
      for (mu=0;mu<D;mu++) phin+=phi[hop[i][mu]];
      
      phi2=phi[i]*phi[i];
      mom2=mom[i]*mom[i];
      H+=mom2-2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
  }

  return H;
}



int main(void)
{ 
  int i;
  
  for (i=0;i<N_s;i++)
  {
     move_phi(eps/2);
     move_mom(eps);
     move_phi(eps/2);
  }

 
  ham=hamiltonian();
  
  printf("hamiltonian = %lf\n",ham)

}

