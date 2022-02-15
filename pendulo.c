#include<stdio.h>
#include<stdlib.h>
#include<math.h>


double V(double v0, double theta){
  double res = v0*(1-cos(theta));
  return(res);
}
	 
void main(int argc, char *argv[]){
  void LU ();
  int i,j,k, N;
  k=0;
  double a[2000],c[2000],e[2000],x[2000],xnew[2000],lnew, error, norma, norma0, p, h, w;
  //caso 2*pi (por arriba)
  //w=6.0;
  double l, Z;
  Z=M_PI+0.5;          /// controla dimensiones grafica
  //  l=-5967.330925;            /// valor aproximado del eigenvalor
  sscanf(argv[1],"%lf",&l);
  h=0.005;         /// paso de aproximacion finita
  N=ceil((2*Z)/h)-1;/// numero de iteraciones
  norma0=0;

  double V0 = 1;
  
  FILE*arch=fopen("amplitud.text", "w");
  FILE*arch1=fopen("densidad.text","w");
  
  for(i=0;i<N;i++){
    //a[i], c[i], e[i] conforman la matriz A
    a[i]=2/(h*h) + V(V0,-Z+(i+1)*h);
    c[i]=-1/(h*h);
    e[i]=-1/(h*h);
    //x[i]=f(hx*(i+1)); //Condicion inicial
  }
  
  
  for(i=0;i<N;i++){
    x[i]=1/sqrt(N);
    a[i]=a[i]-l; //l es la constante que afecta a la identidad
    
  }//cierra for
  
  do{
    
    k++;
    //Llamar LU
    LU (a,x,c,e,xnew,N);
    
    //Norma
    norma=0;
    for(i=0;i<N;i++){
      norma +=(xnew[i]*xnew[i]);
    }
    
    norma=sqrt(norma);
    printf("%lf\n", norma);
    //Normalizando
    for(i=0;i<N;i++){
      xnew[i] = xnew[i]/norma;
    }
    
    //Precision
    error=norma-norma0;
    
    for(i=0;i<N;i++){
      x[i]=xnew[i];
    }
    norma0=norma;
  }
  while(error>1E-12);
  //Por abajo, cambia el signo a l+1
  printf("Eigen Valor %lf\n", (l+1/norma));
  printf("Iteraciones %d\n", k );
  //  fprintf(arch,"%lf, %lf\n",Z,0.0);
  //  fprintf(arch1,"%lf, %lf\n",Z,0.0);
  for(i=0;i<N;i++){
    fprintf(arch,"%lf, %lf\n",Z-(i+1)*h, x[i]);
    fprintf(arch1,"%lf, %lf\n",Z-(i+1)*h, x[i]*x[i]);
    //fprintf(arch,"\n");
  }
  //  fprintf(arch,"%lf, %lf\n",Z,0.0);
  //fprintf(arch1,"%lf, %lf\n",Z,0.0);
  fclose(arch);
  fclose(arch1);
}
void LU(double a[], double b[], double c[], double e[], double x[], int N)
{
  int i;
  i = 0;
  double w[N];
  double u[N];
  double y[N];
  
  w[0] = a[0];
  y[0] = b[0]/w[0];
  for(i = 1;i <= N-1;i++)
    {
      u[i-1] = e[i-1]/w[i-1];
      w[i] = a[i] - u[i-1]*c[i];
      y[i] = (b[i] - c[i]*y[i-1])/w[i];
    }//cierra for1
  //u[N-1] = e[N-1]/w[N-1];
  
  x[N-1] = y[N-1];
  for(i = N-2;i >= 0;i--)
    {
      x[i] = y[i] - u[i]*x[i+1];
    }//cierra for2
}//cierra void LU
