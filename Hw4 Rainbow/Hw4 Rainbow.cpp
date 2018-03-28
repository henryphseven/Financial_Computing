#include <afxwin.h>
#include "ddxddv.h"
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <ctime>
const int M=5;
const double PI=3.14159265359;
const int L=500;
const int MAX=100000;
double Sum;
int i,j,k,n;
//Generate standard normal random variable
double Normal()
{
     double Y,X[2];
	 for(int j=0;j<2;j++)
	 {
	   X[j]=double(rand())/RAND_MAX;
	 }
     Y=sqrt((-2)*log(X[0]))*cos(2*PI*X[1]);
	 return Y;
}
//Solve Cholesky decomposition
void Cholesky(double A[][M],double C[][M],int N)
{
	for(i=0;i<M;i++)for(j=0;j<M;j++)A[i][j]=0;
   	A[0][0]=sqrt(C[0][0]);
   	for(j=1;j<N;j++)A[0][j]=C[0][j]/A[0][0];   
   	for(i=1;i<N-1;i++)
   	{
     Sum=0; 
   	 for(k=0;k<i;k++)Sum+=pow(A[k][i],2);
   	 A[i][i]=sqrt(C[i][i]-Sum);
   	 for(j=i+1;j<N;j++)
   	 {
   	  Sum=0;
	  for(k=0;k<i;k++)Sum+=A[k][i]*A[k][j];
  	  A[i][j]=(C[i][j]-Sum)/A[i][i];
     }
    }
   	Sum=0;
   	for(k=0;k<N-1;k++){Sum+=pow(A[k][N-1],2);}
   	A[N-1][N-1]=sqrt(C[N-1][N-1]-Sum);
   	return;
}
void swap(double *x,double *y)
{
 	 double temp=*x;
 	 *x=*y;
 	 *y=temp;
 	 return;
}
//Use Gauss-Jordan elimination to calculate inverse(A)*B
void GJ(double A[][M],double B[][M],double X[][M],int N)
{
 	 double pivot; 
     double G[M][2*M];
 	 int p;
 	 for(i=0;i<M;i++)for(j=0;j<2*M;j++)G[i][j]=0;
     for(i=0;i<N;i++)for(j=0;j<N;j++)G[i][j]=A[i][j];
     for(i=0;i<N;i++)for(j=0;j<N;j++)G[i][j+M]=B[i][j];
     for(j=0;j<N;j++)
     {   
	 	 pivot=0;
	 	 for(i=j;i<N;i++)
	 	 {
		  if(fabs(G[i][j])>pivot){pivot=G[i][j]; p=i;}
	     }
	     if(p!=j)for(k=0;k<2*M;k++)swap(&G[p][k],&G[j][k]);
	     for(k=0;k<2*M;k++)G[j][k]=G[j][k]/pivot;
	     double multiplier;
	     for(i=0;i<N;i++)
		 {
		   multiplier=G[i][j];
		   for(k=0;k<2*M;k++)
	       {
	        if(i!=j)G[i][k]=G[i][k]-multiplier*G[j][k];
           }
		 }
	 }
	 for(i=0;i<M;i++)for(j=0;j<M;j++)X[i][j]=G[i][j+M];
	 return;
} 
class MyFrameWindow : public CFrameWnd
{
	public:

	 double S[M],q[M],C[M][M],K,r,T;
     int N,m;
     double STDERR, Upper_Limit, Lower_Limit;
     int With,Without,Indirect,Direct;
	 
	MyFrameWindow ( );

	afx_msg void OnPaint( );
	afx_msg void OnInput( );

	DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(MyFrameWindow, CFrameWnd)
	ON_WM_PAINT()
	ON_COMMAND(ID_INPUT, OnInput)
END_MESSAGE_MAP()

class MyApp : public CWinApp
{
	public:

	MyApp ( )  
		: CWinApp ( "Rainbow Option Price" )
	{}

	BOOL InitInstance( )
	{
		CFrameWnd *MyFrame = new MyFrameWindow;
		m_pMainWnd = MyFrame;
		MyFrame->LoadFrame ( IDR_MAINFRAME );
		MyFrame->ShowWindow ( SW_SHOW );
		return TRUE;
	}
} DDVDDXAPP;

class MyDialog : public CDialog
{
	public:

     double S[M],q[M],C[M][M],K,r,T;
     int N,m;
     int With,Without,Indirect,Direct;

	MyDialog ( CWnd* parent ) : CDialog ( IDD_INPUT, parent ) {}
	void DoDataExchange( CDataExchange* pDX );
};

void MyDialog::DoDataExchange( CDataExchange* pDX )
{
	CDialog::DoDataExchange ( pDX );
    DDX_Text ( pDX, IDC_NN, N );
	DDV_MinMaxInt( pDX, N, 1, M );
	DDX_Text ( pDX, IDC_r, r );
	DDV_MinMaxDouble( pDX, r, 0, 1 );
	DDX_Text ( pDX, IDC_T, T );
	DDX_Text ( pDX, IDC_K, K );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, 1000 );
	DDX_Text ( pDX, IDC_S1, S[0] );
	DDX_Text ( pDX, IDC_S2, S[1] );
	DDX_Text ( pDX, IDC_S3, S[2] );
	DDX_Text ( pDX, IDC_S4, S[3] );
	DDX_Text ( pDX, IDC_S5, S[4] );
	DDX_Text ( pDX, IDC_q1, q[0] );
	DDV_MinMaxDouble( pDX, q[0], 0, 1 );
	DDX_Text ( pDX, IDC_q2, q[1] );
	DDV_MinMaxDouble( pDX, q[1], 0, 1 );
	DDX_Text ( pDX, IDC_q3, q[2] );
	DDV_MinMaxDouble( pDX, q[2], 0, 1 );
	DDX_Text ( pDX, IDC_q4, q[3] );
	DDV_MinMaxDouble( pDX, q[3], 0, 1 );
	DDX_Text ( pDX, IDC_q5, q[4] );
	DDV_MinMaxDouble( pDX, q[4], 0, 1 );
	DDX_Text ( pDX, IDC_C11, C[0][0] );
	DDX_Text ( pDX, IDC_C12, C[0][1] );
	DDX_Text ( pDX, IDC_C13, C[0][2] );
	DDX_Text ( pDX, IDC_C14, C[0][3] );
	DDX_Text ( pDX, IDC_C15, C[0][4] );
	DDX_Text ( pDX, IDC_C22, C[1][1] );
	DDX_Text ( pDX, IDC_C23, C[1][2] );
	DDX_Text ( pDX, IDC_C24, C[1][3] );
	DDX_Text ( pDX, IDC_C25, C[1][4] );
	DDX_Text ( pDX, IDC_C33, C[2][2] );
	DDX_Text ( pDX, IDC_C34, C[2][3] );
	DDX_Text ( pDX, IDC_C35, C[2][4] );
	DDX_Text ( pDX, IDC_C44, C[3][3] );
	DDX_Text ( pDX, IDC_C45, C[3][4] );
	DDX_Text ( pDX, IDC_C55, C[4][4] );
	DDX_Check ( pDX, IDC_With, With );
	DDX_Check ( pDX, IDC_Without, Without );
	DDX_Check ( pDX, IDC_Indirect, Indirect );
	DDX_Check ( pDX, IDC_Direct, Direct );
}

MyFrameWindow::MyFrameWindow ( )
{
     N=5,K=120,r=0.1,T=0.5,m=10;
	 for(i=0;i<N;i++){S[i]=100+10*i; q[i]=i/100;}
	 for(i=0;i<N;i++)
	 {
	  for(j=0;j<N;j++)
	  {
	   if(i==j) C[i][j]=0.09;
	   else C[i][j]=0;
	  }
	 }
	 With=0,Without=0,Indirect=0,Direct=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[200];
	CPaintDC dc( this );

	for(i=0;i<M;i++)
	{
 	 for(j=0;j<i;j++)C[i][j]=C[j][i];
    }

	//Calculate A matrix
    double A[M][M];
    Cholesky(A,C,N);
    
	double Z[M],R[M];
	double ST,MaxST;
	double Price,MPrice,SquareSum;
	time_t Start,End;
	double Speed;
	
	double D[M][M],B[M][M],X[M][M];
	double *z[MAX],*w[MAX];
	
	int l=(L*2)*m;
	
	int x_axis=20,y_axis=0;
	
	srand( (unsigned)time( NULL ) );
	
	if(Indirect)
	{
	 _stprintf(info,"%d*%d:",2*L,m);
	 y_axis+=20;
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
	 _stprintf(info,"MDMM");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"Mean");
	 dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"SE");
	 dc.TextOut(x_axis+200,y_axis,info);
	 _stprintf(info,"95%% CI");
	 dc.TextOut(x_axis+300,y_axis,info);
	 _stprintf(info,"Time");
	 dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	 y_axis+=20;
	 
	 if(With)
	 {
 	  for(i=0;i<l;i++)
	  {
	   z[i]=new double[M];
	   w[i]=new double[M]; 
      }

	  Start=time(NULL);
	  Price=0,SquareSum=0;
      for(n=0;n<m;n++)
      {
	   do{
       //Generate Z vector 
	   for(i=0;i<L;i++)
	   {  
        for(j=0;j<N;j++)
	  	{
	   	 z[i][j]=Normal();
	   	 z[i+L][j]=-z[i][j];
  	    }
       } 
     
       //Calculate Var-Cov Matrix
       for(i=0;i<M;i++)for(j=0;j<M;j++)D[i][j]=0;
     
       for(i=0;i<N;i++)
       for(j=0;j<N;j++)
       {
	    for(k=0;k<2*L;k++)
        D[i][j]+=(z[k][i]*z[k][j])/(2*L);
	   }
	 
 	   //Calculate B matrix
	   Cholesky(B,D,N);

	   GJ(B,A,X,N);

	   //Generate N vector
	   for(k=0;k<2*L;k++)
	   {
        for(j=0;j<N;j++)
        {
         w[k][j]=0;
         for(i=0;i<N;i++)w[k][j]+=z[k][i]*X[i][j];
	    }
       }
	 
	   Sum=0;
	   for(k=0;k<2*L;k++)
	   {
  	    //Calculate payoff at time T
  	    MaxST=0;
        for(i=0;i<N;i++)
        {
	     ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*w[k][i]);
	     MaxST=max(MaxST,ST);
	    }
	    Sum+=exp(-r*T)*max(MaxST-K,0);
	   }
	 
	   MPrice=Sum/(2*L);
	   }while(MPrice==0);
	   Price+=MPrice;
	   SquareSum+=pow(MPrice/sqrt((double)(m-1)),2);
      }
	  End=time(NULL);
	  Speed=difftime(End,Start);

	  Price/=m;
      STDERR=sqrt(SquareSum/m-pow(Price,2)/(m-1));
	  Lower_Limit=Price-2*STDERR;
      Upper_Limit=Price+2*STDERR;
      
	  _stprintf(info,"o");
	  dc.TextOut(x_axis,y_axis,info);
	  _stprintf(info,"%.4lf", Price);
	  dc.TextOut(x_axis+100,y_axis,info);
	  _stprintf(info,"%.4lf", STDERR);
	  dc.TextOut(x_axis+200,y_axis,info);
	  _stprintf(info,"[%.4lf,%.4lf]", Lower_Limit, Upper_Limit);
	  dc.TextOut(x_axis+300,y_axis,info);
	  _stprintf(info,"%d", (int)(Speed));
	  dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	  y_axis+=20;
	  
 	  for(i=0;i<l;i++)
	  {
	   delete[] z[i];
	   z[i]=NULL;
	   delete[] w[i];
	   w[i]=NULL;
      }
	 }
	 
	 if(Without)
	 {
      Start=time(NULL);
	  Price=0,SquareSum=0; 
      for(n=0;n<m;n++)
      {
       do{
	   Sum=0;
	   for(j=0;j<L;j++)
	   {
        //Generate Z vector and N vector    
        for(i=0;i<N;i++)Z[i]=Normal();
        for(i=0;i<N;i++)
        {
	     R[i]=0;
  	     for(k=0;k<N;k++)R[i]+=Z[k]*A[k][i];
	    }

	    //Calculate payoff at time T
	    MaxST=0;
        for(i=0;i<N;i++)
	    {
	     ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*R[i]);
	     MaxST=max(MaxST,ST);
	    }
	    Sum+=exp(-r*T)*max(MaxST-K,0);

	    //Antithetic variate approach   
        for(i=0;i<N;i++)
        {
	     R[i]=0;
  	     for(k=0;k<N;k++)R[i]+=(-Z[k])*A[k][i];
	    }
	    MaxST=0;
        for(i=0;i<N;i++)
	    {
	     ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*R[i]);
	     MaxST=max(MaxST,ST);
	    }
	    Sum+=exp(-r*T)*max(MaxST-K,0);
	   }
	   MPrice=Sum/(2*L);	
	   }while((_isnan(MPrice)!=0)||(_finite(MPrice)==0)); 

	   Price+=MPrice;
	   SquareSum+=pow(MPrice/sqrt((double)(m-1)),2);
	  }
	  Price/=m;
	
	  End=time(NULL);
	  Speed=difftime(End,Start);

      STDERR=sqrt(SquareSum/m-pow(Price,2)/(m-1));
	  Lower_Limit=Price-2*STDERR;
      Upper_Limit=Price+2*STDERR;
      
	  _stprintf(info,"x");
	  dc.TextOut(x_axis,y_axis,info);
	  _stprintf(info,"%.4lf", Price);
	  dc.TextOut(x_axis+100,y_axis,info);
	  _stprintf(info,"%.4lf", STDERR);
	  dc.TextOut(x_axis+200,y_axis,info);
	  _stprintf(info,"[%.4lf,%.4lf]", Lower_Limit, Upper_Limit);
	  dc.TextOut(x_axis+300,y_axis,info);
	  _stprintf(info,"%d", (int)(Speed));
	  dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	  y_axis+=20;
	 }
	}
	
	if(Direct)
	{
	 _stprintf(info,"%d:",(2*L)*m);
	 y_axis+=20;
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
	 _stprintf(info,"MDMM");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"Mean");
	 dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"SE");
	 dc.TextOut(x_axis+200,y_axis,info);
	 _stprintf(info,"95%% CI");
	 dc.TextOut(x_axis+300,y_axis,info);
	 _stprintf(info,"Time");
	 dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	 y_axis+=20;
	 
	 if(With)
	 {
 	  for(i=0;i<l;i++)
	  {
	   z[i]=new double[M];
	   w[i]=new double[M]; 
      }

	  Start=time(NULL);
	
	  do{
	  //Generate Z vector 
      for(i=0;i<(l/2);i++)
      {  
       for(j=0;j<N;j++)
	   {
	    z[i][j]=Normal();
	    z[i+(l/2)][j]=-z[i][j];
       }
      } 
     
      //Calculate Var-Cov Matrix
      for(i=0;i<M;i++)for(j=0;j<M;j++)D[i][j]=0;
    
      for(i=0;i<N;i++)
      for(j=0;j<N;j++)
      {
       for(k=0;k<l;k++)
       D[i][j]+=(z[k][i]*z[k][j])/l;
      }
	 
      //Calculate B matrix
      Cholesky(B,D,N);

      GJ(B,A,X,N);

      //Generate N vector
      for(k=0;k<l;k++)
      {
       for(j=0;j<N;j++)
       {
        w[k][j]=0;
        for(i=0;i<N;i++)w[k][j]+=z[k][i]*X[i][j];
       }
      }
	 
      Price=0,SquareSum=0;
	  for(k=0;k<l;k++)
	  { 
	   //Calculate payoff at time T
	   MaxST=0;
       for(i=0;i<N;i++)
	   {
	    ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*w[k][i]);
	    MaxST=max(MaxST,ST);
	   }
	  
	   MPrice=exp(-r*T)*max(MaxST-K,0);
	   Price+=MPrice;
	   SquareSum+=pow(MPrice/sqrt((double)(l-1)),2);
      }
	  Price/=l;
	  }while(Price==0);	
	  End=time(NULL);
	  Speed=difftime(End,Start);

      STDERR=sqrt(SquareSum/l-pow(Price,2)/(l-1));
	  Lower_Limit=Price-2*STDERR;
      Upper_Limit=Price+2*STDERR;
      
	  _stprintf(info,"o");
	  dc.TextOut(x_axis,y_axis,info);
	  _stprintf(info,"%.4lf", Price);
	  dc.TextOut(x_axis+100,y_axis,info);
	  _stprintf(info,"%.4lf", STDERR);
	  dc.TextOut(x_axis+200,y_axis,info);
	  _stprintf(info,"[%.4lf,%.4lf]", Lower_Limit, Upper_Limit);
	  dc.TextOut(x_axis+300,y_axis,info);
	  _stprintf(info,"%d", (int)(Speed));
	  dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	  y_axis+=20;
	  
 	  for(i=0;i<l;i++)
	  {
	   delete[] z[i];
	   z[i]=NULL;
	   delete[] w[i];
	   w[i]=NULL;
      }
	 }
	 
	 if(Without)
	 {
	  Start=time(NULL);
      do{
      Price=0,SquareSum=0;
	  for(j=0;j<(l/2);j++)
 	  {
  	   //Generate Z vector and N vector    
       for(i=0;i<N;i++)Z[i]=Normal();
       for(i=0;i<N;i++)
       {
	    R[i]=0;
  	    for(k=0;k<N;k++)R[i]+=Z[k]*A[k][i];
	   }
	   //Calculate payoff at time T
	   MaxST=0;
       for(i=0;i<N;i++)
	   {
	    ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*R[i]);
	    MaxST=max(MaxST,ST);
	   }
	   MPrice=exp(-r*T)*max(MaxST-K,0);
	   Price+=MPrice;
	   SquareSum+=pow(MPrice/sqrt((double)(l-1)),2);
	  
	   //Antithetic variate approach   
       for(i=0;i<N;i++)
       {
	    R[i]=0;
  	    for(k=0;k<N;k++)R[i]+=(-Z[k])*A[k][i];
	   }
	   MaxST=0;
       for(i=0;i<N;i++)
	   {
	    ST=S[i]*exp((r-q[i]-C[i][i]/2)*T+sqrt(T)*R[i]);
	    MaxST=max(MaxST,ST);
	   }
	   MPrice=exp(-r*T)*max(MaxST-K,0);
	   Price+=MPrice;
	   SquareSum+=pow(MPrice/sqrt((double)(l-1)),2);
      }
	  Price/=l;
	  }while((_isnan(Price)!=0)||(_finite(Price)==0));
	
	  End=time(NULL);
	  Speed=difftime(End,Start);

      STDERR=sqrt(SquareSum/l-pow(Price,2)/(l-1));
	  Lower_Limit=Price-2*STDERR;
      Upper_Limit=Price+2*STDERR;

	  _stprintf(info,"x");
	  dc.TextOut(x_axis,y_axis,info);
	  _stprintf(info,"%.4lf", Price);
	  dc.TextOut(x_axis+100,y_axis,info);
	  _stprintf(info,"%.4lf", STDERR);
	  dc.TextOut(x_axis+200,y_axis,info);
	  _stprintf(info,"[%.4lf,%.4lf]", Lower_Limit, Upper_Limit);
	  dc.TextOut(x_axis+300,y_axis,info);
	  _stprintf(info,"%d", (int)(Speed));
	  dc.TextOut(x_axis+500,y_axis,info);	 	 	 
	  y_axis+=20;
	 }
	}
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.N = N;
	dlg.r = r;
	dlg.T = T;
	dlg.K = K;
	dlg.m = m;
	for(i=0;i<M;i++){dlg.S[i]=S[i]; dlg.q[i]=q[i];}
	for(i=0;i<M;i++)
	{
	 for(j=0;j<i+1;j++)dlg.C[j][i]=C[j][i];
	}
	dlg.With = With;
	dlg.Without = Without;
	dlg.Indirect = Indirect;
	dlg.Direct = Direct;

	if ( dlg.DoModal ( ) == IDOK )
	{
	N = dlg.N;
	r = dlg.r;
	T = dlg.T;
	K = dlg.K;
	m = dlg.m;
	for(i=0;i<M;i++){S[i]=dlg.S[i]; q[i]=dlg.q[i];}
	for(i=0;i<M;i++)
	{
	 for(j=0;j<i+1;j++)C[j][i]=dlg.C[j][i];
	}
	With = dlg.With;
	Without = dlg.Without;
	Indirect = dlg.Indirect;
	Direct = dlg.Direct;
	Invalidate();
	}
}
