#include <afxwin.h>
#include "ddxddv.h"
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
double Sum;
const double PI=3.14159265359;
const int MAX=500;
const int L=10000;
const int sMAX=1000;
//standard normal distribution's pdf
double Standard_Normal_Distribution(double d)
{
   int flag=0;   //Flag =1 if d<0
   if(d<0)
     {
      flag=1;
      d=fabs(d);
     }
   double rr=0.2316419;
   double a1=0.31938153;
   double a2=-0.356563782;
   double a3=1.781477937;
   double a4=-1.821255978;
   double a5=1.330274429;
   double k=1/(1+d*rr);
   double value=1-exp(d*d/(-2))*(a1*k+a2*pow(k,2)+a3*pow(k,3)+a4*pow(k,4)
                  +a5*pow(k,5))/sqrt(2*PI);
   if(flag) return 1-value;
    else return value;
}
//European lookback put's B-S formula
double LookbackPut(double S, double T, double Smax, double r, double q, double Sigma)
{
  double OptionValue;
  double D_1 = (log(Smax/S)+( q-r+(Sigma*Sigma)*0.5)*T)/(Sigma*sqrt(T));
  double D_2 = (log(Smax/S)+( r-q-(Sigma*Sigma)*0.5)*T)/(Sigma*sqrt(T));
  double y = 2.0*log(Smax/S)*( r-q-(Sigma*Sigma)*0.5)/(Sigma*Sigma);
  OptionValue=
    Smax*(Standard_Normal_Distribution(D_1) - (Sigma*Sigma)/(2.0*(r-q))*exp(y)* Standard_Normal_Distribution(-D_2))/exp(r*T) +
    S*Standard_Normal_Distribution(-D_1+ Sigma*sqrt(T)) *
    ((Sigma*Sigma)/(2.0*(r-q)))/exp(q*T) -
    S*Standard_Normal_Distribution(D_1- Sigma*sqrt(T))/exp(q*T);
  return OptionValue;
}
//Generate a standard normal random variable
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
class MyFrameWindow : public CFrameWnd
{
	public:

	 double S,r,T,Sigma,q,Smax_0;
	 int m,n,s;
	 int Black_Scholes,Monte_Carlo,Binomial_Tree,Continuous,Discrete;

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
		: CWinApp ( "Lookback Put Price" )
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

	 double S,r,T,Sigma,q,Smax_0;
	 int m,n,s;
	 int Black_Scholes,Monte_Carlo,Binomial_Tree,Continuous,Discrete;

	MyDialog ( CWnd* parent ) : CDialog ( IDD_INPUT, parent ) {}
	void DoDataExchange( CDataExchange* pDX );
};

void MyDialog::DoDataExchange( CDataExchange* pDX )
{
	CDialog::DoDataExchange ( pDX );
    DDX_Text ( pDX, IDC_S, S );
	DDX_Text ( pDX, IDC_Smax_0, Smax_0 );
	DDX_Text ( pDX, IDC_r, r );
	DDV_MinMaxDouble( pDX, r, 0, 1 );
	DDX_Text ( pDX, IDC_T, T );
	DDX_Text ( pDX, IDC_Sigma, Sigma );
	DDV_MinMaxDouble( pDX, Sigma, 0, 1 );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, MAX );
	DDX_Text ( pDX, IDC_n, n );
	DDV_MinMaxInt( pDX, n, 1, MAX );
	DDX_Text ( pDX, IDC_n2, s );
	DDV_MinMaxInt( pDX, s, 1, n );
	DDX_Check ( pDX, IDC_Black_Scholes, Black_Scholes );
	DDX_Check ( pDX, IDC_Monte_Carlo, Monte_Carlo );
	DDX_Check ( pDX, IDC_Binomial_Tree, Binomial_Tree );
	DDX_Check ( pDX, IDC_Continuous, Continuous );
	DDX_Check ( pDX, IDC_Discrete, Discrete );
}

MyFrameWindow::MyFrameWindow ( )
{
	 S=50,r=0.1,T=0.25,Sigma=0.4,q=0,m=10,n=100,s=50; 
	 Smax_0=S;
 	 Black_Scholes=0,Monte_Carlo=0,Binomial_Tree=0,Continuous=0,Discrete=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );

		int i,j,k,l;
		
		double E_Put,Square_E_Put,St,Smax;
		double STDERR,Upper_Limit,Lower_Limit; 
	
		int div;
		double u, d, pu, pd, delta_t;
   		int flag,flag2,pU;
		double E_Payoffu,E_Payoffd,A_Payoffu,A_Payoffd;

	int x_axis=20,y_axis=0;
	
	if(Black_Scholes)
	{
	 E_Put=LookbackPut(S,T,Smax_0,r,q,Sigma);
     _stprintf(info,"Black-Scholes (European): %.4lf",E_Put);
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
    }
    
	if(Monte_Carlo)
	{
     srand( (unsigned)time( NULL ) );
     
   	 delta_t=T/s;
	 E_Put=0,Square_E_Put=0;
     for(i=0;i<m;i++)
     {
	  Sum=0;
	  for(j=0;j<L;j++)
	  {
	   St=S,Smax=Smax_0;
	   for(k=1;k<=s;k++)
	   {
	    St=St*exp((r-q-Sigma*Sigma/2)*delta_t+Sigma*sqrt(delta_t)*Normal());
	    if(St>Smax)Smax=St;
	   }
	   Sum+=exp(-r*T)*max(Smax-St,0);
	  }
	  E_Put+=Sum/L;
	  Square_E_Put+=pow(Sum/L,2);
	 }
     
     E_Put/=m;
     STDERR=sqrt((double)((Square_E_Put-m*pow(E_Put,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit=E_Put-2*STDERR;
     Upper_Limit=E_Put+2*STDERR;	

     _stprintf(info,"Monte Carlo (European): %.4lf with SE= %.4lf and 95%% CI= [%.4lf,%.4lf]"
     ,E_Put,STDERR,Lower_Limit,Upper_Limit);
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
    }
     
	if(Binomial_Tree)
	{
     _stprintf(info,"Binomial Tree:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Sampling");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"American");
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;
     
   	 div=n/s;
	 delta_t=T/n;
  	 u=exp(Sigma*sqrt(delta_t));
  	 d=exp(-Sigma*sqrt(delta_t));
  	 pu=(exp((r-q)*delta_t)-d)/(u-d);
  	 pd=1-pu;

  	 //S is the current price at (i,j), and p is # of Smax
  	 struct node
  	 {
	  		double S;
	  		double *Smax;
	  		double *E_Put;
	  		double *A_Put;
	  		int p;
	 } *tree[MAX+1];
	 
 	 for(i=0;i<=n;i++)
	 tree[i]=new struct node[MAX+1];
	 
	 for(i=0;i<=n;i++)
  	 {
	  for(j=0;j<=i;j++)
	  {
	   tree[i][j].Smax=new double[2*MAX+1];
	   tree[i][j].E_Put=new double[2*MAX+1];
	   tree[i][j].A_Put=new double[2*MAX+1];
	  }
     }
     
	 //set initial value
	 tree[0][0].S=S;
	 tree[0][0].Smax[0]=Smax_0;
	 tree[0][0].p=1;
	 
	 //calculate St at each node
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       //in case that exactly the same values differ
       if(i-j>j)tree[i][j].S=S*pow(u,(i-j)-j);
       else if((i-j)==j)tree[i][j].S=S;
       else tree[i][j].S=S*pow(d,j-(i-j));
	  }
   	 }
   	 
     if(Continuous)
     {
   	  //inherit Smax
  	  for(i=1;i<=n;i=i+1)
  	  {	 
	   if(tree[i-1][0].Smax[0]>tree[i][0].S)
	   tree[i][0].Smax[0]=tree[i-1][0].Smax[0];
	   else tree[i][0].Smax[0]=tree[i][0].S;
	   tree[i][0].p=1;
	  
	   for(j=1;j<i;j++)
	   {
	    tree[i][j].p=0;
	    flag2=0;
	   
	    //inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
	     if(tree[i-1][j-1].Smax[k]>tree[i][j].S)
	     {
		  tree[i][j].Smax[tree[i][j].p]=tree[i-1][j-1].Smax[k];
	      tree[i][j].p++;
	     }
	     else flag2=1;
        }
	   
	    pU=tree[i][j].p;

	    //inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
	     if(tree[i-1][j].Smax[k]>tree[i][j].S)
	     {
          flag=0;
          for(l=0;l<pU;l++)
          {
           if(tree[i-1][j].Smax[k]==tree[i][j].Smax[l]) flag=1;
		  }
	      if(flag==0)
		  {
	       tree[i][j].Smax[tree[i][j].p]=tree[i-1][j].Smax[k];
	       tree[i][j].p++;
	      }
         }
         else flag2=1;
	    }
	   
	    if(flag2==1)
        {
         tree[i][j].Smax[tree[i][j].p]=tree[i][j].S;
         tree[i][j].p++;
        }
       } 
	  	 
	   if(tree[i-1][i-1].Smax[0]>tree[i][i].S)
	   tree[i][i].Smax[0]=tree[i-1][i-1].Smax[0];
	   else tree[i][i].Smax[0]=tree[i][i].S;
	   tree[i][i].p=1;
   	  }
     
	  //decide each Smax's payoff for each terminal node
	  for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	  {
	   tree[n][j].E_Put[k]=tree[n][j].Smax[k]-tree[n][j].S;
	   tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	  }
	 
	  //backward induction
  	  for(i=n-1;i>=0;i=i-1)
  	  {
   	   for(j=0;j<=i;j++)
   	   {
        for(k=0;k<tree[i][j].p;k++)
        {
         //decide Payoffu
         for(l=0;l<tree[i+1][j].p;l++)
         {
	      if(tree[i][j].Smax[k]==tree[i+1][j].Smax[l])
	      {
	       E_Payoffu=tree[i+1][j].E_Put[l];
	       A_Payoffu=tree[i+1][j].A_Put[l];
	       break;
          }
		 }
		 if(l==tree[i+1][j].p)
		 {
  		  for(l=0;l<tree[i+1][j].p;l++)
          {
	       if(tree[i+1][j].S==tree[i+1][j].Smax[l])
	       {
	        E_Payoffu=tree[i+1][j].E_Put[l];
	        A_Payoffu=tree[i+1][j].A_Put[l];
	        break;
           }
          }
		 }

         //decide Payoffd
         for(l=0;l<tree[i+1][j+1].p;l++)
         {
	      if(tree[i][j].Smax[k]==tree[i+1][j+1].Smax[l])
	      {
	       E_Payoffd=tree[i+1][j+1].E_Put[l];
	       A_Payoffd=tree[i+1][j+1].A_Put[l];
	       break;
          }
		 }
		 if(l==tree[i+1][j+1].p)
		 {
  		  for(l=0;l<tree[i+1][j+1].p;l++)
          {
	       if(tree[i+1][j+1].S==tree[i+1][j+1].Smax[l])
	       {
	        E_Payoffd=tree[i+1][j+1].E_Put[l];
	        A_Payoffd=tree[i+1][j+1].A_Put[l];
	        break;
           }
          }
		 }
        
         //calculate put price
         tree[i][j].E_Put[k]=exp(-r*delta_t)*(E_Payoffu*pu+E_Payoffd*pd);
         tree[i][j].A_Put[k]=
		 max(exp(-r*delta_t)*(A_Payoffu*pu+A_Payoffd*pd),
		 tree[i][j].Smax[k]-tree[i][j].S);
	    }
       }  
      }
   
      _stprintf(info,"Continuous");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",tree[0][0].E_Put[0]);
      dc.TextOut(x_axis+100,y_axis,info);
      _stprintf(info,"%.4lf",tree[0][0].A_Put[0]);
      dc.TextOut(x_axis+200,y_axis,info);
      y_axis+=20;
     }
     
     if(Discrete)
     {
 	  //inherit Smax
  	  for(i=1;i<=n;i=i+1)
  	  {
	   //update						 
	   if(i%div==0){
	   if(tree[i-1][0].Smax[0]>tree[i][0].S)
	   tree[i][0].Smax[0]=tree[i-1][0].Smax[0];
	   else tree[i][0].Smax[0]=tree[i][0].S;
	   tree[i][0].p=1;
	  
	   for(j=1;j<i;j++)
	   {
	    tree[i][j].p=0;
	    flag2=0;
	   
	    //inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
	     if(tree[i-1][j-1].Smax[k]>tree[i][j].S)
	     {
		  tree[i][j].Smax[tree[i][j].p]=tree[i-1][j-1].Smax[k];
	      tree[i][j].p++;
	     }
	     else flag2=1;
	    }
	   
	    pU=tree[i][j].p;

	    //inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
	     if(tree[i-1][j].Smax[k]>tree[i][j].S)
	     {
          flag=0;
          for(l=0;l<pU;l++)
          {
           if(tree[i-1][j].Smax[k]==tree[i][j].Smax[l])flag=1;
          }
		  if(flag==0)
		  {
	       tree[i][j].Smax[tree[i][j].p]=tree[i-1][j].Smax[k];
	       tree[i][j].p++;
	      }
         }
         else flag2=1;
	    }
	   
	    if(flag2==1)
        {
         tree[i][j].Smax[tree[i][j].p]=tree[i][j].S;
         tree[i][j].p++;
        }
       } 
	  	 
	   if(tree[i-1][i-1].Smax[0]>tree[i][i].S)
	   tree[i][i].Smax[0]=tree[i-1][i-1].Smax[0];
	   else tree[i][i].Smax[0]=tree[i][i].S;
	   tree[i][i].p=1;
      }
	  
	  //wait
	  else
	  {
	  tree[i][0].Smax[0]=tree[i-1][0].Smax[0];
	  tree[i][0].p=1;
	  
	  for(j=1;j<i;j++)
	  {
	   tree[i][j].p=0;
	   
	   //inherit from the upper
	   for(k=0;k<tree[i-1][j-1].p;k++)
	   {
		 tree[i][j].Smax[tree[i][j].p]=tree[i-1][j-1].Smax[k];
	     tree[i][j].p++;
	   }
	   
	   pU=tree[i][j].p;

	   //inherit from the lower
	   for(k=0;k<tree[i-1][j].p;k++)
	   {
         flag=0;
         for(l=0;l<pU;l++)
         {
          if(tree[i-1][j].Smax[k]==tree[i][j].Smax[l])flag=1;
		 }
		 if(flag==0)
		 {
	      tree[i][j].Smax[tree[i][j].p]=tree[i-1][j].Smax[k];
	      tree[i][j].p++;
	     } 
	    }
	   }
	  
	   tree[i][i].Smax[0]=tree[i-1][i-1].Smax[0];
	   tree[i][i].p=1;
	   }
   	  }
     
	  //decide each Smax's payoff for each terminal node
	  for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	  {
	   tree[n][j].E_Put[k]=tree[n][j].Smax[k]-tree[n][j].S;
	   tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	  }
	 
	  //backward induction
  	  for(i=n-1;i>=0;i=i-1)
  	  {
   	   for(j=0;j<=i;j++)
   	   {
        for(k=0;k<tree[i][j].p;k++)
        {
         //decide Payoffu
         for(l=0;l<tree[i+1][j].p;l++)
         {
	      if(tree[i][j].Smax[k]==tree[i+1][j].Smax[l])
	      {
	       E_Payoffu=tree[i+1][j].E_Put[l];
	       A_Payoffu=tree[i+1][j].A_Put[l];
	       break;
          }
		 }
		 if(l==tree[i+1][j].p)
		 {
  		  for(l=0;l<tree[i+1][j].p;l++)
          {
	       if(tree[i+1][j].S==tree[i+1][j].Smax[l])
	       {
	        E_Payoffu=tree[i+1][j].E_Put[l];
	        A_Payoffu=tree[i+1][j].A_Put[l];
	        break;
           }
          }
		 }

         //decide Payoffd
         for(l=0;l<tree[i+1][j+1].p;l++)
         {
	      if(tree[i][j].Smax[k]==tree[i+1][j+1].Smax[l])
	      {
	       E_Payoffd=tree[i+1][j+1].E_Put[l];
	       A_Payoffd=tree[i+1][j+1].A_Put[l];
	       break;
          }
		 }
		 if(l==tree[i+1][j+1].p)
		 {
  		  for(l=0;l<tree[i+1][j+1].p;l++)
          {
	       if(tree[i+1][j+1].S==tree[i+1][j+1].Smax[l])
	       {
	        E_Payoffd=tree[i+1][j+1].E_Put[l];
	        A_Payoffd=tree[i+1][j+1].A_Put[l];
	        break;
           }
          }
		 }
        
         //calculate put price
         tree[i][j].E_Put[k]=exp(-r*delta_t)*(E_Payoffu*pu+E_Payoffd*pd);
		 tree[i][j].A_Put[k]
		 =max(exp(-r*delta_t)*(A_Payoffu*pu+A_Payoffd*pd),
		 tree[i][j].Smax[k]-tree[i][j].S);
	    }
       }
      }

      _stprintf(info,"Discrete");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",tree[0][0].E_Put[0]);
      dc.TextOut(x_axis+100,y_axis,info);
      _stprintf(info,"%.4lf",tree[0][0].A_Put[0]);
      dc.TextOut(x_axis+200,y_axis,info);
      y_axis+=20;
     }
     
   	 //release occupied memory
	 for(i=0;i<=n;i++)
  	 {
	  for(j=0;j<=i;j++)
	  {
 	   delete[] tree[i][j].Smax;
	   tree[i][j].Smax=NULL;
	   delete[] tree[i][j].E_Put;
	   tree[i][j].E_Put=NULL;
	   delete[] tree[i][j].A_Put;
	   tree[i][j].A_Put=NULL;
	  }
     }
     
	 for(i=0;i<=n;i++)
	 {
  	  delete[] tree[i];
  	  tree[i]=NULL;
	 }
    }	
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S = S;
	dlg.Smax_0 = Smax_0;
	dlg.r = r;
	dlg.T = T;
	dlg.Sigma = Sigma;
	dlg.q = q;
	dlg.m = m;
	dlg.n = n;
	dlg.s = s;
	dlg.Black_Scholes = Black_Scholes;
	dlg.Monte_Carlo = Monte_Carlo;
	dlg.Binomial_Tree = Binomial_Tree;
	dlg.Continuous = Continuous;
	dlg.Discrete = Discrete;

	if ( dlg.DoModal ( ) == IDOK )
	{
	S=dlg.S;
	Smax_0=dlg.Smax_0;
	r=dlg.r;
	T=dlg.T;
	Sigma=dlg.Sigma;
	q=dlg.q;
	m=dlg.m;
	n=dlg.n;
	s=dlg.s;
	Black_Scholes = dlg.Black_Scholes;
	Monte_Carlo = dlg.Monte_Carlo;
	Binomial_Tree = dlg.Binomial_Tree;
	Continuous = dlg.Continuous;
	Discrete = dlg.Discrete;
	Invalidate();
	}
}
