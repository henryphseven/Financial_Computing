#include <afxwin.h>
#include "ddxddv.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
const double PI=3.14159265359;
const int L=5000;

class MyFrameWindow : public CFrameWnd
{
	public:

	 int Black_Scholes,Monte_Carlo,Binomial_Tree,Call,Put;
     double S, K, r, T, Sigma, q, STDERR, Upper_Limit, Lower_Limit;
	 int n,m;
     double Value;

	MyFrameWindow ( );

	afx_msg void OnPaint( );
	afx_msg void OnInput( );
	afx_msg double Standard_Normal_Distribution(double );
	afx_msg double Normal();
	afx_msg double BS_European_Call();
	afx_msg double BS_European_Put();
	afx_msg double Simulation_European_Call();
	afx_msg double Simulation_European_Put();
	afx_msg double Binomial_Tree_European_Call();
	afx_msg double Binomial_Tree_European_Put();
	afx_msg double Binomial_Tree_American_Call();
    afx_msg double Binomial_Tree_American_Put();
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
		: CWinApp ( "Plain Vanilla Option Price" )
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

     int Black_Scholes,Monte_Carlo,Binomial_Tree,Call,Put;
     double S, K, r, T, Sigma, q;
     int n,m;

	MyDialog ( CWnd* parent ) : CDialog ( IDD_INPUT, parent ) {}
	void DoDataExchange( CDataExchange* pDX );
};

void MyDialog::DoDataExchange( CDataExchange* pDX )
{
	CDialog::DoDataExchange ( pDX );
    DDX_Text ( pDX, IDC_S, S );
	DDX_Text ( pDX, IDC_r, r );
	DDV_MinMaxDouble( pDX, r, 0, 1 );
	DDX_Text ( pDX, IDC_T, T );
	DDX_Text ( pDX, IDC_K, K );
	DDX_Text ( pDX, IDC_Sigma, Sigma );
	DDV_MinMaxDouble( pDX, Sigma, 0, 1 );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_n, n );
	DDV_MinMaxInt( pDX, n, 1, 10000 );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, 10000 );
	DDX_Check ( pDX, IDC_Black_Scholes, Black_Scholes );
	DDX_Check ( pDX, IDC_Monte_Carlo, Monte_Carlo );
	DDX_Check ( pDX, IDC_Binomial_Tree, Binomial_Tree );
	DDX_Check ( pDX, IDC_Call, Call );
	DDX_Check ( pDX, IDC_Put, Put );
}

MyFrameWindow::MyFrameWindow ( )
{
	S=50, K=50, r=0.1, T=0.5, Sigma=0.3, q=0, n=1,m=10;
	Black_Scholes=0,Monte_Carlo=0,Binomial_Tree=0,Call=0,Put=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );
	
	int x_axis=20,y_axis=0;
	
	if(Black_Scholes)
	{
     _stprintf(info,"Black-Scholes:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Option Type");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
     y_axis+=20;
     
     if(Call)
     {
      Value=BS_European_Call();
      
      _stprintf(info,"Call");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      y_axis+=20;
     }
     
     if(Put)
     {
      Value=BS_European_Put();
      
      _stprintf(info,"Put");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      y_axis+=20;
     }
    }
    
	if(Monte_Carlo)
	{
     _stprintf(info,"Monte Carlo:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Option Type");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"SE");
     dc.TextOut(x_axis+200,y_axis,info);
     _stprintf(info,"95%% CI");
     dc.TextOut(x_axis+300,y_axis,info);
     y_axis+=20;
     
     if(Call)
     {
      Value=Simulation_European_Call();
      Lower_Limit=Value-2*STDERR;
      Upper_Limit=Value+2*STDERR;
      
      _stprintf(info,"Call");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      _stprintf(info,"%.4lf",STDERR);
      dc.TextOut(x_axis+200,y_axis,info);
      _stprintf(info,"[%.4lf,%.4lf]",Lower_Limit,Upper_Limit);
      dc.TextOut(x_axis+300,y_axis,info);
      y_axis+=20;
     }
     
     if(Put)
     {
      Value=Simulation_European_Put();
      Lower_Limit=Value-2*STDERR;
      Upper_Limit=Value+2*STDERR;
      
      _stprintf(info,"Put");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      _stprintf(info,"%.4lf",STDERR);
      dc.TextOut(x_axis+200,y_axis,info);
      _stprintf(info,"[%.4lf,%.4lf]",Lower_Limit,Upper_Limit);
      dc.TextOut(x_axis+300,y_axis,info);
      y_axis+=20;
     }
    }
     
	if(Binomial_Tree)
	{
     _stprintf(info,"Binomial Tree:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Option Type");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"American");
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;
     
     if(Call)
     {
      _stprintf(info,"Call");
      dc.TextOut(x_axis,y_axis,info);
      Value=Binomial_Tree_European_Call();
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      Value=Binomial_Tree_American_Call();
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+200,y_axis,info);
      y_axis+=20;
     }
     
     if(Put)
     {
      _stprintf(info,"Put");
      dc.TextOut(x_axis,y_axis,info);
      Value=Binomial_Tree_European_Put();
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+100,y_axis,info);
      Value=Binomial_Tree_American_Put();
      _stprintf(info,"%.4lf",Value);
      dc.TextOut(x_axis+200,y_axis,info);
      y_axis+=20;
     }
    }	
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S = S;
	dlg.r = r;
	dlg.T = T;
    dlg.K = K;
	dlg.Sigma = Sigma;
	dlg.q = q;
	dlg.n = n;
	dlg.m = m;
	dlg.Black_Scholes = Black_Scholes;
	dlg.Monte_Carlo = Monte_Carlo;
	dlg.Binomial_Tree = Binomial_Tree;
	dlg.Call = Call;
	dlg.Put = Put;

	if ( dlg.DoModal ( ) == IDOK )
	{
	S=dlg.S;
	r=dlg.r;
	T=dlg.T;
    K=dlg.K;
	Sigma=dlg.Sigma;
	q=dlg.q;
	n=dlg.n;
	m=dlg.m;
	Black_Scholes = dlg.Black_Scholes;
	Monte_Carlo = dlg.Monte_Carlo;
	Binomial_Tree = dlg.Binomial_Tree;
	Call = dlg.Call;
	Put = dlg.Put;
	Invalidate();
	}
}

//Calculate cumulative probability of standard normal distribution
double MyFrameWindow::Standard_Normal_Distribution(double d)   
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

//Generate a standard normal random variable
double MyFrameWindow::Normal()
{
     double Y,X[2];
     int j;
	 for(j=0;j<2;j++)
	 {
	   X[j]=double(rand())/RAND_MAX;
	 }
     Y=sqrt((-2)*log(X[0]))*cos(2*PI*X[1]);
	 return Y;
}

//The first method: Using Black-Scholes formula
double MyFrameWindow::BS_European_Call()
{
   double d1=(log(S/K)+(r-q+Sigma*Sigma/2)*T)/(Sigma*sqrt(T));
   double d2=d1-Sigma*sqrt(T);
   double European_CallValue=S*exp(-q*T)*Standard_Normal_Distribution(d1)
   -K*exp(-r*T)*Standard_Normal_Distribution(d2);
   return European_CallValue;
}

double MyFrameWindow::BS_European_Put()
{
   double d1=(log(S/K)+(r-q+Sigma*Sigma/2)*T)/(Sigma*sqrt(T));
   double d2=d1-Sigma*sqrt(T);
   double European_PutValue=K*exp(-r*T)*Standard_Normal_Distribution(-d2)
   -S*exp(-q*T)*Standard_Normal_Distribution(-d1);
   return European_PutValue;
}

//The second method: Using Monte Carlo simulation
double MyFrameWindow::Simulation_European_Call()
{
   srand( (unsigned)time( NULL ) );
   double N, ST, V, SumV, MV, European_CallValue, SquareSum;
   double AN, AST, AV, ASumV;
   do{
   European_CallValue=0, SquareSum=0;
   for(int j=0;j<m;j++)
   {
   SumV=0, ASumV=0, MV=0;
   for(int i=0;i<L;i++)
   {
    N=Normal();
    AN=-N;
	ST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*N);
	AST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*AN);
    V=exp(-r*T)*max(ST-K,0);
    AV=exp(-r*T)*max(AST-K,0);
    SumV+=V;
    ASumV+=AV;
   }
   MV=(SumV+ASumV)/(2*L);
   European_CallValue+=MV;
   SquareSum+=pow(MV/sqrt((double)(m-1)),2);
   }
   European_CallValue/=m;
   } while((_isnan(European_CallValue)!=0)||(_finite(European_CallValue)==0));
   STDERR=sqrt(SquareSum/m-pow(European_CallValue,2)/(m-1));
   return European_CallValue;
}

double MyFrameWindow::Simulation_European_Put()
{
   srand( (unsigned)time( NULL ) );
   double N, ST, V, SumV, MV, European_PutValue, SquareSum;
   double AN, AST, AV, ASumV;
   do{
   European_PutValue=0, SquareSum=0;
   for(int j=0;j<m;j++)
   {
   SumV=0, ASumV=0, MV=0;
   for(int i=0;i<L;i++)
   {
    N=Normal();
    AN=-N;
	ST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*N);
	AST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*AN);
    V=exp(-r*T)*max(K-ST,0);
    AV=exp(-r*T)*max(K-AST,0);
    SumV+=V;
    ASumV+=AV;
   }
   MV=(SumV+ASumV)/(2*L);
   European_PutValue+=MV;
   SquareSum+=pow(MV/sqrt((double)(m-1)),2);
   }
   European_PutValue/=m;
   } while((_isnan(European_PutValue)!=0)||(_finite(European_PutValue)==0));
   STDERR=sqrt(SquareSum/m-pow(European_PutValue,2)/(m-1));
   return European_PutValue;
}

//The third method: Using binomial tree 
double MyFrameWindow::Binomial_Tree_European_Call()
{
  //Calculate parameters
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  //Calculate possible stock prices at T and store them in Array[1001]
  double Array[1001];
  double CurrentS=S*pow(U,n);
  int i;
  for(i=0;i<=n;i=i+1)
  {
   Array[i]=max(CurrentS-K,0);
   CurrentS=CurrentS*D*D;
  }
  //Backward Induction
  for(i=n-1;i>=0;i=i-1)
  {
   for(int j=0;j<=i;j++)
   {
    Array[j]=exp(-r*DeltaT)*(Array[j]*Pu+Array[j+1]*Pd);
   }  
  }
  double European_CallValue=Array[0];
  return European_CallValue;
}

double MyFrameWindow::Binomial_Tree_European_Put()
{
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  double Array[1001];
  double CurrentS=S*pow(U,n);
  int i;
  for(i=0;i<=n;i=i+1)
  {
   Array[i]=max(K-CurrentS,0);
   CurrentS=CurrentS*D*D;
  }
  for(i=n-1;i>=0;i=i-1)
  {
   for(int j=0;j<=i;j++)
   {
    Array[j]=exp(-r*DeltaT)*(Array[j]*Pu+Array[j+1]*Pd);
   }  
  }
  double European_PutValue=Array[0];
  return European_PutValue;       
}

double MyFrameWindow::Binomial_Tree_American_Call()
{
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  double Array[1001];
  double CurrentS=S*pow(U,n);
  int i;
  for(i=0;i<=n;i=i+1)
  {
   Array[i]=max(CurrentS-K,0);
   CurrentS=CurrentS*D*D;
  }
  for(i=n-1;i>=0;i=i-1) 
  {
   CurrentS=S*pow(U,i);
   for(int j=0;j<=i;j++)
   {
           //Compare the values at each node and take the larger 
           Array[j]=max(exp(-r*DeltaT)*(Array[j]*Pu+Array[j+1]*Pd), CurrentS-K);
           CurrentS=CurrentS*D*D;
   }
  }
  double American_CallValue=Array[0];
  return American_CallValue;
}

double MyFrameWindow::Binomial_Tree_American_Put()
{
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  double Array[1001];
  double CurrentS=S*pow(U,n);
  int i;
  for(i=0;i<=n;i=i+1)
  {
   Array[i]=max(K-CurrentS,0);
   CurrentS=CurrentS*D*D;
  }
  for(i=n-1;i>=0;i=i-1) 
  {
   CurrentS=S*pow(U,i);
   for(int j=0;j<=i;j++)
   {
           Array[j]=max(exp(-r*DeltaT)*(Array[j]*Pu+Array[j+1]*Pd), K-CurrentS);
           CurrentS=CurrentS*D*D;
   }
  }
  double American_PutValue=Array[0];
  return American_PutValue;
}
