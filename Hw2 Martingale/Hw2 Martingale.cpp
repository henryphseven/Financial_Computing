#include <afxwin.h>
#include "ddxddv.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
const double PI=3.14159265359;
const int L=5000;
//Calculate cumulative probability of standard normal distribution
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
//Generate a standard normal random variable
double Normal()
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
//Black-Scholes formula
double BS_European_Call(double S,double K,double r,double T,double Sigma,double q)
{
   double d1=(log(S/K)+(r-q+Sigma*Sigma/2)*T)/(Sigma*sqrt(T));
   double d2=d1-Sigma*sqrt(T);
   double European_CallValue=S*exp(-q*T)*Standard_Normal_Distribution(d1)
   -K*exp(-r*T)*Standard_Normal_Distribution(d2);
   return European_CallValue;
}
//The first method: Using formula
double Formula(double S,double K[],double r,double T,double Sigma,double q)
{ 
   return BS_European_Call(S, K[1], r, T, Sigma, q)-BS_European_Call(S, K[2], r, T, Sigma, q)
   -BS_European_Call(S, K[3], r, T, Sigma, q)+BS_European_Call(S, K[4], r, T, Sigma, q);
}
//Payoff structure
double Payoff(double S, double K[])
{
       double P;
       if((K[1]<S)&&(S<K[2]))
       {P=S-K[1];}
       else if((K[2]<S)&&(S<K[3]))
       {P=K[2]-K[1];}
       else if((K[3]<S)&&(S<K[4]))
       {P=K[4]-S;}
       else {P=0;}
       return P;
}
//The second method: Using Monte Carlo simulation
double Simulation(double S,double K[],double r,double T,double Sigma,double q,int m,double &STDERR)
{
   srand( (unsigned)time( NULL ) );
   double N, ST, V, SumV, MV, Value=0, SquareSum=0;
   double AN, AST, AV, ASumV;
   for(int j=0;j<m;j++)
   {
   SumV=0, ASumV=0, MV=0;
   for(int i=0;i<L;i++)
   {
    N=Normal();
    AN=-N;
	ST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*N);
	AST=S*exp((r-q-Sigma*Sigma/2)*T+Sigma*sqrt(T)*AN);
    V=Payoff(ST, K);
    AV=Payoff(AST, K);
    SumV+=V;
    ASumV+=AV;
   }
   MV=(SumV+ASumV)/(2*L);
   Value+=MV;
   SquareSum+=pow(MV/sqrt((double)(m-1)),2);
   }
   Value=exp(-r*T)*(Value/m);
   STDERR=sqrt(SquareSum/m-pow(Value,2)/(m-1));
   return Value;
}
class MyFrameWindow : public CFrameWnd
{
	public:

     double S, r, T, Sigma, q, STDERR, Upper_Limit, Lower_Limit;
     double K[5];
     int m,Black_Scholes,Monte_Carlo;
     double Value;

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
		: CWinApp ( "Derivative Price" )
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

    double S, r, T, Sigma, q;
    double K[5];
    int m,Black_Scholes,Monte_Carlo;

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
	DDX_Text ( pDX, IDC_Sigma, Sigma );
	DDV_MinMaxDouble( pDX, Sigma, 0, 1 );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, 10000 );
	DDX_Check ( pDX, IDC_Black_Scholes, Black_Scholes );
	DDX_Check ( pDX, IDC_Monte_Carlo, Monte_Carlo );	
	DDX_Text ( pDX, IDC_K1, K[1] );
	DDX_Text ( pDX, IDC_K2, K[2] );
	DDX_Text ( pDX, IDC_K3, K[3] );
}

MyFrameWindow::MyFrameWindow ( )
{
    S=20, r=0.1, T=0.5, Sigma=0.3, q=0;
	K[1]=10, K[2]=20, K[3]=30;
	m=10,Black_Scholes=0,Monte_Carlo=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );
	
	int x_axis=20,y_axis=20;
	
	_stprintf(info,"<Martingale Pricing Method>");
	dc.TextOut(x_axis,y_axis,info);
	y_axis+=20;
	
	K[4]=K[3]+(K[2]-K[1]);
	
	if(Black_Scholes)
	{
     Value=Formula(S, K, r, T, Sigma, q);
     _stprintf(info,"Black-Scholes: %.4lf",Value);
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
    }
    
    if(Monte_Carlo)
    {
     Value=Simulation(S, K, r, T, Sigma, q, m, STDERR);
     Lower_Limit=Value-2*STDERR;
     Upper_Limit=Value+2*STDERR;
     _stprintf(info,"Monte Carlo: %.4lf with SE= %.4lf and 95%% CI= [%.4lf,%.4lf]"
     ,Value,STDERR,Lower_Limit,Upper_Limit);
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
    }
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S = S;
	dlg.r = r;
	dlg.T = T;
	dlg.Sigma = Sigma;
	dlg.q = q;
	dlg.m = m;
	dlg.Black_Scholes=Black_Scholes;
    dlg.Monte_Carlo=Monte_Carlo;
	dlg.K[1] = K[1];
	dlg.K[2] = K[2];
	dlg.K[3] = K[3];

	if ( dlg.DoModal ( ) == IDOK )
	{
	S=dlg.S;
	r=dlg.r;
	T=dlg.T;
	Sigma=dlg.Sigma;
	q=dlg.q;
	m=dlg.m;
	Black_Scholes=dlg.Black_Scholes;
    Monte_Carlo=dlg.Monte_Carlo;
	K[1]=dlg.K[1];
	K[2]=dlg.K[2];
	K[3]=dlg.K[3];
	Invalidate();
	}
}
