#include <afxwin.h>
#include "ddxddv.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
const double PI=3.14159265359;
const double MIN=0.00000001;
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
//Calculate probability density of standard normal distribution
double NormalDensity(double y)
{
  return 1/sqrt(2*PI)*exp(-y*y/2);
}
//The first method: Using Black-Scholes formula
double BS_European_Call(double S,double K,double r,double T,double Sigma,double q)
{
   double d1=(log(S/K)+(r-q+Sigma*Sigma/2)*T)/(Sigma*sqrt(T));
   double d2=d1-Sigma*sqrt(T);
   double European_CallValue=S*exp(-q*T)*Standard_Normal_Distribution(d1)
   -K*exp(-r*T)*Standard_Normal_Distribution(d2);
   return European_CallValue;
}
double BS_European_Put(double S,double K,double r,double T,double Sigma,double q)
{
   double d1=(log(S/K)+(r-q+Sigma*Sigma/2)*T)/(Sigma*sqrt(T));
   double d2=d1-Sigma*sqrt(T);
   double European_PutValue=K*exp(-r*T)*Standard_Normal_Distribution(-d2)
   -S*exp(-q*T)*Standard_Normal_Distribution(-d1);
   return European_PutValue;
}
//The third method: Using binomial tree 
double Binomial_Tree_European_Call(double S,double K,double r,double T,double Sigma,double q,int n)
{
  //Calculate parameters
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  //Calculate possible stock prices at T and store them in Array[10001]
  double Array[10001];
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
double Binomial_Tree_American_Call(double S,double K,double r,double T,double Sigma,double q,int n)
{
  double U, D, Pu, Pd, DeltaT;
  DeltaT=T/n;
  U=exp(Sigma*sqrt(DeltaT));
  D=exp(-Sigma*sqrt(DeltaT));
  Pu=(exp((r-q)*DeltaT)-D)/(U-D);
  Pd=1-Pu;
  double Array[10001];
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
double Binomial_Tree_European_Put(double S,double K,double r,double T,double Sigma,double q,int n)
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
double Binomial_Tree_American_Put(double S,double K,double r,double T,double Sigma,double q,int n)
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
//Differential of binomial tree 
double diff_Binomial_Tree_European_Call(double S,double K,double r,double T,double Sigma,double q,int n)
{
  return (Binomial_Tree_European_Call(S, K, r, T, Sigma+MIN, q, n)-Binomial_Tree_European_Call(S, K, r, T, Sigma, q, n))/MIN;
}
double diff_Binomial_Tree_American_Call(double S,double K,double r,double T,double Sigma,double q,int n)
{
  return (Binomial_Tree_American_Call(S, K, r, T, Sigma+MIN, q, n)-Binomial_Tree_American_Call(S, K, r, T, Sigma, q, n))/MIN;
}
double diff_Binomial_Tree_European_Put(double S,double K,double r,double T,double Sigma,double q,int n)
{
  return (Binomial_Tree_European_Put(S, K, r, T, Sigma+MIN, q, n)-Binomial_Tree_European_Put(S, K, r, T, Sigma, q, n))/MIN;
}
double diff_Binomial_Tree_American_Put(double S,double K,double r,double T,double Sigma,double q,int n)
{
  return (Binomial_Tree_American_Put(S, K, r, T, Sigma+MIN, q, n)-Binomial_Tree_American_Put(S, K, r, T, Sigma, q, n))/MIN;
}
class MyFrameWindow : public CFrameWnd
{
	public:

     double S, K, r, T, Sigma[2], q, P, tol;
     int n;
	 double fa,fb,fc,a,b,c,error_bound;
	 double a0, b0, error_bound0;
	 time_t Bisection_start, Bisection_end;
	 time_t Newton_start, Newton_end;
	 double Bisection_speed, Newton_speed;
	 int i; double sign;
	 int Call,Put,Bisection,Newton,Black_Scholes,Binomial_Tree;
	 int European;

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
		: CWinApp ( "Implied Volatility" )
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

     double S, K, r, T, q, P, tol;
     int n;
     int Call,Bisection,Newton,Black_Scholes,Binomial_Tree;
     int European;

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
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_P, P );
	DDX_Text ( pDX, IDC_tol, tol );
	DDX_Text ( pDX, IDC_n, n );
	DDV_MinMaxInt( pDX, n, 1, 10000 );
	DDX_Radio ( pDX, IDC_Put, Call );
	DDX_Check ( pDX, IDC_Bisection, Bisection );
	DDX_Check ( pDX, IDC_Newton, Newton );
	DDX_Check ( pDX, IDC_Black_Scholes, Black_Scholes );
	DDX_Check ( pDX, IDC_Binomial_Tree, Binomial_Tree );
	DDX_Radio ( pDX, IDC_Put2, European );
}

MyFrameWindow::MyFrameWindow ( )
{
    n=1;
	S=250; r=0.1; T=0.25; K=250; q=0.1; 
	P=11.15; tol=0.00000001;
	Call=1,Bisection=0,Newton=1,Black_Scholes=1,Binomial_Tree=0;
	European=1; 
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );
	
	int x_axis=20,y_axis=0,temp_x;
	
	if(Call)
	{
     _stprintf(info,"<Call>");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     if(Black_Scholes)
     {
      _stprintf(info,"Black-Scholes:");
      y_axis+=20;
      dc.TextOut(x_axis,y_axis,info);
      y_axis+=20;
      
      _stprintf(info,"Method");
      dc.TextOut(x_axis,y_axis,info);
      _stprintf(info,"European");
      dc.TextOut(x_axis+100,y_axis,info);
      _stprintf(info,"Time");
      dc.TextOut(x_axis+200,y_axis,info);    
      y_axis+=20;
      
      /*Bolzano's Theorem*/
      i=1;
      sign=1;
      while((sign>0)&&(i!=10))
      {
               sign=(BS_European_Call(S, K, r, T, ((double)i/10), q)-P)
               *(BS_European_Call(S, K, r, T, ((double)(i+1)/10), q)-P);
               i++;
      }
      if(i==10)
      {
                a=MIN, b=0.1;
      }
      else 
      {
            a=((double)(i-1)/10), b=((double)i/10);
      }
      a0=a;
      b0=b;
      error_bound0=(b - a)/2.0; 
   	   
      if(Bisection)
      {
       /*Bisection Method*/
       Bisection_start=time(NULL);
       a=a0;
       b=b0;
       error_bound=error_bound0;
       c = (a+b)/2.0;
       fa = BS_European_Call(S, K, r, T, a, q)-P; 
       fb = BS_European_Call(S, K, r, T, b, q)-P; 
       fc = BS_European_Call(S, K, r, T, c, q)-P; 	   
       while (error_bound > tol)
       {
		     if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = BS_European_Call(S, K, r, T, c, q)-P; 
		     error_bound = 0.5*(b - a);
        
       }
       Sigma[0]=c;
       Bisection_end=time(NULL);
       Bisection_speed=difftime(Bisection_end, Bisection_start);
       
       _stprintf(info,"Bisection");
       dc.TextOut(x_axis,y_axis,info);
       _stprintf(info,"%.4lf%%",Sigma[0]*100);
       dc.TextOut(x_axis+100,y_axis,info);
       _stprintf(info,"%d",(int)(Bisection_speed));
       dc.TextOut(x_axis+200,y_axis,info);    
       y_axis+=20;
      }
 	  
      if(Newton)
      {
       /*Newton's Method*/
	   Newton_start=time(NULL);
       a=a0;
       b=b0;
       error_bound=error_bound0;
	   double d1;
	   a=b;
	   while (error_bound > tol)
       {
        b=a;
		d1=(log(S/K)+(r-q+b*b/2)*T)/(b*sqrt(T));
		a=b-((BS_European_Call(S, K, r, T, b, q)-P)/(S*sqrt(T)*NormalDensity(d1)*exp(-q*T)));
        error_bound = 0.5*fabs(b - a);
       }
       Sigma[1]=a;
	   Newton_end=time(NULL);
       Newton_speed=difftime(Newton_end, Newton_start);
	   
       _stprintf(info,"Newton");
       dc.TextOut(x_axis,y_axis,info);
       _stprintf(info,"%.4lf%%",Sigma[1]*100);
       dc.TextOut(x_axis+100,y_axis,info);
       _stprintf(info,"%d",(int)(Newton_speed));
       dc.TextOut(x_axis+200,y_axis,info);    
       y_axis+=20;
      }
     }
     
     if(Binomial_Tree)
     {
      _stprintf(info,"Binomial Tree:");
      y_axis+=20;
      dc.TextOut(x_axis,y_axis,info);
      y_axis+=20;
      
      _stprintf(info,"Method");
      dc.TextOut(x_axis,y_axis,info);
      
      temp_x=x_axis+100;
      
      if(European)
      {
       _stprintf(info,"European");
       dc.TextOut(temp_x,y_axis,info);
       temp_x+=100;
       _stprintf(info,"Time");
       dc.TextOut(temp_x,y_axis,info);  
      }
	   
	  else
      {
       _stprintf(info,"American");
       dc.TextOut(temp_x,y_axis,info);
       temp_x+=100;
       _stprintf(info,"Time");
       dc.TextOut(temp_x,y_axis,info); 
      }	  
	  
      y_axis+=20;
   	   
      if(Bisection)
      {
       _stprintf(info,"Bisection");
       dc.TextOut(x_axis,y_axis,info);
       
       temp_x=x_axis+100;
       
       if(European)
       {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_European_Call(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_European_Call(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
        error_bound0=(b - a)/2.0;  
  	 
        /*Bisection Method*/
        Bisection_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
        c = (a+b)/2.0;
        fa = Binomial_Tree_European_Call(S, K, r, T, a, q, n)-P; 
        fb = Binomial_Tree_European_Call(S, K, r, T, b, q, n)-P; 
        fc = Binomial_Tree_European_Call(S, K, r, T, c, q, n)-P; 	   
        while (error_bound > tol)
        {
		     if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = Binomial_Tree_European_Call(S, K, r, T, c, q, n)-P; 
		     error_bound = 0.5*(b - a);
        }
        Sigma[0]=c;
        Bisection_end=time(NULL);
        Bisection_speed=difftime(Bisection_end, Bisection_start);
   
        _stprintf(info,"%.4lf%%",Sigma[0]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Bisection_speed));
        dc.TextOut(temp_x,y_axis,info);
       }
	   
	   else
       {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_American_Call(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_American_Call(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
        error_bound0=(b - a)/2.0; 
	   
        /*Bisection Method*/
        Bisection_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
        c = (a+b)/2.0;
        fa = Binomial_Tree_American_Call(S, K, r, T, a, q, n)-P; 
        fb = Binomial_Tree_American_Call(S, K, r, T, b, q, n)-P; 
        fc = Binomial_Tree_American_Call(S, K, r, T, c, q, n)-P; 	   
        while (error_bound > tol)
        {
	         if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = Binomial_Tree_American_Call(S, K, r, T, c, q, n)-P; 
		     error_bound = 0.5*(b - a);
        }  
        Sigma[0]=c;
        Bisection_end=time(NULL);
        Bisection_speed=difftime(Bisection_end, Bisection_start);
        
        _stprintf(info,"%.4lf%%",Sigma[0]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Bisection_speed));
        dc.TextOut(temp_x,y_axis,info);
       }  
       
       y_axis+=20;
      }
 	  
      if(Newton)
      {
       _stprintf(info,"Newton");
       dc.TextOut(x_axis,y_axis,info);
       
       temp_x=x_axis+100;
       
       if(European)
       {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_European_Call(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_European_Call(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
        error_bound0=(b - a)/2.0;  
		
	    /*Newton's Method*/
	    Newton_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
	    a=b;
	    while (error_bound > tol)
        {
          b=a;
		  a=b-((Binomial_Tree_European_Call(S, K, r, T, b, q, n)-P)/diff_Binomial_Tree_European_Call(S, K, r, T, b, q, n));
          error_bound = 0.5*fabs(b - a);
        }
        Sigma[1]=a;
	    Newton_end=time(NULL);
        Newton_speed=difftime(Newton_end, Newton_start);

        _stprintf(info,"%.4lf%%",Sigma[1]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Newton_speed));
        dc.TextOut(temp_x,y_axis,info); 
       }
	   
	   else
       {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_American_Call(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_American_Call(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
        error_bound0=(b - a)/2.0; 
   	   
        /*Newton's Method*/
        Newton_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
        a=b;
        while (error_bound > tol)
        {
          b=a;
		  a=b-((Binomial_Tree_American_Call(S, K, r, T, b, q, n)-P)/diff_Binomial_Tree_American_Call(S, K, r, T, b, q, n));
          error_bound = 0.5*fabs(b - a);
        }
        Sigma[1]=a;
        Newton_end=time(NULL);
        Newton_speed=difftime(Newton_end, Newton_start);
        
        _stprintf(info,"%.4lf%%",Sigma[1]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Newton_speed));
        dc.TextOut(temp_x,y_axis,info);	   
       }  
       
       y_axis+=20;
      } //Newton
     } //Binomial Tree
    } //Call
	
	else
	{
      _stprintf(info,"<Put>");
      y_axis+=20;
      dc.TextOut(x_axis,y_axis,info);
      y_axis+=20;
     
      if(Black_Scholes)
      {
       _stprintf(info,"Black-Scholes:");
       y_axis+=20;
       dc.TextOut(x_axis,y_axis,info);
       y_axis+=20;
      
       _stprintf(info,"Method");
       dc.TextOut(x_axis,y_axis,info);
       _stprintf(info,"European");
       dc.TextOut(x_axis+100,y_axis,info);
       _stprintf(info,"Time");
       dc.TextOut(x_axis+200,y_axis,info);    
       y_axis+=20;
      
       /*Bolzano's Theorem*/
       i=1;
       sign=1;
       while((sign>0)&&(i!=10))
       {
               sign=(BS_European_Put(S, K, r, T, ((double)i/10), q)-P)
               *(BS_European_Put(S, K, r, T, ((double)(i+1)/10), q)-P);
               i++;
       }
       if(i==10)
       {
                a=MIN, b=0.1;
       }
       else 
       {
            a=((double)(i-1)/10), b=((double)i/10);
       }
       a0=a;
       b0=b;
   	   error_bound0=(b - a)/2.0; 
   	   
       if(Bisection)
       {
        /*Bisection Method*/
        Bisection_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
	    c = (a+b)/2.0;
        fa = BS_European_Put(S, K, r, T, a, q)-P; 
	    fb = BS_European_Put(S, K, r, T, b, q)-P; 
	    fc = BS_European_Put(S, K, r, T, c, q)-P; 	   
	    while (error_bound > tol)
        {
		     if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = BS_European_Put(S, K, r, T, c, q)-P; 
		     error_bound = 0.5*(b - a);
        }
        Sigma[0]=c;
        Bisection_end=time(NULL);
        Bisection_speed=difftime(Bisection_end, Bisection_start);
       
	    _stprintf(info,"Bisection");
        dc.TextOut(x_axis,y_axis,info);
        _stprintf(info,"%.4lf%%",Sigma[0]*100);
        dc.TextOut(x_axis+100,y_axis,info);
        _stprintf(info,"%d",(int)(Bisection_speed));
        dc.TextOut(x_axis+200,y_axis,info);    
        y_axis+=20;
 	   }
 	  
 	   if(Newton)
 	   {
        /*Newton's Method*/
		Newton_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
		double d1;
		a=b;
		while (error_bound > tol)
        {
        b=a;
		d1=(log(S/K)+(r-q+b*b/2)*T)/(b*sqrt(T));
		a=b-((BS_European_Put(S, K, r, T, b, q)-P)/(S*sqrt(T)*NormalDensity(d1)*exp(-q*T)));
        error_bound = 0.5*fabs(b - a);
        }
        Sigma[1]=a;
		Newton_end=time(NULL);
        Newton_speed=difftime(Newton_end, Newton_start);
	   
	    _stprintf(info,"Newton");
        dc.TextOut(x_axis,y_axis,info);
        _stprintf(info,"%.4lf%%",Sigma[1]*100);
        dc.TextOut(x_axis+100,y_axis,info);
        _stprintf(info,"%d",(int)(Newton_speed));
        dc.TextOut(x_axis+200,y_axis,info);    
        y_axis+=20;
       }
      }
     
      if(Binomial_Tree)
      {
       _stprintf(info,"Binomial Tree:");
       y_axis+=20;
       dc.TextOut(x_axis,y_axis,info);
       y_axis+=20;
      
       _stprintf(info,"Method");
       dc.TextOut(x_axis,y_axis,info);
      
       temp_x=x_axis+100;
      
       if(European)
       {
        _stprintf(info,"European");
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"Time");
        dc.TextOut(temp_x,y_axis,info);   
	   }
	  
	   else
	   {
        _stprintf(info,"American");
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"Time");
        dc.TextOut(temp_x,y_axis,info); 
	   }	  
	   
	   y_axis+=20;
   	   
       if(Bisection)
       {
  	    _stprintf(info,"Bisection");
        dc.TextOut(x_axis,y_axis,info);
       
   	    temp_x=x_axis+100;
       
        if(European)
        {
         /*Bolzano's Theorem*/
         i=1;
         sign=1;
         while((sign>0)&&(i!=10))
         {
               sign=(Binomial_Tree_European_Put(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_European_Put(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
         }
         if(i==10)
         {
                a=MIN, b=0.1;
         }
         else 
         {
            a=((double)(i-1)/10), b=((double)i/10);
         }
         a0=a;
         b0=b;
   	     error_bound0=(b - a)/2.0;  
  	 
	     /*Bisection Method*/
         Bisection_start=time(NULL);
         a=a0;
         b=b0;
         error_bound=error_bound0;
	     c = (a+b)/2.0;
         fa = Binomial_Tree_European_Put(S, K, r, T, a, q, n)-P; 
	     fb = Binomial_Tree_European_Put(S, K, r, T, b, q, n)-P; 
	     fc = Binomial_Tree_European_Put(S, K, r, T, c, q, n)-P; 	   
	     while (error_bound > tol)
         {
		     if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = Binomial_Tree_European_Put(S, K, r, T, c, q, n)-P; 
		     error_bound = 0.5*(b - a);
         }
         Sigma[0]=c;
         Bisection_end=time(NULL);
         Bisection_speed=difftime(Bisection_end, Bisection_start);
   
         _stprintf(info,"%.4lf%%",Sigma[0]*100);
         dc.TextOut(temp_x,y_axis,info);
         temp_x+=100;
         _stprintf(info,"%d",(int)(Bisection_speed));
         dc.TextOut(temp_x,y_axis,info);
	    }
	   
		else
	    {
         /*Bolzano's Theorem*/
         i=1;
         sign=1;
         while((sign>0)&&(i!=10))
         {
               sign=(Binomial_Tree_American_Put(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_American_Put(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
         }
         if(i==10)
         {
                a=MIN, b=0.1;
         }
         else 
         {
            a=((double)(i-1)/10), b=((double)i/10);
         }
         a0=a;
         b0=b;
 	     error_bound0=(b - a)/2.0; 
	   
         /*Bisection Method*/
         Bisection_start=time(NULL);
         a=a0;
         b=b0;
         error_bound=error_bound0;
   	     c = (a+b)/2.0;
         fa = Binomial_Tree_American_Put(S, K, r, T, a, q, n)-P; 
	     fb = Binomial_Tree_American_Put(S, K, r, T, b, q, n)-P; 
	     fc = Binomial_Tree_American_Put(S, K, r, T, c, q, n)-P; 	   
	     while (error_bound > tol)
         {
		     if (fa * fc < 0.0) {
			 b = c;
			 fb = fc;
			 c = (a+c)/2.0;
             }
             else 
             {
			 a = c;
			 fa = fc;
			 c = (b+c)/2.0;
             }
		     fc = Binomial_Tree_American_Put(S, K, r, T, c, q, n)-P; 
		     error_bound = 0.5*(b - a);
         }
         Sigma[0]=c;
         Bisection_end=time(NULL);
         Bisection_speed=difftime(Bisection_end, Bisection_start);
        
         _stprintf(info,"%.4lf%%",Sigma[0]*100);
         dc.TextOut(temp_x,y_axis,info);
         temp_x+=100;
         _stprintf(info,"%d",(int)(Bisection_speed));
         dc.TextOut(temp_x,y_axis,info);
	    }  
        
		y_axis+=20;
 	  }
 	  
 	  if(Newton)
 	  {
	   _stprintf(info,"Newton");
       dc.TextOut(x_axis,y_axis,info);
       
       temp_x=x_axis+100;
       
	   if(European)
	   {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_European_Put(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_European_Put(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
        error_bound0=(b - a)/2.0;  
		
		/*Newton's Method*/
		Newton_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
		a=b;
		while (error_bound > tol)
        {
        b=a;
		a=b-((Binomial_Tree_European_Put(S, K, r, T, b, q, n)-P)/diff_Binomial_Tree_European_Put(S, K, r, T, b, q, n));
        error_bound = 0.5*fabs(b - a);
        }
        Sigma[1]=a;
		Newton_end=time(NULL);
        Newton_speed=difftime(Newton_end, Newton_start);

        _stprintf(info,"%.4lf%%",Sigma[1]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Newton_speed));
        dc.TextOut(temp_x,y_axis,info);
	   }
	   
	   else
	   {
        /*Bolzano's Theorem*/
        i=1;
        sign=1;
        while((sign>0)&&(i!=10))
        {
               sign=(Binomial_Tree_American_Put(S, K, r, T, ((double)i/10), q, n)-P)
               *(Binomial_Tree_American_Put(S, K, r, T, ((double)(i+1)/10), q, n)-P);
               i++;
        }
        if(i==10)
        {
                a=MIN, b=0.1;
        }
        else 
        {
            a=((double)(i-1)/10), b=((double)i/10);
        }
        a0=a;
        b0=b;
   	    error_bound0=(b - a)/2.0; 
   	   
        /*Newton's Method*/
		Newton_start=time(NULL);
        a=a0;
        b=b0;
        error_bound=error_bound0;
		a=b;
		while (error_bound > tol)
        {
        b=a;
		a=b-((Binomial_Tree_American_Put(S, K, r, T, b, q, n)-P)/diff_Binomial_Tree_American_Put(S, K, r, T, b, q, n));
        error_bound = 0.5*fabs(b - a);
        }
        Sigma[1]=a;
		Newton_end=time(NULL);
        Newton_speed=difftime(Newton_end, Newton_start);
        
        _stprintf(info,"%.4lf%%",Sigma[1]*100);
        dc.TextOut(temp_x,y_axis,info);
        temp_x+=100;
        _stprintf(info,"%d",(int)(Newton_speed));
        dc.TextOut(temp_x,y_axis,info);	   
	   }  
	   
       y_axis+=20;
      } //Newton
     } //Binomial Tree
	} //Put
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S = S;
	dlg.K = K;
	dlg.r = r;
	dlg.T = T;
	dlg.q = q;
	dlg.P = P;
	dlg.tol = tol;
	dlg.n = n;
	dlg.Call = Call;
	dlg.Bisection = Bisection;
	dlg.Newton = Newton;
	dlg.Black_Scholes = Black_Scholes;
	dlg.Binomial_Tree = Binomial_Tree;
	dlg.European = European;

	if ( dlg.DoModal ( ) == IDOK )
	{
	S=dlg.S;
	K=dlg.K;
	r=dlg.r;
	T=dlg.T;
	q=dlg.q;
	P=dlg.P;
	tol=dlg.tol;
	n=dlg.n;
	Call = dlg.Call;
	Bisection = dlg.Bisection;
	Newton = dlg.Newton;
	Black_Scholes = dlg.Black_Scholes;
	Binomial_Tree = dlg.Binomial_Tree;
	European = dlg.European;
	Invalidate();
	}
}
