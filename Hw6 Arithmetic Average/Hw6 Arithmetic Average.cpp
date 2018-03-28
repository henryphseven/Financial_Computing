#include <afxwin.h>
#include "ddxddv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
int i,j,k,l;
double Sum;
const int MAX=1000;
const double M_PI=3.14159265359;
//Generate a standard normal random variable
double Normal()
{
     double Y,X[2];
	 for(int j=0;j<2;j++)
	 {
	   X[j]=double(rand())/RAND_MAX;
	 }
     Y=sqrt((-2)*log(X[0]))*cos(2*M_PI*X[1]);
	 return Y;
}
//Standard normal distribution's pdf
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
                  +a5*pow(k,5))/sqrt(2*M_PI);
   if(flag) return 1-value;
    else return value;
}
//Geometric average call
double GeometricCall(double S, double T, double X, double r, double Sigma, double q)
{
  double OptionValue;
  double D_Start=(r-q-(Sigma*Sigma)/6.0)/2.0*T;
  double D=(log(S/X)+
	      (r-q+(Sigma*Sigma)/6.0)/2.0*T)/
		  (Sigma*sqrt(T/3.0));
  OptionValue=exp(D_Start)*S*Standard_Normal_Distribution(D)-
	 X*Standard_Normal_Distribution(D-Sigma*sqrt(T/3.0));
  double Discounted=exp(r*T);
  OptionValue/=Discounted; 
  return OptionValue;
}
//Geometric average put
double GeometricPut(double S, double T, double K, double r, double Sigma, double q)
{
   double mu_GA=(r-q-Sigma*Sigma/6)/2;
   double Sigma_GA=Sigma/sqrt(3.0);
   double d1=(log(S/K)+(mu_GA+Sigma_GA*Sigma_GA/2)*T)/(Sigma_GA*sqrt(T));
   double d2=d1-Sigma_GA*sqrt(T);
   return exp(-r*T)*(K*Standard_Normal_Distribution(-d2)
   -S*exp(mu_GA*T)*Standard_Normal_Distribution(-d1));
}
//Sequential search
int SequentialSearch(double a,double b[],int P) //P: Top index
{
    int i;
	for(i=1;i<P;i++)
	{
		if(a>b[i]) break;
	}
    return i-1;
}
//Binary search
int BinarySearch(double a,double b[],int P) //P: Top index
{
    int bu=0,bd=P;
    while(bd-bu>1)
    {
     if(a>b[(bu+bd)/2]) bd=(bu+bd)/2;
     else bu=(bu+bd)/2;
    }
    return bu;
}
//Interpolation
double Interpolation(double a0,double a2,double b0,double b1,double b2)
{
    if(b0==b2) return a0;
    else return (b1-b0)*(a2-a0)/(b2-b0)+a0;
}
//Ordering: Decreasing 
void bubble_sort(double a[],int p) //p: Number of numbers
{
   int i,j;
   double temp;
   int flag=0;

   for(i=1;flag==0;i++)
   {
      flag=1;			// 將flag設為1
      for(j=0;j<(p-i);j++)
         if(a[j]<a[j+1])
         {
            temp=a[j];	    	// 對換陣列內的值
            a[j]=a[j+1];
            a[j+1]=temp;
            flag=0;				// 對調後將flag設為0
         }
   }
   return;
}
void space(double a,char *combination)
{
	int i,decimal;

	i=0;
	while(a/pow(10.0,i)>1)
	{
		i++;
	}

	decimal=i;

	switch(decimal)
	{
	 case 0: strcat(combination,"  "); break;
     case 1: strcat(combination,"  "); break;
	 case 2: strcat(combination,"   "); break;
	 case 3: strcat(combination,"    "); break;
	 case 4: strcat(combination,"     "); break;
	 case 5: strcat(combination,"      "); break;
	 case 6: strcat(combination,"       "); break;
	 case 7: strcat(combination,"        "); break;
	 case 8: strcat(combination,"         "); break;
	}
}

class MyFrameWindow : public CFrameWnd
{
	public:

 	  //Input
 	  double S,K,r,T,Sigma,q;
 	  double Save_0,T_0;
 	  int n,p;
	  int Call,Put,Show_All;
 	  int Pascal,Insert_K,Log_Space,Binary;
 	  int Monte_Carlo,m,With,Without,Binomial_Tree;

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
		: CWinApp ( "Arithmetic Average Option Price" )
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

 	  //Input
 	  double S,K,r,T,Sigma,q;
 	  double Save_0,T_0;
 	  int n,p;
	  int Call,Put,Show_All;
 	  int Pascal,Insert_K,Log_Space,Binary;
 	  int Monte_Carlo,m,With,Without,Binomial_Tree;

	MyDialog ( CWnd* parent ) : CDialog ( IDD_INPUT, parent ) {}
	void DoDataExchange( CDataExchange* pDX );
};

void MyDialog::DoDataExchange( CDataExchange* pDX )
{
	CDialog::DoDataExchange ( pDX );

	DDX_Text ( pDX, IDC_S, S );
	DDX_Text ( pDX, IDC_K, K );
	DDX_Text ( pDX, IDC_r, r );
	DDV_MinMaxDouble( pDX, r, 0, 1 );
	DDX_Text ( pDX, IDC_T, T );
	DDX_Text ( pDX, IDC_Sigma, Sigma );
	DDV_MinMaxDouble( pDX, Sigma, 0, 1 );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_n, n );
	DDV_MinMaxInt( pDX, n, 1, MAX );
	DDX_Text ( pDX, IDC_p, p );
	DDV_MinMaxInt( pDX, p, 1, MAX );
	DDX_Text ( pDX, IDC_Save_0, Save_0 );
	DDX_Text ( pDX, IDC_T_0, T_0 );
	DDX_Check( pDX, IDC_Call, Call );
	DDX_Check( pDX, IDC_Put, Put );
	DDX_Check( pDX, IDC_Show_All, Show_All );
	DDX_Check( pDX, IDC_Pascal, Pascal );
	DDX_Check( pDX, IDC_Insert_K, Insert_K );
	DDX_Check( pDX, IDC_Log_Space, Log_Space );
	DDX_Check( pDX, IDC_Binary, Binary );	
	DDX_Check( pDX, IDC_Monte_Carlo, Monte_Carlo );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, 10000 );
	DDX_Check( pDX, IDC_With, With );
	DDX_Check( pDX, IDC_Without, Without );
	DDX_Check( pDX, IDC_Binomial_Tree, Binomial_Tree );
}

MyFrameWindow::MyFrameWindow ( )
{
 	  //Input
 	  S=50,K=50,r=0.1,T=1,Sigma=0.4,q=0;
 	  Save_0=S,T_0=0;
 	  n=20,p=4;
	  Call=0,Put=0,Show_All=0;
 	  Pascal=0,Insert_K=0,Log_Space=0,Binary=0;
 	  Monte_Carlo=0,m=10,With=0,Without=0,Binomial_Tree=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );

	  int L=100,P=p-1;

      double Save,St,delta_t=T/n;
      double E_Call,Square_E_Call;
      double E_Put,Square_E_Put;
      int n_0=(int)(T_0/delta_t);
      double Sum_G,Save_G;
	  double E_Call_G,Square_E_Call_G;
      double E_Put_G,Square_E_Put_G;
      double STDERR,Upper_Limit,Lower_Limit,STDERR_G,Upper_Limit_G,Lower_Limit_G;
	  
      double u,d,pu,pd;
      double Au,Ad,Cu,Cd;
      time_t Start,Middle,End_S,End;
      double Speed;
      int position_u,position_d;
      
      int flag;

	  char combination[128];

	  int x_axis=20,y_axis=0;
	  
      if(Call) { //Call: Start
      _stprintf(info,"<Call>");
	  y_axis+=20;
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	
	 if(Monte_Carlo) {		
//Monte Carlo

     srand( (unsigned)time( NULL ) );
	 y_axis+=20;
     _stprintf(info,"Monte Carlo:");
	 dc.TextOut(x_axis,y_axis,info); 
	 y_axis+=20;
     
     do{
	 E_Call=0,Square_E_Call=0;
	 E_Call_G=0,Square_E_Call_G=0;
     for(i=0;i<m;i++)
     {
	  Sum=0,Sum_G=0;
	  for(j=0;j<L;j++)
	  {
	   St=S;
	   Save=Save_0*(n_0+1);
	   Save_G=log(St);
	   for(k=1;k<=n;k++)
	   {
	    St=St*exp((r-q-Sigma*Sigma/2)*delta_t+Sigma*sqrt(delta_t)*Normal());
	    Save+=St;
	    Save_G+=log(St);
	   }
	   Sum+=(exp(-r*T)*max(Save/(n_0+1+n)-K,0))/L;
	   Sum_G+=(exp(-r*T)*
	   (max(Save/(n_0+1+n)-K,0)-max(exp(Save_G/(n+1))-K,0))
	   +GeometricCall(S,T,K,r,Sigma,q))/L;
	  }
	  E_Call+=Sum;
	  Square_E_Call+=pow(Sum,2);
	  E_Call_G+=Sum_G;
	  Square_E_Call_G+=pow(Sum_G,2);
	 }
     E_Call/=m;
     }while((_isnan(E_Call)!=0)||(_finite(E_Call)==0)); 

     STDERR=sqrt((double)((Square_E_Call-m*pow(E_Call,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit=E_Call-2*STDERR;
     Upper_Limit=E_Call+2*STDERR;	
	 
  	 E_Call_G/=m;
     STDERR_G=sqrt((double)((Square_E_Call_G-m*pow(E_Call_G,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit_G=E_Call_G-2*STDERR;
     Upper_Limit_G=E_Call_G+2*STDERR;

	 strcpy(combination,"Control Variate");
	 strcat(combination,"  ");
	 space(max(E_Call,E_Call_G),combination);
	 strcat(combination,"European");
	 space(max(E_Call,E_Call_G),combination);
	 strcat(combination,"    ");
	 space(max(STDERR,STDERR_G),combination);
	 strcat(combination,"SE");
	 space(max(STDERR,STDERR_G),combination);
	 strcat(combination,"        ");
	 space(max(Lower_Limit,Lower_Limit_G),combination);
	 space(max(Upper_Limit,Upper_Limit_G),combination);
	 strcat(combination,"95%  CI");

	 _stprintf(info,"%s",combination);
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
	 if(With)
	 {
      _stprintf(info,"           o                 %.4lf       %.4lf      [%.4lf,%.4lf]"
	  ,E_Call_G,STDERR_G,Lower_Limit_G,Upper_Limit_G);
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	 }
	 
	 if(Without)
	 {
	  _stprintf(info,"           x                 %.4lf       %.4lf      [%.4lf,%.4lf]"
	  ,E_Call,STDERR,Lower_Limit,Upper_Limit);
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	 }

	 } //Monte Carlo

	 if(Binomial_Tree) {
//Binomial tree

     y_axis+=20;
	 _stprintf(info,"Binomial Tree:");
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
  	 u=exp(Sigma*sqrt(delta_t));
  	 d=exp(-Sigma*sqrt(delta_t));
  	 pu=(exp((r-q)*delta_t)-d)/(u-d);
  	 pd=1-pu;
  	 
  	 struct node
  	 {
	  		double S;
	  		double *A;
	  		double *E_Call;
	  		double *A_Call;
	  		int p;
	 } *tree[MAX+1];
	 
 	 for(i=0;i<=n;i++)
	 tree[i]=new struct node[MAX+1];
	 
	 for(i=0;i<=n;i++)
  	 {
	  for(j=0;j<=i;j++)
	  {
	   tree[i][j].A=new double[MAX+2];
	   tree[i][j].E_Call=new double[MAX+2];
	   tree[i][j].A_Call=new double[MAX+2];
	  }
     }
	 
	 //Calculate St at each node
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       //In case that exactly the same values differ
       if(i-j>j) tree[i][j].S=S*pow(u,(i-j)-j);
       else if((i-j)==j) tree[i][j].S=S;
       else tree[i][j].S=S*pow(d,j-(i-j));
	  }
   	 }

//Original
   	 
   	 Start=time(NULL);
   	 //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
       }
       tree[i][j].p=P+1;
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Call[k]=max(tree[n][j].A[k]-K,0);
      tree[n][j].A_Call[k]=tree[n][j].E_Call[k];
	 }
	 
	 Middle=time(NULL);
	 
	 //Backward induction: Using sequential search 
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {		
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=SequentialSearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=SequentialSearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     }
     End_S=time(NULL);
     Speed=difftime(End_S,Start);

	 strcpy(combination,"Pascal  Insert K  Log-Space  Binary  ");
	 space(tree[0][0].E_Call[0],combination);
	 strcat(combination,"European");
	 space(tree[0][0].E_Call[0],combination);
	 strcat(combination," ");
	 space(tree[0][0].A_Call[0],combination);
	 strcat(combination,"American");
	 space(tree[0][0].A_Call[0],combination);
	 space(Speed,combination);
	 strcat(combination,"Time");

	 _stprintf(info,"%s",combination);
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"     x            x              x              x           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
   	 if(Show_All) {
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start)-difftime(End_S,Middle);
     
	 _stprintf(info,"     x            x              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
//Pascal

     Start=time(NULL);
	 tree[0][0].A[0]=Save_0;
     tree[0][0].p=1;
     
   	 //Inherit Save
  	 for(i=1;i<=n;i=i+1)
  	 {	 
	  tree[i][0].A[0]=(tree[i-1][0].A[0]*(n_0+i)+tree[i][0].S)/(n_0+i+1);
	  tree[i][0].p=1;
	  
	  for(j=1;j<i;j++)
	  {
       if(tree[i-1][j-1].p+tree[i-1][j].p>P+1)
   	   {
        tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
        /((n_0+1)+i);
        tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
        /((n_0+1)+i);
        for(k=1;k<=P-1;k++)
        {
         tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
        }
        tree[i][j].p=P+1;
       }
       
       else
       {
        tree[i][j].p=0;
	   
	    //Inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j-1].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
        }

	    //Inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
	    }
	    
        bubble_sort(tree[i][j].A,tree[i][j].p);
       }
      } 
	  
	  tree[i][i].A[0]=(tree[i-1][i-1].A[0]*(n_0+i)+tree[i][i].S)/(n_0+i+1);
	  tree[i][i].p=1;
   	 }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Call[k]=max(tree[n][j].A[k]-K,0);
      tree[n][j].A_Call[k]=tree[n][j].E_Call[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     o            x              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;

//Inserting the strike price   

     Start=time(NULL);
	 //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
       }
       tree[i][j].p=P+1;
       if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
       {
        flag=0;
        for(l=0;l<tree[i][j].p;l++)
        {
         if(K==tree[i][j].A[l]) flag=1;
        }
        if(flag==0)
        {
         tree[i][j].A[tree[i][j].p]=K;
         tree[i][j].p++;
         bubble_sort(tree[i][j].A,tree[i][j].p);
        }
	   }
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Call[k]=max(tree[n][j].A[k]-K,0);
      tree[n][j].A_Call[k]=tree[n][j].E_Call[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     x            o              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
//Using log-space partition

     Start=time(NULL);
	 //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
       }
       tree[i][j].p=P+1;
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Call[k]=max(tree[n][j].A[k]-K,0);
      tree[n][j].A_Call[k]=tree[n][j].E_Call[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     x            x              o              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed)); 
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 } //Show_All: End
	 
//Selected combination

     if((Pascal)||(Insert_K)||(Log_Space)||(Binary)) {
	 
	 Start=time(NULL);

	 if(Pascal) {
     tree[0][0].A[0]=Save_0;
     tree[0][0].p=1;
     
   	 //Inherit Save
  	 for(i=1;i<=n;i=i+1)
  	 {	 
	  tree[i][0].A[0]=(tree[i-1][0].A[0]*(n_0+i)+tree[i][0].S)/(n_0+i+1);
	  tree[i][0].p=1;
	  
	  for(j=1;j<i;j++)
	  {
       if(tree[i-1][j-1].p+tree[i-1][j].p>P+1)
   	   {
        tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
        /((n_0+1)+i);
        tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
        /((n_0+1)+i);
        for(k=1;k<=P-1;k++)
        {
         if(Log_Space) tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
		 else tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
		}
        tree[i][j].p=P+1;
        
        if(Insert_K) {
        if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
        {
         flag=0;
         for(l=0;l<tree[i][j].p;l++)
         {
          if(K==tree[i][j].A[l]) flag=1;
         }
         if(flag==0)
         {
          tree[i][j].A[tree[i][j].p]=K;
          tree[i][j].p++;
          bubble_sort(tree[i][j].A,tree[i][j].p);
         }
	    }
	    }
       }
       
       else
       {
        tree[i][j].p=0;
	   
	    //Inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j-1].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
        }

	    //Inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
	    }
	    
        if(Insert_K) {
        if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
        {
         flag=0;
         for(l=0;l<tree[i][j].p;l++)
         {
          if(K==tree[i][j].A[l]) flag=1;
         }
         if(flag==0)
         {
          tree[i][j].A[tree[i][j].p]=K;
          tree[i][j].p++;
         }
	    }
	    }
	    
        bubble_sort(tree[i][j].A,tree[i][j].p);
       }
      } 
	  
	  tree[i][i].A[0]=(tree[i-1][i-1].A[0]*(n_0+i)+tree[i][i].S)/(n_0+i+1);
	  tree[i][i].p=1;
   	 }
	 }	  
     
     else {
     //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        if(Log_Space) tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
		else tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
	   }
       tree[i][j].p=P+1;
       
       if(Insert_K) {
       if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
       {
        flag=0;
        for(l=0;l<tree[i][j].p;l++)
        {
         if(K==tree[i][j].A[l]) flag=1;
        }
        if(flag==0)
        {
         tree[i][j].A[tree[i][j].p]=K;
         tree[i][j].p++;
         bubble_sort(tree[i][j].A,tree[i][j].p);
        }
	   } //if
	   } //Insert_K
      } //j
     } //i
	 } //else
	 	      
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Call[k]=max(tree[n][j].A[k]-K,0);
      tree[n][j].A_Call[k]=tree[n][j].E_Call[k];
	 }
	 
	 //Backward induction
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {		
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        if(Binary) position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        else position_u=SequentialSearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        if(Binary) position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1);  
        else position_d=SequentialSearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1);             
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Call[position_u],tree[i+1][j].E_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Call[position_d],tree[i+1][j+1].E_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Call[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Call[position_u],tree[i+1][j].A_Call[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Call[position_d],tree[i+1][j+1].A_Call[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Call[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),tree[i][j].A[k]-K);
       }
      }
     } 
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 strcpy(combination,"");
	 
	 if(Pascal) strcat(combination,"     o      ");
	 else strcat(combination,"     x      ");

	 if(Insert_K) strcat(combination,"      o        ");
	 else strcat(combination,"      x        ");

	 if(Log_Space) strcat(combination,"      o         ");
	 else strcat(combination,"      x         ");

	 if(Binary) strcat(combination,"     o          ");
	 else strcat(combination,"     x          ");

	 _stprintf(info,"%s %.4lf         %.4lf        %2d",combination,tree[0][0].E_Call[0],tree[0][0].A_Call[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 } //Selected combination

//Release occupied memory

	 for(i=0;i<=n;i++)
  	 {
	  for(j=0;j<=i;j++)
	  {
 	   delete[] tree[i][j].A;
	   tree[i][j].A=NULL;
	   delete[] tree[i][j].E_Call;
	   tree[i][j].E_Call=NULL;
	   delete[] tree[i][j].A_Call;
	   tree[i][j].A_Call=NULL;
	  }
     }
     
	 for(i=0;i<=n;i++)
	 {
  	  delete[] tree[i];
  	  tree[i]=NULL;
	 }
	 
	 } //Binomial Tree
	 } //Call: End
	 	 	 
	 if(Put) { //Put: Start
      _stprintf(info,"<Put>");
	  y_axis+=20;
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	 
	 if(Monte_Carlo) {
//Monte Carlo

     srand( (unsigned)time( NULL ) );
	 y_axis+=20;
     _stprintf(info,"Monte Carlo:");
	 dc.TextOut(x_axis,y_axis,info); 
	 y_axis+=20;
     
     do{
	 E_Put=0,Square_E_Put=0;
	 E_Put_G=0,Square_E_Put_G=0;
     for(i=0;i<m;i++)
     {
	  Sum=0,Sum_G=0;
	  for(j=0;j<L;j++)
	  {
	   St=S;
	   Save=Save_0*(n_0+1);
	   Save_G=log(St);
	   for(k=1;k<=n;k++)
	   {
	    St=St*exp((r-q-Sigma*Sigma/2)*delta_t+Sigma*sqrt(delta_t)*Normal());
	    Save+=St;
	    Save_G+=log(St);
	   }
	   Sum+=(exp(-r*T)*max(K-Save/(n_0+1+n),0))/L;
	   Sum_G+=(exp(-r*T)*
	   (max(K-Save/(n_0+1+n),0)-max(K-exp(Save_G/(n+1)),0))
	   +GeometricPut(S,T,K,r,Sigma,q))/L;
	  }
	  E_Put+=Sum;
	  Square_E_Put+=pow(Sum,2);
	  E_Put_G+=Sum_G;
	  Square_E_Put_G+=pow(Sum_G,2);
	 }
     E_Put/=m;
     }while((_isnan(E_Put)!=0)||(_finite(E_Put)==0)); 

     STDERR=sqrt((double)((Square_E_Put-m*pow(E_Put,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit=E_Put-2*STDERR;
     Upper_Limit=E_Put+2*STDERR;	
	 
  	 E_Put_G/=m;
     STDERR_G=sqrt((double)((Square_E_Put_G-m*pow(E_Put_G,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit_G=E_Put_G-2*STDERR;
     Upper_Limit_G=E_Put_G+2*STDERR;

	 strcpy(combination,"Control Variate");
	 strcat(combination,"  ");
	 space(max(E_Put,E_Put_G),combination);
	 strcat(combination,"European");
	 space(max(E_Put,E_Put_G),combination);
	 strcat(combination,"    ");
	 space(max(STDERR,STDERR_G),combination);
	 strcat(combination,"SE");
	 space(max(STDERR,STDERR_G),combination);
	 strcat(combination,"        ");
	 space(max(Lower_Limit,Lower_Limit_G),combination);
	 space(max(Upper_Limit,Upper_Limit_G),combination);
	 strcat(combination,"95%  CI");

	 _stprintf(info,"%s",combination);
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
	 if(With)
	 {
      _stprintf(info,"           o                 %.4lf       %.4lf      [%.4lf,%.4lf]"
	  ,E_Put_G,STDERR_G,Lower_Limit_G,Upper_Limit_G);
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	 }
	 
	 if(Without)
	 {
	  _stprintf(info,"           x                 %.4lf       %.4lf      [%.4lf,%.4lf]"
	  ,E_Put,STDERR,Lower_Limit,Upper_Limit);
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	 }

	 } //Monte Carlo

	 if(Binomial_Tree) {
//Binomial tree

     y_axis+=20;
	 _stprintf(info,"Binomial Tree:");
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
  	 u=exp(Sigma*sqrt(delta_t));
  	 d=exp(-Sigma*sqrt(delta_t));
  	 pu=(exp((r-q)*delta_t)-d)/(u-d);
  	 pd=1-pu;
  	 
  	 struct node
  	 {
	  		double S;
	  		double *A;
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
	   tree[i][j].A=new double[MAX+2];
	   tree[i][j].E_Put=new double[MAX+2];
	   tree[i][j].A_Put=new double[MAX+2];
	  }
     }
	 
	 //Calculate St at each node
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       //In case that exactly the same values differ
       if(i-j>j) tree[i][j].S=S*pow(u,(i-j)-j);
       else if((i-j)==j) tree[i][j].S=S;
       else tree[i][j].S=S*pow(d,j-(i-j));
	  }
   	 }

//Original
   	 
   	 Start=time(NULL);
     //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
       }
       tree[i][j].p=P+1;
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Put[k]=max(K-tree[n][j].A[k],0);
      tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	 }
	 
	 Middle=time(NULL);
	 
	 //Backward induction: Using sequential search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {		
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=SequentialSearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=SequentialSearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     }
     End_S=time(NULL);
     Speed=difftime(End_S,Start);

	 strcpy(combination,"Pascal  Insert K  Log-Space  Binary  ");
	 space(tree[0][0].E_Put[0],combination);
	 strcat(combination,"European");
	 space(tree[0][0].E_Put[0],combination);
	 strcat(combination," ");
	 space(tree[0][0].A_Put[0],combination);
	 strcat(combination,"American");
	 space(tree[0][0].A_Put[0],combination);
	 space(Speed,combination);
	 strcat(combination,"Time");

	 _stprintf(info,"%s",combination);
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"     x            x              x              x           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
   	 if(Show_All) {
	 //Backward induction: Using binary search 
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start)-difftime(End_S,Middle);
     
	 _stprintf(info,"     x            x              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
     
//Pascal

	 Start=time(NULL);	
     tree[0][0].A[0]=Save_0;
     tree[0][0].p=1;
     
   	 //Inherit Save
  	 for(i=1;i<=n;i=i+1)
  	 {	 
	  tree[i][0].A[0]=(tree[i-1][0].A[0]*(n_0+i)+tree[i][0].S)/(n_0+i+1);
	  tree[i][0].p=1;
	  
	  for(j=1;j<i;j++)
	  {
       if(tree[i-1][j-1].p+tree[i-1][j].p>P+1)
   	   {
        tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
        /((n_0+1)+i);
        tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
        /((n_0+1)+i);
        for(k=1;k<=P-1;k++)
        {
         tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
        }
        tree[i][j].p=P+1;
       }
       
       else
       {
        tree[i][j].p=0;
	   
	    //Inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j-1].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
        }

	    //Inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
	    }
	    
        bubble_sort(tree[i][j].A,tree[i][j].p);
       }
      } 
	  
	  tree[i][i].A[0]=(tree[i-1][i-1].A[0]*(n_0+i)+tree[i][i].S)/(n_0+i+1);
	  tree[i][i].p=1;
   	 }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Put[k]=max(K-tree[n][j].A[k],0);
      tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     o            x              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;

//Inserting the strike price   

     Start=time(NULL);
	 //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
       }
       tree[i][j].p=P+1;
       if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
       {
        flag=0;
        for(l=0;l<tree[i][j].p;l++)
        {
         if(K==tree[i][j].A[l]) flag=1;
        }
        if(flag==0)
        {
         tree[i][j].A[tree[i][j].p]=K;
         tree[i][j].p++;
         bubble_sort(tree[i][j].A,tree[i][j].p);
        }
	   }
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Put[k]=max(K-tree[n][j].A[k],0);
      tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     x            o              x              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 
//Using log-space partition

     Start=time(NULL);
	 //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
       }
       tree[i][j].p=P+1;
      }
     }
          
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Put[k]=max(K-tree[n][j].A[k],0);
      tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	 }
	 
	 //Backward induction: Using binary search
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        
        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
               
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     }
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 _stprintf(info,"     x            x              o              o           "
	 "%.4lf         %.4lf        %2d"
	 ,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed)); 
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 } //Show_All: End
	 
//Selected combination

     if((Pascal)||(Insert_K)||(Log_Space)||(Binary)) {
	 										 
	 Start=time(NULL); 
	 if(Pascal) {
     tree[0][0].A[0]=Save_0;
     tree[0][0].p=1;
     
   	 //Inherit Save
  	 for(i=1;i<=n;i=i+1)
  	 {	 
	  tree[i][0].A[0]=(tree[i-1][0].A[0]*(n_0+i)+tree[i][0].S)/(n_0+i+1);
	  tree[i][0].p=1;
	  
	  for(j=1;j<i;j++)
	  {
       if(tree[i-1][j-1].p+tree[i-1][j].p>P+1)
   	   {
        tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
        /((n_0+1)+i);
        tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
        /((n_0+1)+i);
        for(k=1;k<=P-1;k++)
        {
         if(Log_Space) tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
		 else tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];
		}
        tree[i][j].p=P+1;
        
        if(Insert_K) {
        if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
        {
         flag=0;
         for(l=0;l<tree[i][j].p;l++)
         {
          if(K==tree[i][j].A[l]) flag=1;
         }
         if(flag==0)
         {
          tree[i][j].A[tree[i][j].p]=K;
          tree[i][j].p++;
          bubble_sort(tree[i][j].A,tree[i][j].p);
         }
	    }
		}
       }
       
       else
       {
        tree[i][j].p=0;
	   
	    //Inherit from the upper
	    for(k=0;k<tree[i-1][j-1].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j-1].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
        }

	    //Inherit from the lower
	    for(k=0;k<tree[i-1][j].p;k++)
	    {
         tree[i][j].A[tree[i][j].p]=(tree[i-1][j].A[k]*(n_0+i)+tree[i][j].S)/(n_0+i+1);
         tree[i][j].p++;
	    }
	    
        if(Insert_K) {
        if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
        {
         flag=0;
         for(l=0;l<tree[i][j].p;l++)
         {
          if(K==tree[i][j].A[l]) flag=1;
         }
         if(flag==0)
         {
          tree[i][j].A[tree[i][j].p]=K;
          tree[i][j].p++;
         }
	    }
	    }
	    
        bubble_sort(tree[i][j].A,tree[i][j].p);
       }
      } 
	  
	  tree[i][i].A[0]=(tree[i-1][i-1].A[0]*(n_0+i)+tree[i][i].S)/(n_0+i+1);
	  tree[i][i].p=1;
   	 }
	 }	  
     
	 else {
     //Calculate average
	 for(i=0;i<=n;i=i+1)
  	 {
	  for(j=0;j<=i;j++)
	  {
       tree[i][j].A[0]=(Save_0*(n_0+1)-S+S*(1-pow(u,i-j+1))/(1-u)+S*pow(u,i-j)*d*(1-pow(d,j))/(1-d))
       /((n_0+1)+i);
       tree[i][j].A[P]=(Save_0*(n_0+1)-S+S*(1-pow(d,j+1))/(1-d)+S*pow(d,j)*u*(1-pow(u,i-j))/(1-u))
       /((n_0+1)+i);
       for(k=1;k<=P-1;k++)
       {
        if(Log_Space) tree[i][j].A[k]=exp(((double)(k)/P)*log(tree[i][j].A[P])+(1-(double)(k)/P)*log(tree[i][j].A[0]));
		else tree[i][j].A[k]=((double)(k)/P)*tree[i][j].A[P]+(1-(double)(k)/P)*tree[i][j].A[0];        
	   }
       tree[i][j].p=P+1;
       
       if(Insert_K) {
       if((tree[i][j].A[P]<=K)&&(K<=tree[i][j].A[0]))
       {
        flag=0;
        for(l=0;l<tree[i][j].p;l++)
        {
         if(K==tree[i][j].A[l]) flag=1;
        }
        if(flag==0)
        {
         tree[i][j].A[tree[i][j].p]=K;
         tree[i][j].p++;
         bubble_sort(tree[i][j].A,tree[i][j].p);
        }
	   }
	   }
      }
     }
	 }
	      
     //Decide terminal node's payoff
   	 for(j=0;j<=n;j++)for(k=0;k<tree[n][j].p;k++)
	 {
      tree[n][j].E_Put[k]=max(K-tree[n][j].A[k],0);
      tree[n][j].A_Put[k]=tree[n][j].E_Put[k];
	 }
	 
	 //Backward induction
	 for(i=n-1;i>=0;i=i-1)
	 {
      for(j=0;j<=i;j++)
   	  {
       for(k=0;k<tree[i][j].p;k++)
       {		
        Au=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j].S)/(n_0+i+2);
        if(Binary) position_u=BinarySearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);
        else position_u=SequentialSearch(Au,tree[i+1][j].A,tree[i+1][j].p-1);

        Ad=((n_0+i+1)*tree[i][j].A[k]+tree[i+1][j+1].S)/(n_0+i+2);
        if(Binary) position_d=BinarySearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1);              
		else position_d=SequentialSearch(Ad,tree[i+1][j+1].A,tree[i+1][j+1].p-1); 
      
        //European     
        Cu=Interpolation(tree[i+1][j].E_Put[position_u],tree[i+1][j].E_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].E_Put[position_d],tree[i+1][j+1].E_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].E_Put[k]=(pu*Cu+pd*Cd)*exp(-r*delta_t);
        
        //American
        Cu=Interpolation(tree[i+1][j].A_Put[position_u],tree[i+1][j].A_Put[position_u+1],tree[i+1][j].A[position_u],Au,tree[i+1][j].A[position_u+1]);
        Cd=Interpolation(tree[i+1][j+1].A_Put[position_d],tree[i+1][j+1].A_Put[position_d+1],tree[i+1][j+1].A[position_d],Ad,tree[i+1][j+1].A[position_d+1]);
        tree[i][j].A_Put[k]=max((pu*Cu+pd*Cd)*exp(-r*delta_t),K-tree[i][j].A[k]);
       }
      }
     } 
     End=time(NULL);
     Speed=difftime(End,Start);
     
	 strcpy(combination,"");

	 if(Pascal) strcat(combination,"     o      ");	 
	 else strcat(combination,"     x      ");

	 if(Insert_K) strcat(combination,"      o        ");
	 else strcat(combination,"      x        ");

	 if(Log_Space) strcat(combination,"      o         ");
	 else strcat(combination,"      x         ");

	 if(Binary) strcat(combination,"     o          ");
	 else strcat(combination,"     x          ");

	 _stprintf(info,"%s %.4lf         %.4lf        %2d",combination,tree[0][0].E_Put[0],tree[0][0].A_Put[0],(int)(Speed));
	 dc.TextOut(x_axis,y_axis,info);
	 y_axis+=20;
	 } //Selected combination
   
//Release occupied memory

	 for(i=0;i<=n;i++)
  	 {
	  for(j=0;j<=i;j++)
	  {
 	   delete[] tree[i][j].A;
	   tree[i][j].A=NULL;
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
	 
	 } //Binomial Tree
	 } //Put: End
}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S = S;
	dlg.K=K;
	dlg.r=r;
	dlg.T=T;
	dlg.Sigma=Sigma;
	dlg.q=q;
	dlg.n=n;
	dlg.p=p;
	dlg.Save_0=Save_0;
	dlg.T_0=T_0;
	dlg.Call=Call;
	dlg.Put=Put;
	dlg.Show_All=Show_All;
	dlg.Pascal=Pascal;
	dlg.Insert_K=Insert_K;
	dlg.Log_Space=Log_Space;
	dlg.Binary=Binary;
	dlg.Monte_Carlo=Monte_Carlo;
	dlg.m=m;
	dlg.With=With;
	dlg.Without=Without;
	dlg.Binomial_Tree=Binomial_Tree;

	if ( dlg.DoModal ( ) == IDOK )
	{
		
		S = dlg.S;
		K=dlg.K;
		r=dlg.r;
		T=dlg.T;
		Sigma=dlg.Sigma;
		q=dlg.q;
		n=dlg.n;
		p=dlg.p;
		Save_0=dlg.Save_0;
		T_0=dlg.T_0;
		Call=dlg.Call;
		Put=dlg.Put;
		Show_All=dlg.Show_All;
		Pascal=dlg.Pascal;
		Insert_K=dlg.Insert_K;
		Log_Space=dlg.Log_Space;
		Binary=dlg.Binary;
		Monte_Carlo=dlg.Monte_Carlo;
		m=dlg.m;
		With=dlg.With;
		Without=dlg.Without;
		Binomial_Tree=dlg.Binomial_Tree;
		Invalidate();
	}
}
	
