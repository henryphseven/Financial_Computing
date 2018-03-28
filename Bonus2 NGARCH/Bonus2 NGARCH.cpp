#include <afxwin.h>
#include "ddxddv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
const int DAYS=200,HALF_RANGE=500,N=10,P=100;
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
//Calculate n!/(ju!jm!jd!)
double Trinomial(int n,int ju,int jd)
{
     int k;
 	 double Sum[4];
 
 	 for(k=0;k<4;k++) Sum[k]=0;
 
 	 for(k=1;k<=n;k++) Sum[0]+=log((double)(k));
 	 for(k=1;k<=ju;k++) Sum[1]+=log((double)(k));
 	 for(k=1;k<=jd;k++) Sum[2]+=log((double)(k));
 	 for(k=1;k<=(n-(ju+jd));k++) Sum[3]+=log((double)(k));
 
 	 return exp(Sum[0]-(Sum[1]+Sum[2]+Sum[3]));
}
//Calculate path probability for j=0,+-1,+-2,...,+-n
void CalculateProb(double P[],double h,double eta,double gamma,double r,double q,int n)
{
     int j;
	 double Pu,Pm,Pd;
     int ju,jm,jd;

     Pu=h/(2*pow(eta,2)*pow(gamma,2))+(r-q-h/2.0)*sqrt(1.0/n)/(2*eta*gamma);
     Pm=1-h/(pow(eta,2)*pow(gamma,2));
     Pd=h/(2*pow(eta,2)*pow(gamma,2))-(r-q-h/2.0)*sqrt(1.0/n)/(2*eta*gamma);

	 for(j=n;j>=-n;j--)
	 {
	  P[n-j]=0;
	  for(ju=0;ju<=n;ju++)
	  {
	   jd=ju-j;
	   jm=n-(ju+jd);
	   if((jm>=0)&&(jd>=0))
	    P[n-j]+=Trinomial(n,ju,jd)*pow(Pu,ju)*pow(Pm,jm)*pow(Pd,jd);
      }
     }
     
     return;
}
//Ordering: Decreasing 
void bubble_sort(int a[],int p) //p: Number of numbers
{
   int i,j;
   int temp;
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
//Get eta
int get_eta(double h0,double h,int n)
{
   double gamma=sqrt(h*n);
   
   if(h0==h) return 1;
   else return int(sqrt(h0)/gamma)+1;
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
class MyFrameWindow : public CFrameWnd
{
	public:

      //Traditional input
 	  double S,K,r,q;
 	  int days,n,p,m;
 	  
 	  //GARCH input
 	  double h,B0,B1,B2,C,lambda;
	  int Monte_Carlo,Tree;

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
		: CWinApp ( "Option Price" )
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

      //Traditional input
 	  double S,K,r,q;
 	  int days,n,p,m;
 	  
 	  //GARCH input
 	  double h,B0,B1,B2,C,lambda;
	  int Monte_Carlo,Tree;

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
	DDX_Text ( pDX, IDC_days, days );
	DDV_MinMaxInt( pDX, days, 1, DAYS );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_n, n );
	DDV_MinMaxInt( pDX, n, 1, N );
	DDX_Text ( pDX, IDC_p, p );
	DDV_MinMaxInt( pDX, p, 1, P );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 10, 10000 );

	DDX_Text ( pDX, IDC_h, h );
	DDX_Text ( pDX, IDC_B0, B0 );
	DDX_Text ( pDX, IDC_B1, B1 );
	DDX_Text ( pDX, IDC_B2, B2 );
	DDX_Text ( pDX, IDC_C, C );
	DDX_Text ( pDX, IDC_lambda, lambda );
	DDX_Check( pDX, IDC_Monte_Carlo, Monte_Carlo );
	DDX_Check( pDX, IDC_Tree, Tree );
}

MyFrameWindow::MyFrameWindow ( )
{
      //Traditional input
 	  S=1000,K=1000,r=0,q=0;
 	  days=3,n=1,p=3,m=10;
 	  
 	  //GARCH input
 	  h=0.0001096,B0=6.575*pow(10.0,-6),B1=0.90,B2=0.04,C=0,lambda=0;
	  Monte_Carlo=0,Tree=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );

	  int i,j,k,l;
      
 	  int L=1000;
 	  double delta_t=1.0/365,T=delta_t*days;

 	  double Cstar=C+lambda;

 	  double St,ht,epsilon;
      double E_Call,Sum_E_Call,Square_E_Call,E_Put,Sum_E_Put,Square_E_Put;
      double STDERR,Upper_Limit,Lower_Limit,STDERR_P,Upper_Limit_P,Lower_Limit_P;
      
      double gamma=sqrt(h*n);
	  double h0,h1,epsilon1;
	  int eta,position;
	  int jj,kk;
	  double Prob[2*N+1],I;

	  int x_axis=20,y_axis=20;
      
      _stprintf(info,"<NGARCH Model>");
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;

	  if(Monte_Carlo) {
   
//Monte Carlo

     srand( (unsigned)time( NULL ) );
     _stprintf(info,"Monte Carlo:");
	 y_axis+=20;
	 dc.TextOut(x_axis,y_axis,info); 
	 y_axis+=20;
     
     do{
	 E_Call=0,Square_E_Call=0;
	 E_Put=0,Square_E_Put=0;
     for(i=0;i<m;i++)
     {
	  Sum_E_Call=0,Sum_E_Put=0;
	  for(j=0;j<L;j++)
	  {
	   St=S; ht=h;
	   for(k=1;k<=days;k++)
	   {
	    epsilon=Normal();
        St=St*exp((r-q-ht/2)+sqrt(ht)*epsilon);
	    ht=B0+B1*ht+B2*ht*pow(epsilon-Cstar,2);
	   }
	   Sum_E_Call+=(exp(-r*T)*max(St-K,0))/L;
	   Sum_E_Put+=(exp(-r*T)*max(K-St,0))/L;
	  }
	  E_Call+=Sum_E_Call;
	  E_Put+=Sum_E_Put;
	  Square_E_Call+=pow(Sum_E_Call,2);
	  Square_E_Put+=pow(Sum_E_Put,2);
	 }
     E_Call/=m;
     E_Put/=m;
     }while((_isnan(E_Call)!=0)||(_finite(E_Call)==0));
     
     STDERR=sqrt((double)((Square_E_Call-m*pow(E_Call,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit=E_Call-2*STDERR;
     Upper_Limit=E_Call+2*STDERR;
     
     STDERR_P=sqrt((double)((Square_E_Put-m*pow(E_Put,2))/(m-1)))/sqrt((double)(m));
     Lower_Limit_P=E_Put-2*STDERR_P;
     Upper_Limit_P=E_Put+2*STDERR_P;

	 _stprintf(info,"Option Type");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"SE");
     dc.TextOut(x_axis+200,y_axis,info);
	 _stprintf(info,"95%% CI");
     dc.TextOut(x_axis+300,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"Call");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"%.4lf",E_Call);
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"%.4lf",STDERR);
     dc.TextOut(x_axis+200,y_axis,info);
	 _stprintf(info,"[%.4lf,%.4lf]",Lower_Limit,Upper_Limit);
     dc.TextOut(x_axis+300,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"Put");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"%.4lf",E_Put);
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"%.4lf",STDERR_P);
     dc.TextOut(x_axis+200,y_axis,info);
	 _stprintf(info,"[%.4lf,%.4lf]",Lower_Limit_P,Upper_Limit_P);
     dc.TextOut(x_axis+300,y_axis,info);
	 y_axis+=20;

	  }

	  if(Tree) {
   
//Multinomial tree

	 _stprintf(info,"Multinomial Tree:");
	 y_axis+=20;
	 dc.TextOut(x_axis,y_axis,info); 
	 y_axis+=20;
     
     struct node 
     {
	   double y;
	   double h[P],E_Call[P],A_Call[P],E_Put[P],A_Put[P]; 
       int h_num;
     } *tree[DAYS+1];
      
     for(i=0;i<=days;i++)
     tree[i]=new struct node[2*HALF_RANGE+1];   
	  
     struct record
     {
      int num;
      int *pos;
	 } *daily;
	 
	 daily=new struct record[DAYS+1];
	 
     for(i=0;i<=days;i++)
     daily[i].pos=new int[2*HALF_RANGE+1];
	 
	 daily[0].pos[0]=HALF_RANGE;
	 daily[0].num=1;
	 
	 tree[0][HALF_RANGE].y=log(S);
	 tree[0][HALF_RANGE].h[0]=h;
	 tree[0][HALF_RANGE].h[p-1]=h;
	 
     try {
	 //Build the tree
	 for(i=0;i<days;i++)
	 {
	  daily[i+1].num=0;
	  for(k=0;k<daily[i].num;k++)
	  {
	   for(l=0;l<p;l++)
	   {
	    tree[i][daily[i].pos[k]].h[l]
        =l*(tree[i][daily[i].pos[k]].h[p-1]-tree[i][daily[i].pos[k]].h[0])/(p-1)
         +tree[i][daily[i].pos[k]].h[0];
	    
	    h0=tree[i][daily[i].pos[k]].h[l];
	    eta=get_eta(h0,h,n);
	    
		//stretch 2*n+1
	    for(j=n;j>=-n;j--)
		{
   		 epsilon1=(j*eta*gamma-(r-q-h0/2.0))/sqrt(h0);
         h1=B0+B1*h0+B2*h0*pow(epsilon1-Cstar,2);
         
         kk=daily[i].pos[k]+j*eta;
         if((kk>2*HALF_RANGE)||(kk<0)) throw kk;

         if(daily[i+1].num==0)
   		 {
		  daily[i+1].pos[daily[i+1].num]=kk;
          daily[i+1].num++;
          
          tree[i+1][kk].y=tree[i][daily[i].pos[k]].y+j*eta*gamma;
		  tree[i+1][kk].h[0]=h1;
          tree[i+1][kk].h_num=1;
		 }
		 else 
		 {	
		  for(jj=0;jj<daily[i+1].num;jj++)
          {
		   if(daily[i+1].pos[jj]==kk)
		   {
		    if(tree[i+1][kk].h_num==1)
		    {
		     if(h1>tree[i+1][kk].h[0])
		     {
		      tree[i+1][kk].h[p-1]=tree[i+1][kk].h[0];
		      tree[i+1][kk].h[0]=h1;
			 }
			 else 
             {
		      tree[i+1][kk].h[p-1]=h1;
			 }
			 tree[i+1][kk].h_num=2;
		    }
		    else if(tree[i+1][kk].h_num==2)
		    {
   		     if(h1>tree[i+1][kk].h[0])
   		     {
			  tree[i+1][kk].h[0]=h1;
			 }
   		     else if(h1<tree[i+1][kk].h[p-1])
   		     {
			  tree[i+1][kk].h[p-1]=h1;
			 }
		    }
	        break;
		   }
		  }
		  if(jj==daily[i+1].num)
		  {
		   daily[i+1].pos[daily[i+1].num]=kk;
           daily[i+1].num++;
          
           tree[i+1][kk].y=tree[i][daily[i].pos[k]].y+j*eta*gamma;
		   tree[i+1][kk].h[0]=h1;
           tree[i+1][kk].h_num=1;
		  }		
		 } 
	    } //jth path
	   } //lth h
	  } //kth node
	  bubble_sort(daily[i+1].pos,daily[i+1].num);
	 } //ith day
	 
     for(k=0;k<daily[days].num;k++)
     {
      for(l=0;l<p;l++)
      {
  	   tree[days][daily[days].pos[k]].h[l]
       =l*(tree[days][daily[days].pos[k]].h[p-1]-tree[days][daily[days].pos[k]].h[0])/(p-1)
        +tree[days][daily[days].pos[k]].h[0];
      }
     }

	 //Decide terminal node's payoff
   	 for(k=0;k<daily[days].num;k++)for(l=0;l<p;l++)
	 {
      tree[days][daily[days].pos[k]].E_Call[l]
      =max(exp(tree[days][daily[days].pos[k]].y)-K,0);
      
      tree[days][daily[days].pos[k]].A_Call[l]
      =tree[days][daily[days].pos[k]].E_Call[l];
      
      tree[days][daily[days].pos[k]].E_Put[l]
      =max(K-exp(tree[days][daily[days].pos[k]].y),0);
      
      tree[days][daily[days].pos[k]].A_Put[l]
      =tree[days][daily[days].pos[k]].E_Put[l];
	 }
	 
	 //Backward induction
	 for(i=days-1;i>=0;i--)
	 { 
      for(k=0;k<daily[i].num;k++)
	  { 
	   for(l=0;l<p;l++)		
	   {
		tree[i][daily[i].pos[k]].E_Call[l]=0;
		tree[i][daily[i].pos[k]].A_Call[l]=0;
		tree[i][daily[i].pos[k]].E_Put[l]=0;
		tree[i][daily[i].pos[k]].A_Put[l]=0;
		
		h0=tree[i][daily[i].pos[k]].h[l];
		eta=get_eta(h0,h,n);
		CalculateProb(Prob,h0,eta,gamma,r,q,n);
		
		for(j=n;j>=-n;j--)
		{
	 	 epsilon1=(j*eta*gamma-(r-q-h0/2.0))/sqrt(h0);
         h1=B0+B1*h0+B2*h0*pow(epsilon1-Cstar,2);
		 kk=daily[i].pos[k]+j*eta;
		 position=BinarySearch(h1,tree[i+1][kk].h,p-1);
		 
         I=Interpolation(tree[i+1][kk].E_Call[position],tree[i+1][kk].E_Call[position+1],tree[i+1][kk].h[position],h1,tree[i+1][kk].h[position+1]); 
         tree[i][daily[i].pos[k]].E_Call[l]+=Prob[n-j]*I;
         
         I=Interpolation(tree[i+1][kk].A_Call[position],tree[i+1][kk].A_Call[position+1],tree[i+1][kk].h[position],h1,tree[i+1][kk].h[position+1]); 
         tree[i][daily[i].pos[k]].A_Call[l]+=Prob[n-j]*I;
         
         I=Interpolation(tree[i+1][kk].E_Put[position],tree[i+1][kk].E_Put[position+1],tree[i+1][kk].h[position],h1,tree[i+1][kk].h[position+1]); 
         tree[i][daily[i].pos[k]].E_Put[l]+=Prob[n-j]*I;
         
         I=Interpolation(tree[i+1][kk].A_Put[position],tree[i+1][kk].A_Put[position+1],tree[i+1][kk].h[position],h1,tree[i+1][kk].h[position+1]); 
         tree[i][daily[i].pos[k]].A_Put[l]+=Prob[n-j]*I;
		}
		
		tree[i][daily[i].pos[k]].E_Call[l]*=exp(-r*delta_t);
		
		tree[i][daily[i].pos[k]].A_Call[l]*=exp(-r*delta_t);
		tree[i][daily[i].pos[k]].A_Call[l]=max(tree[i][daily[i].pos[k]].A_Call[l],exp(tree[i][daily[i].pos[k]].y)-K);
		
		tree[i][daily[i].pos[k]].E_Put[l]*=exp(-r*delta_t);
		
		tree[i][daily[i].pos[k]].A_Put[l]*=exp(-r*delta_t);
		tree[i][daily[i].pos[k]].A_Put[l]=max(tree[i][daily[i].pos[k]].A_Put[l],K-exp(tree[i][daily[i].pos[k]].y));
	   }
      }
     }

	 _stprintf(info,"Option Type");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"European");
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"American");
     dc.TextOut(x_axis+200,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"Call");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"%.4lf",tree[0][daily[0].pos[0]].E_Call[0]);
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"%.4lf",tree[0][daily[0].pos[0]].A_Call[0]);
     dc.TextOut(x_axis+200,y_axis,info);
	 y_axis+=20;

	 _stprintf(info,"Put");
	 dc.TextOut(x_axis,y_axis,info);
	 _stprintf(info,"%.4lf",tree[0][daily[0].pos[0]].E_Put[0]);
     dc.TextOut(x_axis+100,y_axis,info);
	 _stprintf(info,"%.4lf",tree[0][daily[0].pos[0]].A_Put[0]);
     dc.TextOut(x_axis+200,y_axis,info);
	 y_axis+=20;

	 } //try
	 
	 catch(int kk)
	 {
	     if(kk>2*HALF_RANGE) 
		 {   
		  _stprintf(info,"The upper bound is %d, but the node may jump up to %d.",2*HALF_RANGE,kk);
	      dc.TextOut(x_axis,y_axis,info);
	      y_axis+=20;
          _stprintf(info,"Therefore, you should reset a higher upper bound.");
	      dc.TextOut(x_axis,y_axis,info);
	      y_axis+=20;
	     }
         else if(kk<0)
         {
		  _stprintf(info,"The lower bound is 0, but the node may jump down to %d.",kk);
	      dc.TextOut(x_axis,y_axis,info);
	      y_axis+=20;
          _stprintf(info,"Therefore, you should reset a higher starting point.");
	      dc.TextOut(x_axis,y_axis,info);
	      y_axis+=20;
	     }
	 }
	 
//Release occupied memory
     	 
     for(i=0;i<=days;i++)
     {
      delete[] daily[i].pos;
      daily[i].pos=NULL;
	 }
	 
	 delete[] daily;
	 daily=NULL;
	 
	 for(i=0;i<=days;i++)
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
	dlg.K = K;
	dlg.r = r;
	dlg.days = days;
    dlg.q = q;
	dlg.n = n;
	dlg.m = m;
	dlg.p = p;

	dlg.h = h;
	dlg.B0 = B0;
	dlg.B1 = B1;
	dlg.B2 = B2;
    dlg.C = C;
	dlg.lambda = lambda;
	dlg.p = p;
	dlg.Monte_Carlo=Monte_Carlo;
	dlg.Tree=Tree;

	if ( dlg.DoModal ( ) == IDOK )
	{
	 S=dlg.S;
	 K=dlg.K;
	 r=dlg.r;
	 days=dlg.days;
	 q=dlg.q;
	 n=dlg.n;
	 m=dlg.m;
	 p=dlg.p;

	 h=dlg.h;
	 B0=dlg.B0;
	 B1=dlg.B1;
	 B2=dlg.B2;
	 C=dlg.C;
	 lambda=dlg.lambda;
	 Monte_Carlo=dlg.Monte_Carlo;
	 Tree=dlg.Tree;

	 Invalidate();
	}
}
