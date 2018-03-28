#include <afxwin.h>
#include "ddxddv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
const int M=1000;
double Sum;
int index(double a)
{
 	if(a==1) return 1;
 	else return 0; 
}
void swap(double *x,double *y)
{
 	 double temp=*x;
 	 *x=*y;
 	 *y=temp;
 	 return;
}
//Use Gauss-Jordan elimination to solve Tg=a_g
void GJ(double *T[],double a_g[],double g[],int N)
{
 	 int i,j,k;
      
     double pivot;
	 double *G[M-1];
     for(i=0;i<M-1;i++) G[i]=new double[M]; 
     
 	 int p;
 	 for(i=0;i<M-1;i++)for(j=0;j<M;j++) G[i][j]=0;
     for(i=0;i<N;i++)for(j=0;j<N;j++) G[i][j]=*(*(T+i)+j);
     for(i=0;i<N;i++) G[i][M-1]=a_g[i];
     for(j=0;j<N;j++)
     {   
	 	 pivot=0;
	 	 for(i=j;i<N;i++)
	 	 {
		  if(fabs(G[i][j])>pivot) {pivot=G[i][j]; p=i;}
	     }
	     if(p!=j) for(k=0;k<M;k++) swap(&G[p][k],&G[j][k]);
	     for(k=0;k<M;k++) G[j][k]=G[j][k]/pivot;
	     double multiplier;
	     for(i=0;i<N;i++)
		 {
		   multiplier=G[i][j];
		   for(k=0;k<M;k++)
	       {
	        if(i!=j)G[i][k]=G[i][k]-multiplier*G[j][k];
           }
		 }
	 }
	 for(i=0;i<N;i++) g[i]=G[i][M-1];
	 
     for(i=0;i<M-1;i++)
     {
      delete[] G[i]; 
	  G[i]=NULL;
     }
     
	 return;
} 
//Multiply a matrix T with a vector a_g to get g=T*a_g 
void Multiplication(double *T[],double a_g[],double g[],int N)
{
     int i,j;
     
     for(i=0;i<N;i++) g[i]=0;
     for(i=0;i<N;i++)for(j=0;j<N;j++) g[i]+=*(*(T+i)+j)*a_g[j];
     
     return;
}

class MyFrameWindow : public CFrameWnd
{
	public:

     //Input
     double S,K,r,T,Sigma,q,multiple;
     int m,n; //m: partition of the stock price, n: partition of time
     int Implicit,Explicit;

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

     //Input
     double S,K,r,T,Sigma,q,multiple;
     int m,n; //m: partition of the stock price, n: partition of time
     int Implicit,Explicit;

	MyDialog ( CWnd* parent ) : CDialog ( IDD_INPUT, parent ) {}
	void DoDataExchange( CDataExchange* pDX );
};

void MyDialog::DoDataExchange( CDataExchange* pDX )
{
	CDialog::DoDataExchange ( pDX );

    DDX_Text ( pDX, IDC_S, S );
	DDX_Text ( pDX, IDC_multiple, multiple );
    DDX_Text ( pDX, IDC_K, K );
	DDX_Text ( pDX, IDC_r, r );
	DDV_MinMaxDouble( pDX, r, 0, 1 );
	DDX_Text ( pDX, IDC_T, T );
	DDX_Text ( pDX, IDC_Sigma, Sigma );
	DDV_MinMaxDouble( pDX, Sigma, 0, 1 );
	DDX_Text ( pDX, IDC_q, q );
	DDV_MinMaxDouble( pDX, q, 0, 1 );
	DDX_Text ( pDX, IDC_m, m );
	DDV_MinMaxInt( pDX, m, 1, M );
	DDX_Text ( pDX, IDC_n, n );
	DDX_Check ( pDX, IDC_Implicit, Implicit );
	DDX_Check ( pDX, IDC_Explicit, Explicit );
}

MyFrameWindow::MyFrameWindow ( )
{
     //Input
     S=50,K=50,r=0.1,T=5.0/12,Sigma=0.4,q=0,multiple=2;
     m=20,n=10; //m: partition of the stock price, n: partition of time
     Implicit=0,Explicit=0;
}

void MyFrameWindow::OnPaint( )
{
	char info[128];
	CPaintDC dc( this );
	
	 int i,j,k,l;
     
     double Smax=multiple*S,delta_S=Smax/m,delta_t=T/n;

     struct node
     {
      double E_Call,A_Call,E_Put,A_Put;
     } *f[2];
     
     for(i=0;i<2;i++) f[i]=new struct node[M+1];
     
     int sign,a_sign;
	 double a[M],b[M],c[M]; 
	 double g[M-1],a_g[M-1];
	 double *C[M-1];
     for(i=0;i<M-1;i++) C[i]=new double[M-1]; 
     
     int x_axis=20,y_axis=20;
     
     _stprintf(info,"<Finite Difference Method>");
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     if(Implicit) { //Implicit: Start
     
     _stprintf(info,"Implicit:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
   
     sign=index(pow(-1.0,n));
     
	 for(j=0;j<=m;j++) 
     {
      f[sign][j].E_Call=max(j*delta_S-K,0);
      f[sign][j].A_Call=f[sign][j].E_Call;
      
      f[sign][j].E_Put=max(K-j*delta_S,0);
      f[sign][j].A_Put=f[sign][j].E_Put;
     }

     for(i=n-1;i>=0;i--)
     {
      sign=index(pow(-1.0,i));
      a_sign=index(pow(-1.0,i+1));
      
  	  f[sign][0].E_Call=0; f[sign][0].A_Call=f[sign][0].E_Call;
	  f[sign][m].E_Call=max(Smax-K,0); f[sign][m].A_Call=f[sign][m].E_Call;
      
      f[sign][0].E_Put=K; f[sign][0].A_Put=f[sign][0].E_Put;
	  f[sign][m].E_Put=0; f[sign][m].A_Put=f[sign][m].E_Put;
	  
 	  for(j=1;j<=m-1;j++)
 	  {
	   a[j]=(1.0/2)*(r-q)*j*delta_t-(1.0/2)*Sigma*Sigma*j*j*delta_t;
	   b[j]=1+Sigma*Sigma*j*j*delta_t+r*delta_t;
	   c[j]=-(1.0/2)*(r-q)*j*delta_t-(1.0/2)*Sigma*Sigma*j*j*delta_t;
	  }
	  
	  for(k=0;k<m-1;k++)for(l=0;l<m-1;l++) C[k][l]=0;
	  C[0][0]=b[m-1];
	  C[0][1]=a[m-1];
      for(k=1;k<m-2;k++)
	  {
       C[k][k-1]=c[(m-1)-k];
       C[k][k]=b[(m-1)-k];
       C[k][k+1]=a[(m-1)-k];
      }
      C[m-2][m-3]=c[1];
      C[m-2][m-2]=b[1];
	  
	  //European call
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].E_Call;
	  a_g[0]-=c[m-1]*f[sign][m].E_Call;
	  a_g[m-2]-=a[1]*f[sign][0].E_Call;
	  
	  GJ(C,a_g,g,m-1);
	  
	  for(j=0;j<m-1;j++) 
       f[sign][(m-1)-j].E_Call=max(g[j],0);
       
      //American call
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].A_Call;
	  a_g[0]-=c[m-1]*f[sign][m].A_Call;
	  a_g[m-2]-=a[1]*f[sign][0].A_Call;
	  
	  GJ(C,a_g,g,m-1);
	  
	  for(j=0;j<m-1;j++) 
       f[sign][(m-1)-j].A_Call=max(max(g[j],((m-1)-j)*delta_S-K),0);
      
      //European put
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].E_Put;
	  a_g[0]-=c[m-1]*f[sign][m].E_Put;
	  a_g[m-2]-=a[1]*f[sign][0].E_Put;
	  
	  GJ(C,a_g,g,m-1);
	  
	  for(j=0;j<m-1;j++) 
       f[sign][(m-1)-j].E_Put=max(g[j],0);
       
      //American put
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].A_Put;
	  a_g[0]-=c[m-1]*f[sign][m].A_Put;
	  a_g[m-2]-=a[1]*f[sign][0].A_Put;
	  
	  GJ(C,a_g,g,m-1);
	  
	  for(j=0;j<m-1;j++) 
       f[sign][(m-1)-j].A_Put=max(max(g[j],K-((m-1)-j)*delta_S),0);
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
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].E_Call);
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].A_Call);
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Put");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].E_Put);
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].A_Put);
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;

	 } //Implicit: End  

	 if(Explicit) { //Explicit: Start
     
     _stprintf(info,"Explicit:");
     y_axis+=20;
     dc.TextOut(x_axis,y_axis,info);
     y_axis+=20;
     
     try {
     if(0.5*pow(Sigma*S/delta_S,2)*delta_t>0.5) throw 1;

	 sign=index(pow(-1.0,n));
     
	 for(j=0;j<=m;j++) 
     {
      f[sign][j].E_Call=max(j*delta_S-K,0);
      f[sign][j].A_Call=f[sign][j].E_Call;
      
      f[sign][j].E_Put=max(K-j*delta_S,0);
      f[sign][j].A_Put=f[sign][j].E_Put;
     }
     
     for(i=n-1;i>=0;i--)
     {
      sign=index(pow(-1.0,i));
      a_sign=index(pow(-1.0,i+1));
      
  	  f[sign][0].E_Call=0; f[sign][0].A_Call=f[sign][0].E_Call;
	  f[sign][m].E_Call=max(Smax-K,0); f[sign][m].A_Call=f[sign][m].E_Call;
      
      f[sign][0].E_Put=K; f[sign][0].A_Put=f[sign][0].E_Put;
	  f[sign][m].E_Put=0; f[sign][m].A_Put=f[sign][m].E_Put;
	  
 	  for(j=1;j<=m-1;j++)
 	  {
	   a[j]=(-(1.0/2)*(r-q)*j*delta_t+(1.0/2)*Sigma*Sigma*j*j*delta_t)/(1+r*delta_t);
	   b[j]=(1-Sigma*Sigma*j*j*delta_t)/(1+r*delta_t);
	   c[j]=((1.0/2)*(r-q)*j*delta_t+(1.0/2)*Sigma*Sigma*j*j*delta_t)/(1+r*delta_t);
	  }
	  
	  for(k=0;k<m-1;k++)for(l=0;l<m-1;l++) C[k][l]=0;
	  C[0][0]=b[m-1];
	  C[0][1]=a[m-1];
      for(k=1;k<m-2;k++)
	  {
       C[k][k-1]=c[(m-1)-k];
       C[k][k]=b[(m-1)-k];
       C[k][k+1]=a[(m-1)-k];
      }
      C[m-2][m-3]=c[1];
      C[m-2][m-2]=b[1];
	  
	  //European call
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].E_Call;
	  
	  Multiplication(C,a_g,g,m-1);
	  
      for(j=1;j<m-2;j++) f[sign][(m-1)-j].E_Call=max(g[j],0);
      f[sign][m-1].E_Call=max(g[0]+c[m-1]*f[a_sign][m].E_Call,0);
      f[sign][1].E_Call=max(g[m-2]+a[1]*f[a_sign][0].E_Call,0);
       
      //American call
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].A_Call;
	  
	  Multiplication(C,a_g,g,m-1);
	  
      for(j=1;j<m-2;j++) f[sign][(m-1)-j].A_Call=max(max(g[j],((m-1)-j)*delta_S-K),0);
      f[sign][m-1].A_Call=max(max(g[0]+c[m-1]*f[a_sign][m].A_Call,(m-1)*delta_S-K),0);
      f[sign][1].A_Call=max(max(g[m-2]+a[1]*f[a_sign][0].A_Call,delta_S-K),0);
      
      //European put
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].E_Put;
	  
	  Multiplication(C,a_g,g,m-1);
	  
      for(j=1;j<m-2;j++) f[sign][(m-1)-j].E_Put=max(g[j],0);
      f[sign][m-1].E_Put=max(g[0]+c[m-1]*f[a_sign][m].E_Put,0);
      f[sign][1].E_Put=max(g[m-2]+a[1]*f[a_sign][0].E_Put,0);
       
      //American put
	  for(j=0;j<m-1;j++) a_g[j]=f[a_sign][(m-1)-j].A_Put;
	  
	  Multiplication(C,a_g,g,m-1);
	  
      for(j=1;j<m-2;j++) f[sign][(m-1)-j].A_Put=max(max(g[j],K-((m-1)-j)*delta_S),0);
      f[sign][m-1].A_Put=max(max(g[0]+c[m-1]*f[a_sign][m].A_Put,K-(m-1)*delta_S),0);
      f[sign][1].A_Put=max(max(g[m-2]+a[1]*f[a_sign][0].A_Put,K-delta_S),0);
      
      //Check convergence
      for(j=1;j<m-2;j++)
      {
       if((f[sign][j].E_Call>Smax)||(f[sign][j].A_Call>Smax)||(f[sign][j].E_Put>K)||(f[sign][j].A_Put>K))
	   throw 2;
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
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].E_Call);
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].A_Call);
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;
     
     _stprintf(info,"Put");
     dc.TextOut(x_axis,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].E_Put);
     dc.TextOut(x_axis+100,y_axis,info);
     _stprintf(info,"%.4lf",f[sign][int(m/multiple)].A_Put);
     dc.TextOut(x_axis+200,y_axis,info);
     y_axis+=20;
	 } //try
	 
	 catch(int diverge)
	 {
      _stprintf(info,"The explicit solution will diverge!"); 
	  dc.TextOut(x_axis,y_axis,info);
	  y_axis+=20;
	  _stprintf(info,"Error code: %d",diverge);
	  dc.TextOut(x_axis,y_axis,info);
	 }
	 
	 } //Explicit: End
	 
//Release occupied memory

     for(i=0;i<2;i++) 
	 {
      delete[] f[i];
      f[i]=NULL;
     }
     
	 for(i=0;i<M-1;i++)
     {
      delete[] C[i]; 
	  C[i]=NULL;
     }

}

void MyFrameWindow::OnInput( )
{
	MyDialog dlg ( this );

	dlg.S=S;
	dlg.multiple=multiple;
	dlg.K=K;
	dlg.r=r;
	dlg.T=T;
	dlg.Sigma=Sigma;
	dlg.q=q;
	dlg.m=m;
	dlg.n=n;
	dlg.Implicit=Implicit;
	dlg.Explicit=Explicit;

	if ( dlg.DoModal ( ) == IDOK )
	{
		S=dlg.S;
		multiple=dlg.multiple;
		K=dlg.K;
		r=dlg.r;
		T=dlg.T;
        Sigma=dlg.Sigma;
		q=dlg.q;
		m=dlg.m;
		n=dlg.n;
		Implicit=dlg.Implicit;
		Explicit=dlg.Explicit;
		Invalidate();
	}
}
	
