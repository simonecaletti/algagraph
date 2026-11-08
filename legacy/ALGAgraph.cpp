#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TF2.h>
#include <TF3.h>
#include <TLine.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixDBase.h>
#include <TAxis.h>
#include "Riostream.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include <TDecompLU.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixT.h>
using namespace std;

TApplication *app;


class alga{


private:
  string  m_str;
  bool    m_quadrics;
  int     m_n;
  void    m_matrixForm();
  TMatrixDSym  *m_A,*m_B;
 
public: 
  //le funzioni
  alga(string="x^2+y^2-1");
  ~alga();
  TVectorD eigenVal();
  TMatrixD eigenVect();
  vector <double> Invar();
  void normalForm(vector<double> inv, TVectorD eigenVal);
  void type(vector<double>);
  void draw();
  void drawNF(vector<double> inv, TVectorD eigenVal);
  void rewrite(string eq);
  int n();
};

alga::alga(string eq){                     //costruttore alga

  m_A = NULL;
  m_B = NULL;

  m_str = eq;
  //quadrica o conica?
  string z="z";
  m_quadrics=false;
  m_n=3;

  // if (eq.find("z")!=string::npos)
  for(unsigned int i=0; i<eq.length();i++)
    if(eq[i]==z){
     m_quadrics=true;
     m_n=4;
    }
  m_matrixForm();
}

int alga::n(){
  int n=m_n-1; 
 return n;
}

void alga::m_matrixForm(){        //metodo private matrixForm

  string plus="+",minus="-";
  vector <string> mono; 
  int lung=m_str.size(),k=0;
  for(int i=1;i<lung;i++){
    if(m_str[i]==plus || m_str[i]==minus){  
      mono.push_back(m_str.substr(k,i-k)); 
      k=i;
    }
  }
  mono.push_back(m_str.substr(k,m_str.size()));

  //setta il numero di coeff (lungh di a[]) e altra roba
  int c;
  vector<string> test;
  if(m_quadrics==true){    
    c=10;
    test.push_back("x^2");
    test.push_back("y^2");
    test.push_back("z^2");
    test.push_back("x*y");
    test.push_back("x*z");
    test.push_back("y*z");
    test.push_back("x");
    test.push_back("y");
    test.push_back("z");
  }
  if(m_quadrics==false){
    c=6;
    test.push_back("x^2");
    test.push_back("x*y");
    test.push_back("y^2");
    test.push_back("x");
    test.push_back("y");
  }

  double *a = new double[c];
  for(int i=0;i<c;i++)     
    a[i]=0;

  size_t pos=0;
  string per="*";             
  string str;                
  for (int i=0;i<c-1;i++){
    unsigned int j;
    bool found=false;
    for(j=0;j<mono.size();j++){
      if (mono[j].find(test[i])!=string::npos){
        pos = mono[j].find(test[i]);
	if (pos!=0){
	  if (mono[j][pos]=='*') pos = pos-1;
	  str  = mono[j].substr(0,pos);
          if (str=="+"){
            a[i]=1;
	  } else if (str=="-"){
            a[i]=-1;
	  } else {
	    a[i] = atof(str.c_str());
	  }
	} else{ 
          a[i]=1;
	}
        found=true;
	break;
      }
    }
    if (found){
      mono.erase(mono.begin()+j);
    }
  }
  if (mono.size()>0){
    a[c-1]=atof(mono[0].c_str());
  } 

  //assengno ora i valori alla matrice
  if(m_quadrics==true){
    m_A = new TMatrixDSym(m_n);
    (*m_A)(0,0)=a[0];
    (*m_A)(0,1)=(a[3])/2;
    (*m_A)(1,1)=a[1];
    (*m_A)(0,2)=(a[4])/2;
    (*m_A)(1,2)=(a[5])/2;
    (*m_A)(2,2)=a[2];
    (*m_A)(1,0)=(a[3])/2;
    (*m_A)(2,0)=(a[4])/2;
    (*m_A)(2,1)=(a[5])/2;
    (*m_A)(3,0)=(a[6])/2;
    (*m_A)(0,3)=(a[6])/2;
    (*m_A)(3,1)=(a[7])/2;
    (*m_A)(1,3)=(a[7])/2;
    (*m_A)(3,2)=(a[8])/2;
    (*m_A)(2,3)=(a[8])/2;
    (*m_A)(3,3)= a[9] ;
    
    m_B = new TMatrixDSym(m_n-1);
    (*m_B)(0,0)=a[0];
    (*m_B)(0,1)=(a[3])/2;
    (*m_B)(1,0)=(a[3])/2;
    (*m_B)(1,1)=a[1];
    (*m_B)(1,2)=(a[5])/2;
    (*m_B)(2,1)=(a[5])/2;
    (*m_B)(2,2)=a[2];
    (*m_B)(0,2)=(a[4])/2;
    (*m_B)(2,0)=(a[4])/2;
  }
  if(m_quadrics==false){
    m_A = new TMatrixDSym(m_n);
    (*m_A)(0,0)=a[0];
    (*m_A)(0,1)=(a[1])/2;
    (*m_A)(1,1)=a[2];
    (*m_A)(0,2)=(a[3])/2;
    (*m_A)(1,2)=(a[4])/2;
    (*m_A)(2,2)=a[5];
    (*m_A)(1,0)=(a[1])/2;
    (*m_A)(2,0)=(a[3])/2;
    (*m_A)(2,1)=(a[4])/2;
    
   m_B = new TMatrixDSym(m_n-1);
    (*m_B)(0,0)=a[0];
    (*m_B)(0,1)=(a[1])/2;
    (*m_B)(1,1)=a[2];
    (*m_B)(1,0)=(a[1])/2;
  }
  return;
}

TVectorD alga::eigenVal(){     //metodo autovalori
  const TMatrixDSymEigen eigen(*m_B);                //quadriche & coniche
  const TVectorD eigenVal = eigen.GetEigenValues();
  return eigenVal;
}

TMatrixD alga::eigenVect(){                //metodo autovettori
  const TMatrixDSymEigen eigen2(*m_B);                     //quadriche & coniche
  const TMatrixD  P=eigen2.GetEigenVectors();
  return P;
}

vector <double> alga::Invar(){ //metodo Invarianti
    vector <double> inv;
 if(m_quadrics==true){       //quadriche
    double I4, I3, I2=0, I1=0;
    TMatrixD Q(2,2);
    I4=(*m_A).Determinant();
    inv.push_back(I4);
    I3=(*m_B).Determinant();
    inv.push_back(I3);                     
    Q(0,0)=(*m_B)(0,0);
    Q(0,1)=(*m_B)(0,1);
    Q(1,0)=(*m_B)(0,1);
    Q(1,1)=(*m_B)(1,1);
    I2 += (Q).Determinant();
    Q(0,0)=(*m_B)(0,0);
    Q(0,1)=(*m_B)(0,2);
    Q(1,0)=(*m_B)(0,2);
    Q(1,1)=(*m_B)(2,2);
    I2 += (Q).Determinant();
    Q(0,0)=(*m_B)(1,1);
    Q(0,1)=(*m_B)(1,2);
    Q(1,0)=(*m_B)(1,2);
    Q(1,1)=(*m_B)(2,2);
    I2 += (Q).Determinant();
    inv.push_back(I2);
    for(int i=0;i<m_n-1;i++)
      I1 += (*m_B)(i,i);
    inv.push_back(I1);
 }
  else{                     //coniche
    double I3, I2, I1=0; 
    I3=(*m_A).Determinant();
    inv.push_back(I3);
    I2=(*m_B).Determinant();
    inv.push_back(I2);
    for(int i=0;i<m_n-1;i++)
      I1 += (*m_B)(i,i);
    inv.push_back(I1);
  }
  return inv;
}

void alga::normalForm(vector<double> inv, TVectorD eigenVal){  //metodo Forma Canonica
  if(m_quadrics==true){ 
    vector <double> autoval;          //quadriche
    for(int i=0;i<eigenVal.GetNoElements();i++)
      if(eigenVal[i]!=0)
	autoval.push_back(eigenVal[i]);

    if(inv[0]!=0){  
      vector <string> stampa;
      if(autoval.size()==2) { 
	if(autoval[1]*sqrt(-(autoval[1]*autoval[0])/inv[0])>0)
	  stampa.push_back("*x^2+");
	else
	  stampa.push_back("*x^2");

	stampa.push_back("*y^2");

	cout<<"The equation in normal form is:"<<endl;    
	cout<<autoval[0]*sqrt(-(autoval[0]*autoval[1])/inv[0])<<stampa[0];
	cout<<autoval[1]*sqrt(-(autoval[1]*autoval[0])/inv[0])<<stampa[1];
	cout<<"+2*z = 0"<<endl;
      }
      if(autoval.size()==3){    
	cout<<"The equation in normal form is:"<<endl;    
	if(autoval[1]>0)
	  stampa.push_back("*x^2+");
	else
	  stampa.push_back("*x^2");	
	if(autoval[2]>0)
	  stampa.push_back("*y^2+");
	else
	  stampa.push_back("*y^2");
	if(inv[0]/inv[1]>0)
	  stampa.push_back("*z^2+");
	else
	  stampa.push_back("*z^2");
     
	for(int i=0;i<3;i++)     
	  cout<<autoval[i]<<stampa[i];

	cout<<inv[0]/inv[1]<<" = 0"<<endl;
      }
    }
    else
      {
	cout<<"You have a degenerate case, I do not restore normal form."<<endl;
      }    
  }
  else{                          //coniche
    if(inv[0]!=0){  
      if(eigenVal[0]==0) {    
	cout<<"The equation in normal form is:"<<endl;    
	cout<<eigenVal[1]<<"*x^2+"<<2*sqrt(abs(inv[0]/eigenVal[1]))<<"*y = 0"<<endl;
      }
      if(eigenVal[1]==0){    
	cout<<"The equation in normal form is:"<<endl;    
	cout<<eigenVal[0]<<"*x^2+"<<2*sqrt(abs(inv[0]/eigenVal[0]))<<"*y = 0"<<endl;
      }
      else{
	vector <string> stampa;
	if(eigenVal[1]>0)
	  stampa.push_back("*x^2+");
	else
	  stampa.push_back("*x^2");
	if((inv[0]/inv[1])>0)
	  stampa.push_back("*y^2+");
	else
	  stampa.push_back("*y^2");

	  cout<<"The equation in normal form is:"<<endl;
	  cout<<eigenVal[0]<<stampa[0]<<eigenVal[1]<<stampa[1]<<inv[0]/inv[1]<<" = 0"<<endl;
	
      }
    }
    else
      {
	cout<<"You have a degenerate case, I do not restore normal form."<<endl;
      }
  }
   
}

void alga::type(vector<double> inv){        //metodo Studio della conica/quadrica
  if(m_quadrics==true){             //quadrica
    if(inv[0]!=0){  
      if(inv[1]==0){
	if(inv[0]<0)
	  cout<<"It's an elliptical paraboloid."<<endl;
	if(inv[0]>0)
	  cout<<"It is an hyperbolic paraboloid (saddle paraboloid)."<<endl;
      }
      else{
	if(inv[0]>0){
	  if(inv[2]>0 && inv[3]*inv[1]>0)
	    cout<<"It's ah ellipsoid in complex points."<<endl;
	  else
	    cout<<"It's an hyperboloid of one sheet."<<endl;
	}
	if(inv[0]<0){
	  if(inv[2]>0 && inv[3]*inv[1]>0)
	    cout<<"It's an ellipsoid in real points."<<endl;
	  else
	    cout<<"It's an elliptical hyperboloid (of two sheets)."<<endl;
	}
      }
    }
  }
  else{                        //conica
    if(inv[0]==0){
      if(inv[1]<0)
	cout<<"Are two distinct real straight lines."<<endl;
      if(inv[1]==0)
	cout<<"Pair of parallel/coincident/complex conjugate straight lines."<<endl;
      if(inv[1]>0)
	cout<<"Pair of conjugate imaginary straight lines."<<endl;
    }
    else{
      if(inv[1]<0){
	if(inv[2]==0)
	  cout<<"It is an equilateral hyperbola."<<endl;
	else
	  cout<<"It is a generic hyperbola."<<endl;
      }
      if(inv[1]==0){
	cout<<"It's a generic parabola."<<endl;     
      }
      if(inv[1]>0){
	if(inv[2]*inv[0]<0)
	  cout<<"It is an ellipse in real points."<<endl;
	if(inv[2]*inv[0]>0)
	  cout<<"It is an ellipse in complex points"<<endl;
      }
    }
  }
}

void alga::draw(){                 //metodo di stampa
  if(m_quadrics==true){      //quadriche
    double xmin, xmax, ymin, ymax,zmin,zmax;
    cout<< "Define the range in which to graph the quadric: (xmin,xmax,ymin,ymax,zmin,zmax)." <<endl;
    cin>> xmin>>xmax>>ymin>>ymax>>zmin>>zmax;
    TF3 f("f",m_str.c_str(),xmin,xmax,ymin,ymax,zmin,zmax);
    //stampa la quadrica
    TCanvas c1;
    c1.SetTitle("Quadric");
    double val[1];
    f.SetContour(1,val);
    f.Draw("CONT1");
    app->Run(true);
    c1.Close();
  } 
  else {                  //coniche
    double xmin, xmax, ymin, ymax;
    cout<< "Define the range in which to graph the conic: (xmin,xmax,ymin,ymax)." <<endl;
    cin>> xmin>>xmax>>ymin>>ymax;
    TF2 f("f",m_str.c_str(),xmin,xmax,ymin,ymax);
    //stampa la conica tagliata a z=0
    TCanvas c1;
    c1.SetTitle("Conic");
    f.SetNpx(1000);
    f.SetNpy(1000);
    double val[1];
    f.SetContour(1,val);
    f.Draw("CONT1");
    TLine x(xmin,0,xmax,0);
    x.Draw("SAME");
    TLine y(0,ymin,0,ymax);
    y.Draw("SAME");
    app->Run(true);
    c1.Close();
  }
  return;
}
 
void alga::drawNF(vector<double> inv, TVectorD eigenVal){   //metodo draw FC
  if(m_quadrics==true){            //quadriche
    vector <double> autoval;
    for(int i=0;i<3;i++)
      if(eigenVal[i]!=0)
	autoval.push_back(eigenVal[i]);

    if(autoval.size()==2 && inv[0]!=0){          //paraboloide
      double xmin, xmax, ymin, ymax,zmin,zmax;
      cout<< "Define the range in which to graph the quadric: (xmin,xmax,ymin,ymax,zmin,zmax)." <<endl;
      cin>> xmin>>xmax>>ymin>>ymax>>zmin>>zmax;
      TF3 fc("fc","[0]*x^2+[1]*y^2+2*z^2",xmin,xmax,ymin,ymax,zmin,zmax);
      TCanvas c2;
      c2.SetTitle("Normal Form");
      fc.SetParameter(0,autoval[0]*sqrt(-(autoval[0]*autoval[1])/inv[0]));
      fc.SetParameter(1,autoval[1]*sqrt(-(autoval[1]*autoval[0])/inv[0]));
      double val[1];
      fc.SetContour(1,val);
      fc.Draw("CONT1");
      app->Run(true);
      c2.Close();
    }
    if(autoval.size()==3 && inv[0]!=0)
      {
	double xmin, xmax, ymin, ymax,zmin,zmax;
	cout<< "Define the range in which to graph the quadric: (xmin,xmax,ymin,ymax,zmin,zmax)." <<endl;
	cin>> xmin>>xmax>>ymin>>ymax>>zmin>>zmax;
	TF3 fc("fc","[0]*x^2+[1]*y^2+[2]*z^2+[3]",xmin,xmax,ymin,ymax,zmin,zmax);	
	TCanvas c2;
	c2.SetTitle("Normal Form");
	fc.SetParameter(0,autoval[0]);
	fc.SetParameter(1,autoval[1]);
	fc.SetParameter(2,autoval[2]);
	fc.SetParameter(3,inv[0]/inv[1]);
	double val[1];
	fc.SetContour(1,val);
	fc.Draw("CONT1");
	app->Run(true);
	c2.Close();
      }
  }
  if(m_quadrics==false){           //coniche
    if(inv[1]!=0 && inv[0]!=0){
      double xmin,xmax,ymin,ymax;
      cout<< "Define the range in which to graph the conic: (xmin,xmax,ymin,ymax)." <<endl;
      cin>> xmin>>xmax>>ymin>>ymax;
      TF2 fc("fc","[0]*x^2+[1]*y^2+[2]",xmin,xmax,ymin,ymax);
      TCanvas c2;
      c2.SetTitle("Normal Form");
      fc.SetParameter(0,eigenVal[0]);
      fc.SetParameter(1,eigenVal[1]);
      fc.SetParameter(2,inv[0]/inv[1]);
      fc.SetNpx(1000);
      fc.SetNpy(1000);
      double val1[1];
      fc.SetContour(1,val1);
      fc.Draw("CONT1");
      TLine x(xmin,0,xmax,0);
      x.Draw("SAME");
      TLine y(0,ymin,0,ymax);
      y.Draw("SAME");
      app->Run(true);
      c2.Close();
    }
    if(inv[1]==0 && inv[0]!=0)
      {
	double xmin,xmax,ymin,ymax;
	cout<< "Define the range in which to graph the conic: (xmin,xmax,ymin,ymax)." <<endl;
	cin>> xmin>>xmax>>ymin>>ymax;
	TF2 fc("fc","[0]*x^2+[1]*y",xmin,xmax,ymin,ymax);
	if(eigenVal[0]==0){
	  TCanvas c2;
	  c2.SetTitle("Normal Form");
	  fc.SetParameter(0,eigenVal[1]);
	  fc.SetParameter(1,sqrt(abs(inv[0]/eigenVal[1])));
	  fc.SetNpx(1000);
	  fc.SetNpy(1000);
	  double val1[1];
	  fc.SetContour(1,val1);
	  fc.Draw("CONT1");
	  TLine x(xmin,0,xmax,0);
	  x.Draw("SAME");
	  TLine y(0,ymin,0,ymax);
	  y.Draw("SAME");
	  app->Run(true);
	  c2.Close();
	}
	
	if(eigenVal[1]==0){
	  TCanvas c2;
	  c2.SetTitle("Normal Form");
	  fc.SetParameter(0,eigenVal[0]);
	  fc.SetParameter(1,(-1)*sqrt(abs(inv[0]/eigenVal[0])));
	  fc.SetNpx(1000);
	  fc.SetNpy(1000);
	  double val1[1];
	  fc.SetContour(1,val1);
	  fc.Draw("CONT1");
	  TLine x(xmin,0,xmax,0);
	  x.Draw("SAME");
	  TLine y(0,ymin,0,ymax);
	  y.Draw("SAME");
	  app->Run(true);
	  c2.Close();
	}
	
      }
  }
}

void alga::rewrite(string eq){

  m_A=NULL;
  m_B=NULL;

  m_str = eq;
  //quadrica o conica?
  string z="z";
  m_quadrics=false;
  m_n=3;
  for(unsigned int i=0; i<eq.length();i++)
    if(eq[i]==z){
     m_quadrics=true;
     m_n=4;
    }   
  m_matrixForm();
}

alga::~alga(){             //distruttore classe alga
  if (m_A)  delete m_A;
  if (m_B)  delete m_B;
}

int main(){

  app = new TApplication("app",0,NULL);

  int times=1;
  //inserisco l'equazione
  vector <string> eq;
  string pass;
  cout<<"\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"\t\t\tALGAgraph ver 1.0:\n\twrite 'help' after entering the first equation for more info."<<endl;
  cout<<"\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  cout<<"\nALGAgraph[0]->Enter the eq. :"; 
  cin >> pass; 
  eq.push_back(pass);

  alga q(pass);
  int Neq=0;
  int pos=Neq;
  //calcolo a prescindere dei valori utili
  vector <double> inv= q.Invar();
  vector <TVectorD> eVal;
  vector <TMatrixD> P;
  eVal.push_back(q.eigenVal());
  P.push_back(q.eigenVect());

  //inserire il comando desiderato
  while(1){
    cout<<"ALGAgraph["<<times<<"]->";
    string command;
    cin>>command;
    bool find_command=false;    

    if(command=="eigenVal"){
      cout<<"EigenValues:"<<endl;
      for(int i=0; i<q.n()/*eVal.GetNoElements()*/; i++)
	cout<<eVal[pos][i]<<endl;
      find_command=true;
    }
    if(command=="eigenVect"){
      cout<<"EigenVectors:"<<endl;
      for(int i=0;i<q.n();i++){
    	for(int j=0;j<q.n();j++) {  
    	  cout<< P[pos](i,j)<<" " ;
    	}
    	cout<<endl;
      }
      find_command=true;
    }
    if(command=="invar"){
      cout<<"Invariants:"<<endl;
      for(int i=0; i<inv.size(); i++)
    	cout<<"I"<<inv.size()-i<<": "<<inv[i]<<endl;
      find_command=true;
    }
    if(command=="normalForm"){   
      q.normalForm(inv,eVal[pos]);
      find_command=true;
    }
    if(command=="draw"){
      q.draw();
      find_command=true;
    }
    if(command=="drawNF"){
      q.drawNF(inv,eVal[pos]);
      find_command=true;
    }
    if(command=="type"){
      q.type(inv);
      find_command=true;
    }
    if(command=="new"){
      string pass;
      cout<<"Enter the eq. : ";
      cin>>pass;
      q.rewrite(pass);
      inv  = q.Invar();
      eVal.push_back(q.eigenVal());
      P.push_back(q.eigenVect());
      eq.push_back(pass);
      Neq++;
      pos=Neq;
      find_command=true;
    }
    if(command==".q"){
      cout<<"Thanks for using ALGAgraph."<<endl;
      find_command=true;     
      break;
    }
    //comandi di "estetica"
    if(command=="list"){
      for(int i=0;i<eq.size();i++)
	cout<<i+1<<"°-> "<<eq[i]<<endl;
      find_command=true;    
    }
    if(command=="select"){
      int sel;
      cin>>sel;
      if(sel>eq.size())
	cout<<"The selected equation does not exist"<<endl;
      else{
	pass=eq[sel-1];
	q.rewrite(pass);
	inv= q.Invar();
	pos=sel;
      }
      find_command=true;   
    }
    if(command=="clear"){
      eq.clear();
      string pass="x^2+y^2-1";
      q.rewrite(pass);
      inv= q.Invar();
      eVal.clear();
      P.clear();
      Neq=0;
      pos=Neq;
      find_command=true;    
    }
    if(command=="help"){
      cout<<"->type:\nReturns output the type of conic/ quadric ( ellipse , hyperbola ,elliptical paraboloid , etc. ).\nIn the case of a quadric it does not distinguish the different degenerate quadrics."<<endl;
      cout<<"->eigenVal:\nIt returns the eigenvalues ​​relative to the conic/quadric."<<endl;
      cout<<"->eigenVect:\nIt returns the matrix P composed of the eigenvectors put in columns."<<endl;
      cout<<"->invar:\nIt returns invariants related to the matrices associated to the conic."<<endl;
      cout<<"->normalForm:\nIt returns a normal form (one possible) of the conical/quadric (if not degenerate)."<<endl;
      cout<<"->draw:\nDraws on a canvas, with different functions, the conic/quadric inserted."<<endl;
      cout<<"->drawNF:\nDraws as the previous command the conical/quadric inserted but in canonical form."<<endl;
      cout<<"->new:\nIt allows the user to enter a new equation to be studied by saving in memory the previous one."<<endl;
      cout<<"->list:\nReturns the sorted list of all equations entered from power program."<<endl;
      cout<<"->select:\nSelect between the equations stored in memory that the user chooses.\nIt is therefore possible to reuse the equation in all functions."<<endl;
      cout<<"->clear:\nDeletes all user-entered functions from the program memory."<<endl;
      cout<<"->.q:\nQuit ALGAgraph but leaves open the terminal window used."<<endl;
      find_command=true;
    }
    if(find_command==false){
      cout<<"Command not recognized."<<endl;
    }
    times++;
  }

  

  return 0;
}
