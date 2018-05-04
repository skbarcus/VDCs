#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <math.h>
//#define theta1 21.;

Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;

Double_t Eprime(Double_t *Ep, Double_t *par)
{
  //Ep = E / ( 1 + (2*E/M)*pow(sin(angle[0]/2.),2.) );
  //return Ep;
  
  //theta = TMath::ASin(pow( (M/(2*E))*(E/Ep[0] - 1) ,0.5))*2 *1/deg2rad;
  theta = TMath::ACos(M/E+1-M/Ep[0]) * 1/deg2rad;
  return theta;
}

// define the parameteric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

void VDC_Cross_Talk_Simulation() 
{
  Int_t n = 1000;
  double t0 = 0;
  double dt = 10;
  double p0[4] = {10,20,1,2};
  double p1[4] = {-10,20,1,2};

   // draw original line
   TPolyLine3D *l0 = new TPolyLine3D(n);
   for (int i = 0; i <n;++i) {
      double t = t0+ dt*i/n;
      double x,y,z;
      line(t,p0,x,y,z);
      l0->SetPoint(i,x,y,z);
   }
   l0->SetLineColor(kBlue);
   l0->Draw("same");

   // draw original line
   TPolyLine3D *l1 = new TPolyLine3D(n);
   for (int i = 0; i <n;++i) {
      double t = t0+ dt*i/n;
      double x,y,z;
      line(t,p1,x,y,z);
      l1->SetPoint(i,x,y,z);
   }
   l1->SetLineColor(kBlue);
   l1->Draw("same");
}
