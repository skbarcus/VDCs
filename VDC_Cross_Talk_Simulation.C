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
Int_t debug = 0;
Int_t nevts = 10000;
Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;
Double_t Y = 0.;
Double_t dl2p1 = 0.;   //Smallest distance from particle trajectory to wire1.
Double_t dl2p2 = 0.;   //Smallest distance from particle trajectory to wire2.
Double_t dl2c1 = 0.;   //Smallest distance from particle trajectory to circle1.
Double_t dl2c2 = 0.;   //Smallest distance from particle trajectory to circle2.
Double_t dc12p1 = 0.;  //Distance from circle1 to point1 or trajectory to point if trajectory passes through circle1.
Double_t dc22p2 = 0.;  //Distance from circle2 to point1 or trajectory to point if trajectory passes through circle2.
Double_t dl2c1pp = 0;  //Temporary variable to determine which intersect on circle1 is closest to trajecotry.
Double_t dl2c1pn = 0;
Double_t dl2c1np = 0;
Double_t dl2c1nn = 0;
Double_t dl2c2pp = 0;  //Temporary variable to determine which intersect on circle2 is closest to trajecotry.
Double_t dl2c2pn = 0;
Double_t dl2c2np = 0;
Double_t dl2c2nn = 0;
Double_t dc1 = 0.;     //Distance from particle trajectory to circle1 following the straight field lines.
Double_t dc2 = 0.;     //Distance from particle trajectory to circle2 following the straight field lines.
Double_t xp1 = -2.1215e-3;     //X coordinate of wire1's point.
Double_t xp2 = 2.1215e-3;     //X coordinate of wire2's point.
Double_t yp1 = 0.;     //Y coordinate of wire1's point.
Double_t yp2 = 0.;     //Y coordinate of wire2's point.
Double_t px[2] = {xp1,xp2};
Double_t py[2] = {yp1,yp2};
Double_t xc1p = 0.;    //X coordinate of perpendicular line to trajectory intersecting circle1 positive solution.
Double_t xc2p = 0.;    //X coordinate of perpendicular line to trajectory intersecting circle2 positive solution.
Double_t yc1pp = 0.;   //Y coordinate of perpendicular line to trajectory intersecting circle1 positive solution.
Double_t yc2pp = 0.;   //Y coordinate of perpendicular line to trajectory intersecting circle1 positive solution.
Double_t xc1n = 0.;    //X coordinate of perpendicular line to trajectory intersecting circle1 negative solution.
Double_t xc2n = 0.;    //X coordinate of perpendicular line to trajectory intersecting circle2 negative solution.
Double_t yc1pn = 0.;   //Y coordinate of perpendicular line to trajectory intersecting circle1 negative solution.
Double_t yc2pn = 0.;   //Y coordinate of perpendicular line to trajectory intersecting circle1 negative solution.
Double_t yc1np = 0.;
Double_t yc2np = 0.;
Double_t yc1nn = 0.;
Double_t yc2nn = 0.;
Double_t xc1 = 0.;     //X coordinate of closest point on circle 1 intersecting line perpendicular to trajectory.
Double_t yc1 = 0.;     //Y coordinate of closest point on circle 1 intersecting line perpendicular to trajectory.
Double_t xc2 = 0.;     //X coordinate of closest point on circle 2 intersecting line perpendicular to trajectory.
Double_t yc2 = 0.;     //Y coordinate of closest point on circle 2 intersecting line perpendicular to trajectory.
Double_t xl1 = 0.;     //X coordinate of perpendicular line (circle1) to trajectory intersecting the trajectory.
Double_t yl1 = 0.;     //Y coordinate of perpendicular line (circle1) to trajectory intersecting the trajectory.
Double_t xl2 = 0.;     //X coordinate of perpendicular line (circle2) to trajectory intersecting the trajectory.
Double_t yl2 = 0.;     //Y coordinate of perpendicular line (circle2) to trajectory intersecting the trajectory.
Double_t timeoutc1 = 0.;    //Time for particle to drift from trajectory to circle1 (=0 if trajectory is inside of circle1).
Double_t timeoutc2 = 0.;    //Time for particle to drift from trajectory to circle2 (=0 if trajectory is inside of circle2).
Double_t timeinc1 = 0.;     //Time for particle to drift to p1 inside circle1.
Double_t timeinc2 = 0.;     //Time for particle to drift to p2 inside circle2.
Double_t time1 = 0.;   //Total drift time to reach p1 (wire 1).
Double_t time2 = 0.;   //Total drift time to reach p2 (wire 2).
Double_t vout = 50.4e-6/(1.e-9);    //Drift velocity of particle in the linear field region (outside the circles).
Double_t vin = 2.*(50.4e-6/(1.e-9));     //Drift elocity in the quasiradial region (inside the circles). 
Double_t m = 0.;       //Slope of particle trajectory.
Double_t b = 0.;       //Y-intercept of particle trajectory.
Double_t posnegm = 0;  //Randomly assign slope to be positive or negative.
Double_t r = 8.e-3;     //Radius of circle (quasiradial field region around wires).
Double_t maxdth = 55.;  //Maximum detector theta angle of cosmic ray (degrees).
Double_t mindth = 45.;  //Minimum detector theta angle of cosmic ray (degrees).
Double_t maxb = 13.e-3;
Double_t minb = -13.e-3;

Double_t line(Double_t *X, Double_t *par)
{
  Y = m*X[0]+b;
  return Y;
}

Double_t circle1p(Double_t *X, Double_t *par)
{
  Y = pow(r*r-(X[0]-xp1)*(X[0]-xp1),0.5) + yp1;
  return Y;
}

Double_t circle1n(Double_t *X, Double_t *par)
{
  Y = -pow(r*r-(X[0]-xp1)*(X[0]-xp1),0.5) + yp1;
  return Y;
}

Double_t circle2p(Double_t *X, Double_t *par)
{
  Y = pow(r*r-(X[0]-xp2)*(X[0]-xp2),0.5) + yp2;
  return Y;
}

Double_t circle2n(Double_t *X, Double_t *par)
{
  Y = -pow(r*r-(X[0]-xp2)*(X[0]-xp2),0.5) + yp2;
  return Y;
}

void VDC_Cross_Talk_Simulation() 
{
  //Define new random number generator (could use TRandom2,3 instead).
  TRandom *r1=new TRandom();
  r1->SetSeed(0);            //New random seed each time.
  
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Difine histograms for wires' timing spectrum.
  TH1* h1 = new TH1D("h1", "Time Spectrum Wire 1", 100, 0.0, 0.2e-6);
  TH1* h2 = new TH1D("h2", "Time Spectrum Wire 2", 100, 0.0, 0.2e-6);
  TH2* h3 = new TH2D("h3", "Time Spectrum of Wire 2 vs. Wire 1", 1000, -0.01e-6, 0.4e-6, 1000, -0.01e-6, 0.4e-6);
  
  
  for(Int_t i=0;i<nevts;i++)
    {
      posnegm = r1->Uniform(-1.,1.);
      if(posnegm>=0)
	{
	  posnegm = 1;
	}
      else
	{
	  posnegm = -1;
	}
      m = r1->Uniform(1/tan(maxdth*deg2rad),1/tan(mindth*deg2rad))*posnegm;  //Select a random slope uniformly between two values.
      b = r1->Uniform(minb,maxb);  //Select a random y-intercept uniformly between two values. 
      //m = -0.75;
      //b = 3;
      
      //Calculate intersect of line perpendicular to trajectory passing through circle center with the particle trajectory.
      xl1 = (1/(m+1/m))*(yp1+xp1/m-b);
      yl1 = m*xl1+b;

      xl2 = (1/(m+1/m))*(yp2+xp2/m-b);
      yl2 = m*xl2+b;
      
      //Calculate shortest distances from particle trajectory (line) to each wire (points).
      dl2p1 = fabs(m*xp1-yp1+b)/pow(m*m+1,0.5);
      dl2p2 = fabs(m*xp2-yp2+b)/pow(m*m+1,0.5);
      
      //Calculate intersect of line perpendicular to trajectory passing through circle center with the circle (quasiradial region).
      //Wire 1:
      xc1p = (2.*xp1+pow(4.*xp1*xp1-4.*xp1*xp1+(4.*r*r)/(1.+1./(m*m)),0.5))/2.;
      xc1n = (2.*xp1-pow(4.*xp1*xp1-4.*xp1*xp1+(4.*r*r)/(1.+1./(m*m)),0.5))/2.;
      
      yc1pp = (2.*yp1+pow(4.*yp1*yp1-4.*(yp1*yp1+(xc1p-xp1)*(xc1p-xp1)-r*r),0.5))/2.;
      yc1pn = (2.*yp1-pow(4.*yp1*yp1-4.*(yp1*yp1+(xc1p-xp1)*(xc1p-xp1)-r*r),0.5))/2.;
      yc1np = (2.*yp1+pow(4.*yp1*yp1-4.*(yp1*yp1+(xc1n-xp1)*(xc1n-xp1)-r*r),0.5))/2.;
      yc1nn = (2.*yp1-pow(4.*yp1*yp1-4.*(yp1*yp1+(xc1n-xp1)*(xc1n-xp1)-r*r),0.5))/2.;
      
      dl2c1pp = pow((xc1p-xl1)*(xc1p-xl1) + (yc1pp-yl1)*(yc1pp-yl1),0.5);
      dl2c1pn = pow((xc1p-xl1)*(xc1p-xl1) + (yc1pn-yl1)*(yc1pn-yl1),0.5);
      dl2c1np = pow((xc1n-xl1)*(xc1n-xl1) + (yc1np-yl1)*(yc1np-yl1),0.5);
      dl2c1nn = pow((xc1n-xl1)*(xc1n-xl1) + (yc1nn-yl1)*(yc1nn-yl1),0.5);
      
      if(dl2c1pp < dl2c1pn && dl2c1pp < dl2c1np && dl2c1pp < dl2c1nn)
	{
	  xc1 = xc1p;
	  yc1 = yc1pp;
	}
      
      if(dl2c1pn < dl2c1pp && dl2c1pn < dl2c1np && dl2c1pn < dl2c1nn)
	{
	  xc1 = xc1p;
	  yc1 = yc1pn;
	}
      
      if(dl2c1np < dl2c1pn && dl2c1np < dl2c1pp && dl2c1np < dl2c1nn)
	{
	  xc1 = xc1n;
	  yc1 = yc1np;
	}
      
      if(dl2c1nn < dl2c1pn && dl2c1nn < dl2c1np && dl2c1nn < dl2c1pp)
	{
	  xc1 = xc1n;
	  yc1 = yc1nn;
	}

      //Wire 2:
      xc2p = (2.*xp2+pow(4.*xp2*xp2-4.*xp2*xp2+(4.*r*r)/(1.+1./(m*m)),0.5))/2.;
      xc2n = (2.*xp2-pow(4.*xp2*xp2-4.*xp2*xp2+(4.*r*r)/(1.+1./(m*m)),0.5))/2.;
      
      yc2pp = (2.*yp2+pow(4.*yp2*yp2-4.*(yp2*yp2+(xc2p-xp2)*(xc2p-xp2)-r*r),0.5))/2.;
      yc2pn = (2.*yp2-pow(4.*yp2*yp2-4.*(yp2*yp2+(xc2p-xp2)*(xc2p-xp2)-r*r),0.5))/2.;
      yc2np = (2.*yp2+pow(4.*yp2*yp2-4.*(yp2*yp2+(xc2n-xp2)*(xc2n-xp2)-r*r),0.5))/2.;
      yc2nn = (2.*yp2-pow(4.*yp2*yp2-4.*(yp2*yp2+(xc2n-xp2)*(xc2n-xp2)-r*r),0.5))/2.;
      
      dl2c2pp = pow((xc2p-xl2)*(xc2p-xl2) + (yc2pp-yl2)*(yc2pp-yl2),0.5);
      dl2c2pn = pow((xc2p-xl2)*(xc2p-xl2) + (yc2pn-yl2)*(yc2pn-yl2),0.5);
      dl2c2np = pow((xc2n-xl2)*(xc2n-xl2) + (yc2np-yl2)*(yc2np-yl2),0.5);
      dl2c2nn = pow((xc2n-xl2)*(xc2n-xl2) + (yc2nn-yl2)*(yc2nn-yl2),0.5);
      
      if(dl2c2pp < dl2c2pn && dl2c2pp < dl2c2np && dl2c2pp < dl2c2nn)
	{
	  xc2 = xc2p;
	  yc2 = yc2pp;
	}
      
      if(dl2c2pn < dl2c2pp && dl2c2pn < dl2c2np && dl2c2pn < dl2c2nn)
	{
	  xc2 = xc2p;
	  yc2 = yc2pn;
	}
      
      if(dl2c2np < dl2c2pn && dl2c2np < dl2c2pp && dl2c2np < dl2c2nn)
	{
	  xc2 = xc2n;
	  yc2 = yc2np;
	}
      
      if(dl2c2nn < dl2c2pn && dl2c2nn < dl2c2np && dl2c2nn < dl2c2pp)
	{
	  xc2 = xc2n;
	  yc2 = yc2nn;
	}
      
      //Calculate shortest distance from particle trajectory (line) to each quasiradial field area (circle) around the wires.
      dl2c1 = pow( (xc1-xl1)*(xc1-xl1)+(yc1-yl1)*(yc1-yl1)  ,0.5);

      dl2c2 = pow( (xc2-xl2)*(xc2-xl2)+(yc2-yl2)*(yc2-yl2)  ,0.5);
      
      //Calculate distance from particle trajectory to the circle following the straight field lines.
      dc1 = fabs(yc1-m*xc1-b);

      dc2 = fabs(yc2-m*xc2-b);
      
      //Calculate distance from circle to point or trajectory to point if trajectory passes through circle1.
      //Wire 1:
      if(r > pow((xp1-xl1)*(xp1-xl1)+(yp1-yl1)*(yp1-yl1),0.5))
	{
	  dc12p1 = pow((xl1-xp1)*(xl1-xp1)+(yl1-yp1)*(yl1-yp1),0.5); 
	}
      else
	{
	  dc12p1 = pow((xc1-xp1)*(xc1-xp1)+(yc1-yp1)*(yc1-yp1),0.5);
	}

      //Wire 2:
      if(r > pow((xp2-xl2)*(xp2-xl2)+(yp2-yl2)*(yp2-yl2),0.5))
	{
	  dc22p2 = pow((xl2-xp2)*(xl2-xp2)+(yl2-yp2)*(yl2-yp2),0.5); 
	}
      else
	{
	  dc22p2 = pow((xc2-xp2)*(xc2-xp2)+(yc2-yp2)*(yc2-yp2),0.5);
	}

      //Calculate the time required to reach the wires.
      //Wire 1:
      
      if(r < pow((xp1-xl1)*(xp1-xl1)+(yp1-yl1)*(yp1-yl1),0.5))
	{
	  //timeoutc1 = dl2c1/vout;//Not following field lines (straight path to wire).
	  timeoutc1 = dc1/vout;
	}
      else
	{
	  timeoutc1 = 0;
	}

      if(r < pow((xp1-xl1)*(xp1-xl1)+(yp1-yl1)*(yp1-yl1),0.5))
	{
	  timeinc1 = dc12p1/vin;
	}
      else
	{
	  timeinc1 = dl2p1/vin;
	}

      time1 = timeoutc1 + timeinc1;

      //Wire 2:

      if(r < pow((xp2-xl2)*(xp2-xl2)+(yp2-yl2)*(yp2-yl2),0.5))
	{
	  //timeoutc2 = dl2c2/vout;//Not following field lines (straight path to wire).
	  timeoutc2 = dc2/vout;
	}
      else
	{
	  timeoutc2 = 0;
	}

      if(r < pow((xp2-xl2)*(xp2-xl2)+(yp2-yl2)*(yp2-yl2),0.5))
	{
	  timeinc2 = dc22p2/vin;
	}
      else
	{
	  timeinc2 = dl2p2/vin;
	}

      time2 = timeoutc2 + timeinc2;

      //Print individual variables if debugging.      
      if(debug==1)
	{
	  cout<<"m = "<<m<<"   b = "<<b<<endl;
	  cout<<"dl2p1 = "<<dl2p1<<"   dl2p2 = "<<dl2p2<<endl;
	  cout<<"xc1p = "<<xc1p<<"   xc1n = "<<xc1n<<"   yc1pp = "<<yc1pp<<"   yc1pn = "<<yc1pn<<"   yc1np = "<<yc1np<<"   yc1nn = "<<yc1nn<<endl;
	  cout<<"xl1 = "<<xl1<<"   yl2 = "<<yl1<<endl;
	  cout<<"dl2c1 = "<<dl2c1<<endl;
	  cout<<"xc1 = "<<xc1<<"   yc1 = "<<yc1<<endl;
	  cout<<"dl2c1pp = "<<dl2c1pp<<"   dl2c1pn = "<<dl2c1pn<<"   dl2c1np = "<<dl2c1np<<"   dl2c1nn = "<<dl2c1nn<<endl;
	  cout<<"dc1 = "<<dc1<<endl;
	  cout<<"dc12p1 = "<<dc12p1<<endl;
	  cout<<"timeoutc1 = "<<timeoutc1<<"   timeinc1 = "<<timeinc1<<"   time1 = "<<time1<<endl;
	}

      TF1 *ftrajectory = new TF1("ftrajectory", line, -3.e-2, 3.e-2, 1);
      //Fill histograms with wire time data.
      h1->Fill(time1);
      h2->Fill(time2);
      h3->Fill(time1,time2);
      
    }
  
  cout<<"m = "<<m<<"   b = "<<b<<endl;
  cout<<"dl2p1 = "<<dl2p1<<"   dl2p2 = "<<dl2p2<<endl;
  cout<<"xc1p = "<<xc1p<<"   xc1n = "<<xc1n<<"   yc1pp = "<<yc1pp<<"   yc1pn = "<<yc1pn<<"   yc1np = "<<yc1np<<"   yc1nn = "<<yc1nn<<endl;
  cout<<"xl1 = "<<xl1<<"   yl2 = "<<yl1<<endl;
  cout<<"dl2c1 = "<<dl2c1<<endl;
  cout<<"xc1 = "<<xc1<<"   yc1 = "<<yc1<<endl;
  cout<<"dl2c1pp = "<<dl2c1pp<<"   dl2c1pn = "<<dl2c1pn<<"   dl2c1np = "<<dl2c1np<<"   dl2c1nn = "<<dl2c1nn<<endl;
  cout<<"dc1 = "<<dc1<<endl;
  cout<<"dc12p1 = "<<dc12p1<<endl;
  cout<<"timeoutc1 = "<<timeoutc1<<"   timeinc1 = "<<timeinc1<<"   time1 = "<<time1<<endl;
	  
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  
  //c1->GetXaxis()->SetRange(-10,10);
  //TList *l = c1->GetListOfPrimitives();
  //cout<<l<<endl;

  //c1->GetListOfPrimitives()->ls();
  //c1->GetListOfPrimitives()->Print();

  ftrajectory->SetNpx(1000);
  ftrajectory->Draw();
  ftrajectory->SetMaximum(13.e-3);
  ftrajectory->SetMinimum(-13.e-3);
  Int_t q = 5;
  //Draw lines perpendicular to the particle trajectory passing through each circle's center.
  TF1 *fperp1 = new TF1("fperp1", Form("-x/%f+%f+%f/%f",m,yp1,xp1,m), -1.e-2, 1.e-2);//Form("-x/%f+%f+%f/%f",m,yp1,xp1,m)
  fperp1->Draw("same");
  fperp1->SetLineColor(6);
  TF1 *fperp2 = new TF1("fperp2", Form("-x/%f+%f+%f/%f",m,yp2,xp2,m), -1.e-2, 1.e-2);//Form("-x/%f+%f+%f/%f",m,yp1,xp1,m)
  fperp2->Draw("same");
  fperp2->SetLineColor(6);

  TGraph *gr = new TGraph(2,px,py);
  gr->Draw("p");
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(1);

  TF1 *fc1p = new TF1("fc1p", circle1p, -r+xp1, r+xp1, 1);
  TF1 *fc1n = new TF1("fc1n", circle1n, -r+xp1, r+xp1, 1);
  TF1 *fc2p = new TF1("fc2p", circle2p, -r+xp2, r+xp2, 1);
  TF1 *fc2n = new TF1("fc2n", circle2n, -r+xp2, r+xp2, 1);
  
  fc1p->Draw("same");
  fc1p->SetNpx(10000);
  fc1p->SetLineColor(1);
  fc1n->Draw("same");
  fc1n->SetNpx(10000);
  fc1n->SetLineColor(1);
  fc2p->Draw("same");
  fc2p->SetNpx(10000);
  fc2p->SetLineColor(1);
  fc2n->Draw("same");
  fc2n->SetNpx(10000);
  fc2n->SetLineColor(1);

  //Lines from circle intersect (or trajectory if inside circle radius) to circle center.
  //Wire 1 :
  if(r > pow((xp1-xl1)*(xp1-xl1)+(yp1-yl1)*(yp1-yl1),0.5))
    {
      TLine *line1 = new TLine(xp1, yp1, xl1, yl1);
      line1->SetLineColor(kBlue);
      line1->SetLineWidth(2);
      line1->Draw("SAME");
    }
  else
    {
      TLine *line1 = new TLine(xp1, yp1, xc1, yc1);
      line1->SetLineColor(kBlue);
      line1->SetLineWidth(2);
      line1->Draw("SAME");
    }

  //Wire 2:
  if(r > pow((xp2-xl2)*(xp2-xl2)+(yp2-yl2)*(yp2-yl2),0.5))
    {
      TLine *line2 = new TLine(xp2, yp2, xl2, yl2);
      line2->SetLineColor(kBlue);
      line2->SetLineWidth(2);
      line2->Draw("SAME");
    }
  else
    {
      TLine *line2 = new TLine(xp2, yp2, xc2, yc2);
      line2->SetLineColor(kBlue);
      line2->SetLineWidth(2);
      line2->Draw("SAME");
    }

  //Line representing distance from trajectory to circle following straight field lines.
  //Wire 1:
  if(r < pow((xp1-xl1)*(xp1-xl1)+(yp1-yl1)*(yp1-yl1),0.5))
    {
      TLine *line3 = new TLine(xc1, yc1, xc1, m*xc1+b);
      line3->SetLineColor(3);
      line3->SetLineWidth(2);
      line3->Draw("SAME");
    }

  //Wire 2:
  if(r < pow((xp2-xl2)*(xp2-xl2)+(yp2-yl2)*(yp2-yl2),0.5))
    {
      TLine *line4 = new TLine(xc2, yc2, xc2, m*xc2+b);
      line4->SetLineColor(3);
      line4->SetLineWidth(2);
      line4->Draw("SAME");
    }

  //Lines representing the cell around the two wires.
  //TLine *line5 = new TLine(-10, -10, -10, 10);
  //line5->SetLineColor(kBlue);
  //line5->SetLineWidth(2);
  //line5->Draw("SAME");;

  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  h1->Draw("");
  h1->GetXaxis()->SetTitle("Time (s)");
  h1->GetYaxis()->SetTitle("Counts");

  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  h2->Draw("");
  h2->GetXaxis()->SetTitle("Time (s)");
  h2->GetYaxis()->SetTitle("Counts");

  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();
  h3->Draw("colz");
  h3->GetXaxis()->SetTitle("Wire 1 Time (s)");
  h3->GetYaxis()->SetTitle("Wire 2 Time (s)");

  st->Stop();
  cout<<"   CPU time = "<<st->CpuTime()<<"   Real time = "<<st->RealTime()<<endl;
}
