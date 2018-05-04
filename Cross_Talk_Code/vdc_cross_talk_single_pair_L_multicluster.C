#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
//#include <THaRun.h>

#include <TStopwatch.h>
#include <TMath.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

void vdc_cross_talk_single_pair_L_multicluster(int run){

  //gROOT->Reset();

  const TString rootfilePath = "./rootfiles/";
  std::ostringstream str;
  str << rootfilePath<<"triton_L_"<<run;
  TString basename = str.str().c_str();
  TString rootfile = basename + "-1.6.root";

  TChain* T = new TChain("T");

  T->SetBranchStatus("*",0);               //Turns off all branches (faster). 
  T->SetBranchStatus("L.vdc.u1.dist",1);   //Turns on single branch.
  T->SetBranchStatus("L.vdc.u1.time",1);
  T->SetBranchStatus("L.vdc.u1.wire",1);
  T->SetBranchStatus("L.vdc.u1.nclust",1);
  T->SetBranchStatus("L.vdc.u1.clsiz",1);
  T->SetBranchStatus("L.vdc.u1.clpivot",1);
  T->SetBranchStatus("L.vdc.u1.clpos",1);
  T->SetBranchStatus("L.vdc.u1.nhit",1);
  T->SetBranchStatus("L.tr.d_th",1);
  T->SetBranchStatus("L.tr.d_ph",1);

  Long_t split=0;
  char* file = 0;
  
  //====adding splits rootfiles =======================//
  
  Long_t q=0;
  while ( !gSystem->AccessPathName(rootfile.Data()) ) 
    {
      T->Add(rootfile.Data());
      cout << "ROOT file " << rootfile << " added to TChain." << endl;
      q++;
      //rootfile = basename + Form("_%d",q) + ".root";
      rootfile = basename + "-1.6_" + q + ".root";
    }
  
  //==finish adding splits rootfiles=====================//
    gStyle->SetOptStat(0);
    
    Int_t wire1, wire2;
    
    cout << "\nreplay: Please enter wire 1 number (lower number): ";
    cin >> wire1;
    
    cout << "\nreplay: Please enter wire 2 number (higher number): ";
    cin >> wire2;
    //======================================================//
    
    //Define a new stopwatch.
    TStopwatch *st=new TStopwatch();
    st->Start(kTRUE);

    Int_t nevts = T->GetEntries();
    cout<<"Number of Events = "<<nevts<<endl;

    Int_t in_cl = 0.;            //Flag to determine if target wire was part of a cluster. (0=not in cluster. 1=in cluster)
    Int_t ncl = 0.;             //How many clusters were found in a given event.
    Double_t ct_y = 0.;            //Event was likely cross talk.
    Double_t ct_n = 0.;            //Event was likely real.
    Double_t total_counts = 0.;
    Int_t debug = 0;
    Double_t d_th_max = 55.; //Maximum detector theta.
    Double_t d_th_min = 40.; //Minimum detector theta.

    Double_t wire[400] = {0.};
    Double_t time[400] = {0.};
    Double_t dist[400] = {0.};
    Double_t nclust = 0.;
    Double_t clsiz[25] = {0.};
    Double_t clpivot[25] = {0.};
    Double_t clpos[25] = {0.};
    Double_t nhit = 0.;
    Double_t d_th[400] = {0.};
    Double_t d_ph[400] = {0.};

    char plane[] = {'x','x','\0'};
    plane[0] = 'u';
    plane[1] = '1';

  //Open file to write output to.
  std::ofstream output (Form("Cross_Talk_L_%i_%s_wire%i_wire%i.txt",run,plane,wire1,wire2), std::ofstream::out);
  output<<"Wires      Total Counts    CT Counts    CT (%)"<<endl;

    T->SetBranchAddress("L.vdc.u1.wire",wire);
    T->SetBranchAddress("L.vdc.u1.time",time);
    T->SetBranchAddress("L.vdc.u1.dist",dist);
    T->SetBranchAddress("L.vdc.u1.nclust",&nclust);
    T->SetBranchAddress("L.vdc.u1.clsiz",clsiz);
    T->SetBranchAddress("L.vdc.u1.clpivot",clpivot);
    T->SetBranchAddress("L.vdc.u1.clpos",clpos);
    T->SetBranchAddress("L.vdc.u1.nhit",&nhit);
    T->SetBranchAddress("L.tr.d_th",d_th);  
    T->SetBranchAddress("L.tr.d_ph",d_ph);

    TH2D *h1 = new TH2D("h1",Form("Crosstalk Run %i",run) , 100, -0.1e-6, 0.4e-6, 100, -0.1e-6, 0.4e-6);

    //Positive slope lower cross talk bound.
    TF1 *line1 = new TF1("line1","1.0*x-0.03e-6", -0.1e-6, 0.4e-6);
    //Positive slope upper cross talk bound.
    TF1 *line2 = new TF1("line2","1.0*x+0.03e-6", -0.15e-6, 0.4e-6);
    //Negative slope cross talk bound.
    TF1 *line3 = new TF1("line3","-1.0*x+0.15e-6", -0.e-6, 0.14e-6);

    Double_t xwire, ywire, wire_total=0., wire_total_temp=0.;
    Double_t fired_total = 0.;
    Double_t both_fired = 0.;
 
    Int_t event = 2728;
    /*
    cout<<"Event # = "<<event<<endl;
    cout<<"**********************************************************************"<<endl;
    T->GetEntry(event);
    for(Int_t j=0;j<100;j++)
      {
	cout<<wire[j]<<"   ";
      }
    cout<<endl;
    cout<<"**********************************************************************"<<endl;
    */

    //90150 247336, 90151 284684
    for(Int_t j = 0; j<nevts;j++)
      {	
	T->GetEntry(j);
	
	if(j%1000==0)
	  {
	    cout<<j<<endl;
	  }

	for(Int_t k=0;k<25;k++)
	  {
	    
	    if(clsiz[k]==8 && (clpivot[k]-4)<=wire1 && wire1<=(clpivot[k]+4))
	      {
		in_cl = 1;
	      }
	    if((clsiz[k]==7 || clsiz[k]==6) && (clpivot[k]-3)<=wire1 && wire1<=(clpivot[k]+3))
	      {
		in_cl = 1;
	      }
	    if((clsiz[k]==5 || clsiz[k]==4) && (clpivot[k]-2)<=wire1 && wire1<=(clpivot[k]+2))
	      {
		in_cl = 1;
	      }
	    ncl = nclust;
	    if(in_cl==1)
	      {
		for(Int_t m=0;m<400;m++)
		  {
		    if(wire[m]==wire1 && wire[m+1]==wire2 && TMath::ATan(d_th[k])*180/3.141592653 <d_th_max && TMath::ATan(d_th[k])*180/3.141592653 > d_th_min)
		      {
			
			if(debug==1)
			  {
			    cout<<"Event # = "<<j<<endl;
			    cout<<"Numbers of clusters = "<<nclust<<endl;
			    cout<<"Cluster sizes:"<<endl;
			    for(Int_t n = 0;n<20;n++)
			      {
				cout<<clsiz[n]<<"   ";
			      }
			    cout<<endl;
			    
			    cout<<"Cluster pivot wires:"<<endl;
			    for(Int_t n = 0;n<20;n++)
			      {
				cout<<clpivot[n]<<"   ";
			      }
			    cout<<endl;
			    
			    cout<<"Number of active wires: "<<endl;
			    for(Int_t n = 0;n<50;n++)
			      {
				cout<<wire[n]<<"   ";
			      }
			    cout<<endl;
			    cout<<"Timing Difference between wires (time wire1 -time wire2) = "<<(time[m]-time[m+1])/1e-9<<" ns"<<endl;
			    cout<<"**********************************************************************"<<endl;
			  }

			h1->Fill(time[m],time[m+1]);
			
			total_counts++;

			//Check if the wires firing is in the cross talk region.
			if(time[m+1]>line1->Eval(time[m]) && time[m+1]<line2->Eval(time[m]) && time[m+1]>line3->Eval(time[m]))
			  {
			    ct_y++;
			  }
			else
			  {
			    ct_n++;			   
			  }
		      }
		  }
	      }
	     in_cl = 0.;
	  }

	//Reset wire array to zero or data from event j-1 left in array.
	for(Int_t k=0;k<400;k++)
	  {
	    wire[k] = 0.;
	    time[k] = 0.;
	    dist[k] = 0.;
	    d_th[k] = 0.;
	    d_ph[k] = 0.;
	  }
	for(Int_t k=0;k<25;k++)
	  {
	    clsiz[k] = 0.;
	    clpivot[k] = 0.;
	    clpos[k] = 0.;
	  }

	    nhit = 0.;
	    nclust = 0.;
      }
    
    cout<<"Total counts = "<<total_counts<<endl;
    cout<<"Total number of times the target wires fired in a good cluster is "<<ct_y+ct_n<<". Cross talk occured "<<ct_y<<" times, and did not occur "<<ct_n<<" times. The percentage of events suspected to be due to cross talk is "<<(ct_y/(ct_y+ct_n))*100.<<"%."<<endl;

    output<<wire1<<" "<<wire2<<"    "<<ct_y+ct_n<<"             "<<ct_y<<"            "<<(ct_y/(ct_y+ct_n))*100.<<endl;
    output.close();

    TCanvas* c1=new TCanvas("c1");
    c1->SetGrid();

    h1->GetXaxis()->SetTitle(Form("Time of Wire %i (s)",wire1));
    h1->GetYaxis()->SetTitle(Form("Time of wire %i (s)",wire2));
    h1->Draw("colz");

    line1->Draw("SAME");
    line2->Draw("SAME");
    line3->Draw("SAME");

    c1->SaveAs(Form("./%i/%s/%i_%s_wire%i_wire%i_multicluster.png",run,plane,run,plane,wire1,wire2));
    c1->SaveAs(Form("./%i/%s/%i_%s_wire%i_wire%i_multicluster.C",run,plane,run,plane,wire1,wire2));

    st->Stop();
    cout<<"   CPU time = "<<st->CpuTime()<<"   Real time = "<<st->RealTime()<<endl;
}  

int main()
{
  vdc_cross_talk_single_pair_L_multicluster(361);
}
