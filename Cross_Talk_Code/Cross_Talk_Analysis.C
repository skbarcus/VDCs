#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>

void Cross_Talk_Analysis() {

  gROOT->Reset();
  
  //Open files.
  FILE *f1 = fopen("/home/skbarcus/Tritium/Analysis/VDCs/Cross_Talk_Code/Cross_Talk_361_u1.txt","r");
  FILE *f2 = fopen("/home/skbarcus/Tritium/Analysis/VDCs/Cross_Talk_Code/Cross_Talk_361_u2.txt","r");
  FILE *f3 = fopen("/home/skbarcus/Tritium/Analysis/VDCs/Cross_Talk_Code/Cross_Talk_361_v1.txt","r");
  FILE *f4 = fopen("/home/skbarcus/Tritium/Analysis/VDCs/Cross_Talk_Code/Cross_Talk_361_v2.txt","r");
  
  Int_t debug = 0;
  Float_t v,w,x,y,z;          //Temporary variables to hold data from data file.
  Int_t ncols;          //Set how many columns of data we have in the data file.
  Int_t nlines_u1 = 0;     //Counts number of lines in the data file. 
  Int_t nlines_u2 = 0;
  Int_t nlines_v1 = 0;
  Int_t nlines_v2 = 0;
  Int_t skip = 1;       //Gives number of lines to skip at top of data file. 
  char* str[300];       //Variable to read lines of the data file.
  Double_t wire1_u1[400];    //Arrays to hold X and Y data.
  Double_t wire2_u1[400];
  Double_t tot_counts_u1[400];
  Double_t CT_counts_u1[400];
  Double_t CT_percent_u1[400];
  Double_t CT_error_u1[400];
  Double_t wire1_u2[400];    //Arrays to hold X and Y data.
  Double_t wire2_u2[400];
  Double_t tot_counts_u2[400];
  Double_t CT_counts_u2[400];
  Double_t CT_percent_u2[400];
  Double_t CT_error_u2[400];
  Double_t wire1_v1[400];    //Arrays to hold X and Y data.
  Double_t wire2_v1[400];
  Double_t tot_counts_v1[400];
  Double_t CT_counts_v1[400];
  Double_t CT_percent_v1[400];
  Double_t CT_error_v1[400];
  Double_t wire1_v2[400];    //Arrays to hold X and Y data.
  Double_t wire2_v2[400];
  Double_t tot_counts_v2[400];
  Double_t CT_counts_v2[400];
  Double_t CT_percent_v2[400];
  Double_t CT_error_v2[400];
  Double_t avg_CT_u1 = 0.;
  Double_t avg_CT_u2 = 0.;
  Double_t avg_CT_v1 = 0.;
  Double_t avg_CT_v2 = 0.;
  Double_t avg_CT_CT_u1 = 0.;
  Double_t avg_CT_CT_u2 = 0.;
  Double_t avg_CT_CT_v1 = 0.;
  Double_t avg_CT_CT_v2 = 0.;
  Double_t CT_wires_u1 = 0.;
  Double_t CT_wires_u2 = 0.;
  Double_t CT_wires_v1 = 0.;
  Double_t CT_wires_v2 = 0.;

  //Read in blank disc data.
  while (1) {
    //Skips the first 'skip' lines of the file. 
    if (nlines_u1 < skip)
      {
	fgets(str,100,f1);
	nlines_u1++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(f1,"%f %f %f %f %f", &v, &w, &x, &y, &z);
	if (ncols < 0) break;    
	wire1_u1[nlines_u1-skip] = v;
	wire2_u1[nlines_u1-skip] = w;
	tot_counts_u1[nlines_u1-skip] = x;
	CT_counts_u1[nlines_u1-skip] = y;
	if(x==0)
	  {
	    CT_percent_u1[nlines_u1-skip] = 0.;
	    CT_error_u1[nlines_u1-skip] = 0.;
	  }
	else
	  {
	    CT_percent_u1[nlines_u1-skip] = z;
	    CT_error_u1[nlines_u1-skip] = (pow(y,0.5)/x)*100.;
	  }

	//Count the number of lines in the file. 
	nlines_u1++;
      }
  }

  while (1) {
    //Skips the first 'skip' lines of the file. 
    if (nlines_u2 < skip)
      {
	fgets(str,100,f2);
	nlines_u2++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(f2,"%f %f %f %f %f", &v, &w, &x, &y, &z);
	if (ncols < 0) break;    
	wire1_u2[nlines_u2-skip] = v;
	wire2_u2[nlines_u2-skip] = w;
	tot_counts_u2[nlines_u2-skip] = x;
	CT_counts_u2[nlines_u2-skip] = y;
	if(x==0)
	  {
	    CT_percent_u2[nlines_u2-skip] = 0.;
	    CT_error_u2[nlines_u2-skip] = 0.;
	  }
	else
	  {
	    CT_percent_u2[nlines_u2-skip] = z;
	    CT_error_u2[nlines_u2-skip] = (pow(y,0.5)/x)*100.;
	  }

	//Count the number of lines in the file. 
	nlines_u2++;
      }
  }

  while (1) {
    //Skips the first 'skip' lines of the file. 
    if (nlines_v1 < skip)
      {
	fgets(str,100,f3);
	nlines_v1++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(f3,"%f %f %f %f %f", &v, &w, &x, &y, &z);
	if (ncols < 0) break;    
	wire1_v1[nlines_v1-skip] = v;
	wire2_v1[nlines_v1-skip] = w;
	tot_counts_v1[nlines_v1-skip] = x;
	CT_counts_v1[nlines_v1-skip] = y;
	if(x==0)
	  {
	    CT_percent_v1[nlines_v1-skip] = 0.;
	    CT_error_v1[nlines_v1-skip] = 0.;
	  }
	else
	  {
	    CT_percent_v1[nlines_v1-skip] = z;
	    CT_error_v1[nlines_v1-skip] = (pow(y,0.5)/x)*100.;
	  }

	//Count the number of lines in the file. 
	nlines_v1++;
      }
  }

  while (1) {
    //Skips the first 'skip' lines of the file. 
    if (nlines_v2 < skip)
      {
	fgets(str,100,f4);
	nlines_v2++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(f4,"%f %f %f %f %f", &v, &w, &x, &y, &z);
	if (ncols < 0) break;    
	wire1_v2[nlines_v2-skip] = v;
	wire2_v2[nlines_v2-skip] = w;
	tot_counts_v2[nlines_v2-skip] = x;
	CT_counts_v2[nlines_v2-skip] = y;
	if(x==0)
	  {
	    CT_percent_v2[nlines_v2-skip] = 0.;
	    CT_error_v2[nlines_v2-skip] = 0.;
	  }
	else
	  {
	    CT_percent_v2[nlines_v2-skip] = z;
	    CT_error_v2[nlines_v2-skip] = (pow(y,0.5)/x)*100.;
	  }

	//Count the number of lines in the file. 
	nlines_v2++;
      }
  }
  
  //Print the data read from the file. 
  if(debug==1)
    {
      for(int i=0; i<400; i++)
	{
	  cout<<"wire1_u1["<<i<<"] = "<<wire1_u1[i]<<"   wire2_u1["<<i<<"] = "<<wire2_u1[i]<<"   tot_counts_u1["<<i<<"] = "<<tot_counts_u1[i]<<"   CT_counts_u1["<<i<<"] = "<<CT_counts_u1[i]<<"   CT_percent_u1["<<i<<"] = "<<CT_percent_u1[i]<<"%   CT_error_u1["<<i<<"] = "<<CT_error_u1[i]<<" %"<<endl;
	}

for(int i=0; i<400; i++)
	{
	  cout<<"wire1_u2["<<i<<"] = "<<wire1_u2[i]<<"   wire2_u2["<<i<<"] = "<<wire2_u2[i]<<"   tot_counts_u2["<<i<<"] = "<<tot_counts_u2[i]<<"   CT_counts_u2["<<i<<"] = "<<CT_counts_u2[i]<<"   CT_percent_u2["<<i<<"] = "<<CT_percent_u2[i]<<"%   CT_error_u2["<<i<<"] = "<<CT_error_u2[i]<<" %"<<endl;
	}

for(int i=0; i<400; i++)
	{
	  cout<<"wire1_v1["<<i<<"] = "<<wire1_v1[i]<<"   wire2_v1["<<i<<"] = "<<wire2_v1[i]<<"   tot_counts_v1["<<i<<"] = "<<tot_counts_v1[i]<<"   CT_counts_v1["<<i<<"] = "<<CT_counts_v1[i]<<"   CT_percent_v1["<<i<<"] = "<<CT_percent_v1[i]<<"%   CT_error_v1["<<i<<"] = "<<CT_error_v1[i]<<" %"<<endl;
	}

for(int i=0; i<400; i++)
	{
	  cout<<"wire1_v2["<<i<<"] = "<<wire1_v2[i]<<"   wire2_v2["<<i<<"] = "<<wire2_v2[i]<<"   tot_counts_v2["<<i<<"] = "<<tot_counts_v2[i]<<"   CT_counts_v2["<<i<<"] = "<<CT_counts_v2[i]<<"   CT_percent_v2["<<i<<"] = "<<CT_percent_v2[i]<<"%   CT_error_v2["<<i<<"] = "<<CT_error_v2[i]<<" %"<<endl;
	}
    }

  //Print number of lines with data.
  printf(" found %d points for U1\n",nlines_u1 - skip);
  printf(" found %d points for U2\n",nlines_u2 - skip);
  printf(" found %d points for V1\n",nlines_v1 - skip);
  printf(" found %d points for V2\n",nlines_v2 - skip);

  //Create some histograms.
  TH1F *hu1 = new TH1F("hu1","U1 Number of Cross Talk Hits;Number of Cross Talk Hits;Number of Wires",50,0,20);
  TH1F *hu2 = new TH1F("hu2","U2 Number of Cross Talk Hits;Number of Cross Talk Hits;Number of Wires",50,0,20);
  TH1F *hv1 = new TH1F("hv1","V1 Number of Cross Talk Hits;Number of Cross Talk Hits;Number of Wires",50,0,20);
  TH1F *hv2 = new TH1F("hv2","V2 Number of Cross Talk Hits;Number of Cross Talk Hits;Number of Wires",50,0,20);

  //Calculate average cross talk across planes two ways. 1st is a normal average and second is average of only the wires which had cross talk.
  for(Int_t i=0;i<368;i++)
    {
      avg_CT_u1 = avg_CT_u1 + CT_percent_u1[i];
      avg_CT_u2 = avg_CT_u2 + CT_percent_u2[i];
      avg_CT_v1 = avg_CT_v1 + CT_percent_v1[i];
      avg_CT_v2 = avg_CT_v2 + CT_percent_v2[i];
      
      if(CT_counts_u1[i]!=0)
	{
	  avg_CT_CT_u1 = avg_CT_CT_u1 + CT_percent_u1[i];
	  CT_wires_u1++;
	}
      if(CT_counts_u2[i]!=0)
	{
	  avg_CT_CT_u2 = avg_CT_CT_u2 + CT_percent_u2[i];
	  CT_wires_u2++;
	}
      if(CT_counts_v1[i]!=0)
	{
	  avg_CT_CT_v1 = avg_CT_CT_v1 + CT_percent_v1[i];
	  CT_wires_v1++;
	}
      if(CT_counts_v2[i]!=0)
	{
	  avg_CT_CT_v2 = avg_CT_CT_v2 + CT_percent_v2[i];
	  CT_wires_v2++;
	}
      //cout<<"CT_counts_u1 = "<<CT_counts_u1[i]<<"   avg_CT_u1 = "<<avg_CT_u1<<"   avg_CT_CT_u1 = "<<avg_CT_CT_u1<<endl;
      //cout<<"CT_wires_u1 = "<<CT_wires_u1<<endl;

      hu1->Fill(CT_counts_u1[i]);
      hu2->Fill(CT_counts_u2[i]);
      hv1->Fill(CT_counts_v1[i]);
      hv2->Fill(CT_counts_v2[i]);
    }

  avg_CT_u1 =  avg_CT_u1/368.;
  avg_CT_u2 =  avg_CT_u2/368.;
  avg_CT_v1 =  avg_CT_v1/368.;
  avg_CT_v2 =  avg_CT_v2/368.;
  avg_CT_CT_u1 =  avg_CT_CT_u1/CT_wires_u1;
  avg_CT_CT_u2 =  avg_CT_CT_u2/CT_wires_u2;
  avg_CT_CT_v1 =  avg_CT_CT_v1/CT_wires_v1;
  avg_CT_CT_v2 =  avg_CT_CT_v2/CT_wires_v2;

  cout<<"Average cross talk per wire pair U1 = "<<avg_CT_u1<<endl;
  cout<<"Average cross talk per wire pair U2 = "<<avg_CT_u2<<endl;
  cout<<"Average cross talk per wire pair V1 = "<<avg_CT_v1<<endl;
  cout<<"Average cross talk per wire pair V2 = "<<avg_CT_v2<<endl;
  cout<<"Average cross talk per wire pair U1 with CT = "<<avg_CT_CT_u1<<endl;
  cout<<"Average cross talk per wire pair U2 with CT = "<<avg_CT_CT_u2<<endl;
  cout<<"Average cross talk per wire pair V1 with CT = "<<avg_CT_CT_v1<<endl;
  cout<<"Average cross talk per wire pair V2 with CT = "<<avg_CT_CT_v2<<endl;
  
  //Close data files. 
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  
  //Make a new canvas to plot data.
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph1 = new TGraphErrors(nlines_u1-skip,wire1_u1,CT_percent_u1,0,CT_error_u1);
  //Draw the new TGraph called graph on the canvas. 
  graph1->Draw("AP");
  
  //Set X axis
  graph1->GetXaxis()->SetLimits(0,400);
  //Set Y axis Min and Max (not sure why different from X).
  graph1->SetMinimum(0);
  graph1->SetMaximum(5);
  graph1->SetLineWidth(1);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->SetMarkerStyle(3);
  graph1->SetTitle("Run 90150 Plane U1 Percentage of Hits Flagged as Cross Talk per Wire Pairs; Number of First Wire in Wire Pair; Percentage of Hits Flagged as Cross Talk (%)");
  //graph_expected.SetFillColor(kYellow);
  //graph_expected.DrawClone("E3AL"); // E3 draws the band

  // Draw the Legend
  //TLegend leg1(0.9,.7,.56,.9,"Legend Title");
  //leg1.SetFillColor(0);
  //leg1.AddEntry(graph1,"Curve Name");
  //leg1.DrawClone("Same");                     //Lets you draw multiple curves to the same canvas.

  //Make a new canvas to plot data.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph2 = new TGraphErrors(nlines_u2-skip,wire1_u2,CT_percent_u2,0,CT_error_u2);
  //Draw the new TGraph called graph on the canvas. 
  graph2->Draw("AP");
  
  //Set X axis
  graph2->GetXaxis()->SetLimits(0,400);
  //Set Y axis Min and Max (not sure why different from X).
  graph2->SetMinimum(0);
  graph2->SetMaximum(5);
  graph2->SetLineWidth(1);
  graph2->SetLineColor(1);
  graph2->SetFillColor(0);
  graph2->SetMarkerStyle(3);
  graph2->SetTitle("Run 90150 Plane U2 Percentage of Hits Flagged as Cross Talk per Wire Pairs; Number of First Wire in Wire Pair; Percentage of Hits Flagged as Cross Talk (%)");


  //Make a new canvas to plot data.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph3 = new TGraphErrors(nlines_v1-skip,wire1_v1,CT_percent_v1,0,CT_error_v1);
  //Draw the new TGraph called graph on the canvas. 
  graph3->Draw("AP");
  
  //Set X axis
  graph3->GetXaxis()->SetLimits(0,400);
  //Set Y axis Min and Max (not sure why different from X).
  graph3->SetMinimum(0);
  graph3->SetMaximum(5);
  graph3->SetLineWidth(1);
  graph3->SetLineColor(1);
  graph3->SetFillColor(0);
  graph3->SetMarkerStyle(3);
  graph3->SetTitle("Run 90150 Plane V1 Percentage of Hits Flagged as Cross Talk per Wire Pairs; Number of First Wire in Wire Pair; Percentage of Hits Flagged as Cross Talk (%)");


  //Make a new canvas to plot data.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();

  //Make a new TGraph to plot on the canvas. Parameters are (number of points, x values array, y values array).
  graph4 = new TGraphErrors(nlines_v2-skip,wire1_v2,CT_percent_v2,0,CT_error_v2);
  //Draw the new TGraph called graph on the canvas. 
  graph4->Draw("AP");
  
  //Set X axis
  graph4->GetXaxis()->SetLimits(0,400);
  //Set Y axis Min and Max (not sure why different from X).
  graph4->SetMinimum(0);
  graph4->SetMaximum(5);
  graph4->SetLineWidth(1);
  graph4->SetLineColor(1);
  graph4->SetFillColor(0);
  graph4->SetMarkerStyle(3);
  graph4->SetTitle("Run 90150 Plane V2 Percentage of Hits Flagged as Cross Talk per Wire Pairs; Number of First Wire in Wire Pair; Percentage of Hits Flagged as Cross Talk (%)");

  //Make a new canvas to plot data.
  TCanvas* c5=new TCanvas("c5");
  c5->SetGrid();
  hu1->SetFillColor(4);
  hu1->Draw("b");

  //Make a new canvas to plot data.
  TCanvas* c6=new TCanvas("c6");
  c6->SetGrid();
  hu2->SetFillColor(4);
  hu2->Draw("b");

  //Make a new canvas to plot data.
  TCanvas* c7=new TCanvas("c7");
  c7->SetGrid();
  hv1->SetFillColor(4);
  hv1->Draw("b");

  //Make a new canvas to plot data.
  TCanvas* c8=new TCanvas("c8");
  c8->SetGrid();
  hv2->SetFillColor(4);
  hv2->Draw("b");
}

