#include <string>
#include <vector>
using namespace std;

void CJ2ROOT()
{
  ifstream infile("two_events.txt");
  
  TFile *f = TFile::Open("hvector.root","RECREATE");

  int event, nHits;
  std::vector<int> pad;
  std::vector<int> tdc;
  std::vector<double> adc;

   // Create a TTree
   TTree *t = new TTree("tvec","Tree with vectors");
   t->Branch("event",&event);
   t->Branch("nHits",&nHits);
   t->Branch("pad",&pad);
   t->Branch("tdc",&tdc);
   t->Branch("adc",&adc);

   int fevent, fpad, ftdc, nHits;
   double fadc;

   int event_counter = 1;
   int fHits = 0;
   
   while (!infile.eof())
     {
       if (infile >> fevent >> fpad >> ftdc >> fadc)
	{
	  if (fevent != event_counter) //if there is a change on event, store and clean the vectors from the previous event
	    {
	      event = fevent;
	      nHits = fHits;
	      t->Fill();
	      
	      pad.clear();
	      tdc.clear();
	      adc.clear();
	      fHits = 0;
	      
	      event_counter++;
	    }
 
	  pad.push_back(fpad);
	  tdc.push_back(ftdc);
	  adc.push_back(fadc);
	  fHits++;
	  if(fHits%1000==0) cout<<fHits<<endl; //counter
	}
     }

   //We need to repeat these lines in order to store the last event
   event = fevent;
   nHits = fHits;
   t->Fill();

   
   f->Write();
   delete f;
   
     

}
