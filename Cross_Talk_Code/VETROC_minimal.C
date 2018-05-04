// Analysis of FADC data
// S. Barcus, R. Michaels, May 2015

// These event limits used if maxevent=-1 (1st arg)
#define FIRSTEVENT 40000
#define LASTEVENT  120000

#define FIRSTWORD  0xb0b0b0b4
#define MAXROC     50
//#define MAXRAW     6000
#define MAXRAW     12000


#define UNKNOWN -1
#define TRUE  1
#define FALSE 0

#define VETROC_DATA_BLOCK_HEADER      0x00000000
#define VETROC_DATA_BLOCK_TRAILER     0x08000000
#define VETROC_DATA_EVENT_HEADER      0x10000000
#define VETROC_DATA_TRIGGER_TIME      0x18000000
#define VETROC_DATA_INVALID           0x70000000
#define VETROC_DATA_FILLER            0x78000000
#define VETROC_DUMMY_DATA             0xf800fafa
#define VETROC_DATA_TDCEVT            0x40000000

#include <iostream>
#include <string>
#include <vector>
#include "THaCodaFile.h"
#include "THaEtClient.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include <vector>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "./UTEvent_minimal.h"

using namespace std;

void usage();
void decode(int* data);
void clear();
//void vetrocDataDecode(unsigned int data);

// Global data 

Int_t *frawdata;
Int_t *irn, *rocpos, *roclen;
Int_t myroc,mychan;
Int_t evlen, evtype, evnum;
Int_t ntdcdata;
Int_t nposdata;
Int_t *tdcdata;
Int_t *rowdata;
Int_t *coldata;
Int_t *chnlhits;
Int_t chnl;
const Int_t N_PMT = 128;      //Number of PMTs in array.
char h1title[100];
char c1title[100];
Int_t tdcChData[N_PMT][MAXRAW];  //Make a 2d array with elements equal to the number of PMTs and then store the individual channel data in those elements. Element 1 is the number of PMTs and element 2 is the data for that PMTs channel. 
Int_t ntdcChDataTotal;         //Total number of hits in all channels for one event. 
Int_t *ntdcChData;             //Total number of hits in each individual PMT channel.
Int_t *ntdcChDataEvtDisp;             //Total number of hits in each individual PMT channel for the event display.

//Int_t Data[2000][N_PMT][100];   //Create 3D array to store the event in which the data is taking place, the channel that saw hits, and the TDC times of those hits.
Int_t evt_number = -1;
Int_t nhits = 0;
Int_t hits_per_evt[N_PMT] = {0};

Int_t *numdataevents;
Int_t ndataevents;

//Int_t *trgtime0;       //filled with trig time 0
Int_t *firsttdcevt;    //filled with first tdc value after trig time 0
Int_t ntrgtime0;       //counts number of trig time 0 events
Int_t nfirsttdcevt;    //counts number of first tdc events seen
Int_t checkfirst;      //check to be sure that the only tdc event recorded is the first tdc event after trg time 0
//Int_t event;

Int_t *evthdr;           //Stores number of event headers seen
Int_t nevthdr;           //Counts number of event headers seen
Int_t *afterevthdr;      //Stores number of words seen after an event header until the next header
Int_t nafterevthdr;      //Counts number of words seen after an event header until the next header 
Int_t *sumnafterevthdr;  //Array which stores the number of times a certain number of words followed an evt hdr
//Int_t nsumnafterevthdr;
Int_t nevthdrreset;

Int_t ntdcevt;             //counts number of words seen after and event header that are tdc events
Int_t *tdcchdif;           //stores time difference between consecutive tdc events in the same channel
Int_t ntdcchdif;           //counts the number of consecutive hit times compared
Int_t firstevt = 0;        //checks if first evt (0) or not (1)
Int_t nwordsfirstevt = 0;  //counts number of words in the first event
Int_t oldref = 1;          //0 = Don't use old reference channel method from TEDF tests with pulser. 1 = Do use old method.

Int_t userefch = 1;            //0-> don't use a reference channel. 1-> use ch 0 as the reference to all others (removes start time jitter).
Int_t refchoffset;             //stores the value of the first TDC time in the reference channel
//Int_t firstreftime = 0;        //stores if this is the first TDC time in ref ch or not for each event
Int_t refch = 112;               //set which channel should be used as the reference channel
Int_t detector_geometry = 0;   //0 = prototype detector PMT geometry. 1 = GRINCH detector PMT geometry.

Int_t mydebug = 1;             //0 = Turns off some of the histos used for analysis. 1 = Turns on the histos for debugging.

unsigned long long trgtime0;
unsigned long long trgtime1;
unsigned long long trgtime;

Int_t total_events = 0;
Int_t total_hits = 0;
Int_t tooManyHitsCounter = 0;

UTEventMin *event;
TTree *eventtree;

int main(int argc, char* argv[])
{
  Int_t load_status = gSystem->Load("./UTEvent_minimal_h.so");
  cout << "Status: " << load_status << endl;


  THaCodaData *coda;      
  char cname[100],ctitle[100];

  Int_t maxevent = 1000;

  Int_t debug = 1;             
  Int_t choice1 = 1;  /* source of data */
  myroc=9;
  mychan=4;

//  cout << "AAA"<<endl;
  
  if (argc < 2) {
     usage();
     return 1;
  }
 
  irn = new Int_t[MAXROC];
  rocpos = new Int_t[MAXROC];
  roclen = new Int_t[MAXROC];
  tdcdata = new Int_t[MAXRAW];
  rowdata = new Int_t[MAXRAW];
  coldata = new Int_t[MAXRAW];
  chnlhits = new Int_t[MAXRAW];
  numdataevents = new Int_t[MAXRAW];
  firsttdcevt = new Int_t[MAXRAW];
//  cout << "BBB"<<endl;
  evthdr = new Int_t[10000000];
  afterevthdr = new Int_t[10000000];
  tdcchdif = new Int_t[MAXRAW];
  sumnafterevthdr = new Int_t[MAXRAW];
  ntdcChData = new Int_t[N_PMT];
  ntdcChDataEvtDisp = new Int_t[N_PMT];
  ntdcdata = 0;
  nposdata = 0;
  ndataevents = 0;
  ntrgtime0 = 0;
  nfirsttdcevt = 0;
  checkfirst = 0;
  nevthdr = 0;
//  cout << "CCC"<<endl;
  nafterevthdr = 0;
  nevthdrreset = 0;
  ntdcevt = 0;
  ntdcchdif = 0;
  ntdcChDataTotal = 0;
  ntdcChData[N_PMT] = {0};
  ntdcChDataEvtDisp[N_PMT] = {0};

  if (argc > 1) 
    {
      maxevent = atoi(argv[1]);
    }
  if (argc > 2) 
    {
      mychan = atoi(argv[2]);
    }
 
//  cout << "VETROC analysis "<<endl;
//  cout << "Events to process "<<maxevent<<endl;

  frawdata = new Int_t[MAXRAW]; 



  // Initialize root and output
  TROOT evanana("evanroot","Hall A VETROC analysis, 1st version");
  TFile hfile("minimal.root","RECREATE","VETROC data");

  //Build TTree to store TDC data
  eventtree = new TTree("eventtree","TDC data");
  eventtree->Branch("event",&event);

  Int_t mychan=4;
  sprintf(ctitle,"Raw samples on channel %d",mychan);

  if (choice1 == 1) 
    {  // CODA File
      
      // CODA file "run.dat" may be a link to CODA datafile on disk
      TString filename("run_minimal.dat");
      
      coda = new THaCodaFile();
      if (coda->codaOpen(filename) != 0) 
	{
	  cout << "ERROR:  Cannot open CODA data" << endl;
	  goto end1;
	}
      
    
    } 
  else 
    {         // Online ET connection
      
      int mymode = 1;
      TString mycomputer("adaq1");
      TString mysession("Compton");
      
      coda = new THaEtClient();
      if (coda->codaOpen(mycomputer, mysession, mymode) != 0) 
	{
	  cout << "ERROR:  Cannot open ET connection" << endl;
	  goto end1;
	}
      
    }

  // Loop over events
  int status,evlo,evhi;
  
  clear();
  
  evlo = 0;
  evhi = maxevent;

  for (int iev = 0; iev < evhi; iev++) 
    {//the loop over the event bounds
      
      if (iev > 0 && ((iev%1000)==0) ) printf("%d events\n",iev);

      clear();

      status = coda->codaRead();  

      if (status != 0) 
	{  // EOF or no data
	  
	  if ( status == -1) 
	    {
	      if (choice1 == 1) 
		{
		  cout << "End of CODA file. Bye bye." << endl;
		  evhi=iev;
		}
	      if (choice1 == 2) cout << "CODA/ET not running. Bye bye." << endl;
	    } 
	  else 
	    {
	      cout << "ERROR: codaRread status = " << hex << status << endl;
	    }
	  goto end1;
  

	} 
      else 
	{   // have data ...
	  decode( coda->getEvBuffer() );

//	  cout<<"Finished decode call, iev = " << iev << endl;
	  
	}
    }
  
 end1:
  
  coda->codaClose();

  
  eventtree->Write();     //Write 7 REM
  hfile.Write();
  hfile.Close();

  return 0;

}; //end of main function

void usage() 
{  
  cout << "Usage:  'evanana [maxevents] [chan]' " << endl;
}; 

void clear() 
{
  ntdcdata = 0;
  nposdata = 0;
  ndataevents = 0;
  ntrgtime0 = 0;
  nfirsttdcevt = 0;
  checkfirst = 0;
  nafterevthdr = 0;
  nevthdrreset = 0;
  ntdcevt = 0;
  ntdcchdif = 0;
  ntdcChDataTotal = 0;
  total_hits = 0;
  nhits = 0;
  memset (frawdata, 0, MAXRAW);
  memset (tdcdata, 0, MAXRAW);
  memset (rowdata, 0, MAXRAW);
  memset (coldata, 0, MAXRAW);
  memset (chnlhits, 0, MAXRAW);
  memset (numdataevents, 0, MAXRAW);
  memset (firsttdcevt, 0, MAXRAW);
  memset (sumnafterevthdr, 0, MAXRAW);
  memset (tdcChData,0,MAXRAW);
  memset (ntdcChData,0,MAXRAW);
  memset (ntdcChDataEvtDisp,0,MAXRAW);
}

void decode (int* data) 
{

  //REM stuff
  const Int_t maxHitsPerEvent = 2000;
  Int_t channelID[maxHitsPerEvent] = {0};
  Int_t tdcTime[maxHitsPerEvent] = {0};
  Bool_t fallingFlag[maxHitsPerEvent] = {false};
  //REM stuff
 
  Int_t found_dat;
  evlen = data[0] + 1;
  evtype = data[1]>>16;
  evnum = data[4];
  static int dodump = 1;  // dump the raw data
  static int debug = 1;   // debug the decoding
  if (dodump) 
    { // Event dump
//      cout << "\n\n Event number " << dec << evnum;
//      cout << " length " << evlen << " type " << evtype << endl;
      int ipt = 0;
      for (int j = 0; j < (evlen/5); j++) 
	{
//	  cout << dec << "\n evbuffer[" << ipt << "] = ";
	  for (int k=j; k<j+5; k++) 
	    {
//	      cout << hex << data[ipt++] << " ";
	    }
//	  cout << endl;
	}
      if (ipt < evlen) 
	{
//	  cout << dec << "\n evbuffer[" << ipt << "] = ";
	  for (int k=ipt; k<evlen; k++) 
	    {
//	      cout << hex << data[ipt++] << " ";
	    }
//	  cout << endl;
	}
    }

  // Decoding 
  if (evtype != -999) 
    {

      // First find pointers to ROCs.  Useful for multi-crate analysis.
      
      // n1 = pointer to first word of ROC
      int pos = data[2]+3;
      int nroc = 0;
      while( pos+1 < evlen && nroc < MAXROC ) 
	{
	  int len  = data[pos];
	  int iroc = (data[pos+1]&0xff0000)>>16;
	  if(iroc>=MAXROC || nroc >= MAXROC) 
	    {
	      cout << "Decoder error:  ";
	      cout << "  ROC num " <<dec<<iroc;
	      cout << "  nroc "<<nroc<<endl;
	      return;
	    }
	  // Save position and length of each found ROC data block
	  rocpos[iroc]  = pos;
	  roclen[iroc]  = len;
	  irn[nroc++] = iroc;
	  pos += len+1;
	}
      Int_t found = 0;
      for (int j=0; j<nroc; j++) 
	{
	  if(myroc == irn[j]) found=1;
	}
      if (!found) 
	{
	  cout << "ROC "<<dec<<myroc<<" not in datastream !!"<<endl;
	  return;
	}

      if (debug)
	{
//	  cout << "Roc info "<<nroc<<endl;
	  for (int j=0; j<nroc; j++) 
	    {
	      Int_t iroc = irn[j];
//	      cout << "irn "<<dec<<iroc<<"   pos "<<rocpos[iroc];
//	      cout <<"   len "<<roclen[iroc]<<endl;
	    }
	}    

      // Go over the data in myroc

      found_dat = 0;
 
      for (int j = rocpos[myroc]; j < rocpos[myroc]+roclen[myroc]; j++) 
	{
	  
	  if (found_dat) break;
	  
//	  if (debug) printf("data[%d] = 0x%x = (dec) %d \n",j,data[j],data[j]);
//	  if (debug) printf("Looking for firstword = 0x%x \n",FIRSTWORD);
	  
	  Int_t icnt = 0;
	  Int_t numwords = 0;
//	  if ((data[j] != 0xffffffff) && ((data[j] & FIRSTWORD) ==  FIRSTWORD))  -- VERY SUSCEPTIBLE TO FALSE POSITIVES!!!! -- REM -- 2016-12-12
          if (data[j] == FIRSTWORD)
	    {
//	     cout << "FOUND FIRSTWORD, j = " << j << endl;
	      found_dat = 1;
              numwords = data[j+1];
//	      printf("datablock: 0x%x numwords = 0x%x = %d",data[j],data[j+1],data[j+1]);
	      if (numwords > MAXRAW) 
		{
		  printf("error:  numwords %d exceeds MAXRAW = %d.  Truncating.\n",numwords,MAXRAW);
		  numwords = MAXRAW-1;
		}
//	      if (debug) printf("found event header = 0x%x  numwords = 0x%x = %d  roclen %d \n",data[j],numwords,numwords,roclen[myroc]);
	     
	      for (int k=j+2; k<j+numwords; k++) 
		{
//		  printf("data[%d] = 0x%x = %d\n",k,data[k],data[k]);
		  if(userefch == 1)
		    {
		      //Count TDC Events
		      if((data[k] & 0x78000000) == VETROC_DATA_EVENT_HEADER) 
			{
			  evt_number++;
			  event->EventID = evt_number;
			  event->nHits = 0;
//			  cout<<"$$$ Total Events = "<<evt_number<<"  event->EventID = "<<event->EventID<<endl;
			  for(Int_t i = 0; (data[k+3+i] & 0x78000000) == VETROC_DATA_TDCEVT; i++)
			    {
                            event->nHits++;
                            }
//			cout << "EventID = " << evt_number << "      nHits = " << event->nHits << endl;
			}

			//REM here is the event loop 
		      if((data[k] & 0x78000000) == VETROC_DATA_EVENT_HEADER)//Check to see if data is a TDC event header.
			{
			  for(Int_t i = 0; (data[k+3+i] & 0x78000000) == VETROC_DATA_TDCEVT; i++)//Loop over all TDC events after this event header.
			    {
                              if(i>=maxHitsPerEvent)
                                {
                                tooManyHitsCounter++;
                                cout << tooManyHitsCounter << ". TOO MANY HITS! MISSING EVENT HEADER?  >" << i << " hits." << endl;
                                cerr << tooManyHitsCounter << ". TOO MANY HITS! MISSING EVENT HEADER?  >" << i << " hits." << endl;
                                break;
                                }
                              //REM
                              channelID[i] = ((data[k+3+i]>>16)&0x7f);
                              tdcTime[i] = ((data[k+3+i])&0xffff);
                              //cout << "i = " << i << "    channel = " << channelID[i] << "     tdcTime = " << tdcTime[i] << endl;
                              Bool_t tempBool = false;
                              if( ((data[k+3+i]>>26)&0x1) == 0x1 )	//0 = rising, 1 = falling
                              {
                                tempBool = true;
                              }
                              fallingFlag[i] = tempBool;
                              //REM
			    }
                          event->channelID = channelID;
                          event->tdcTime = tdcTime;
                          event->fallingFlag = fallingFlag;

			  eventtree->Fill();   //Fill 3
			  //Reset TDCVal array to all zeros after each event.
			  //Reset the data storage arrays.
			  for(int m=0; m<128; m++)
			    {
			      ntdcChDataEvtDisp[m] = 0;
			    }
			}
		    }
		  
		}//End loop over data[k].
	    }
	}
    }
}


