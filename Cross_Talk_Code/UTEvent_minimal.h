
#include"TROOT.h"


class UTEventMin 
{
public:
  Int_t EventID;
  Int_t nHits;
  Int_t *channelID;	//[nHits]
  Int_t *tdcTime;	//[nHits]
  Bool_t *fallingFlag;	//[nHits]
  UTEventMin(){EventID=0; nHits=0;}
  ~UTEventMin(){EventID=0; nHits=0;}
};

