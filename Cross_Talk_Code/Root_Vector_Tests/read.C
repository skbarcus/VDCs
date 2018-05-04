void read()
{
   TFile *f = TFile::Open("hvector.root","READ");
   if (!f) { return; }
   
   TTree *t; f->GetObject("tvec",t);

   int nHits;
   std::vector<float> *adc = 0;
   tvec->SetBranchAddress("nHits",nHits);
   
   TBranch *bvpx = 0;
   t->SetBranchAddress("adc",&adc);
   
   // for (Int_t i = 0; i < 25000; i++) {
   //    Long64_t tentry = t->LoadTree(i);



   tvec->GetEntry(0);
   
      for (UInt_t j = 0; j < adc->size(); ++j) {
	cout<<adc->at(j)<<endl;
      }
      
}
