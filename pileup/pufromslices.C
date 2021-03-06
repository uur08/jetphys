// A handy script for fetching the number of entries in a pthat sliced MC sample
// The user should replace fileformat and dirname according to present needs
// Aouthor Hannu Siikonen hannu.siikonen@cern.ch

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <map>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>

#include "../settings.h"

void pufromslices() {
  TString dirname="/eos/cms/store/group/phys_smp/Multijet/13TeV/MC/P825ns80X_Moriond17";
  std::regex fileformat("QCD_Pt_([0-9]*)to([0-9]*|Inf)_TuneCUETP8M_13TeV_pythia8.root");
  std::cmatch match;
  TFile *output = new TFile("pu.root","RECREATE");
  TH1D *summary = new TH1D("pileupmc","",600,0,60);
  vector<int> pthatmin =
    {30,50,80,120,170,300,470,600,800,1000,1400,1800,2400,3200};


  TSystemDirectory dir(dirname.Data(), dirname.Data());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        TFile *f = new TFile((dirname+"/"+fname).Data());
        TTree *t = (TTree*) f->Get("ak4/ProcessedTree");
        std::regex_match(fname.Data(), match, fileformat);
        int number = stoi(match[1]);
        const char *histname = Form("pileupmc%d",number);
        TH1D *hist = new TH1D(histname,"",600,0,60);
        t->Draw(Form("EvtHdr_.mTrPu>>%s",histname));
        int pos = 0;
        while (pthatmin[pos]!=number)
          ++pos;
        summary->Add(hist,_jp_pthatsigmas[pos]/_jp_pthatnevts[pos]);
      }
    }
  }
  output->cd();
  summary->Write();
  //output->Write();
}
