// Purpose:  Define basic histograms for jet physics studies
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 9, 2015
#ifndef __basicHistos_h__
#define __basicHistos_h__

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectory.h"

#include <string>
#include <map>
#include <vector>

class basicHistos {

 public:

  // phase space
  std::string trigname;
  std::string cotrig;
  double ymin;
  double ymax;
  double pttrg;
  double ptmin;
  double ptmax;
  bool ismc;
  
  /*
  // raw spectrum
  TH1D *hpt;
  TH1D *hpt_pre;
  TH1D *hcopt;
  TH1D *hpt_noid;
  TH1D *hpt_nojetid;
  TH1D *hpt_noevtid;
  TH1D *hptevt;
  TH1D *hpttmp;

  // delete-m jackknife
  std::vector<TH1D*> hpt_jk;
  TH2D *h2jk;
  
  //INCLUSIVE JETS
  TH1D *hpt_tmp;
  TH1D *hpt_evtcount;
  TH1D *hpt_evt;
  TH1D *hpt_jet;
  
  // 1 GeV ins
  TH1D *hpt0;
  
  // leading and non-leading jets
  TH1D *hpt1;
  TH1D *hpt2;
  TH1D *hpt3;
  
  // High chf jets pt spectrum
  TH1D *hpt4;
  //
  */
  
  // dijet mass
  TH1D *hdjmass;
  TH1D *hdjmass_half;
  TH1D *hdjmass0;
  
  TH1D *djmass_matched; // Matched reco jets with gen jets under dR criteria
  TH1D *hdj_leading;
  TH1D *hdj_subleading;
  
  TProfile *pdjmass_ptratio;
  TProfile *pdjmass0_ptratio;
  
  // GEN-LEVEL dijet mass
  
  TH1D *gen_hdjmass;
  TH1D *gen_hdjmass_half;
  TH1D *gen_hdjmass0;
  
  TH1D *gen_djmass_matched; // Matched gen jets with reco under dR criteria
  TH1D *gen_hdj_leading;
  TH1D *gen_hdj_subleading;
  
  /*
  // basic properties
  TProfile *ppt;
  TProfile *pmass;
  TH1D *hmass;
  TProfile *pjec;
  TProfile *pjec2;
  TProfile *punc;
  TH1D *hnpvgood;
  TH1D *hrho;
  */
  
  /*
  TProfile *pjec_l1;
  TProfile *pjec_l2l3;
  TProfile *pjec_res;
  */
  /*
  TH1D *htrpu;
  TH1D *hitpu;
  TH1D *hootpuearly;
  TH1D *hootpulate;
  TH2D *h2itvsoot;
  */
  
  /*
  // pile-up information
  TProfile *pa;
  TProfile *ptrpu;
  TProfile *pnpv;
  TProfile *pnpvall;
  TProfile *prho;
  TProfile *pnpvvsrho;
  TProfile *prhovsnpv;
  TProfile *prhovsnpvall;
  TH2D *h2rhovsnpv;
  //
  TProfile *prhovstrpu;
  TProfile *pnpvvstrpu;
  TProfile *pnpvallvstrpu;
  TProfile *pitpuvstrpu;
  TH1D *htrpu2;
  TH1D *hjet_vstrpu;
  */
  TH1D *hlumi_vstrpu;
  
  // luminosity
  TH1D *hlumi;
  TH1D *hlumi2;
  std::map<int, std::map<int, float> > lums;
  double lumsum;
  double lumsum2;
  
  // inclusive efficiencies
  TProfile *peff;
  TProfile *pideff;
  TProfile *pvtxeff;
  TProfile *pdqmeff;
  
  /*
  // control plots of components (JEC)
  TProfile *pncand;
  TProfile *pnch;
  TProfile *pnne;
  TProfile *pnnh;
  TProfile *pnce;
  TProfile *pnmu;
  TProfile *pnhh;
  TProfile *pnhe;
  TProfile *pchf;
  TProfile *pchf_vs_cm;
  
  // CHF studies  
  TProfile *pchf_pt;
  TProfile *pchf_eta;
  TProfile *pchf_phi;
    
  TProfile *pnhf_pt;
  TProfile *pnhf_eta;
  TProfile *pnhf_phi;
    
  TProfile *pncand_pt;
  TProfile *pncand_eta;
  TProfile *pncand_phi;
  // 
	
  TProfile *pnef;
  TProfile *pnhf;
  TProfile *pcef;
  TProfile *pmuf;
  TProfile *phhf;
  TProfile *phef;
  TProfile *pbeta;
  TProfile *pbetastar;
  TH1D *hncand;
  TH1D *hnch;
  TH1D *hnne;
  TH1D *hnnh;
  TH1D *hnce;
  TH1D *hnmu;
  TH1D *hnhh;
  TH1D *hnhe;
  TH1D *hchf;

  TH1D *hchf_leading;
  TH1D *hchf_subleading;
  TH1D *hchf_thirdjet;
  TH1D *hnef;
  TH1D *hnhf;
  TH1D *hcef;
  TH1D *hmuf;
  TH1D *hhhf;
  TH1D *hhef;
  TH1D *hbeta;
  TH1D *hbetastar;
  
  // control plots of components (JEC tag-and-probe)
  TProfile *pncandtp;
  TProfile *pnchtp;
  TProfile *pnnetp;
  TProfile *pnnhtp;
  TProfile *pncetp;
  TProfile *pnmutp;
  TProfile *pnhhtp;
  TProfile *pnhetp;
  TProfile *pchftp;
  TProfile *pneftp;
  TProfile *pnhftp;
  TProfile *pceftp;
  TProfile *pmuftp;
  TProfile *phhftp;
  TProfile *pheftp;
  TProfile *pbetatp;
  TProfile *pbetastartp;
  //
  TProfile *pchftp2;
  TProfile *pneftp2;
  TProfile *pnhftp2;
  TProfile *pceftp2;
  TProfile *pmuftp2;
  TProfile *phhftp2;
  TProfile *pheftp2;
  //
  TH1D *hncandtp;
  TH1D *hnchtp;
  TH1D *hnnetp;
  TH1D *hnnhtp;
  TH1D *hncetp;
  TH1D *hnmutp;
  TH1D *hnhhtp;
  TH1D *hnhetp;
  TH1D *hchftp;
  TH1D *hneftp;
  TH1D *hnhftp;
  TH1D *hceftp;
  TH1D *hmuftp;
  TH1D *hhhftp;
  TH1D *hheftp;
  TH1D *hbetatp;
  TH1D *hbetastartp;
  //
  TProfile *pncandtp_vsnpv;
  TProfile *pnchtp_vsnpv;
  TProfile *pnnetp_vsnpv;
  TProfile *pnnhtp_vsnpv;
  TProfile *pncetp_vsnpv;
  TProfile *pnmutp_vsnpv;
  TProfile *pnhhtp_vsnpv;
  TProfile *pnhetp_vsnpv;
  TProfile *pchftp_vsnpv;
  TProfile *pneftp_vsnpv;
  TProfile *pnhftp_vsnpv;
  TProfile *pceftp_vsnpv;
  TProfile *pmuftp_vsnpv;
  TProfile *phhftp_vsnpv;
  TProfile *pheftp_vsnpv;
  TProfile *pbetatp_vsnpv;
  TProfile *pbetastartp_vsnpv;
  //
  TProfile *pchftp_vstrpu;
  TProfile *pneftp_vstrpu;
  TProfile *pnhftp_vstrpu;
  TProfile *pceftp_vstrpu;
  TProfile *pmuftp_vstrpu;
  //TProfile *phhftp_vstrpu;
  //TProfile *pheftp_vstrpu;
  TProfile *pbetatp_vstrpu;
  TProfile *pbetastartp_vstrpu;
  
  // control plots for topology (JEC)
  TH1D *hselpt;
  TH1D *hy;
  TH1D *hy2;
  TH1D *heta;
  TH1D *heta2;
  TH1D *hphi;
  TH1D *hdphi;
  TH1D *hdpt;
  TProfile *pdpt;
  TH1D *hjet;
  TH1D *hmet;
  TH1D *hmetphi;
  // control plots for vertex
  TH1D *hpvndof;
  TH1D *hpvx;
  TH1D *hpvy;
  TH1D *hpvz;
  TH1D *hpvr;
  TH1D *hpvrho;
  // closure plots for JEC
  TH1D *hmpf;
  TH1D *hmpf1;
  TH1D *hmpf2;
  TH1D *hmpfx;
  TH1D *hmpfy;
  TH1D *hmpfz;
  TProfile *pmpf;
  TProfile *pmpf1;
  TProfile *pmpf2;
  TProfile *pmpfx;
  TProfile *pmpfy;
  TProfile *pmpfz;
  */
  
  /*
  // Control plots of resolutions:
  // dijet asymmetry binned in dijet pT and 3rd jet pT
  TH3D *hdjasymm;
  TH3D *hdjmpf;
  TH3D *hdjasymmtp;
  TH3D *hdjmpftp;

  TH1D *hr21;
  TH1D *hr31;
  TH1D *hr32;
  TProfile *pr21;
  TProfile *pr31;
  TProfile *pr32;
  TProfile *px21;
  TProfile *px31;
  TProfile *px32;

  TH1D *hyeta;
  TH1D *hyeta2;
  TH2D *hbetabetastar;
  TH2D *hetaphi;
  */
  
  /*
  //////// CHF second peak investigation //////
  TH2D *h_low_chf_etaphi;
  TH2D *h_high_chf_etaphi;
  */
  // MC checks
  /*
  TH1D *hpt_jt30;
  TH1D *hpt_jt60;
  //TH1D *hpt_jt80;
  TH1D *hpt_jt110;
  //TH1D *hpt_jt150;
  TH1D *hpt_jt190;
  TH1D *hpt_jt240;
  TH1D *hpt_jt300;
  TH1D *hpt_jt370;
  //
  TH1D *hpt0_jt30;
  TH1D *hpt0_jt60;
  //TH1D *hpt0_jt80;
  TH1D *hpt0_jt110;
  //TH1D *hpt0_jt150;
  TH1D *hpt0_jt190;
  TH1D *hpt0_jt240;
  TH1D *hpt0_jt300;
  TH1D *hpt0_jt370;
  */
  /*
  // unfolding studies (Mikael)
  TH2D *mT;
  TH2D *mTf;
  TH2D *mTuw;
  TH2D *mTfuw;
  TH1D *mx;
  TH1D *mxf;
  TH1D *mxuw;
  TH1D *mxfuw;
  TH1D *my;
  TH1D *myuw;
  TH1D *myf;
  TH1D *myfuw;
  */
  
  // unfolding studies dijet mass
  TH2D *matrix_gen_reco;
  // 2D delta mass vs mass plots for resolution studies
  TH2D *h2jetres;
  // Profile plot to monitor mean value of mass resolution
  TProfile *pdjmass_res;
  
  // Acceptance and background studies for inclusive jets
  TProfile *paccept_vsPt;
  TProfile *pbg_vsPt;
  
   
  TH1D *hpthat;
  TH1D *hpthatnlo;
  /*
  TH1D *hpt_noid_g;
  TH1D *hpt_nojetid_g;
  TH1D *hpt_noevtid_g;
  TH1D *hpt_r;
  TH1D *hpt_g;
  TH1D *hpt_gtw;
  TH1D *hpt_gg;
  TH1D *hpt_gg0;
  TH1D *hpt_g0;
  TH1D *hpt_g0tw;
  TH1D *hpt_g0_tmp;
  //

  TProfile *ppt_r;
  TProfile *ppt_g;
  */
  
  /*
  // Response closure
  TProfile3D *p3rvsnpv;
  TProfile3D *p3rvsnpvW;
  TProfile2D *p2rvsnpv;
  TH2D *h2r_r;
  TH2D *h2r_g;
  TProfile *p2r_r;
  TProfile *p2r_g;
  TProfile *p2r_ruw;
  TProfile *p2r_guw;

  //Rapidity closure
  TH2D *h2dy_r;
  TH2D *h2dy_g;
  TProfile *p2dy_r;
  TProfile *p2dy_g;
  TProfile *p2dy_ruw;
  TProfile *p2dy_guw;
  TProfile2D *pdy_r;
  TProfile2D *pdy_g;
  */
  
  basicHistos(TDirectory *dir, std::string trigname="", std::string cotrig="",
	      double ymin = 0., double ymax = 2.0,
	      double pttrg = 10., double ptmin = 10., double ptmax = 50.,
	      bool ismc = false);
  ~basicHistos();

// private:
  TDirectory *dir;
};

#endif // __basicHistos_h__
