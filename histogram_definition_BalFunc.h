#include "call_libraries.h"
#include "input_variables.h"

// jet nchtrk bin
const int nch_bin = 11;
//double nch_binedge[nch_bin+1] = {0., 25., 35., 45., 55., 65., 75., 85., 97., 110. 200.};
//double nch_binedge[nch_bin+1] = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200.};
//double nch_binedge[nch_bin+1] = {0., 10., 20., 30., 40., 50., 60., 72., 85., 200.};
double nch_binedge[nch_bin+1] = {0., 10., 20., 30., 40., 50., 60., 70., 85.,110.,150., 200.};
//double nch_binedge[nch_bin+1] = {0., 20., 30., 40., 50., 60., 70., 80., 92., 200.};
TH1D* hnch_bin_hist = new TH1D("hnch_bin_hist", "", nch_bin, nch_binedge);

// jet pt bin
const int pt_bin = 3;
double pt_binedge[pt_bin+1] = {119.9, 150., 200., 500.1};
//double pt_binedge[pt_bin+1] = {119.9, 160., 220., 500.1};
TH1D* hpt_bin_hist = new TH1D("hpt_bin_hist", "", pt_bin, pt_binedge);

// centrality bin
const int NCentbin_bal = 5;
double CentbinEdge_bal[NCentbin_bal+1] = {0., 20., 60., 100., 160., 200};
TH1D* hCentBin_bal = new TH1D("hCentBin_bal", "", NCentbin_bal, CentbinEdge_bal);

// event histo
// for reco
TH1D* hnjets_afterCut = new TH1D("hnjets_afterCut", "", 20, 0, 20);
//for gen
TH1D* hnjets_Gen_afterCut = new TH1D("hnjets_afterCut_gen", "", 20, 0, 20);

// jet hist
//reco
TH1D* hJet_ntrk = new TH1D("hJet_ntrk", "", 300, 0, 300);
TH1D* hJet_nchtrk = new TH1D("hJet_nchtrk", "", 150, 0, 150);
TH2D* hJet_nchtrk_pT = new TH2D("hJet_nchtrk_pT", "", 150, 0, 150, 80, 100, 500);
TH1D* hJet_nchtrk_[nch_bin];
TH1D* hJet_pt = new TH1D("hJet_pt", "", 200, 0., 1000);
//gen
TH1D* hJet_Gen_ntrk = new TH1D("hJet_Gen_ntrk", "", 300, 0, 300);
TH1D* hJet_Gen_nchtrk = new TH1D("hJet_Gen_nchtrk", "", 150, 0, 150);
TH2D* hJet_Gen_nchtrk_pT = new TH2D("hJet_Gen_nchtrk_pT", "", 150, 0, 150, 80, 100, 500);
TH1D* hJet_Gen_nchtrk_[nch_bin];
TH1D* hJet_Gen_pt = new TH1D("hJet_Gen_pt", "", 200, 0., 1000);

// track histo
// all tracks
TH1D* halltrk_pt = new TH1D("halltrk_pt", "", 200, 0., 2000.);
TH1D* halltrk_eta = new TH1D("halltrk_eta", "", 48, -2.4, 2.4);
TH1D* halltrk_phi = new TH1D("halltrk_phi", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* halltrk_pt_jetaxis = new TH1D("halltrk_pt_jetaxis", "", 500, 0., 50.);
TH1D* halltrk_eta_jetaxis = new TH1D("halltrk_eta_jetaxis", "", 150, 0., 15.);
TH1D* halltrk_phi_jetaxis = new TH1D("halltrk_phi_jetaxis", "", 64, -TMath::Pi(), TMath::Pi());

// for gen
TH1D* halltrk_Gen_pt = new TH1D("halltrk_Gen_pt", "", 200, 0., 2000.);
TH1D* halltrk_Gen_eta = new TH1D("halltrk_Gen_eta", "", 48, -2.4, 2.4);
TH1D* halltrk_Gen_phi = new TH1D("halltrk_Gen_phi", "", 64, -TMath::Pi(), TMath::Pi());

//charge tracks
// for reco
TH1D* hchtrk_pt = new TH1D("hchtrk_pt", "", 200, 0., 2000.);
TH1D* hchtrk_eta = new TH1D("hchtrk_eta", "", 48, -2.4, 2.4);
TH1D* hchtrk_phi = new TH1D("hchtrk_phi", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_pt_jetaxis = new TH1D("hchtrk_pt_jetaxis", "", 500, 0., 50.);
TH1D* hchtrk_eta_jetaxis = new TH1D("hchtrk_eta_jetaxis", "", 150, 0., 15.);
TH1D* hchtrk_phi_jetaxis = new TH1D("hchtrk_phi_jetaxis", "", 64, -TMath::Pi(), TMath::Pi());

// for gen
TH1D* hchtrk_Gen_pt = new TH1D("hchtrk_Gen_pt", "", 200, 0., 2000.);
TH1D* hchtrk_Gen_eta = new TH1D("hchtrk_Gen_eta", "", 48, -2.4, 2.4);
TH1D* hchtrk_Gen_phi = new TH1D("hchtrk_Gen_phi", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_Gen_pt_jetaxis = new TH1D("hchtrk_Gen_pt_jetaxis", "", 500, 0., 50.);
TH1D* hchtrk_Gen_eta_jetaxis = new TH1D("hchtrk_Gen_eta_jetaxis", "", 150, 0., 15.);
TH1D* hchtrk_Gen_phi_jetaxis = new TH1D("hchtrk_Gen_phi_jetaxis", "", 64, -TMath::Pi(), TMath::Pi());

// +ve charge tracks
//for reco
TH1D* hchtrk_pt_p = new TH1D("hchtrk_pt_p", "", 200, 0., 2000.);
TH1D* hchtrk_eta_p = new TH1D("hchtrk_eta_p", "", 48, -2.4, 2.4);
TH1D* hchtrk_phi_p = new TH1D("hchtrk_phi_p", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_pt_jetaxis_p = new TH1D("hchtrk_pt_jetaxis_p", "", 500, 0., 50.);
TH1D* hchtrk_eta_jetaxis_p = new TH1D("hchtrk_eta_jetaxis_p", "", 150, 0., 15.);
TH1D* hchtrk_phi_jetaxis_p = new TH1D("hchtrk_phi_jetaxis_p", "", 64, -TMath::Pi(), TMath::Pi());

//for gen
TH1D* hchtrk_Gen_pt_p = new TH1D("hchtrk_Gen_pt_p", "", 200, 0., 2000.);
TH1D* hchtrk_Gen_eta_p = new TH1D("hchtrk_Gen_eta_p", "", 48, -2.4, 2.4);
TH1D* hchtrk_Gen_phi_p = new TH1D("hchtrk_Gen_phi_p", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_Gen_pt_jetaxis_p = new TH1D("hchtrk_Gen_pt_jetaxis_p", "", 500, 0., 50.);
TH1D* hchtrk_Gen_eta_jetaxis_p = new TH1D("hchtrk_Gen_eta_jetaxis_p", "", 150, 0., 15.);
TH1D* hchtrk_Gen_phi_jetaxis_p = new TH1D("hchtrk_Gen_phi_jetaxis_p", "", 64, -TMath::Pi(), TMath::Pi());

// -ve charge tracks
//for reco
TH1D* hchtrk_pt_m = new TH1D("hchtrk_pt_m", "", 200, 0., 2000.);
TH1D* hchtrk_eta_m = new TH1D("hchtrk_eta_m", "", 48, -2.4, 2.4);
TH1D* hchtrk_phi_m = new TH1D("hchtrk_phi_m", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_pt_jetaxis_m = new TH1D("hchtrk_pt_jetaxis_m", "", 500, 0., 50.);
TH1D* hchtrk_eta_jetaxis_m = new TH1D("hchtrk_eta_jetaxis_m", "", 150, 0., 15.);
TH1D* hchtrk_phi_jetaxis_m = new TH1D("hchtrk_phi_jetaxis_m", "", 64, -TMath::Pi(), TMath::Pi());

//for gen
TH1D* hchtrk_Gen_pt_m = new TH1D("hchtrk_Gen_pt_m", "", 200, 0., 2000.);
TH1D* hchtrk_Gen_eta_m = new TH1D("hchtrk_Gen_eta_m", "", 48, -2.4, 2.4);
TH1D* hchtrk_Gen_phi_m = new TH1D("hchtrk_Gen_phi_m", "", 64, -TMath::Pi(), TMath::Pi());
TH1D* hchtrk_Gen_pt_jetaxis_m = new TH1D("hchtrk_Gen_pt_jetaxis_m", "", 500, 0., 50.);
TH1D* hchtrk_Gen_eta_jetaxis_m = new TH1D("hchtrk_Gen_eta_jetaxis_m", "", 150, 0., 15.);
TH1D* hchtrk_Gen_phi_jetaxis_m = new TH1D("hchtrk_Gen_phi_jetaxis_m", "", 64, -TMath::Pi(), TMath::Pi());

// for Ntrig
int    bins4D_newtrk_ntrk[4]   =   { 1  , pt_bin,         8,  NCentbin_bal};
double xmin4D_newtrk_ntrk[4]   =   { 1. , 0.,             0., 0.};
double xmax4D_newtrk_ntrk[4]   =   { 2. , (double)pt_bin, 8., (double)NCentbin_bal};
//for reco
THnSparseD * hnewchtrk_Ntrig = new THnSparseD("hnewchtrk_Ntrig", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);
THnSparseD * hnewchtrk_Ntrig_p = new THnSparseD("hnewchtrk_Ntrig_p", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);
THnSparseD * hnewchtrk_Ntrig_m = new THnSparseD("hnewchtrk_Ntrig_m", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);

//for gen
THnSparseD * hnewchtrk_Gen_Ntrig = new THnSparseD("hnewchtrk_Gen_Ntrig", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);
THnSparseD * hnewchtrk_Gen_Ntrig_p = new THnSparseD("hnewchtrk_Gen_Ntrig_p", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);
THnSparseD * hnewchtrk_Gen_Ntrig_m = new THnSparseD("hnewchtrk_Gen_Ntrig_m", "", 4, bins4D_newtrk_ntrk, xmin4D_newtrk_ntrk, xmax4D_newtrk_ntrk);

// signal histo

double etaW = (4.0*newtrk_eta_max)/nDEtaBins;
double phiW = 2.0*(TMath::Pi())/nDPhiBins;

double Deta_min =  -(2.0*newtrk_eta_max) - etaW/2.;
double Deta_max = (2.0*newtrk_eta_max) + etaW/2.;

double Dphi_min = -( TMath::Pi() - phiW ) / 2.0;
double Dphi_max = ( TMath::Pi() * 3.0 - phiW) / 2.0;


const int ndim_hsignal = 5;
int bin_hsignal_newchtrk_jt03[ndim_hsignal] =     {nDEtaBins+1,  nDPhiBins-1,  pt_bin,         8,  NCentbin_bal        };
double xmin_hsignal_newchtrk_jt03[ndim_hsignal] = {Deta_min,  Dphi_min,  0.,             0., 0.                  };
double xmax_hsignal_newchtrk_jt03[ndim_hsignal] = {Deta_max,  Dphi_max,  (double)pt_bin, 8., (double)NCentbin_bal};

const int ndim_qinv = 5;
int bin_qinv[ndim_qinv] =     {500,  50 , 200 , 8 , NCentbin_bal        };
double xmin_qinv[ndim_qinv] = {0. ,  0. , 0.  , 0., 0.,                 };
double xmax_qinv[ndim_qinv] = {10.,  50., 200., 8., (double)NCentbin_bal};

//for reco
THnSparseD* hsignal_newchtrk_jt03 = new THnSparseD("hsignal_newchtrk_jt03", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03); 
THnSparseD* hsignal_newchtrk_jt03_pp = new THnSparseD("hsignal_newchtrk_jt03_pp", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_newchtrk_jt03_mm = new THnSparseD("hsignal_newchtrk_jt03_mm", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_newchtrk_jt03_pm = new THnSparseD("hsignal_newchtrk_jt03_pm", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_newchtrk_jt03_mp = new THnSparseD("hsignal_newchtrk_jt03_mp", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);

// for qinv
THnSparseD* hsignal_qinv_ls = new THnSparseD("hsignal_qinv_ls", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);
THnSparseD* hsignal_qinv_us = new THnSparseD("hsignal_qinv_us", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);

//for gen
THnSparseD* hsignal_Gen_newchtrk_jt03 = new THnSparseD("hsignal_Gen_newchtrk_jt03", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03); 
THnSparseD* hsignal_Gen_newchtrk_jt03_pp = new THnSparseD("hsignal_Gen_newchtrk_jt03_pp", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_Gen_newchtrk_jt03_mm = new THnSparseD("hsignal_Gen_newchtrk_jt03_mm", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_Gen_newchtrk_jt03_pm = new THnSparseD("hsignal_Gen_newchtrk_jt03_pm", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);
THnSparseD* hsignal_Gen_newchtrk_jt03_mp = new THnSparseD("hsignal_Gen_newchtrk_jt03_mp", "", ndim_hsignal, bin_hsignal_newchtrk_jt03, xmin_hsignal_newchtrk_jt03, xmax_hsignal_newchtrk_jt03);

// for qinv
THnSparseD* hsignal_gen_qinv_ls = new THnSparseD("hsignal_gen_qinv_ls", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);
THnSparseD* hsignal_gen_qinv_us = new THnSparseD("hsignal_gen_qinv_us", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);

// 1D histogram
// for reco
TH1D* hnewchtrk_Ntrig_1D[NCentbin_bal][pt_bin];
TH1D* hnewchtrk_Ntrig_p_1D[NCentbin_bal][pt_bin];
TH1D* hnewchtrk_Ntrig_m_1D[NCentbin_bal][pt_bin];

// for gen
TH1D* hnewchtrk_Gen_Ntrig_1D[NCentbin_bal][pt_bin];
TH1D* hnewchtrk_Gen_Ntrig_p_1D[NCentbin_bal][pt_bin];
TH1D* hnewchtrk_Gen_Ntrig_m_1D[NCentbin_bal][pt_bin];

// 2D histogram
//for reco
TH2D* hsignal_newchtrk_jt03_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_newchtrk_jt03_pp_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_newchtrk_jt03_mm_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_newchtrk_jt03_pm_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_newchtrk_jt03_mp_2D[NCentbin_bal][pt_bin];

//for gen
TH2D* hsignal_Gen_newchtrk_jt03_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_Gen_newchtrk_jt03_pp_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_Gen_newchtrk_jt03_mm_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_Gen_newchtrk_jt03_pm_2D[NCentbin_bal][pt_bin];
TH2D* hsignal_Gen_newchtrk_jt03_mp_2D[NCentbin_bal][pt_bin];

// mixing histo

const int ndim_hmixing = 5;
int bin_hmixing_newchtrk_jt03[ndim_hmixing] =     {nDEtaBins+1,  nDPhiBins-1,  pt_bin,         8,  NCentbin_bal        };
double xmin_hmixing_newchtrk_jt03[ndim_hmixing] = {Deta_min,  Dphi_min,  0.,             0., 0.                  };
double xmax_hmixing_newchtrk_jt03[ndim_hmixing] = {Deta_max,  Dphi_max,  (double)pt_bin, 8., (double)NCentbin_bal};

//for reco
THnSparseD* hmixing_newchtrk_jt03 = new THnSparseD("hmixing_newchtrk_jt03", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_newchtrk_jt03_pp = new THnSparseD("hmixing_newchtrk_jt03_pp", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_newchtrk_jt03_mm = new THnSparseD("hmixing_newchtrk_jt03_mm", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_newchtrk_jt03_pm = new THnSparseD("hmixing_newchtrk_jt03_pm", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_newchtrk_jt03_mp = new THnSparseD("hmixing_newchtrk_jt03_mp", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);

// for qinv
THnSparseD* hmixing_qinv_ls = new THnSparseD("hmixing_qinv_ls", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);
THnSparseD* hmixing_qinv_us = new THnSparseD("hmixing_qinv_us", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);

//for gen
THnSparseD* hmixing_Gen_newchtrk_jt03 = new THnSparseD("hmixing_Gen_newchtrk_jt03", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_Gen_newchtrk_jt03_pp = new THnSparseD("hmixing_Gen_newchtrk_jt03_pp", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_Gen_newchtrk_jt03_mm = new THnSparseD("hmixing_Gen_newchtrk_jt03_mm", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_Gen_newchtrk_jt03_pm = new THnSparseD("hmixing_Gen_newchtrk_jt03_pm", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);
THnSparseD* hmixing_Gen_newchtrk_jt03_mp = new THnSparseD("hmixing_Gen_newchtrk_jt03_mp", "", ndim_hmixing, bin_hmixing_newchtrk_jt03, xmin_hmixing_newchtrk_jt03, xmax_hmixing_newchtrk_jt03);

// for qinv
THnSparseD* hmixing_gen_qinv_ls = new THnSparseD("hmixing_gen_qinv_ls", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);
THnSparseD* hmixing_gen_qinv_us = new THnSparseD("hmixing_gen_qinv_us", "", ndim_qinv, bin_qinv, xmin_qinv, xmax_qinv);

// 2D histogram
//for reco
TH2D* hmixing_newchtrk_jt03_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_newchtrk_jt03_pp_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_newchtrk_jt03_mm_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_newchtrk_jt03_pm_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_newchtrk_jt03_mp_2D[NCentbin_bal][pt_bin];

//for gen
TH2D* hmixing_Gen_newchtrk_jt03_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_Gen_newchtrk_jt03_pp_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_Gen_newchtrk_jt03_mm_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_Gen_newchtrk_jt03_pm_2D[NCentbin_bal][pt_bin];
TH2D* hmixing_Gen_newchtrk_jt03_mp_2D[NCentbin_bal][pt_bin];

void array_hist_def()
{
  for(int inchbin = 0; inchbin < nch_bin; inchbin++)
    {
      // jet histo
      // for reco
      hJet_nchtrk_[inchbin] = new TH1D(Form("hJet_nchtrk_%d", inchbin), "", 200, 0, 200);

      // for gen
      hJet_Gen_nchtrk_[inchbin] = new TH1D(Form("hJet_Gen_nchtrk_%d", inchbin), "", 200, 0, 200);
    }

  // 2D histogram bin;

  for(int ictbin = 0; ictbin < NCentbin_bal; ictbin++)
    {
      for(int iptbin = 0; iptbin < pt_bin; iptbin++)
	{
	  //for reco
	  // for Ntrig
	  hnewchtrk_Ntrig_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Ntrig_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Ntrig_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);
	  hnewchtrk_Ntrig_p_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Ntrig_p_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Ntrig_p_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);
	  hnewchtrk_Ntrig_m_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Ntrig_m_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Ntrig_m_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);

	  // for signal
	  hsignal_newchtrk_jt03_2D[ictbin][iptbin] = new TH2D(Form("hsignal_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), Form("hsignal_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_newchtrk_jt03_pp_2D[ictbin][iptbin] = new TH2D(Form("hsignal_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), Form("hsignal_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_newchtrk_jt03_mm_2D[ictbin][iptbin] = new TH2D(Form("hsignal_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), Form("hsignal_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_newchtrk_jt03_pm_2D[ictbin][iptbin] = new TH2D(Form("hsignal_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), Form("hsignal_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_newchtrk_jt03_mp_2D[ictbin][iptbin] = new TH2D(Form("hsignal_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), Form("hsignal_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);

	  // for mixing
	   hmixing_newchtrk_jt03_2D[ictbin][iptbin] = new TH2D(Form("hmixing_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), Form("hmixing_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_newchtrk_jt03_pp_2D[ictbin][iptbin] = new TH2D(Form("hmixing_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), Form("hmixing_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_newchtrk_jt03_mm_2D[ictbin][iptbin] = new TH2D(Form("hmixing_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), Form("hmixing_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_newchtrk_jt03_pm_2D[ictbin][iptbin] = new TH2D(Form("hmixing_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), Form("hmixing_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_newchtrk_jt03_mp_2D[ictbin][iptbin] = new TH2D(Form("hmixing_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), Form("hmixing_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  
	  //for Gen
	  hnewchtrk_Gen_Ntrig_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Gen_Ntrig_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Gen_Ntrig_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);
	  hnewchtrk_Gen_Ntrig_p_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Gen_Ntrig_p_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Gen_Ntrig_p_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);
	  hnewchtrk_Gen_Ntrig_m_1D[ictbin][iptbin] = new TH1D(Form("hnewchtrk_Gen_Ntrig_m_1D_%d_%d", ictbin, iptbin), Form("hnewchtrk_Gen_Ntrig_m_1D_%d_%d", ictbin, iptbin), 1, 1., 2.);
	  
	  // for signal
	  hsignal_Gen_newchtrk_jt03_2D[ictbin][iptbin] = new TH2D(Form("hsignal_Gen_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), Form("hsignal_Gen_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin] = new TH2D(Form("hsignal_Gen_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), Form("hsignal_Gen_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin] = new TH2D(Form("hsignal_Gen_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), Form("hsignal_Gen_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin] = new TH2D(Form("hsignal_Gen_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), Form("hsignal_Gen_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hsignal_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin] = new TH2D(Form("hsignal_Gen_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), Form("hsignal_Gen_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);

	  // for mixing
	   hmixing_Gen_newchtrk_jt03_2D[ictbin][iptbin] = new TH2D(Form("hmixing_Gen_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), Form("hmixing_Gen_newchtrk_jt03_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin] = new TH2D(Form("hmixing_Gen_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), Form("hmixing_Gen_newchtrk_jt03_pp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin] = new TH2D(Form("hmixing_Gen_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), Form("hmixing_Gen_newchtrk_jt03_mm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin] = new TH2D(Form("hmixing_Gen_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), Form("hmixing_Gen_newchtrk_jt03_pm_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	  hmixing_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin] = new TH2D(Form("hmixing_Gen_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), Form("hmixing_Gen_newchtrk_jt03_mp_2D_%d_%d", ictbin, iptbin), nDEtaBins+1, Deta_min, Deta_max, nDPhiBins-1, Dphi_min, Dphi_max);
	}
    }
}

void sumw2_bal()
{
  // event histo
  //for reco
  hnjets_afterCut->Sumw2();

  //for gen
  hnjets_Gen_afterCut->Sumw2();
  
  // jet histo
  // for reco
  hJet_ntrk->Sumw2();
  hJet_nchtrk->Sumw2();
  hJet_nchtrk_pT->Sumw2();
  hJet_pt->Sumw2();
  
  for(int inchbin = 0; inchbin < nch_bin; inchbin++)
    {
      hJet_nchtrk_[inchbin]->Sumw2();
    }

  // for gen
  hJet_Gen_ntrk->Sumw2();
  hJet_Gen_nchtrk->Sumw2();
  hJet_Gen_nchtrk_pT->Sumw2();
  hJet_Gen_pt->Sumw2();
  
  for(int inchbin = 0; inchbin < nch_bin; inchbin++)
    {
      hJet_Gen_nchtrk_[inchbin]->Sumw2();
    }
  
  // track histo
  //all tracks
  halltrk_pt->Sumw2();
  halltrk_eta->Sumw2();
  halltrk_phi->Sumw2();
  halltrk_pt_jetaxis->Sumw2();
  halltrk_eta_jetaxis->Sumw2();
  halltrk_phi_jetaxis->Sumw2();

  // for gen
  halltrk_Gen_pt->Sumw2();
  halltrk_Gen_eta->Sumw2();
  halltrk_Gen_phi->Sumw2();
  
  //charge tracks
  // for reco
  hchtrk_pt->Sumw2();
  hchtrk_eta->Sumw2();
  hchtrk_phi->Sumw2();
  hchtrk_pt_jetaxis->Sumw2();
  hchtrk_eta_jetaxis->Sumw2();
  hchtrk_phi_jetaxis->Sumw2();

  // for gen
  hchtrk_Gen_pt->Sumw2();
  hchtrk_Gen_eta->Sumw2();
  hchtrk_Gen_phi->Sumw2();
  hchtrk_Gen_pt_jetaxis->Sumw2();
  hchtrk_Gen_eta_jetaxis->Sumw2();
  hchtrk_Gen_phi_jetaxis->Sumw2();

  // +ve charge tracks
  //for reco
  hchtrk_pt_p->Sumw2();
  hchtrk_eta_p->Sumw2();
  hchtrk_phi_p->Sumw2();
  hchtrk_pt_jetaxis_p->Sumw2();
  hchtrk_eta_jetaxis_p->Sumw2();
  hchtrk_phi_jetaxis_p->Sumw2();

  //for gen
  hchtrk_Gen_pt_p->Sumw2();
  hchtrk_Gen_eta_p->Sumw2();
  hchtrk_Gen_phi_p->Sumw2();
  hchtrk_Gen_pt_jetaxis_p->Sumw2();
  hchtrk_Gen_eta_jetaxis_p->Sumw2();
  hchtrk_Gen_phi_jetaxis_p->Sumw2();

  // -ve charge tracks
  // for reco
  hchtrk_pt_m->Sumw2();
  hchtrk_eta_m->Sumw2();
  hchtrk_phi_m->Sumw2();
  hchtrk_pt_jetaxis_m->Sumw2();
  hchtrk_eta_jetaxis_m->Sumw2();
  hchtrk_phi_jetaxis_m->Sumw2();

  // for gen
  hchtrk_Gen_pt_m->Sumw2();
  hchtrk_Gen_eta_m->Sumw2();
  hchtrk_Gen_phi_m->Sumw2();
  hchtrk_Gen_pt_jetaxis_m->Sumw2();
  hchtrk_Gen_eta_jetaxis_m->Sumw2();
  hchtrk_Gen_phi_jetaxis_m->Sumw2();

  //Ntrig
  //for reco
  hnewchtrk_Ntrig->Sumw2();
  hnewchtrk_Ntrig_p->Sumw2();
  hnewchtrk_Ntrig_m->Sumw2();

  //for gen
  hnewchtrk_Gen_Ntrig->Sumw2();
  hnewchtrk_Gen_Ntrig_p->Sumw2();
  hnewchtrk_Gen_Ntrig_m->Sumw2();
  
  // correlation 
  //signal histo
  //for reco
  hsignal_newchtrk_jt03->Sumw2();
  hsignal_newchtrk_jt03_pp->Sumw2();
  hsignal_newchtrk_jt03_mm->Sumw2();
  hsignal_newchtrk_jt03_pm->Sumw2();
  hsignal_newchtrk_jt03_mp->Sumw2();

  // for qinv
  hsignal_qinv_ls->Sumw2();
  hsignal_qinv_us->Sumw2();
  
  //for gen
  hsignal_Gen_newchtrk_jt03->Sumw2();
  hsignal_Gen_newchtrk_jt03_pp->Sumw2();
  hsignal_Gen_newchtrk_jt03_mm->Sumw2();
  hsignal_Gen_newchtrk_jt03_pm->Sumw2();
  hsignal_Gen_newchtrk_jt03_mp->Sumw2();

  // for qinv
  hsignal_gen_qinv_ls->Sumw2();
  hsignal_gen_qinv_us->Sumw2();
  
  // mixing histo
  //for reco
  hmixing_newchtrk_jt03->Sumw2();
  hmixing_newchtrk_jt03_pp->Sumw2();
  hmixing_newchtrk_jt03_mm->Sumw2();
  hmixing_newchtrk_jt03_pm->Sumw2();
  hmixing_newchtrk_jt03_mp->Sumw2();

  // for qinv
  hmixing_qinv_ls->Sumw2();
  hmixing_qinv_us->Sumw2();
  
  //for gen
  hmixing_Gen_newchtrk_jt03->Sumw2();
  hmixing_Gen_newchtrk_jt03_pp->Sumw2();
  hmixing_Gen_newchtrk_jt03_mm->Sumw2();
  hmixing_Gen_newchtrk_jt03_pm->Sumw2();
  hmixing_Gen_newchtrk_jt03_mp->Sumw2();

  // for qinv
  hmixing_gen_qinv_ls->Sumw2();
  hmixing_gen_qinv_us->Sumw2();
  
  for(int ictbin = 0; ictbin < NCentbin_bal; ictbin++)
    {
      for(int iptbin = 0; iptbin < pt_bin; iptbin++)
        {
	  //Ntrig
	  //for reco
	  hnewchtrk_Ntrig_1D[ictbin][iptbin]->Sumw2();
	  hnewchtrk_Ntrig_p_1D[ictbin][iptbin]->Sumw2();
	  hnewchtrk_Ntrig_m_1D[ictbin][iptbin]->Sumw2();
	  
	  //for gen
	  hnewchtrk_Gen_Ntrig_1D[ictbin][iptbin]->Sumw2();
	  hnewchtrk_Gen_Ntrig_p_1D[ictbin][iptbin]->Sumw2();
	  hnewchtrk_Gen_Ntrig_m_1D[ictbin][iptbin]->Sumw2();
	  
	  // correlation 
	  //signal histo
	  //for reco
	  hsignal_newchtrk_jt03_2D[ictbin][iptbin]->Sumw2();
	  hsignal_newchtrk_jt03_pp_2D[ictbin][iptbin]->Sumw2();
	  hsignal_newchtrk_jt03_mm_2D[ictbin][iptbin]->Sumw2();
	  hsignal_newchtrk_jt03_pm_2D[ictbin][iptbin]->Sumw2();
	  hsignal_newchtrk_jt03_mp_2D[ictbin][iptbin]->Sumw2();
	  
	  //for gen
	  hsignal_Gen_newchtrk_jt03_2D[ictbin][iptbin]->Sumw2();
	  hsignal_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin]->Sumw2();
	  hsignal_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin]->Sumw2();
	  hsignal_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin]->Sumw2();
	  hsignal_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin]->Sumw2();
	  
	  // mixing histo
	  //for reco
	  hmixing_newchtrk_jt03_2D[ictbin][iptbin]->Sumw2();
	  hmixing_newchtrk_jt03_pp_2D[ictbin][iptbin]->Sumw2();
	  hmixing_newchtrk_jt03_mm_2D[ictbin][iptbin]->Sumw2();
	  hmixing_newchtrk_jt03_pm_2D[ictbin][iptbin]->Sumw2();
	  hmixing_newchtrk_jt03_mp_2D[ictbin][iptbin]->Sumw2();
	  
	  //for gen
	  hmixing_Gen_newchtrk_jt03_2D[ictbin][iptbin]->Sumw2();
	  hmixing_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin]->Sumw2();
	  hmixing_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin]->Sumw2();
	  hmixing_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin]->Sumw2();
	  hmixing_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin]->Sumw2();
	  
	}
    }
}

void write_event_hist_bal(bool ismc)
{
  hnjets_afterCut->Write();

  // for gen
  if(ismc)
    {
      hnjets_Gen_afterCut->Write();
    }
}

void write_jet_hist_bal(bool ismc)
{
  // for reco
  hJet_pt->Write();
  hJet_ntrk->Write();
  hJet_nchtrk->Write();
  hJet_nchtrk_pT->Write();
  
  for(int inchbin = 0; inchbin < nch_bin; inchbin++)
    {
      hJet_nchtrk_[inchbin]->Write();
    }

  // for gen
  if(ismc)
    {
      hJet_Gen_pt->Write();
      hJet_Gen_ntrk->Write();
      hJet_Gen_nchtrk->Write();
      hJet_Gen_nchtrk_pT->Write();
      
      for(int inchbin = 0; inchbin < nch_bin; inchbin++)
	{
	  hJet_Gen_nchtrk_[inchbin]->Write();
	}
    }
}
void write_track_hist_bal(bool ismc)
{
  //all tracks

  halltrk_pt->Write();
  halltrk_eta->Write();
  halltrk_phi->Write();
  /*
  halltrk_pt_jetaxis->Write();
  halltrk_eta_jetaxis->Write();
  halltrk_phi_jetaxis->Write();
  */

  // charge tracks
  // for reco
  hchtrk_pt->Write();
  hchtrk_eta->Write();
  hchtrk_phi->Write();
  hchtrk_pt_jetaxis->Write();
  hchtrk_eta_jetaxis->Write();
  hchtrk_phi_jetaxis->Write();

  // +ve charge tracks
  //for reco
  hchtrk_pt_p->Write();
  hchtrk_eta_p->Write();
  hchtrk_phi_p->Write();
  hchtrk_pt_jetaxis_p->Write();
  hchtrk_eta_jetaxis_p->Write();
  hchtrk_phi_jetaxis_p->Write();

  // -ve charge tracks
  //for reco
  hchtrk_pt_m->Write();
  hchtrk_eta_m->Write();
  hchtrk_phi_m->Write();
  hchtrk_pt_jetaxis_m->Write();
  hchtrk_eta_jetaxis_m->Write();
  hchtrk_phi_jetaxis_m->Write();


  if(ismc)
    {
      halltrk_Gen_pt->Write();
      halltrk_Gen_eta->Write();
      halltrk_Gen_phi->Write();
      
      // charge tracks
      // for gen
      hchtrk_Gen_pt->Write();
      hchtrk_Gen_eta->Write();
      hchtrk_Gen_phi->Write();
      hchtrk_Gen_pt_jetaxis->Write();
      hchtrk_Gen_eta_jetaxis->Write();
      hchtrk_Gen_phi_jetaxis->Write();
      
      // +ve charge tracks
      //for gen
      hchtrk_Gen_pt_p->Write();
      hchtrk_Gen_eta_p->Write();
      hchtrk_Gen_phi_p->Write();
      hchtrk_Gen_pt_jetaxis_p->Write();
      hchtrk_Gen_eta_jetaxis_p->Write();
      hchtrk_Gen_phi_jetaxis_p->Write();
      
      // -ve charge tracks
      //for gen
      hchtrk_Gen_pt_m->Write();
      hchtrk_Gen_eta_m->Write();
      hchtrk_Gen_phi_m->Write();
      hchtrk_Gen_pt_jetaxis_m->Write();
      hchtrk_Gen_eta_jetaxis_m->Write();
      hchtrk_Gen_phi_jetaxis_m->Write();
    }
}

void write_corr_hist_bal(bool ismc)
{
  //Ntrig
  //for reco
  hnewchtrk_Ntrig->Write();
  hnewchtrk_Ntrig_p->Write();
  hnewchtrk_Ntrig_m->Write();

  // signal histo
  //for reco
  hsignal_newchtrk_jt03->Write();
  hsignal_newchtrk_jt03_pp->Write();
  hsignal_newchtrk_jt03_mm->Write();
  hsignal_newchtrk_jt03_pm->Write();
  hsignal_newchtrk_jt03_mp->Write();

  // for qinv
  hsignal_qinv_ls->Write();
  hsignal_qinv_us->Write();
  
  // mixing histo
  //for reco
  hmixing_newchtrk_jt03->Write();
  hmixing_newchtrk_jt03_pp->Write();
  hmixing_newchtrk_jt03_mm->Write();
  hmixing_newchtrk_jt03_pm->Write();
  hmixing_newchtrk_jt03_mp->Write();

  // for qinv
  hmixing_qinv_ls->Write();
  hmixing_qinv_us->Write();
  
  //for gen
  if(ismc)
    {
      
      //Ntrig
      hnewchtrk_Gen_Ntrig->Write();
      hnewchtrk_Gen_Ntrig_p->Write();
      hnewchtrk_Gen_Ntrig_m->Write();

      // signal
      hsignal_Gen_newchtrk_jt03->Write();
      hsignal_Gen_newchtrk_jt03_pp->Write();
      hsignal_Gen_newchtrk_jt03_mm->Write();
      hsignal_Gen_newchtrk_jt03_pm->Write();
      hsignal_Gen_newchtrk_jt03_mp->Write();
      
      // for qinv
      hsignal_gen_qinv_ls->Write();
      hsignal_gen_qinv_us->Write();
  
      // mixing
      hmixing_Gen_newchtrk_jt03->Write();
      hmixing_Gen_newchtrk_jt03_pp->Write();
      hmixing_Gen_newchtrk_jt03_mm->Write();
      hmixing_Gen_newchtrk_jt03_pm->Write();
      hmixing_Gen_newchtrk_jt03_mp->Write();

      // for qinv
      hmixing_gen_qinv_ls->Write();
      hmixing_gen_qinv_us->Write();
    }
}

void write_corr_hist_bal_2D(bool ismc)
{
  for(int ictbin = 0; ictbin < NCentbin_bal; ictbin++)
    {
      for(int iptbin = 0; iptbin < pt_bin; iptbin++)
        {
	  //Ntrig
	  //for reco
	  hnewchtrk_Ntrig_1D[ictbin][iptbin]->Write();
	  hnewchtrk_Ntrig_p_1D[ictbin][iptbin]->Write();
	  hnewchtrk_Ntrig_m_1D[ictbin][iptbin]->Write();
	  
	  // correlation 
	  //signal histo
	  //for reco
	  hsignal_newchtrk_jt03_2D[ictbin][iptbin]->Write();
	  hsignal_newchtrk_jt03_pp_2D[ictbin][iptbin]->Write();
	  hsignal_newchtrk_jt03_mm_2D[ictbin][iptbin]->Write();
	  hsignal_newchtrk_jt03_pm_2D[ictbin][iptbin]->Write();
	  hsignal_newchtrk_jt03_mp_2D[ictbin][iptbin]->Write();
	  
	  // mixing histo
	  //for reco
	  hmixing_newchtrk_jt03_2D[ictbin][iptbin]->Write();
	  hmixing_newchtrk_jt03_pp_2D[ictbin][iptbin]->Write();
	  hmixing_newchtrk_jt03_mm_2D[ictbin][iptbin]->Write();
	  hmixing_newchtrk_jt03_pm_2D[ictbin][iptbin]->Write();
	  hmixing_newchtrk_jt03_mp_2D[ictbin][iptbin]->Write();

	  if(ismc)
	    {
	      //for gen
	      hnewchtrk_Gen_Ntrig_1D[ictbin][iptbin]->Write();
	      hnewchtrk_Gen_Ntrig_p_1D[ictbin][iptbin]->Write();
	      hnewchtrk_Gen_Ntrig_m_1D[ictbin][iptbin]->Write();
	      
	      //for gen
	      hsignal_Gen_newchtrk_jt03_2D[ictbin][iptbin]->Write();
	      hsignal_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin]->Write();
	      hsignal_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin]->Write();
	      hsignal_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin]->Write();
	      hsignal_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin]->Write();
	      
	      //for gen
	      hmixing_Gen_newchtrk_jt03_2D[ictbin][iptbin]->Write();
	      hmixing_Gen_newchtrk_jt03_pp_2D[ictbin][iptbin]->Write();
	      hmixing_Gen_newchtrk_jt03_mm_2D[ictbin][iptbin]->Write();
	      hmixing_Gen_newchtrk_jt03_pm_2D[ictbin][iptbin]->Write();
	      hmixing_Gen_newchtrk_jt03_mp_2D[ictbin][iptbin]->Write();
	    }
	}
    }
}
