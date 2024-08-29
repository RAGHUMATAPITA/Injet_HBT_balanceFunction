#include "call_libraries.h"  // call libraries from ROOT and C++

bool is_MC = false; // for MC: is_MC = true; for data: is_MC = false
bool is_JES_JER = false; // for JES and JER histogram, is_JES_JER = true; else false
//~~~~~~~~~~~~~~~~~~~~~~Centrality, Vz, and pT weight~~~~~~~~~~~~~~~~~~~~~~~~~  
// apply cent bin and vertx z weight to match MC reco with data
bool is_CentBin_and_VZ_Weight = false;
TString CentBin_and_VZ_Fun;
TString pp_CentBin_and_VZ_Fun = "pp_hCent_fVzx_Weight_DataOverMC.root";
TString PbPb_CentBin_and_VZ_Fun = "PbPb_AOD_DataMC_hCent_fVzx_Weight_DataOverMC.root";
//TString PbPb_CentBin_and_VZ_Fun = "PbPb_hCent_fVzx_Weight_DataOverMC.root";

// apply pt weight to match MC reco with data
bool is_ptWeight = false;
TString ptWeightFun;
TString pp_ptWeight = "pp_fptFunc_Weight_MCoverData_ak4PFJets_Jet80Trig_Jetpt120_500"; // pt weight for Jet80 Trigger for pp
TString PbPb_ptWeight = "PbPb_fptFunc_Weight_MCoverData_ak4PFJets_Jet80Trig_Jetpt120_500"; // pt weight for Jet80 Trigger for PbPb

//~~~~~~~~~~~~~~~~~~~~~~Some cuts~~~~~~~~~~~~~~~~~~~~~~~~~
//event quantities
float vz_cut_min = -15.0; //-ve vz acceptance
float vz_cut_max = 15.0; //+ve vz acceptance
bool use_cent = false; // for multiplicty, it will false
int centCut = 160; // Maximum hiBin cut
int mult_Cut = 400; // Maximum multiplicity cut

std::vector<TString> event_filter_str; // skimed event filters

float pthat_cut = 50.0; // pthat cut

//Jets quantities
float jet_radius = 0.4;

float jet_pt_min_cut = 120.0; //120 jet min pT cut 
float jet_pt_max_cut = 500.0; // jet max pT cut

float jet_eta_min_cut = -1.6; // jet min eta cut 
float jet_eta_max_cut = 1.6; // jet min eta cut

//~~~~~~~~~~~~~~~~~~~~~~JET Trigger~~~~~~~~~~~~~~~~~~~~~~~~~
bool is_JetTrigger = false;
TString jet_trigger;
TString PbPb_jet_trigger = "HLT_HIPuAK4CaloJet80Eta5p1_v1"; // jet 80 trigger
//TString PbPb_jet_trigger = "HLT_HIPuAK4CaloJet60Eta5p1_v1"; // jet 60 trigger

TString pp_jet_trigger = "HLT_HIAK4CaloJet80_v1";

//~~~~~~~~~~~~~~~~~~~~~~JET Correction~~~~~~~~~~~~~~~~~~~~~~~~~
TString jet_collection;
TString PbPb_jet_collection = "akCs4PFJetAnalyzer"; //PF jets with CS background subtraction
//TString PbPb_jet_collection = "akFlowPuCs4PFJetAnalyzer"; //PF jets with Flow background subtraction
//TString PbPb_jet_collection = "akPu4CaloJetAnalyzer"; // Calo jet collection in forest

TString pp_jet_collection = "ak4PFJetAnalyzer"; // PF jets
//TString pp_jet_collection = "ak3PFJetAnalyzer"; // PF jets

//~~~~~~~~~~~~~~~~~~~~For JEC Correction~~~~~~~~~~~~~~~~~~~~~~~~~
TString JEC_file;
TString JEC_file_data;
TString PbPb_JEC_file = "Autumn18_HI_V8_MC_L2Relative_AK4PF.txt"; // for PbPb ak4PF jets
TString PbPb_JEC_file_data = "Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"; // residual for PbPb ak4PF jets for data

//TString PbPb_JEC_file = "Autumn18_HI_V8_MC_L2Relative_AK4Calo.txt"; // for PbPb ak4Calo jets
//TString PbPb_JEC_file_data = "Autumn18_HI_V8_DATA_L2L3Residual_AK4Calo.txt"; residual for PbPb ak4Calo jets for data

TString pp_JEC_file = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt"; // for pp ak4PF jets
TString pp_JEC_file_data = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt"; // residual for pp ak4PF jets for data 

//TString pp_JEC_file = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK3PF.txt"; // for pp ak3PF jets
//TString pp_JEC_file_data = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK3PF.txt"; // residual for pp ak3PF jets for data 

//~~~~~~~~~~~~~~~~~~~~For JEU Correction~~~~~~~~~~~~~~~~~~~~~~~~~
bool is_JEU = false;
bool is_JEU_up = false;
bool is_JEU_down = false;
string JEU_file_data;
string pp_JEU_file_data = "JEC_files/Spring18_ppRef5TeV_V6_DATA_Uncertainty_AK4PF.txt"; // for ak4PF Jets
string PbPb_JEU_file_data = "JEC_files/Autumn18_HI_V8_DATA_Uncertainty_AK4PF.txt"; // for akCS4PF Jets

//~~~~~~~~~~~~~~~~~~~~For uncertainty JER Correction~~~~~~~~~~~~~~~~~~~~~~~~~
bool is_JER = false;
bool is_JER_nom = false;
bool is_JER_up = false;
bool is_JER_down = false;
string JER_file_data;
string PbPb_JER_file_data = "JEC_files/JER_Uncertainty_Nom_Up_Down.root"; // for akCS4PF Jets
string pp_JER_file_data = "JEC_files/JER_Uncertainty_Nom_Up_Down.root"; // for akCS4PF Jets

//~~~~~~~~~~~~~~~~~~~~For JER Smearing~~~~~~~~~~~~~~~~~~~~~~~~~
bool is_JER_Correction = false;
TString JER_File;
TString JER_File_pp = "PbPb_JERFitFunc_Nominal_JetTrigger80_CentVzpTWeight_pThat50_JetpT40_1000_akCs4PF_CentBin3"; // for ak4PF Jets
//TString JER_File_PbPb = "PbPb_JERFitFunc_Nominal_NoTrigger_NoWeight_pThat15_JetpT40_1000_akCs4PF_CentBin"; // for akCS4PF Jets
TString JER_File_PbPb = "PbPb_JERFitFunc_Nominal_JetTrigger80_CentVzpTWeight_pThat50_JetpT40_1000_akCs4PF_CentBin";
/*
//Skimmed event filters
std::vector<TString> event_filter_str{"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "phfCoincFilter", "pVertexFilterCutdz1p0"}; // event filters to be applied (pPb - 2016)
std::vector<TString> event_filter_str{"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "phfCoincFilter", "pVertexFilterCutdz1p0"}; // event filters to be applied (XeXe - 2017)
*/

//number of events to mix
int bkgFactor = 30;

//============= Track information =========================
//const std::vector<double> trk_pt_bins{0.7, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0, 300.0}; //trk pT bin range for correlations
float trk_minpt_cut = 0.5;
float trk_maxpt_cut = 999999.; 
float trk_eta_cut = 2.4;
float trk_pt_resolution_cut = 0.1;
float trk_dca_xy_cut = 3.0;
float trk_dca_z_cut = 3.0;
float chi2_ndf_nlayer_cut = 0.18;
float calo_matching = 0.5;
int nhits = 11;
//TString trk_eff_file = "eff_table_p-going_HIJING.root"; //track efficiency table

double newtrk_eta_min = 0.;
double newtrk_eta_max = 100.;
double newtrk_jt03_min = 0.;
double newtrk_jt03_max = 999999.0;

// for 2D histogram binning
const int nDEtaBins = 80; // 80
const int nDPhiBins = 80; // 80, 32
