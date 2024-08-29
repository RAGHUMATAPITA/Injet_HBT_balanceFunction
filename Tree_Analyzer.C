// call_libraries.h and input_variables.h is called inside histogram_definition_BalFunc.h
#include "coordinateTools.h" // for coordinate transformation w.r.t jet axis
#include "function_defination_BalFunc.h" //histogram_definition_BalFunc.h is called within function_defination_BalFunc.h
#include "histogram_definition.h" // define histograms
#include "read_tree.h" // read the TChains
#include "particleid.h"  // call for particle id
#include "JetCorrector.h" // reader for JEC
#include "JetUncertainty.h" // reader for JEU
#include "function_defination.h" // function defined here

void Tree_Analyzer(TString input_file, int itxtoutFile, TString out_file, TString colliding_system, int isMC)
{
  clock_t sec_start, sec_end;
  sec_start = clock();

  TDatime* date = new TDatime();

  print_start(); // start timing print

  // calling function to define array histograms
  array_hist_def();
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Strat Initializing PbPb and pp used input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(isMC == 1)
    {
      is_MC = true; 
    }
  else 
    {
      is_MC = false;
    }

  // for Centrality and Vz weight
  TFile* hCentVzFile;
  TF1* fVzweight;
  TF1* fCentweight;
  TH1D* hCentweight;

  // for pt weight 
  TFile* fptFile_pp;
  TFile* fptFile_PbPb[NCentbin];
  TF1* fptWeight_PbPb[NCentbin];
  TF1* fptWeight;

  // for JER correction from JER fit
  TFile* fJERFile_pp;
  TFile* fJERFile_PbPb[NCentbin];
  TF1* fJERWeight_PbPb[NCentbin];
  TF1* fJERWeight;

  // for uncertainity in JER
  TH1D* hJER_Uncertainty;
  TFile* fJER_Uncertainty;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~:Setting for pp ref data and MC:~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(colliding_system == "pp")
    {
      use_cent = false;
      if(is_MC)
	{
	  is_CentBin_and_VZ_Weight = true;
	  is_ptWeight = true;
	  is_JES_JER = false;
	}

      //if(!is_MC)
      if(!is_MC || is_MC)
	{
          is_JEU = false;
          is_JEU_up = false;
          is_JEU_down = false;
	}

      if(is_MC)
        {
          is_JER = false;
          is_JER_nom = false;
          is_JER_up = false;
          is_JER_down = false;
        }

      is_JetTrigger = true;
      jet_trigger = pp_jet_trigger;

      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  CentBin_and_VZ_Fun = pp_CentBin_and_VZ_Fun;
	  hCentVzFile = TFile::Open(Form("ReweightFile/%s", pp_CentBin_and_VZ_Fun.Data()));
	  fVzweight = (TF1*)hCentVzFile->Get("hvzFun");
	}

      if(is_MC && is_ptWeight)
	{
	  ptWeightFun = pp_ptWeight;
	  fptFile_pp = TFile::Open(Form("ReweightFile/%s.root", ptWeightFun.Data()));
	  fptWeight = (TF1*)fptFile_pp->Get("hptFun");
	}

      if(is_MC && is_JER)
	{
	  JER_File = JER_File_pp;
          fJERFile_pp = TFile::Open(Form("ReweightFile/%s.root", JER_File.Data()));
          fJERWeight = (TF1*)fJERFile_pp->Get("JER_Fit");
	}

      jet_collection = pp_jet_collection;
      JEC_file = pp_JEC_file;
      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  JEC_file_data = pp_JEC_file_data;
	  if(is_JEU)	  
	    {
	      JEU_file_data = pp_JEU_file_data;
	    }
	}

      if(is_MC && is_JER)
        {
          JER_file_data = pp_JER_file_data;
        }

      event_filter_str.resize(0);
      event_filter_str.push_back("pBeamScrapingFilter");
      event_filter_str.push_back("pPAprimaryVertexFilter");
      event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
      
    }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~:Setting for PbPb data and MC:~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(colliding_system == "PbPb")
    {
      use_cent = true;
      if(is_MC) 
	{
	  is_CentBin_and_VZ_Weight = true;
	  is_ptWeight = true;
	  is_JES_JER = false;
	}
      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  is_JEU = false;
	  is_JEU_up = false;
	  is_JEU_down = false;
	}

      if(is_MC)
	{
	  is_JER = false;
	  is_JER_nom = false;
	  is_JER_up = false;
	  is_JER_down = false;
	}

      is_JetTrigger = true; // default true
      jet_trigger = PbPb_jet_trigger;

      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  CentBin_and_VZ_Fun = PbPb_CentBin_and_VZ_Fun;
	  hCentVzFile = TFile::Open(Form("ReweightFile/%s", PbPb_CentBin_and_VZ_Fun.Data()));
	  fVzweight = (TF1*)hCentVzFile->Get("hvzFun");
	  //hCentweight = (TH1D*)hCentVzFile->Get("hCent");
	  fCentweight = (TF1*)hCentVzFile->Get("hCentFun");
	}

      if(is_MC && is_ptWeight)
        {
          ptWeightFun = PbPb_ptWeight;
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
            {
              fptFile_PbPb[ictt] = TFile::Open(Form("ReweightFile/%s_%d.root", ptWeightFun.Data(), ictt));
	      fptWeight_PbPb[ictt] = (TF1*)fptFile_PbPb[ictt]->Get("hptFun");
            }
        }

      if(is_MC && is_JER) // not to write again did this
        {
          JER_File = JER_File_PbPb;
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
            {
	      fJERFile_PbPb[ictt] = TFile::Open(Form("ReweightFile/%s%d.root", JER_File.Data(), ictt));
	      fJERWeight_PbPb[ictt] = (TF1*)fJERFile_PbPb[ictt]->Get("JER_Fit");
	    }
        }

      jet_collection = PbPb_jet_collection;
      JEC_file = PbPb_JEC_file;
      //if(!is_MC)
      if(!is_MC || is_MC)
	{
	  JEC_file_data = PbPb_JEC_file_data;
	  if(is_JEU)
            {
              JEU_file_data = PbPb_JEU_file_data;
            }
	}

      if(is_MC && is_JER)
	{
	  JER_file_data = PbPb_JER_file_data;
	}

      event_filter_str.resize(0);
      event_filter_str.push_back("pprimaryVertexFilter");
      event_filter_str.push_back("HBHENoiseFilterResultRun2Loose");
      event_filter_str.push_back("collisionEventSelectionAOD");
      event_filter_str.push_back("phfCoincFilter2Th4");
      event_filter_str.push_back("pclusterCompatibilityFilter");
    }

  //~~~~~~~~~~~~~~~~~~~~~End Initializing PbPb and pp used input variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Event quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"colliding_system is: "<<colliding_system.Data()<<std::endl;
  std::cout<<"use Centrality: "<<std::boolalpha<<use_cent<<std::endl;
  std::cout<<"is it Monte Carlo: "<<std::boolalpha<<is_MC<<std::endl;
  std::cout<<"is_JES_JER: "<<std::boolalpha<<is_JES_JER<<std::endl;
  std::cout<<"skimed event_filter_str size is: "<<event_filter_str.size()<<std::endl;
  std::cout<<"skimed event filters are: ";
  for(unsigned int ifl = 0; ifl < event_filter_str.size(); ifl++)
    {
      std::cout<<event_filter_str[ifl].Data()<<", ";
    }
  std::cout<<std::endl;
  std::cout<<"vertex z cut: "<<vz_cut_min<<" to "<<vz_cut_max<<" cm"<<std::endl;
  if(colliding_system == "PbPb")
    {
      std::cout<<"Maximum hiBin cut: "<<centCut<<std::endl;
    }
  if(is_MC)
    {
      std::cout<<"Minimum pThat cut: "<<pthat_cut<<" GeV"<<std::endl;   
    }

  std::cout<<"is_CentBin_and_VZ_Weight: "<<std::boolalpha<<is_CentBin_and_VZ_Weight<<std::endl;
  if(is_MC && is_CentBin_and_VZ_Weight)
    {
      std::cout<<"Cent bin and Vertex z weight file is: "<<hCentVzFile->GetName()<<std::endl;
      //hCentVzFile->Close();
      //delete hCentVzFile;           
    }
  std::cout<<endl;

  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~Jet quantities~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  std::cout<<endl;
  std::cout<<"is_JetTrigger: "<<std::boolalpha<<is_JetTrigger<<std::endl;
  if(is_JetTrigger)
    {
      std::cout<<"Applied jet trigger: "<<jet_trigger.Data()<<std::endl;
    }
  std::cout<<"Jet Collection: "<<jet_collection.Data()<<std::endl;
  std::cout<<"JEC File used: "<<JEC_file.Data()<<std::endl;
  //if(!is_MC) 
  if(!is_MC || is_MC) 
    {
      std::cout<<"JEC File for data used: "<<JEC_file_data.Data()<<std::endl; 
      std::cout<<"is_JEU: "<<std::boolalpha<<is_JEU<<std::endl;
      if(is_JEU)
	{
	  std::cout<<"JEU File for data used: "<<JEU_file_data.c_str()<<std::endl; 
	}
    }

  if(is_MC && is_JER)
    {
      std::cout<<"is_JER: "<<std::boolalpha<<is_JER<<std::endl;
      std::cout<<"JER File for data used: "<<JER_file_data.c_str()<<std::endl;
    } 
  
  std::cout<<"jet Eta cut: "<<jet_eta_min_cut<<" to "<<jet_eta_max_cut<<std::endl;
  std::cout<<"jet pT cut: "<<jet_pt_min_cut<<" to "<<jet_pt_max_cut<<" GeV"<<std::endl;
  std::cout<<"is_ptWeight : "<<std::boolalpha<<is_ptWeight<<std::endl;
  if(is_MC && is_ptWeight)
    {
      if(colliding_system == "pp")
	{
	  std::cout<<"pt weight file is : "<<fptFile_pp->GetName()<<std::endl;
	  fptFile_pp->Close();
          delete fptFile_pp;
	}
      else if(colliding_system == "PbPb")
	{
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
	    {
	      std::cout<<"pt weight file is : "<<fptFile_PbPb[ictt]->GetName()<<std::endl;
	      fptFile_PbPb[ictt]->Close();
	      delete fptFile_PbPb[ictt];
	    }
	}
    }
  
  if(is_MC && is_JER) // not write again for JER uncertainty, did this
    {
      if(colliding_system == "PbPb")
	{
	  for(int ictt = 0; ictt < NCentbin-1; ictt++)
	    {
	      std::cout<<"JER weight file is : "<<fJERFile_PbPb[ictt]->GetName()<<std::endl;
	      fJERFile_PbPb[ictt]->Close();
	      delete fJERFile_PbPb[ictt];
	    }
	} 
    }

  std::cout<<endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End printing the input informations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::cout<<endl;

  sumw2(); // calling sumw2 for all the histograms
  sumw2_bal(); // calling sumw2 for all the histograms
  TH1::SetDefaultSumw2(); // calling sumw2 for all the 1D histograms 
    
  //  input JEC file
  vector<string> JECFiles;
  //JECFiles.push_back(Form("JEC_files/%s", JEC_file.c_str()));
  JECFiles.push_back(Form("JEC_files/%s", JEC_file.Data()));
  if(!is_MC)
    {
      JECFiles.push_back(Form("JEC_files/%s", JEC_file_data.Data())); // for data only
      std::cout<<"It is runing for data (L2L3 residual file is attached)"<<std::endl;
      std::cout<<endl;
    }

  JetCorrector JEC(JECFiles);
  
  // for JEU in data
  JetUncertainty* JEU = NULL;
  
  //if(!is_MC && is_JEU)
  if((!is_MC || is_MC)  && is_JEU)
    {
      JEU = new JetUncertainty(JEU_file_data);
    }

  // for uncertainity in JER
  if(is_MC && is_JER)
    {
      fJER_Uncertainty = new TFile(Form("%s",JER_file_data.c_str()));
      if(is_JER_nom)
	{
	  hJER_Uncertainty = (TH1D*)fJER_Uncertainty->Get("h_nom");
	}
      else if(is_JER_up)
	{
	  hJER_Uncertainty = (TH1D*)fJER_Uncertainty->Get("h_up");
	}
      else if(is_JER_down)
	{
	  hJER_Uncertainty = (TH1D*)fJER_Uncertainty->Get("h_down");
	}
    }

  //  open input forest/skim file
  fstream openInputFile;
  openInputFile.open(Form("%s",input_file.Data()), ios::in);
  if(!openInputFile.is_open())
    {
      cout << "List of input files not founded!" << endl;
      return;
    }

  // Make a chain and a vector of file names
  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(openInputFile, file_chain))
    {
      if(colliding_system == "pp")
	{
	  if(is_MC) file_name_vector.push_back(file_chain.c_str());
	  else file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
	}
      else if (colliding_system == "PbPb")
	{
	  file_name_vector.push_back(Form("root://cmsxrootd.fnal.gov/%s", file_chain.c_str()));
	  //file_name_vector.push_back(Form("davs://xrootd-vanderbilt.sites.opensciencegrid.org:1094/%s", file_chain.c_str()));
	}
    }
  openInputFile.close();

  // Read the trees to be added in the Chain
  TChain *hlt_tree = new TChain("hltanalysis/HltTree");
  TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *ski_tree = new TChain("skimanalysis/HltTree");
  TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data()));
  TChain *trk_tree = new TChain("ppTrack/trackTree");
  TChain *gen_tree = new TChain();
  if(is_MC)
    {
      gen_tree = new TChain("HiGenParticleAna/hi");
    }
  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *testfile = TFile::Open(*listIterator);
      
      if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered))
	{
	  cout << "File: " << *listIterator << " failed to open" << endl;
	  continue;
	}

      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;

      hlt_tree->Add(*listIterator);
      hea_tree->Add(*listIterator);
      ski_tree->Add(*listIterator);
      jet_tree->Add(*listIterator);
      trk_tree->Add(*listIterator);
      if(is_MC)
	{
	  gen_tree->Add(*listIterator);
	}
    }

  hlt_tree->AddFriend(hea_tree);
  hlt_tree->AddFriend(ski_tree);	
  hlt_tree->AddFriend(jet_tree);
  hlt_tree->AddFriend(trk_tree);
  if(is_MC)
    {
      hlt_tree->AddFriend(gen_tree);
    }
  
  // calling function to read forest/skim tree
  read_tree(hlt_tree, is_MC, jet_trigger.Data(), colliding_system.Data(), event_filter_str); // access the tree informations
  
  int nevents = hlt_tree->GetEntries(); // number of events
 
  cout << "Total number of events in those files: "<< nevents << endl;

  //~~~~~~~~~~~~Define vectors use for Signal and Mixed event correlation                                                                  
  // 1D event vector to store event quantities
  std::vector<float> Evtw_vec_1D;
  std::vector<float> vertexz_vec_1D;
  std::vector<int> hiBin_vec_1D;

  // 2D jet vector
  // for gen
  std::vector<std::vector<TVector3>> jet_vec_2D_gn;
  std::vector<std::vector<double>> jetw_Vec_2D_gn;
  std::vector<std::vector<int>>Matched_refpartonB_Vec_gn;
  // for reco and data
  std::vector<std::vector<TVector3>> jet_vec_2D_rc;
  std::vector<std::vector<double>> jetw_Vec_2D_rc;
  std::vector<std::vector<int>> Matched_refpartonB_Vec_rc;

  // 3D trk vector
  // for gen
  std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec_3D_gn;
  std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec_3D_gn;
  std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec_3D_gn;
  std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec_3D_gn;
  std::vector<std::vector<std::vector<int>>> jet_trk_charge_3D_gn;
  // for reco and data
  std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec_3D_rc;
  std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec_3D_rc;
  std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec_3D_rc;
  std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec_3D_rc;
  std::vector<std::vector<std::vector<int>>> jet_trk_charge_3D_rc;

  int count = 0;

  //for(int i = 0; i < nevents; i++) //event loop start
  for(int i = 0; i < 10000; i++) //event loop start
    {
      hlt_tree->GetEntry(i);
      
      if(i%10000 == 0)
	{
	  std::cout<<i<<"  events running"<<std::endl;
	}
      
      if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue; // apply z vertex cut
      
      hEvents->AddBinContent(1,1);

      // determine centrality here
      int centbin;
      int ctbin = -1;
      
      //std::cout<<"Cent Bin is: "<<hiBin<<"  "<<getHiBinFromhiHF(hiHF, 2)<<std::endl;
      //int hiCentBin = hiBin;

      //int hiCentBin = getHiBinFromhiHF(hiHF, 0); // 0 = nominal, 1 = up, 2 = down

      if(use_cent) // if you use centrality
	{
	  if(is_MC)
            {
	      if(hiBin <= 9) continue; 
              centbin = hiBin - 10; // match MC multiplicity with data multiplicity
            }
          else
            {
              centbin = hiBin;
	      //centbin = hiCentBin;
            }

	  if(centbin >= centCut || centbin < 0) continue;
	  //if(centbin >= centCut || centbin < 100) continue;
	  ctbin = hCentBin->FindBin(centbin) - 1;
	}
      else // if you use multiplicity 
	{
	  centbin = 1;
	  if(centbin >= mult_Cut) continue;
	  ctbin = 1;
	}

      hEvents->AddBinContent(2,1);
      
      if(is_MC && is_ptWeight && colliding_system == "PbPb")
        {
	  if(ctbin == 0)
	    {
	      fptWeight = fptWeight_PbPb[0];
	    }
	  else if(ctbin == 1)
	    {
	      fptWeight = fptWeight_PbPb[1];
	    }
	  else if(ctbin == 2)
	    {
	      fptWeight = fptWeight_PbPb[2];
	    }
	  else if(ctbin == 3)
	    {
	      fptWeight = fptWeight_PbPb[3];
	    }
	  else{std::cout<<"Centrality bins are more than 3"<<std::endl; break;}
	}

      if(is_MC && is_JER && colliding_system == "PbPb") // not to write again did this
        {
	  if(ctbin == 0)
	    {
	      fJERWeight = fJERWeight_PbPb[0];
	    }
	  else if(ctbin == 1)
	    {
	      fJERWeight = fJERWeight_PbPb[1];
	    }
	  else if(ctbin == 2)
	    {
	      fJERWeight = fJERWeight_PbPb[2];
	    }
	  else if(ctbin == 3)
	    {
	      fJERWeight = fJERWeight_PbPb[3];
	    }
	  else{std::cout<<"Centrality bins are more than 3"<<std::endl; break;}
	}

      double ptHatw = 1.;
      double pTHat = 0.;
      if(is_MC)
	{
	  ptHatw = weight;
	  pTHat = pthat;
	  if(pTHat == 0 || pTHat <= pthat_cut) continue; // apply pTHat cut
	}
	 
      hEvents->AddBinContent(3,1);
      
      bool skimmed_evtfilter = false;

      for(int ii = 0; ii < (int) event_filter_str.size(); ii++)
	{
	  if (event_filter_bool[ii] != 1) // condition for the skimmed event filters
	    {
	      skimmed_evtfilter = true;
	      break;
	    }
	}

      if(skimmed_evtfilter) continue; // apply the skimmed event filters

      hEvents->AddBinContent(4,1);

      if(is_JetTrigger)
	{
	  if(jet_trigger_bit != 1) continue; // apply jet trigger
	}
      
      hEvents->AddBinContent(5,1);

      if(nref <= 0) continue; // if there is no jets in an event
      
      hEvents->AddBinContent(6,1);

      // determine reco/data event weight here
      double Evtw = ptHatw; // Even weight
      
      if(is_MC && is_CentBin_and_VZ_Weight)
	{
	  if(colliding_system == "pp")
	    {
	      Evtw = ptHatw*(fVzweight->Eval(vertexz));
	    }
	  else if(colliding_system == "PbPb")
	    {
	      //Evtw = ptHatw*(fVzweight->Eval(vertexz))*(hCentweight->GetBinContent(hCentweight->FindBin(centbin)));
	      Evtw = ptHatw*(fVzweight->Eval(vertexz))*(fCentweight->Eval(centbin));
	    }
	}
      
      hEvents->AddBinContent(7,1);

      // Fill the bascis event quantities histograms
      if(is_MC)
	{
	  hpthat->Fill(pTHat, Evtw);
	  hpthatW->Fill(ptHatw);
	}

      hCent->Fill(centbin, Evtw);
      hZvtx->Fill(vertexz, Evtw);

      //~~~~~~~~~~~~Define 1D jet vectors use to fill 2D jet vector
      // for gen
      std::vector<TVector3> jet_vec_1D_gn;
      std::vector<double> jetw_Vec_1D_gn;
      std::vector<int> refpartonB_Vec_gn;
      // for reco and data
      std::vector<TVector3> jet_vec_1D_rc;
      std::vector<double> jetw_Vec_1D_rc;
      std::vector<int> refpartonB_Vec_rc;
      
      //~~~~~~~~~~~~Define 2D trk vectors use to fill 3D trk vector
      // for gen
      std::vector<std::vector<TVector3>> oldtrk_vec_2D_gn;
      std::vector<std::vector<TVector3>> newtrk_vec_2D_gn;
      std::vector<std::vector<TVector3>> oldchtrk_vec_2D_gn;
      std::vector<std::vector<TVector3>> newchtrk_vec_2D_gn;
      std::vector<std::vector<int>> trk_charge_2D_gn;
      // for reco and data
      std::vector<std::vector<TVector3>> oldtrk_vec_2D_rc;
      std::vector<std::vector<TVector3>> newtrk_vec_2D_rc;
      std::vector<std::vector<TVector3>> oldchtrk_vec_2D_rc;
      std::vector<std::vector<TVector3>> newchtrk_vec_2D_rc;
      std::vector<std::vector<int>> trk_charge_2D_rc;
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //Jet loop start
      
      double ld_RawpT=-999., ld_RawEta=-999., ld_RawPhi=-999.; // leading Raw jet quantities                          
      double sld_RawpT=-999., sld_RawEta=-999, sld_RawPhi=-999; // subleading Raw jet quantities                       
      double ld_CorrpT=-999., ld_CorrEta=-999., ld_CorrPhi=-999.; // leading jet quantities                       
      double sld_CorrpT=-999., sld_CorrEta=-999, sld_CorrPhi=-999; // subleading jet quantities    
      double ld_RawWTApT=-999.,  ld_RawWTAEta=-999., ld_RawWTAPhi=-999.; // leading Raw WTA jet quantities             
      double sld_RawWTApT=-999., sld_RawWTAEta=-999, sld_RawWTAPhi=-999; // subleading Raw WTA jet quantities     
      double ld_CorrWTApT=-999.,  ld_CorrWTAEta=-999., ld_CorrWTAPhi=-999.; // leading Jt WTA jet quantities     
      double sld_CorrWTApT=-999., sld_CorrWTAEta=-999, sld_CorrWTAPhi=-999; // subleading Jt WTA jet quantities

      for (int j = 0; j < nref; j++) //Jet loop start
	{
	  double jt_pt = rawpt[j];
	  double jt_eta = jteta[j];
	  double jt_phi = jtphi[j];
	  //double jt_eta = WTAeta[j];
	  //double jt_phi = WTAphi[j];
	  double jt_WTAeta = WTAeta[j];
	  double jt_WTAphi = WTAphi[j];
	  int jt_Flavour = 1;
	  
	  if(trackMax[j]/rawpt[j] < 0.01)continue; // Cut for jets for very low maxium pT track
	  if(trackMax[j]/rawpt[j] > 0.98)continue; // Cut for jets where all the pT is taken by one track

	  // JEC correction
	  JEC.SetJetPT(jt_pt); 
	  JEC.SetJetEta(jt_eta); 
	  JEC.SetJetPhi(jt_phi);

	  float jet_pt_corr = JEC.GetCorrectedPT();
	  float jet_pt_corr_before = JEC.GetCorrectedPT();

	  // apply Jet pt eta cut
	  if(jet_pt_corr < jet_pt_min_cut || jet_pt_corr > jet_pt_max_cut) continue;
	  if(fabs(jt_eta) >= jet_eta_max_cut) continue;

	  count++; // count total number of jets
	  
	  // JEU correction
	  //if(!is_MC && is_JEU)
	  if((!is_MC || is_MC)  && is_JEU)
	    {
	      JEU->SetJetPT(jet_pt_corr);
	      JEU->SetJetEta(jt_eta);
	      JEU->SetJetPhi(jt_phi);

	      if(is_JEU_down && !is_JEU_up)
		{
		  jet_pt_corr = jet_pt_corr*(1. - (JEU->GetUncertainty().first));
		}
	      else if(is_JEU_up && !is_JEU_down)
		{
		  jet_pt_corr = jet_pt_corr*(1. + (JEU->GetUncertainty().second));
		}
	    }
	  
	  if(is_MC && is_JER)
	    {
	      double resolution_factor = 1.;
	      if(jt_eta <= jet_eta_min_cut || jt_eta >= jet_eta_max_cut)
		{
		  resolution_factor = 1.;
		}
	      else
		{
		  resolution_factor = hJER_Uncertainty->GetBinContent( hJER_Uncertainty->GetXaxis()->FindBin(jt_eta) );
		}

	      double extraResolution = TMath::Sqrt(TMath::Max(resolution_factor*resolution_factor-1.0,0.0)); // found jet resolution
	      double sigma_smear  = extraResolution*fJERWeight->Eval(jet_pt_corr); // some % worst --> from JetMET
	      
	      gRandom->SetSeed(0);
	      gRandom = new TRandom3(0);
	      double JER_smear = gRandom->Gaus(1.,sigma_smear);
	      while( JER_smear < 0 ){ JER_smear = gRandom->Gaus(1.,sigma_smear); }
	      jet_pt_corr = jet_pt_corr*JER_smear;

	      //std::cout<<"jet pt: "<<ctbin<<"  "<<jet_pt_corr_before<<"  "<<jet_pt_corr<<"  "<<fJERWeight->Eval(jet_pt_corr)<<"  "<<jteta[j]<<"  "<<resolution_factor<<"  "<<sigma_smear<<"  "<<JER_smear<<std::endl;
	    } // if MC and JER
	  
	  // find leading sub leading jets
	  find_leading_subleading_Jets(jt_pt, jt_eta, jt_phi, ld_RawpT, ld_RawEta, ld_RawPhi, sld_RawpT, sld_RawEta, sld_RawPhi);
	  find_leading_subleading_Jets(jet_pt_corr, jt_eta, jt_phi, ld_CorrpT, ld_CorrEta, ld_CorrPhi, sld_CorrpT, sld_CorrEta, sld_CorrPhi);
	  find_leading_subleading_Jets(jt_pt, jt_WTAeta, jt_WTAphi, ld_RawWTApT, ld_RawWTAEta, ld_RawWTAPhi, sld_RawWTApT, sld_RawWTAEta, sld_RawWTAPhi);
	  find_leading_subleading_Jets(jet_pt_corr, jt_WTAeta, jt_WTAphi, ld_CorrWTApT, ld_CorrWTAEta, ld_CorrWTAPhi, sld_CorrWTApT, sld_CorrWTAEta, sld_CorrWTAPhi);

	  	  
	  double Jet_RawpT_Eta_Phi_ctbin[4] = {jt_pt, jt_eta, jt_phi, (double)ctbin};
	  double Jet_RawpT_WTAEta_WTAPhi_ctbin[4] = {jt_pt, jt_WTAeta, jt_WTAphi, (double)ctbin};
	  
	  double Jet_CorrpT_Eta_Phi_ctbin[4] = {jet_pt_corr, jt_eta, jt_phi, (double)ctbin};
	  double Jet_CorrpT_WTAEta_WTAPhi_ctbin[4] = {jet_pt_corr, jt_WTAeta, jt_WTAphi, (double)ctbin};

	  // apply pt weight for inclusive jets in MC reco
	  double recoweight_corrpt = Evtw;

	  if(is_MC && is_ptWeight)
	    {
	      recoweight_corrpt = Evtw*(1./fptWeight->Eval(jet_pt_corr));	      
	    }

	  // w/o weight 
	  
	  hJet_RawpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_RawpT_Eta_Phi_ctbin);
	  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin);
	  
	  hJet_CorrpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_CorrpT_Eta_Phi_ctbin);
	  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin);
	  
	  // w/ weight 
	  hJet_RawpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_RawpT_Eta_Phi_ctbin, recoweight_corrpt);
	  hJet_RawpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_RawpT_WTAEta_WTAPhi_ctbin, recoweight_corrpt);
	  
	  hJet_CorrpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_CorrpT_Eta_Phi_ctbin, recoweight_corrpt);
	  hJet_CorrpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_CorrpT_WTAEta_WTAPhi_ctbin, recoweight_corrpt);

	  if(is_MC)
	    {
	      double ref_pt = refpt[j]; // gen pt matched with reco pt, -999 entries are unmatched gen
	      
	      double Jet_refpT_ctbin[2] = {ref_pt, (double)ctbin};

	      hJet_refpT_ctbin_nopTCut_noW->Fill(Jet_refpT_ctbin);
	      hJet_refpT_ctbin_nopTCut_W->Fill(Jet_refpT_ctbin, recoweight_corrpt);
	      
	      if(ref_pt > jet_pt_min_cut && ref_pt < jet_pt_max_cut)
		{
		  hJet_refpT_ctbin_pTCut_noW->Fill(Jet_refpT_ctbin);
		  hJet_refpT_ctbin_pTCut_W->Fill(Jet_refpT_ctbin, recoweight_corrpt);
		}
	      
	      int refpartonB = -99;
	      
	      if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6)
		{
		  refpartonB = fabs(refparton_flavorForB[j]);
		}
	      else if(fabs(refparton_flavorForB[j]) == 21)
		{
		  refpartonB = 7;
		}
	      else
		{
		  refpartonB = 0; 
		}
	      
	      if(refpartonB == -99) continue;

	      jt_Flavour = refpartonB;

	      // matched, unmatched, and total jet pT distributions

	      double Jet_CorrpT_refpartonB_ctbin[3] = {jet_pt_corr, (double)jt_Flavour, (double)ctbin};
	      
	      if((int)ref_pt < 0) // unmatched
		{
		  hJet_UCorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, recoweight_corrpt);
		}
	      else if((int)ref_pt >= 0) //matched
		{
		  hJet_MCorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, recoweight_corrpt);
		}
	      //total 
	      hJet_CorrpT_refpartonB_ctbin_W->Fill(Jet_CorrpT_refpartonB_ctbin, recoweight_corrpt);
	      hJet_CorrpT_refpartonB_ctbin_noW->Fill(Jet_CorrpT_refpartonB_ctbin);

	      if(ref_pt > 0.)
		{
		  hJet_refpT_ctbin_pTCut0_noW->Fill(Jet_refpT_ctbin);
		  hJet_refpT_ctbin_pTCut0_noW->Fill(Jet_refpT_ctbin, recoweight_corrpt);
		}
	      
	      if(is_JES_JER)
		{
		  if(ref_pt <= 0.) continue; // discrad the unmacthed jet (ref_pt == -999.)
		  if(jet_pt_corr  <= 0.) continue;
		  if(jt_pt <= 0.) continue;

		  double JES_ratio_rawpT_vs_refpT = jt_pt/ref_pt;
		  double JER_ratio_rawpT_vs_refpT = (jt_pt - ref_pt)/ref_pt;
		  
		  double JES_ratio_CorrpT_vs_refpT = jet_pt_corr/ref_pt;
		  double JER_ratio_CorrpT_vs_refpT = (jet_pt_corr - ref_pt)/ref_pt;
		  
		  //JES
		  double Jes_RawpT_refpt_refpartonB_ctbin[4] = {JES_ratio_rawpT_vs_refpT, ref_pt, (double)refpartonB, (double)ctbin};
		  double Jes_CorrpT_refpt_refpartonB_ctbin[4] = {JES_ratio_CorrpT_vs_refpT, ref_pt, (double)refpartonB, (double)ctbin};
		  
		  //JER
		  double Jer_RawpT_refpt_refpartonB_ctbin[4] = {JER_ratio_rawpT_vs_refpT, ref_pt, (double)refpartonB, (double)ctbin};
		  double Jer_CorrpT_refpt_refpartonB_ctbin[4] = {JER_ratio_CorrpT_vs_refpT, ref_pt, (double)refpartonB, (double)ctbin};
		  
		  //JES
		  //w/o weight
		  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jes_RawpT_refpt_refpartonB_ctbin);
		  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jes_CorrpT_refpt_refpartonB_ctbin);
		  
		  //w/ weight
		  hJes_rawpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jes_RawpT_refpt_refpartonB_ctbin, recoweight_corrpt);
		  hJes_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jes_CorrpT_refpt_refpartonB_ctbin, recoweight_corrpt);
		  
		  //JER
		  //w/o weight 		      
		  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jer_RawpT_refpt_refpartonB_ctbin);
		  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_noW->Fill(Jer_CorrpT_refpt_refpartonB_ctbin);
		  
		  //w/ weight
		  hJer_rawpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jer_RawpT_refpt_refpartonB_ctbin, recoweight_corrpt);
		  hJer_CorrpT_refpT_pTCut_refpartonB_ctbin_W->Fill(Jer_CorrpT_refpt_refpartonB_ctbin, recoweight_corrpt);
		} // if jes and jer
	    } // is MC
	  
	  TVector3 saved_jet_rc;
          saved_jet_rc.SetPtEtaPhi(jet_pt_corr, jt_eta, jt_phi);
          jet_vec_1D_rc.push_back(saved_jet_rc);
	  
	  jetw_Vec_1D_rc.push_back(recoweight_corrpt); 
          refpartonB_Vec_rc.push_back(jt_Flavour);
	  
	  //~~~~~~~~~~~~Define 1D trk vectors use to fill 2D trk vector
          std::vector<TVector3> oldtrk_vec_1D_rc;
          std::vector<TVector3> newtrk_vec_1D_rc;
          std::vector<TVector3> oldchtrk_vec_1D_rc;
          std::vector<TVector3> newchtrk_vec_1D_rc;
          std::vector<int> trk_charge_vec_1D_rc;

	  //start track loop within jet loop
	  for(int irctrk = 0; irctrk < ntrk; irctrk++)
	    {
	      float trk_pt = (float)trkpt[irctrk];
	      float trk_eta = (float)trketa[irctrk];
	      float trk_phi = (float)trkphi[irctrk];
	      bool trk_hp = (bool)highpur[irctrk];
	      int trk_chg = (int)trkcharge[irctrk];
	      int trk_nhits = (int)trknhits[irctrk];
	      float trk_chi2 = (float)trkchi2[irctrk];
	      int trk_ndf = (int)trkndof[irctrk];
	      int trk_nlayers = (int)trknlayer[irctrk];
	      float trk_pterr = (float)trkpterr[irctrk];
	      float trk_dxy = (float)trkdcaxy[irctrk];
	      float trk_dxyerr = (float)trkdcaxyerr[irctrk];
	      float trk_dz = (float)trkdcaz[irctrk];
	      float trk_dzerr = (float)trkdcazerr[irctrk];
	      int trk_algo = 1;
	      double trk_mva = 0.;
	      if(colliding_system=="PbPb")
		{
		  trk_algo = (int)trkalgo[irctrk];
		  trk_mva = (float)trkmva[irctrk];
		}

	      //std::cout<<"reco track: "<<trk_pt<<" "<<trk_eta<<" "<<trk_phi<<" "<<trk_hp<<" "<<trk_chg<<" "<<trk_nhits<<" "<<trk_algo<<" "<<trk_chi2<<" "<<trk_ndf<<" "<<trk_nlayers<<" "<<trk_pterr<<" "<<trk_dxy<<" "<<trk_dxyerr<<" "<<trk_dz<<" "<<trk_dzerr<<std::endl;
	      
	      double RecoTrack_Jet_DeltaEta = jt_eta - trk_eta;
	      double RecoTrack_Jet_DeltaPhi = TVector2::Phi_mpi_pi(jt_phi - trk_phi);
	      
	      double RecoTrack_Jet_DeltaR = TMath::Sqrt(pow(RecoTrack_Jet_DeltaEta,2) + pow(RecoTrack_Jet_DeltaPhi,2));
	      
	      if(RecoTrack_Jet_DeltaR >= jet_radius) continue; // condition for track to be inside the jet

	      if(colliding_system=="PbPb")
		{
		  if(!trk_hp) continue;
		  if(trk_chg == 0) continue;
		  if(trk_pt < trk_minpt_cut || trk_pt > trk_maxpt_cut) continue;
		  if(fabs(trk_eta) > trk_eta_cut) continue;
		  if(trk_nhits < nhits) continue;
		  if(trk_chi2/trk_ndf/trk_nlayers > chi2_ndf_nlayer_cut) continue;
		  if(trk_algo == 6 && trk_mva < 0.98) continue;
		  if(trk_pterr/trk_pt > trk_pt_resolution_cut) continue;
		  if(trk_dxy/trk_dxyerr > trk_dca_xy_cut) continue;
		  if(trk_dz/trk_dzerr > trk_dca_z_cut) continue;
		}
	      else if(colliding_system=="pp")
		{
		  if(!trk_hp) continue;
		  if(trk_chg == 0) continue;
		  if(trk_pt < trk_minpt_cut || trk_pt > trk_maxpt_cut) continue;
		  if(fabs(trk_eta) > trk_eta_cut) continue;
		  if(trk_pterr/trk_pt > trk_pt_resolution_cut) continue;
		  if(trk_dxy/trk_dxyerr > trk_dca_xy_cut) continue;
		  if(trk_dz/trk_dzerr > trk_dca_z_cut) continue;
		}
	      
	      float new_trk_pt = ptWRTJet(jet_pt_corr, jt_eta, jt_phi, trk_pt, trk_eta, trk_phi);
              float new_trk_eta = etaWRTJet(jet_pt_corr, jt_eta, jt_phi, trk_pt, trk_eta, trk_phi);
              float new_trk_phi = phiWRTJet(jet_pt_corr, jt_eta, jt_phi, trk_pt, trk_eta, trk_phi);

	      /*
	      float new_trk_pt = trk_pt;
	      float new_trk_eta = trk_eta;
	      float new_trk_phi = trk_phi;
	      */
	      
	      //push back to 1D trk vector
              trk_charge_vec_1D_rc.push_back(trk_chg);

	      TVector3 saved_oldtrk_rc, saved_newtrk_rc, saved_oldchtrk_rc, saved_newchtrk_rc;
              saved_oldtrk_rc.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
              oldtrk_vec_1D_rc.push_back(saved_oldtrk_rc);
              saved_newtrk_rc.SetPtEtaPhi(new_trk_pt, new_trk_eta, new_trk_phi);
              newtrk_vec_1D_rc.push_back(saved_newtrk_rc);

	      if(fabs((int)trk_chg) > 0)
                {
                  saved_oldchtrk_rc.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
                  oldchtrk_vec_1D_rc.push_back(saved_oldchtrk_rc);
                  saved_newchtrk_rc.SetPtEtaPhi(new_trk_pt, new_trk_eta, new_trk_phi);
                  newchtrk_vec_1D_rc.push_back(saved_newchtrk_rc);
                }
	      
	    } // track loop end

	  //push back to 2D trk vector
          oldtrk_vec_2D_rc.push_back(oldtrk_vec_1D_rc);
          newtrk_vec_2D_rc.push_back(newtrk_vec_1D_rc);
          oldchtrk_vec_2D_rc.push_back(oldchtrk_vec_1D_rc);
          newchtrk_vec_2D_rc.push_back(newchtrk_vec_1D_rc);
          trk_charge_2D_rc.push_back(trk_charge_vec_1D_rc);

	  //clear 1D trk vector
          oldtrk_vec_1D_rc.clear();
          newtrk_vec_1D_rc.clear();
          oldchtrk_vec_1D_rc.clear();
          newchtrk_vec_1D_rc.clear();
          trk_charge_vec_1D_rc.clear();
	  
	} // Jet loop end

      //push back to 2D jet vector                                                                                                       
      jet_vec_2D_rc.push_back(jet_vec_1D_rc);
      jetw_Vec_2D_rc.push_back(jetw_Vec_1D_rc);
      Matched_refpartonB_Vec_rc.push_back(refpartonB_Vec_rc); 
      
      //clear 1D jet vector
      jet_vec_1D_rc.clear();
      jetw_Vec_1D_rc.clear();
      refpartonB_Vec_rc.clear();
      
      //push back to 3D trk vector
      jet_oldtrk_vec_3D_rc.push_back(oldtrk_vec_2D_rc);
      jet_newtrk_vec_3D_rc.push_back(newtrk_vec_2D_rc);
      jet_oldchtrk_vec_3D_rc.push_back(oldchtrk_vec_2D_rc);
      jet_newchtrk_vec_3D_rc.push_back(newchtrk_vec_2D_rc);
      jet_trk_charge_3D_rc.push_back(trk_charge_2D_rc);

      //clear 2D trk vector                                                                                                              
      oldtrk_vec_2D_rc.clear();
      newtrk_vec_2D_rc.clear();
      oldchtrk_vec_2D_rc.clear();
      newchtrk_vec_2D_rc.clear();
      trk_charge_2D_rc.clear();
	
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Reco/Data jet loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop start:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // gen jet loop
      if(is_MC)
        {
	  double genEvtw = Evtw;
	  
	  double ld_genpT=-999., ld_genEta=-999., ld_genPhi=-999.; // leading Raw jet quantities                                     
	  double sld_genpT=-999., sld_genEta=-999, sld_genPhi=-999; // subleading Raw jet quantities          
	  double ld_genWTApT=-999., ld_genWTAEta=-999., ld_genWTAPhi=-999.; // leading Raw jet quantities                         
	  double sld_genWTApT=-999., sld_genWTAEta=-999, sld_genWTAPhi=-999; // subleading Raw jet quantities          
	  
	  for (int ign = 0; ign < (int)ngen; ign++) //gen Jet loop start                                        
            {
	      double gen_jet_pt = gen_jtpt[ign];
	      double gen_jet_eta = gen_jteta[ign];
	      double gen_jet_phi = gen_jtphi[ign];
	      //double gen_jet_eta = gen_WTAeta[ign];
	      //double gen_jet_phi = gen_WTAphi[ign];
	      double gen_jet_WTAeta = gen_WTAeta[ign];
	      double gen_jet_WTAphi = gen_WTAphi[ign];

	      //apply get jet pt eta cut
	      if(gen_jet_pt < jet_pt_min_cut || gen_jet_pt > jet_pt_max_cut) continue;
	      if(fabs(gen_jet_eta) >= jet_eta_max_cut) continue;
	      
	      find_leading_subleading_Jets(gen_jet_pt, gen_jet_eta, gen_jet_phi, ld_genpT, ld_genEta, ld_genPhi, sld_genpT, sld_genEta, sld_genPhi); // standard axis
	      
              find_leading_subleading_Jets(gen_jet_pt, gen_jet_WTAeta, gen_jet_WTAphi, ld_genWTApT, ld_genWTAEta, ld_genWTAPhi, sld_genWTApT, sld_genWTAEta, sld_genWTAPhi); // WTA axis
	      
	      int nrfgn = 0;
	      int index_ref = -1;
	      int index_gen = -1;

	      for (int irf = 0; irf < (int)nref; irf++) //ref Jet loop start to calculate the parton flavour
		{
		  if(fabs(gen_jtpt[ign] - refpt[irf]) < std::numeric_limits<float>::min()) // match gen jet with ref jet to find ref eta and phi
		    {
		      index_ref = irf;
		      index_gen = ign;
		      nrfgn++;
		    }
		}
	    	    
	      if(nrfgn > 1)
		{
		  std::cout<<"More than one gen Jet matched with ref Jet, Please check"<<std::endl;
		}

	      if(nrfgn == 1)
		{
		  if(index_ref >= 0 && index_gen >= 0)
		    {
		      int refpartonB_gn = -99;
		      if(fabs(refparton_flavorForB[index_ref]) >= 1 && fabs(refparton_flavorForB[index_ref]) <= 6)
			{
			  refpartonB_gn = fabs(refparton_flavorForB[index_ref]);
			}
		      else if(fabs(refparton_flavorForB[index_ref]) == 21)
			{
			  refpartonB_gn = 7;
			}
		      else
			{
			  refpartonB_gn = 0;
			}

		      if(refpartonB_gn == -99) continue;

		      double gen_jet_Matchedpt = gen_jtpt[index_gen];
		      double gen_jet_Matchedeta = gen_jteta[index_gen];
		      double gen_jet_Matchedphi = gen_jtphi[index_gen];
		      //double gen_jet_Matchedeta = gen_WTAeta[index_gen];
		      //double gen_jet_Matchedphi = gen_WTAphi[index_gen];
		      double gen_jet_MatchedWTAeta = gen_WTAeta[index_gen];
		      double gen_jet_MatchedWTAphi = gen_WTAphi[index_gen];

		      double Jet_GenpT_Eta_Phi_ctbin[4] = {gen_jet_Matchedpt, gen_jet_Matchedeta, gen_jet_Matchedphi, (double)ctbin};
		      double Jet_GenpT_WTAEta_WTAPhi_ctbin[4] = {gen_jet_Matchedpt, gen_jet_MatchedWTAeta, gen_jet_MatchedWTAphi, (double)ctbin};
		      
		      //w/o weight
		      hJet_GenpT_Eta_Phi_ctbin_pTCut_noW->Fill(Jet_GenpT_Eta_Phi_ctbin);
		      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_noW->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin);
		      
		      //w/ weight
		      hJet_GenpT_Eta_Phi_ctbin_pTCut_W->Fill(Jet_GenpT_Eta_Phi_ctbin, genEvtw);
		      hJet_GenpT_WTAEta_WTAPhi_ctbin_pTCut_W->Fill(Jet_GenpT_WTAEta_WTAPhi_ctbin, genEvtw);
		      
		      TVector3 saved_jet_gn;
		      saved_jet_gn.SetPtEtaPhi(gen_jet_Matchedpt, gen_jet_Matchedeta, gen_jet_Matchedphi);
		      jet_vec_1D_gn.push_back(saved_jet_gn);
		      
		      jetw_Vec_1D_gn.push_back(genEvtw); 
		      refpartonB_Vec_gn.push_back(refpartonB_gn); 
		      
		      //~~~~~~~~~~~~Define 1D trk vectors use to fill 2D trk vector
		      std::vector<TVector3> oldtrk_vec_1D_gn;
		      std::vector<TVector3> newtrk_vec_1D_gn;
		      std::vector<TVector3> oldchtrk_vec_1D_gn;
		      std::vector<TVector3> newchtrk_vec_1D_gn;
		      std::vector<int> trk_charge_vec_1D_gn;
		      
		      for(int igntrk = 0; igntrk < gen_trkpt->size(); igntrk++)
			{                                                                                                              
			  float gentrk_pt = float((*gen_trkpt)[igntrk]);                                                             
			  float gentrk_eta = float((*gen_trketa)[igntrk]);                                                             
			  float gentrk_phi = float((*gen_trkphi)[igntrk]);
			  int gentrk_chg = int((*gen_trkchg)[igntrk]);                                                             

		      
			  double GenTrack_Jet_DeltaEta = gen_jet_Matchedeta - gentrk_eta;                                              
			  double GenTrack_Jet_DeltaPhi = TVector2::Phi_mpi_pi(gen_jet_Matchedphi - gentrk_phi);                        
                          
			  double GenTrack_Jet_DeltaR = TMath::Sqrt(pow(GenTrack_Jet_DeltaEta,2) + pow(GenTrack_Jet_DeltaPhi,2));       
                          
			  if(GenTrack_Jet_DeltaR >= jet_radius) continue; // condition for track to be inside the jet

			  if(gentrk_chg == 0) continue;
			  if(gentrk_pt < trk_minpt_cut || gentrk_pt > trk_maxpt_cut) continue;
			  if(fabs(gentrk_eta) > trk_eta_cut) continue;

			  /*
			  float new_gentrk_pt = ptWRTJet(gen_jet_Matchedpt, gen_jet_Matchedeta, gen_jet_Matchedphi, gentrk_pt, gentrk_eta, gentrk_phi);
			  float new_gentrk_eta = etaWRTJet(gen_jet_Matchedpt, gen_jet_Matchedeta, gen_jet_Matchedphi, gentrk_pt, gentrk_eta, gentrk_phi);
			  float new_gentrk_phi = phiWRTJet(gen_jet_Matchedpt, gen_jet_Matchedeta, gen_jet_Matchedphi, gentrk_pt, gentrk_eta, gentrk_phi);
			  */

			  float new_gentrk_pt = gentrk_pt;
			  float new_gentrk_eta = gentrk_eta;
			  float new_gentrk_phi = gentrk_phi;
			  
			  //push back to 1D trk vector
			  trk_charge_vec_1D_gn.push_back(gentrk_chg);
	      
			  TVector3 saved_oldtrk_gn, saved_newtrk_gn, saved_oldchtrk_gn, saved_newchtrk_gn;
			  saved_oldtrk_gn.SetPtEtaPhi(gentrk_pt, gentrk_eta, gentrk_phi);
			  oldtrk_vec_1D_gn.push_back(saved_oldtrk_gn);
			  saved_newtrk_gn.SetPtEtaPhi(new_gentrk_pt, new_gentrk_eta, new_gentrk_phi);
			  newtrk_vec_1D_gn.push_back(saved_newtrk_gn);
			  
			  if(fabs((int)gentrk_chg) > 0)
			    {
			      saved_oldchtrk_gn.SetPtEtaPhi(gentrk_pt, gentrk_eta, gentrk_phi);
			      oldchtrk_vec_1D_gn.push_back(saved_oldchtrk_gn);
			      saved_newchtrk_gn.SetPtEtaPhi(new_gentrk_pt, new_gentrk_eta, new_gentrk_phi);
			      newchtrk_vec_1D_gn.push_back(saved_newchtrk_gn);
			    }
			} // gen trk loop end

		      //push back to 2D trk vector
		      oldtrk_vec_2D_gn.push_back(oldtrk_vec_1D_gn);
		      newtrk_vec_2D_gn.push_back(newtrk_vec_1D_gn);
		      oldchtrk_vec_2D_gn.push_back(oldchtrk_vec_1D_gn);
		      newchtrk_vec_2D_gn.push_back(newchtrk_vec_1D_gn);
		      trk_charge_2D_gn.push_back(trk_charge_vec_1D_gn);
		      
		      //clear 1D trk vector
		      oldtrk_vec_1D_gn.clear();
		      newtrk_vec_1D_gn.clear();
		      oldchtrk_vec_1D_gn.clear();
		      newchtrk_vec_1D_gn.clear();
		      trk_charge_vec_1D_gn.clear();
		      
		    } // idx_ref >=0 && idx_gen >=0
		} // if (nrfgn ==1)
	    } // gen jet loop

	  //push back to 2D jet vector
	  jet_vec_2D_gn.push_back(jet_vec_1D_gn);
	  jetw_Vec_2D_gn.push_back(jetw_Vec_1D_gn);
	  Matched_refpartonB_Vec_gn.push_back(refpartonB_Vec_gn);
	  
	  //clear 1D jet vector
	  jet_vec_1D_gn.clear();
	  jetw_Vec_1D_gn.clear();
	  refpartonB_Vec_gn.clear();

	  //push back to 3D trk vector
	  jet_oldtrk_vec_3D_gn.push_back(oldtrk_vec_2D_gn);
	  jet_newtrk_vec_3D_gn.push_back(newtrk_vec_2D_gn);
	  jet_oldchtrk_vec_3D_gn.push_back(oldchtrk_vec_2D_gn);
	  jet_newchtrk_vec_3D_gn.push_back(newchtrk_vec_2D_gn);
	  jet_trk_charge_3D_gn.push_back(trk_charge_2D_gn);
	  
	  //clear 2D trk vector
	  oldtrk_vec_2D_gn.clear();
	  newtrk_vec_2D_gn.clear();
	  oldchtrk_vec_2D_gn.clear();
	  newchtrk_vec_2D_gn.clear();
	  trk_charge_2D_gn.clear();
	  
	} // is MC condition
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:Gen jet loop end:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      // push back event quatities (pthatweight, vertexz, and centbin etc) before closing event loop 
      Evtw_vec_1D.push_back(Evtw);
      vertexz_vec_1D.push_back(vertexz);
      hiBin_vec_1D.push_back(centbin);
      
    } // event loop end

  //calling function for signal correlation
  if(is_MC)
    {
      std::cout<<endl;
      std::cout<<"~~~~~~~Signal correlation for Gen level~~~~~~~~~~"<<std::endl;
      // for gen
      signal_corr(Evtw_vec_1D, hiBin_vec_1D, jet_vec_2D_gn, Matched_refpartonB_Vec_gn, jetw_Vec_2D_gn, jet_oldtrk_vec_3D_gn, jet_newtrk_vec_3D_gn, jet_oldchtrk_vec_3D_gn, jet_newchtrk_vec_3D_gn, jet_trk_charge_3D_gn, false);
    }

  std::cout<<endl;
  std::cout<<"~~~~~~~Signal correlation for Reco/data level~~~~~~~~~~"<<std::endl;
  // for reco and data
  signal_corr(Evtw_vec_1D, hiBin_vec_1D, jet_vec_2D_rc, Matched_refpartonB_Vec_rc, jetw_Vec_2D_rc, jet_oldtrk_vec_3D_rc, jet_newtrk_vec_3D_rc, jet_oldchtrk_vec_3D_rc, jet_newchtrk_vec_3D_rc, jet_trk_charge_3D_rc, true);
  
  //calling function for mixing correlation
  if(is_MC)
    {
      std::cout<<endl;
      std::cout<<"~~~~~~~Mixing correlation for Gen level~~~~~~~~~~"<<std::endl;
      //for gen
      mixing_corr(Evtw_vec_1D, hiBin_vec_1D, jet_vec_2D_gn, Matched_refpartonB_Vec_gn, jetw_Vec_2D_gn, jet_oldtrk_vec_3D_gn, jet_newtrk_vec_3D_gn, jet_oldchtrk_vec_3D_gn, jet_newchtrk_vec_3D_gn, jet_trk_charge_3D_gn, vertexz_vec_1D, false, colliding_system);
    }

  std::cout<<endl;
  std::cout<<"~~~~~~~Mixing correlation for Reco/data level~~~~~~~~~~"<<std::endl;
  // for reco and data
  mixing_corr(Evtw_vec_1D, hiBin_vec_1D, jet_vec_2D_rc, Matched_refpartonB_Vec_rc, jetw_Vec_2D_rc, jet_oldtrk_vec_3D_rc, jet_newtrk_vec_3D_rc, jet_oldchtrk_vec_3D_rc, jet_newchtrk_vec_3D_rc, jet_trk_charge_3D_rc, vertexz_vec_1D, true, colliding_system);
    
  // clear event vector
  Evtw_vec_1D.clear();
  vertexz_vec_1D.clear();
  hiBin_vec_1D.clear();
  
  //clear jet vector
  jet_vec_2D_rc.clear();
  Matched_refpartonB_Vec_rc.clear();

  //clear trk vector
  jet_oldtrk_vec_3D_rc.clear();
  jet_newtrk_vec_3D_rc.clear();
  jet_oldchtrk_vec_3D_rc.clear();
  jet_newchtrk_vec_3D_rc.clear();
  jet_trk_charge_3D_rc.clear();

  //clear gen jet vector
  jet_vec_2D_gn.clear();
  Matched_refpartonB_Vec_gn.clear();
  
  //clear gen trk vector
  jet_oldtrk_vec_3D_gn.clear();
  jet_newtrk_vec_3D_gn.clear();
  jet_oldchtrk_vec_3D_gn.clear();
  jet_newchtrk_vec_3D_gn.clear();
  jet_trk_charge_3D_gn.clear();

  std::cout<<std::endl;
  std::cout<<"Total nref is: "<<count<<std::endl;

  delete hea_tree;
  delete hlt_tree;
  delete ski_tree;
  delete jet_tree;
  delete trk_tree;
  delete gen_tree;

  std::string outfilename = Form("%s/%s_Outfile_pTHat%1.1f_JetpT%1.1f_JetEta%1.1f_%d",out_file.Data(), colliding_system.Data(), pthat_cut, jet_pt_min_cut, jet_eta_max_cut, itxtoutFile);
  
  std::replace(outfilename.begin(), outfilename.end(), '.', 'p'); // replace . to p
  std::replace(outfilename.begin(), outfilename.end(), '-', 'N'); // replace - to N for negative

  TFile* fout = new TFile(Form("%s.root", outfilename.c_str()), "recreate");	 

  fout->mkdir("Event_Hist");
  fout->cd("Event_Hist");
  Write_Event_hist(is_MC);
  write_event_hist_bal(is_MC);
  
  fout->mkdir("Jet_QA_Hist");
  fout->cd("Jet_QA_Hist");
  Write_Jet_QA_hist(is_MC);
  write_jet_hist_bal(is_MC);
  
  if(is_MC)
    {
      fout->mkdir("JES_JER_Hist");
      fout->cd("JES_JER_Hist");
      Write_JES_JER_hist(is_JES_JER);
    }
  fout->mkdir("Track_Hist");
  fout->cd("Track_Hist");
  write_track_hist_bal(is_MC);

  fout->mkdir("Corr_Hist");
  fout->cd("Corr_Hist");
  //Write_Corr_hist(is_MC, is_Gen_Reco_Correlation);
  write_corr_hist_bal(is_MC);

  fout->mkdir("Corr_Hist_2D");
  fout->cd("Corr_Hist_2D");
  write_corr_hist_bal_2D(is_MC);
  
  fout->Write();
  fout->Close();

  sec_end = clock(); // stop time counting                                                                                                            
  cout << "========================================" << endl;
  cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
  cout << "========================================" << endl;

  print_stop(); // Print time, date and hour when it stops
  
} // void Tree_Analyzer() end

// main program
int main(int argc, char **argv)
{
  using namespace std;

  TString inputfile;
  int itxtout;
  TString outfile;
  TString coll_sys;
  int ismc;

  if(argc == 1)
    {
      std::cout<<"You did not pass any argument to the code other than the program name"<<std::endl;
      std::cout<<"You need to pass 6 arguments including the program name"<<std::endl;
    }

  if(argc >=1 && argc <=5)
    {
      std::cout<<"Only "<<argc<<" arguments you have given including the program name"<<std::endl;
      std::cout<<"You need to pass 6 arguments including the program name"<<std::endl;
    }

  if(argc == 6)
    {
      std::cout<<std::endl;
      std::cout<<"You have given "<< argc <<" arguments including the program name;  Your program will run"<<std::endl;
      std::cout<<std::endl;

      inputfile = argv[1];
      itxtout = atoi(argv[2]);
      outfile = argv[3];
      coll_sys = argv[4];
      ismc = atoi(argv[5]);

      Tree_Analyzer(inputfile, itxtout, outfile, coll_sys, ismc);
    }
  return 0;
}
