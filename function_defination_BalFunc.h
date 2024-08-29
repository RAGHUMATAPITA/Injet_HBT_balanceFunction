#include "histogram_definition_BalFunc.h"

double GetQInv(TVector3 trk1, TVector3 trk2)
{
  TLorentzVector trk_vec1(trk1, 0.139577);
  TLorentzVector trk_vec2(trk2, 0.139577);

  TLorentzVector QVec = trk_vec1 - trk_vec2;
  return fabs(QVec.Mag());
}

double GetkT(TVector3 trk1, TVector3 trk2)
{
  TLorentzVector trk_vec1(trk1, 0.139577);
  TLorentzVector trk_vec2(trk2, 0.139577);

  TLorentzVector kTVec = trk_vec1 + trk_vec2;
  return fabs(0.5*kTVec.Pt());
}

void signal_corr(std::vector<float> evtw_vec, std::vector<int> hiBin_vect, std::vector<std::vector<TVector3>> jet_vec, std::vector<std::vector<int>> refpartonB_vec, std::vector<std::vector<double>> jetw_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec, std::vector<std::vector<std::vector<int>>> jet_trk_charge, bool isrc)
{
  if(jet_newtrk_vec.size() != jet_vec.size() || jet_newtrk_vec.size() != evtw_vec.size())
    {
      std::cout<<"event numbers are not same, pleaes check"<<std::endl;
    }
  
  else std::cout<<"Total events is: "<<jet_newtrk_vec.size()<<"  "<<jet_vec.size()<<"  "<<evtw_vec.size()<<std::endl;
  
  std::cout<<endl;
  std::cout<<"~~~~~start signal correlation~~~~~~~~"<<std::endl;
  
  for(int ievt = 0; ievt < (int)jet_newtrk_vec.size(); ievt++) // event loop
    {
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for signal correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      
      double evtw = evtw_vec[ievt];
      int centbin = hiBin_vect[ievt];
      int ctbin = hCentBin_bal->FindBin((double)centbin) - 1;
      
      int njets = 0;
      
      for(int ijet = 0; ijet < (int)jet_newtrk_vec[ievt].size(); ijet++) // jet loop
	{
	  TVector3 jetVec = jet_vec[ievt][ijet];
	  double jet_pt = jetVec.Pt();
	  double jet_eta = jetVec.Eta();
	  double jet_phi = jetVec.Phi();
	  int ntrk_jet = jet_newtrk_vec[ievt][ijet].size();
	  int nchtrk_jet = jet_newchtrk_vec[ievt][ijet].size();
	  int refpartonB = refpartonB_vec[ievt][ijet]; 
	  double jetw = jetw_vec[ievt][ijet];

	  int nchbin  = hnch_bin_hist->FindBin(nchtrk_jet) - 1;
	  
	  //extract the jet pt bin                                                                                                        
	  int ptbin = hpt_bin_hist->FindBin(jet_pt) - 1;
	  
	  if(ptbin >= pt_bin)
	    {
	      std::cout<<"pT bin is exceeding: "<<jet_pt<<std::endl;
	      continue;
	    }

	  if(isrc)
	    {
	      hJet_pt->Fill(jet_pt, jetw);
	      hJet_ntrk->Fill(ntrk_jet, jetw);
	      hJet_nchtrk->Fill(nchtrk_jet, jetw);
	      hJet_nchtrk_pT->Fill(nchtrk_jet, jet_pt, jetw);
	      hJet_nchtrk_[nchbin]->Fill(nchtrk_jet, jetw);
	    }
	  else
	    {
	      hJet_Gen_pt->Fill(jet_pt, jetw);
              hJet_Gen_ntrk->Fill(ntrk_jet, jetw);
              hJet_Gen_nchtrk->Fill(nchtrk_jet, jetw);
              hJet_Gen_nchtrk_pT->Fill(nchtrk_jet, jet_pt, jetw);
              hJet_Gen_nchtrk_[nchbin]->Fill(nchtrk_jet, jetw);
	    }
	  for(int itrk = 0; itrk < (int)jet_newtrk_vec[ievt][ijet].size(); itrk++) // trigger trk loop
	    {
	      TVector3 old_trgtrk_vec = jet_oldtrk_vec[ievt][ijet][itrk];
	      double old_trgtrk_pt = old_trgtrk_vec.Pt();
	      double old_trgtrk_eta = old_trgtrk_vec.Eta();
	      double old_trgtrk_phi = old_trgtrk_vec.Phi();
	      TVector3 new_trgtrk_vec = jet_newtrk_vec[ievt][ijet][itrk];
	      double new_trgtrk_pt = new_trgtrk_vec.Pt();
	      double new_trgtrk_eta = new_trgtrk_vec.Eta();
	      double new_trgtrk_phi = new_trgtrk_vec.Phi();
	      int trgtrk_charge = jet_trk_charge[ievt][ijet][itrk];

	      if(isrc)
		{
		  halltrk_pt->Fill(old_trgtrk_pt);
		  halltrk_eta->Fill(old_trgtrk_eta);
		  halltrk_phi->Fill(old_trgtrk_phi);
		}
	      else
		{
		  halltrk_Gen_pt->Fill(old_trgtrk_pt);
		  halltrk_Gen_eta->Fill(old_trgtrk_eta);
		  halltrk_Gen_phi->Fill(old_trgtrk_phi);
		}
	      
	      if(TMath::Abs(trgtrk_charge) > 0 )
		{
		  if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max)
		    {
		      if(new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max)
			{
			  double newchtrk_Ntrig[4] = {1., (double)ptbin, (double)refpartonB, (double)ctbin};
			  
			  if(isrc)
			    {
			      hchtrk_pt->Fill(old_trgtrk_pt, jetw);
			      hchtrk_eta->Fill(old_trgtrk_eta, jetw);
			      hchtrk_phi->Fill(old_trgtrk_phi, jetw);
			      hchtrk_pt_jetaxis->Fill(new_trgtrk_pt, jetw);
			      hchtrk_eta_jetaxis->Fill(new_trgtrk_eta, jetw);
			      hchtrk_phi_jetaxis->Fill(new_trgtrk_phi, jetw);
			      
			      hnewchtrk_Ntrig->Fill(newchtrk_Ntrig, jetw);
			      hnewchtrk_Ntrig_1D[ctbin][ptbin]->Fill(1., jetw);
			      
			      if(trgtrk_charge > 0)
				{
				  hchtrk_pt_p->Fill(old_trgtrk_pt, jetw);
				  hchtrk_eta_p->Fill(old_trgtrk_eta, jetw);
				  hchtrk_phi_p->Fill(old_trgtrk_phi, jetw);
				  hchtrk_pt_jetaxis_p->Fill(new_trgtrk_pt, jetw);
				  hchtrk_eta_jetaxis_p->Fill(new_trgtrk_eta, jetw);
				  hchtrk_phi_jetaxis_p->Fill(new_trgtrk_phi, jetw);
				  
				  hnewchtrk_Ntrig_p->Fill(newchtrk_Ntrig, jetw);
				  hnewchtrk_Ntrig_p_1D[ctbin][ptbin]->Fill(1., jetw);
				}
			      else if(trgtrk_charge < 0)
				{
				  hchtrk_pt_m->Fill(old_trgtrk_pt, jetw);
				  hchtrk_eta_m->Fill(old_trgtrk_eta, jetw);
				  hchtrk_phi_m->Fill(old_trgtrk_phi, jetw);
				  hchtrk_pt_jetaxis_m->Fill(new_trgtrk_pt, jetw);
				  hchtrk_eta_jetaxis_m->Fill(new_trgtrk_eta, jetw);
				  hchtrk_phi_jetaxis_m->Fill(new_trgtrk_phi, jetw);
				  
				  hnewchtrk_Ntrig_m->Fill(newchtrk_Ntrig, jetw);
				  hnewchtrk_Ntrig_m_1D[ctbin][ptbin]->Fill(1., jetw);
				}
			    }
			  else
			    {
			      hchtrk_Gen_pt->Fill(old_trgtrk_pt, jetw);
			      hchtrk_Gen_eta->Fill(old_trgtrk_eta, jetw);
			      hchtrk_Gen_phi->Fill(old_trgtrk_phi, jetw);
			      hchtrk_Gen_pt_jetaxis->Fill(new_trgtrk_pt, jetw);
			      hchtrk_Gen_eta_jetaxis->Fill(new_trgtrk_eta, jetw);
			      hchtrk_Gen_phi_jetaxis->Fill(new_trgtrk_phi, jetw);
			      
			      hnewchtrk_Gen_Ntrig->Fill(newchtrk_Ntrig, jetw);
			      hnewchtrk_Gen_Ntrig_1D[ctbin][ptbin]->Fill(1., jetw);
			      
			      if(trgtrk_charge > 0)
				{
				  hchtrk_Gen_pt_p->Fill(old_trgtrk_pt, jetw);
				  hchtrk_Gen_eta_p->Fill(old_trgtrk_eta, jetw);
				  hchtrk_Gen_phi_p->Fill(old_trgtrk_phi, jetw);
				  hchtrk_Gen_pt_jetaxis_p->Fill(new_trgtrk_pt, jetw);
				  hchtrk_Gen_eta_jetaxis_p->Fill(new_trgtrk_eta, jetw);
				  hchtrk_Gen_phi_jetaxis_p->Fill(new_trgtrk_phi, jetw);
				  
				  hnewchtrk_Gen_Ntrig_p->Fill(newchtrk_Ntrig, jetw);
				  hnewchtrk_Gen_Ntrig_p_1D[ctbin][ptbin]->Fill(1., jetw);
				}
			      else if(trgtrk_charge < 0)
				{
				  hchtrk_Gen_pt_m->Fill(old_trgtrk_pt, jetw);
				  hchtrk_Gen_eta_m->Fill(old_trgtrk_eta, jetw);
				  hchtrk_Gen_phi_m->Fill(old_trgtrk_phi, jetw);
				  hchtrk_Gen_pt_jetaxis_m->Fill(new_trgtrk_pt, jetw);
				  hchtrk_Gen_eta_jetaxis_m->Fill(new_trgtrk_eta, jetw);
				  hchtrk_Gen_phi_jetaxis_m->Fill(new_trgtrk_phi, jetw);
				  
				  hnewchtrk_Gen_Ntrig_m->Fill(newchtrk_Ntrig, jetw);
				  hnewchtrk_Gen_Ntrig_m_1D[ctbin][ptbin]->Fill(1., jetw);
				}
			    }
			} // trk pt condition
		    } // trk eta condition
		} // charge condition
	      
	      for(int jtrk = 0; jtrk < (int)jet_newtrk_vec[ievt][ijet].size(); jtrk++) // associate trk loop
		{
		  if(itrk == jtrk) continue; // to avoid auto correlation
		  
		  TVector3 old_assotrk_vec = jet_oldtrk_vec[ievt][ijet][jtrk];
		  double old_assotrk_pt = old_assotrk_vec.Pt();
		  double old_assotrk_eta = old_assotrk_vec.Eta();
		  double old_assotrk_phi = old_assotrk_vec.Phi();
		  TVector3 new_assotrk_vec = jet_newtrk_vec[ievt][ijet][jtrk];
		  double new_assotrk_pt = new_assotrk_vec.Pt();
		  double new_assotrk_eta = new_assotrk_vec.Eta();
		  double new_assotrk_phi = new_assotrk_vec.Phi();
		  int assotrk_charge = jet_trk_charge[ievt][ijet][jtrk];

		  double deta_jetaxis = new_trgtrk_eta - new_assotrk_eta;
     		  double dphi_jetaxis = new_trgtrk_phi - new_assotrk_phi;
   		  double dphi2_jetaxis = new_assotrk_phi - new_trgtrk_phi;

		  //double qinv = GetQInv(old_trgtrk_vec, old_assotrk_vec);
		  //double kt = GetkT(old_trgtrk_vec, old_assotrk_vec);

		  double qinv = GetQInv(new_trgtrk_vec, new_assotrk_vec);
		  double kt = GetkT(new_trgtrk_vec, new_assotrk_vec);
		  
		  if(dphi_jetaxis > 1.5*TMath::Pi())
		    {
		      dphi_jetaxis = dphi_jetaxis - 2.0*TMath::Pi();
		    }
		  else if(dphi_jetaxis < -0.5*TMath::Pi())
		    {
		      dphi_jetaxis = dphi_jetaxis + 2.0*TMath::Pi();
		    }
		  
		  if(dphi2_jetaxis > 1.5*TMath::Pi())
		    {
		      dphi2_jetaxis = dphi2_jetaxis - 2.0*TMath::Pi();
		    }
		  else if(dphi2_jetaxis < -0.5*TMath::Pi())
		    {
		      dphi2_jetaxis = dphi2_jetaxis + 2.0*TMath::Pi();
		    }

		  if(TMath::Abs(trgtrk_charge) > 0 && TMath::Abs(assotrk_charge) > 0)
		    {
		      if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max && TMath::Abs(new_assotrk_eta) <= newtrk_eta_max)
			{
			  if((new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max) && (new_assotrk_pt >= newtrk_jt03_min && new_assotrk_pt <= newtrk_jt03_max))
			    {
			      //std::cout<<old_assotrk_pt<<"  "<<new_assotrk_pt<<"  "<<old_assotrk_eta<<"  "<<new_assotrk_eta<<"  "<<old_assotrk_phi<<"  "<<new_assotrk_phi<<std::endl;
			      
			      double signal_axis_newchtrk_jt03_1[5] = {TMath::Abs(deta_jetaxis), dphi_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
			      double signal_axis_newchtrk_jt03_2[5] = {-TMath::Abs(deta_jetaxis), dphi_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
			      double signal_axis_newchtrk_jt03_3[5] = {TMath::Abs(deta_jetaxis), dphi2_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
			      double signal_axis_newchtrk_jt03_4[5] = {-TMath::Abs(deta_jetaxis), dphi2_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
			      // for qinv
			      double signal_qinv[5] = {qinv, kt, (double)nchtrk_jet, (double)refpartonB, (double)ctbin};
			      
			      if(isrc)
				{
				  hsignal_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);		     
				  hsignal_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				  hsignal_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				  hsignal_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				  hsignal_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				  hsignal_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				  hsignal_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				  hsignal_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

				  // for qinv
				  if(trgtrk_charge == assotrk_charge)
				    {
				      hsignal_qinv_ls->Fill(signal_qinv, jetw);
				    }
				  else hsignal_qinv_us->Fill(signal_qinv, jetw);
				  
				  if(trgtrk_charge > 0 && assotrk_charge > 0)
				    {
				      hsignal_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);
				      
				      hsignal_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				      hsignal_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				      hsignal_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				      hsignal_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				  if(trgtrk_charge < 0 && assotrk_charge < 0)
				    {
				      
				      hsignal_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				  if(trgtrk_charge > 0 && assotrk_charge < 0)
				    {
				      hsignal_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				  if(trgtrk_charge < 0 && assotrk_charge > 0)
				    {
				      hsignal_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				}
			      else
				{
				  hsignal_Gen_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);		     
				  hsignal_Gen_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				  hsignal_Gen_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				  hsignal_Gen_newchtrk_jt03->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				  hsignal_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				  hsignal_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
				  hsignal_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				  hsignal_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

				  // for qinv
                                  if(trgtrk_charge == assotrk_charge)
                                    {
                                      hsignal_gen_qinv_ls->Fill(signal_qinv, jetw);
                                    }
                                  else hsignal_gen_qinv_us->Fill(signal_qinv, jetw);

				  if(trgtrk_charge > 0 && assotrk_charge > 0)
				    {
				      hsignal_Gen_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pp->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				      
				    }
				  if(trgtrk_charge < 0 && assotrk_charge < 0)
				    {
				      
				      hsignal_Gen_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mm->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				  if(trgtrk_charge > 0 && assotrk_charge < 0)
				    {
				      hsignal_Gen_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_pm->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				  if(trgtrk_charge < 0 && assotrk_charge > 0)
				    {
				      hsignal_Gen_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_1, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_2, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_3, jetw/4.);
				      hsignal_Gen_newchtrk_jt03_mp->Fill(signal_axis_newchtrk_jt03_4, jetw/4.);

				      hsignal_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                      hsignal_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
				    }
				}
			    } // trk pt condition
			} // trk eta condition
		    } // trk charge condition
		} // asso trk loop
	    } // trg trk loop end
	  njets++;
	} // jet loop end
      if(njets > 0)
	{
	  if(isrc) hnjets_afterCut->Fill(njets);
	  else hnjets_Gen_afterCut->Fill(njets);
	}
    } // event loop end
  std::cout<<endl;
  std::cout<<"~~~~~end signal correlation~~~~~~~~"<<std::endl;
} // function loop

void mixing_corr(std::vector<float> evtw_vec, std::vector<int> hiBin_vect, std::vector<std::vector<TVector3>> jet_vec, std::vector<std::vector<int>> refpartonB_vec, std::vector<std::vector<double>> jetw_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_oldchtrk_vec, std::vector<std::vector<std::vector<TVector3>>> jet_newchtrk_vec, std::vector<std::vector<std::vector<int>>> jet_trk_charge, std::vector<float> vertexz, bool isrc, TString colliding_system)
{
  std::cout<<endl;
  std::cout<<"~~~~~start mixed event correlation~~~~~~~~"<<std::endl;
  
  for(int ievt = 0; ievt < (int)jet_newtrk_vec.size(); ievt++) // 1st event loop
    {
      if(ievt%10000 == 0) std::cout<<ievt<<"  events running for mixing correlation of total events "<< jet_newtrk_vec.size() <<std::endl;
      
      double evtw = evtw_vec[ievt];
      int icentbin = hiBin_vect[ievt];
      int ctbin = hCentBin_bal->FindBin((double)icentbin) - 1;
      float ivertexz = vertexz[ievt];
      
      // mixing algorithm
      int mixstart = ievt+1;
      int mixend   = (int)jet_newtrk_vec.size();
      
      if(mixstart > (0.7*(jet_newtrk_vec.size())))
	{
	  int mixstart = (0.5*(jet_newtrk_vec.size()));
	  int mixend   = (int)jet_newtrk_vec.size();
	}
      
       int nmix = 0;
       /*
       int mixstart = ievt - bkgFactor/2;
       int mixend   = ievt + bkgFactor/2 + 1;

      if(ievt < bkgFactor/2)
	{
	  mixstart = 0;
	  mixend   = bkgFactor + 1;
	}
      else if(ievt > (int)jet_newtrk_vec.size() - bkgFactor/2 - 1)
	{
	  mixstart = (int)jet_newtrk_vec.size() - bkgFactor - 1;
	  mixend   = (int)jet_newtrk_vec.size();
	}
      if( mixend > (int)jet_newtrk_vec.size() )
	mixend = (int)jet_newtrk_vec.size();
       */
       
      for(int jevt = mixstart; jevt < mixend; jevt++) // 2nd event loop
	{
	  if(ievt == jevt) continue;
	  int jcentbin = hiBin_vect[jevt];
	  float jvertexz = vertexz[jevt];
	  
	  if (fabs(jvertexz - ivertexz) > 2.0) continue; // vertex Z condition 2 cm
	  if(colliding_system == "PbPb")
	    {
	      if(ctbin < 2)
		{
		  if(fabs(icentbin - jcentbin) > 10.) continue; // Centrality condition 5%
		}
	      else
		{
		  if(fabs(icentbin - jcentbin) > 20.) continue;
		}
	    }

	  nmix++;
	  if(nmix > bkgFactor) break;
	  
	  for(int ijet = 0; ijet < (int)jet_newtrk_vec[ievt].size(); ijet++) // 1st jet loop
	    {
	      TVector3 jetVec = jet_vec[ievt][ijet];
	      double jet_pt = jetVec.Pt();
	      double jet_eta = jetVec.Eta();
	      double jet_phi = jetVec.Phi();
	      int nchtrk_jet = jet_newchtrk_vec[ievt][ijet].size();
	      int ntrk_jet = jet_newtrk_vec[ievt][ijet].size();
	      int refpartonB = refpartonB_vec[ievt][ijet]; 
	      double jetw = jetw_vec[ievt][ijet];
	      
	      int nchbin  = hnch_bin_hist->FindBin(nchtrk_jet) - 1;
	      
	      //extract the jet pt bin
	      int ptbin = hpt_bin_hist->FindBin(jet_pt) - 1;
	      
	      if(ptbin >= pt_bin)
		{
		  std::cout<<"pT bin is exceeding in Mixing: "<<jet_pt<<std::endl;
		  continue;
		}
	      
	      for(int jjet = 0; jjet < (int)jet_newtrk_vec[jevt].size(); jjet++) // 2nd jet loop
		{
		  TVector3 jetVec2 = jet_vec[jevt][jjet];
		  double jet_pt2 = jetVec.Pt();
		  double jet_eta2 = jetVec.Eta();
		  double jet_phi2 = jetVec.Phi();
		  int nchtrk_jet2 = jet_newchtrk_vec[jevt][jjet].size();
		  int ntrk_jet2 = jet_newtrk_vec[jevt][jjet].size();
		  int refpartonB2 = refpartonB_vec[jevt][jjet];
		  double jetw2 = jetw_vec[jevt][jjet];

		  double Dr_jet = jetVec.DeltaR(jetVec2);
		  double Dpt_jet = fabs(jet_pt - jet_pt2)/jet_pt;
		  
		  //if(Dr_jet >= 0.4 || Dpt_jet >= 0.1) continue;
		  
		  if(Dr_jet >= 1.0) continue; 

		  for(int itrk = 0; itrk < (int)jet_newtrk_vec[ievt][ijet].size(); itrk++) // trigger trk loop
		    {
		      TVector3 old_trgtrk_vec = jet_oldtrk_vec[ievt][ijet][itrk];
		      double old_trgtrk_pt = old_trgtrk_vec.Pt();
		      double old_trgtrk_eta = old_trgtrk_vec.Eta();
		      double old_trgtrk_phi = old_trgtrk_vec.Phi();
		      TVector3 new_trgtrk_vec = jet_newtrk_vec[ievt][ijet][itrk];
		      double new_trgtrk_pt = new_trgtrk_vec.Pt();
		      double new_trgtrk_eta = new_trgtrk_vec.Eta();
		      double new_trgtrk_phi = new_trgtrk_vec.Phi();
		      int trgtrk_charge = jet_trk_charge[ievt][ijet][itrk];
		      		      
		      for(int jtrk = 0; jtrk < (int)jet_newtrk_vec[jevt][jjet].size(); jtrk++) // associate trk loop
			{
			  TVector3 old_assotrk_vec = jet_oldtrk_vec[jevt][jjet][jtrk];
			  double old_assotrk_pt = old_assotrk_vec.Pt();
			  double old_assotrk_eta = old_assotrk_vec.Eta();
			  double old_assotrk_phi = old_assotrk_vec.Phi();
			  TVector3 new_assotrk_vec = jet_newtrk_vec[jevt][jjet][jtrk];
			  double new_assotrk_pt = new_assotrk_vec.Pt();
			  double new_assotrk_eta = new_assotrk_vec.Eta();
			  double new_assotrk_phi = new_assotrk_vec.Phi();
			  int assotrk_charge = jet_trk_charge[jevt][jjet][jtrk];
		
			  double deta_jetaxis = new_trgtrk_eta - new_assotrk_eta;
			  double dphi_jetaxis = new_trgtrk_phi - new_assotrk_phi;
			  double dphi2_jetaxis = new_assotrk_phi - new_trgtrk_phi;
			  
			  //double qinv = GetQInv(old_trgtrk_vec, old_assotrk_vec);
			  //double kt = GetkT(old_trgtrk_vec, old_assotrk_vec);

			  double qinv = GetQInv(new_trgtrk_vec, new_assotrk_vec);
			  double kt = GetkT(new_trgtrk_vec, new_assotrk_vec);
			  
			  if(dphi_jetaxis > 1.5*TMath::Pi())
			    {
			      dphi_jetaxis = dphi_jetaxis - 2.0*TMath::Pi();
			    }
			  else if(dphi_jetaxis < -0.5*TMath::Pi())
			    {
			      dphi_jetaxis = dphi_jetaxis + 2.0*TMath::Pi();
			    }
			  
			  if(dphi2_jetaxis > 1.5*TMath::Pi())
			    {
			      dphi2_jetaxis = dphi2_jetaxis - 2.0*TMath::Pi();
			    }
			  else if(dphi2_jetaxis < -0.5*TMath::Pi())
			    {
			      dphi2_jetaxis = dphi2_jetaxis + 2.0*TMath::Pi();
			    }
			  
			  if(TMath::Abs(trgtrk_charge) > 0 && TMath::Abs(assotrk_charge) > 0)
			    {
			      if(TMath::Abs(new_trgtrk_eta) <= newtrk_eta_max && TMath::Abs(new_assotrk_eta) <= newtrk_eta_max)
				{
				  if((new_trgtrk_pt >= newtrk_jt03_min && new_trgtrk_pt <= newtrk_jt03_max) && (new_assotrk_pt >= newtrk_jt03_min && new_assotrk_pt <= newtrk_jt03_max))
				    {
				      
                                      double mixing_axis_newchtrk_jt03_1[5] = {TMath::Abs(deta_jetaxis), dphi_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
                                      double mixing_axis_newchtrk_jt03_2[5] = {-TMath::Abs(deta_jetaxis), dphi_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
				      double mixing_axis_newchtrk_jt03_3[5] = {TMath::Abs(deta_jetaxis), dphi2_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};
				      double mixing_axis_newchtrk_jt03_4[5] = {-TMath::Abs(deta_jetaxis), dphi2_jetaxis, (double)ptbin, (double)refpartonB, (double)ctbin};

				      // for qinv
				      double mixing_qinv[5] = {qinv, kt, (double)nchtrk_jet, (double)refpartonB, (double)ctbin};
			      
				      if(isrc)
					{
					  hmixing_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					  hmixing_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					  hmixing_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					  hmixing_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					  hmixing_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
					  hmixing_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
					  hmixing_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					  hmixing_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					  //for qinv
					  if(trgtrk_charge == assotrk_charge)
					    {
					      hmixing_qinv_ls->Fill(mixing_qinv, jetw);
					    }
					  else hmixing_qinv_us->Fill(mixing_qinv, jetw);
					  
					  if(trgtrk_charge > 0 && assotrk_charge > 0)
					    {
					      hmixing_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);
					      
					      hmixing_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
					      hmixing_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
					      hmixing_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					      hmixing_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					    }
					  if(trgtrk_charge < 0 && assotrk_charge < 0)
					    {
					      hmixing_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					    }
					  if(trgtrk_charge > 0 && assotrk_charge < 0)
					    {
					      hmixing_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					    }
					  if(trgtrk_charge < 0 && assotrk_charge > 0)
					    {
					      hmixing_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					    }
					}
				      else
					{
					  hmixing_Gen_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					  hmixing_Gen_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					  hmixing_Gen_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					  hmixing_Gen_newchtrk_jt03->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					  hmixing_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                          hmixing_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                          hmixing_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                          hmixing_Gen_newchtrk_jt03_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);

					  //for qinv
					  if(trgtrk_charge == assotrk_charge)
					    {
					      hmixing_gen_qinv_ls->Fill(mixing_qinv, jetw);
					    }
					  else hmixing_gen_qinv_us->Fill(mixing_qinv, jetw);
					  
					  if(trgtrk_charge > 0 && assotrk_charge > 0)
					    {
					      hmixing_Gen_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pp->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					    }
					  if(trgtrk_charge < 0 && assotrk_charge < 0)
					    {
					      hmixing_Gen_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mm->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					    }
					  if(trgtrk_charge > 0 && assotrk_charge < 0)
					    {
					      hmixing_Gen_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_pm->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_pm_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					    }
					  if(trgtrk_charge < 0 && assotrk_charge > 0)
					    {
					      hmixing_Gen_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_1, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_2, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_3, jetw/4.);
					      hmixing_Gen_newchtrk_jt03_mp->Fill(mixing_axis_newchtrk_jt03_4, jetw/4.);

					      hmixing_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
                                              hmixing_Gen_newchtrk_jt03_mp_2D[ctbin][ptbin]->Fill(-TMath::Abs(deta_jetaxis), dphi2_jetaxis, jetw/4.);
					    }
					}
				    } // trk pt condition
				} // trk eta condition
			    } // trk charge condition
			} // asso trk loop
		    } // trg trk loop end
		} // 1st jet loop end
	    } // 2nd jet loop end
	} // 2nd event loop end
    } // 1st event loop end
  std::cout<<endl;
  std::cout<<"~~~~~end mixed event correlation~~~~~~~~"<<std::endl;
} // function loop
      
void print_start()
{
  cout << endl;
  time_t init = time(0);
  char* init_time = ctime(&init); // convert now to string form                                        
  cout << "Starting at : " << init_time << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
  cout << endl;
}

void print_stop()
{
  time_t end = time(0);
  char* end_time = ctime(&end); // convert now to string form                                                   
  cout << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
  cout << "Stopping at : " << end_time << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << endl;
}
