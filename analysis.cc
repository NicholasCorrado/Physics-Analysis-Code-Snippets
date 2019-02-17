#include <stdio.h>
#include <stdlib.h>
#include <iostream>


#include <bitset>

#include <algorithm>

#include <set>

#include "analysis.h"
#include "userinfo.h"
#include "helix/Helix.h"
#include "ddk_skim/findKs.h"
#include "particle/combination.h"
#include "particle/utility.h"
#include "particle/FittedMomentum.h"

#include "tables/brecon.h"
#include "tables/fullrecon.h"
#include "tables/evtcls.h"
#include "fullrecon/frec_util.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "belle.h"
#include "benergy/BeamEnergy.h" 

#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/kmakemother.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h" 

#include BELLETDF_H
#include HEPEVT_H
#include MDST_H
#include MDST_OBS_H
#include EVTCLS_H
#include EVTVTX_H
#include LEVEL4_H
#include TRG_H

#include "trg/TrgBit.h" 

#include "kid/atc_pid.h"
#include "eid/eid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "ip/IpProfile.h"
#include "hamlet/Hamlet.h"
#include "hamlet/Fbtag_MultDimLikelihood0.h"

#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "particle/gammac.h"

#if defined(BELLE_NAMESPACE)
namespace Belle
{
#endif

  //***********************************************************
  void Analysis::begin_run(BelleEvent*, int*) {
    //cout << "begin_run" << endl;

    m_event_run = 0;
    m_run_seq++;

    // Determine if this is data or MC
    Belle_runhead_Manager &rhd_mgr = Belle_runhead_Manager::get_manager();
    std::vector<Belle_runhead>::const_iterator rhd = rhd_mgr.begin();
    if( rhd == rhd_mgr.end() ) {
      std::cout << "Cannot get Belle_runhead!" << std::endl;
      isMC = false;
    } else {
      isMC = ( rhd->ExpMC() == 2 ? true : false );
    }

    BeamEnergy::begin_run();// initialize BeamEnergy class
    m_eBeamCM = BeamEnergy::Ecm();
    m_eBeam = BeamEnergy::E_HER() + BeamEnergy::E_LER();
    if (m_eBeam == 0.0) m_eBeam = BeamEnergy::E_HER_orig() + BeamEnergy::E_LER_orig();
    if (m_eBeam == 0.0) cout << "eBeam = 0.0 -- expect problems!" << endl;
    //m_pBeam   = BeamEnergy::p_beam();
    // m_pBeamM = BeamEnergy::p_cm(m_pBeam);
    m_cmBoost = - BeamEnergy::CMBoost();
    
    m_pExp = BeamEnergy::p_beam();
    m_pExpCM = BeamEnergy::p_cm(m_pBeam);
    // do NOT use! June 7, 2017  VS&NC m_pExp = BeamEnergy::p_beam2();

    //    cout << "-VS-D-m_pBeam, m_eBeamCM = " << m_pBeam << ", " << m_eBeamCM << endl;
    //    cout << "-VS-D-isMC = " << isMC << endl;

    eid::init_data();// initialize electron ID (done in Hamlet as well)
    
    IpProfile::begin_run();// get run-by-run IP data from DB
    IpProfile::dump();
    m_ip = IpProfile::e_position();
  }



  void Analysis::event(BelleEvent* evptr, int* status) {

    //cout << "event" << endl;
    *status = 0; 
    
    m_event_seq++; // total number of processed events
    m_event_run++; // events in the current run

    // as advertised previously! 
    //    if (! IpProfile::usable()) {
    //      return;
    //    }

    Belle_event_Manager& EvtMgr = Belle_event_Manager::get_manager();
    /*
    Gen_hepevt_Manager &GenMgr = Gen_hepevt_Manager::get_manager();
    Mdst_charged_Manager& ChargedMgr = Mdst_charged_Manager::get_manager();
    Mdst_pi0_Manager& Pi0Mgr = Mdst_pi0_Manager::get_manager();
    Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();
    Mdst_vee2_Manager& Vee2Mgr = Mdst_vee2_Manager::get_manager(); 
    Evtcls_gamma_cluster_Manager& Gam_info_Manager =  Evtcls_gamma_cluster_Manager::get_manager();
    */
    Evtcls_hadron_info_Manager &HadInfoMgr = Evtcls_hadron_info_Manager::get_manager();
    Evtcls_hadronic_flag_Manager& HadFlagMgr = Evtcls_hadronic_flag_Manager::get_manager();
    Evtcls_flag_Manager& FlagMgr = Evtcls_flag_Manager::get_manager();

    m_exp   = EvtMgr[0].ExpNo();
    m_run   = EvtMgr[0].RunNo();
    m_event = EvtMgr[0].EvtNo()&0x0FFFFFFF;

    // Get HadInfo
    if ( &HadInfoMgr ) {
      if (HadInfoMgr.count()) 
	{
	  m_evt_Ntrk         = HadInfoMgr[0].Ntrk();
	  m_evt_Ncls         = HadInfoMgr[0].Ncls();
	  m_evt_Psum         = HadInfoMgr[0].Psum();
	  m_evt_Esum         = HadInfoMgr[0].Esum();
	  m_evt_Evis         = HadInfoMgr[0].Evis();
	  m_evt_Pz           = HadInfoMgr[0].Pz();
	  m_evt_HeavyJetMass = HadInfoMgr[0].HeavyJetMass();
	  m_evt_Thrust       = HadInfoMgr[0].Thrust();
	  m_evt_R2           = HadInfoMgr[0].R2();
	}
    }
    
    m_HadronB_flag = -1;
    m_HadronA_flag = -1;
    
    if ( &HadFlagMgr ) {
      if (HadFlagMgr.count()>0) {
	m_HadronB_flag = HadFlagMgr[0].hadronic_flag(2);
	m_HadronA_flag = HadFlagMgr[0].hadronic_flag(0);
      }
    }	
    
    m_Tau_flag = -1;
    if ( &FlagMgr ) {
      if (FlagMgr.count()>0) {
	m_Tau_flag = FlagMgr[0].flag(4);
      }
    }

    //--------------------------------------------------------------------

    // trigger (tsim)

    unsigned trg_bit1_in = 0;
    unsigned trg_bit2_in = 0;

    unsigned trg_bit1_out = 0;
    unsigned trg_bit2_out = 0;

    int trg_type = -1;
    int trg_timing_src = -1;

    // http://totoro.phyast.pitt.edu/twiki/bin/view/PittBelle/BelleTriggerSimulationPittLocal

    // for trigger bit definitions see /belle/belle/b20090127_0910.p41/src/cal/cal-trg/basf_if/TrgTool.icc

    // also see /belle/belle/b20090127_0910.p41/src/cal/cal-trg/trg/TrgBit.h
    // and /belle/belle/b20090127_0910.p41/src/cal/cal-trg/trg/TrgBitDef.h

    TrgBit TrgBits;
    int trg_dimu = -1;
    int trg_dimu_noz = -1;

    int expNo, runNo, evtNo;
    bool McFlag;
    getEventInfo(expNo,runNo,evtNo,McFlag);
    
    if (McFlag == true) {
      if(expNo>30){
	Rectrg_summary3_Manager& TrgMgr = Rectrg_summary3_Manager::get_manager();
	vector<Rectrg_summary3>::iterator itft = TrgMgr.begin();
	int count = TrgMgr.count();
	if(count == 0 ){
	  cout << "VS-W-Warning: trg count == 0" << endl;
	}else{
	  trg_type = TrgBits.getType();
	  trg_timing_src = TrgBits.getTimingSource();
	  trg_bit1_in = (*itft).input(0);
	  trg_bit2_in = (*itft).input(1);
	  trg_bit1_out = (*itft).final(0);
	  trg_bit2_out = (*itft).final(1);
	  trg_dimu = TrgBits.getPSNM(TrgBit::dimu) ;
	  trg_dimu_noz = TrgBits.getPSNM(TrgBit::dimu_noz) ;
	}
      }else if(expNo<30){
	Rectrg_summary_Manager& TrgMgr = Rectrg_summary_Manager::get_manager();
	vector<Rectrg_summary>::iterator itft = TrgMgr.begin();
	int count = TrgMgr.count();
	if(count == 0 ){
	  cout << "VS-W-Warning: trg count == 0" << endl;
	}else{
	  trg_type = TrgBits.getType();
	  trg_timing_src = TrgBits.getTimingSource();
	  trg_bit1_in = (*itft).input(0);
	  trg_bit2_in = (*itft).input(1);
	  trg_bit1_out = (*itft).final(0);
	  trg_bit2_out = (*itft).final(1);
	  trg_dimu = TrgBits.getPSNM(TrgBit::dimu) ;
	  trg_dimu_noz = TrgBits.getPSNM(TrgBit::dimu_noz) ;
	}
      }
    }

    unsigned trigger = trg_bit1_out + trg_bit2_out;

    /*
    if ( trigger == 0 ) {
      cout << "VS-D-Event would have failed to trigger!" << endl;
    }
    else {
      cout << "VS-D-Event would have triggered. trigger = " << trigger << endl;
    }
    */

    //    cout << "VS-D-Trigger info: trg_type = " << trg_type << endl;
    //    cout << "VS-D-Trigger info: trg_timing_src = 0x" << std::hex << trg_timing_src << std::dec << endl;
    // cout << "VS-D-Trigger info: trg_dimu, trg_dimu_noz = " << trg_dimu << ", " << trg_dimu_noz << endl;
    //    cout << "VS-D-Trigger info: trg_bit1_in,  trg_bit2_in  = " << bitset<32>(trg_bit1_in) << ", " << bitset<32>(trg_bit2_in) << endl;
    //    cout << "VS-D-Trigger info: trg_bit1_out, trg_bit2_out = " << bitset<32>(trg_bit1_out) << ", " << bitset<32>(trg_bit2_out) << endl;

    //    binary2char(trg_bit1_out);
    //    binary2char(trg_bit2_out);

    m_trigbits = trigger;

    m_istrig = 0;
    if ( m_trigbits != 0 ) m_istrig = 1;

    //--------------------------------------------------------------------

    clearValuesMC();
    cout << "BEFORE: isr = " << m_isr << endl;
    analyzeHepEvtRecordSignal(); // MC Truth analysis
    cout << "AFTER: isr = " << m_isr << endl;
    reweight(m_mc_gam_isr_boost_e);

    //if (mc_reweight == 1) {
      analyzeBackgroundSignal();   // Minimal MC truth analysis for background events. Looking only at Upsilon(nS) transitions.
      
      selectChargedPions();  // identify charged pion candidates
      selectChargedMuons();  // identify muon candidates
      
      countChargedTracks();  // count high-quality charged tracks identified as pi+, pi-, mu+ and/or mu-
      
      selectPhotons();       // identify photon candidates
      selectRhos();          // identify rho0->pi+pi- candidates
      selectYnS();           // identify Upsilon(nS)->mu+mu- candidates
      selectWbj();           // identify Wbj->rho0 upsilon(nS)->pi+pi- mu+mu- candidates
      selectY5S();           // identify upsilon(5S)->Wbj gamma->rho0 upsilon(nS)->pi+pi- mu+mu- candidates
      
      prepareCandidates();
      //}
      //cout << "testing2... " << m_mc_fsr_YnS_e << endl;
      
    std::sort(candidatesRec.begin(), candidatesRec.end(), bestEnergyY5S());
    //std::sort(candidatesRec.begin(), candidatesRec.end(), bestGamWidth());


    m_mc  = ( isMC ) ? 1 : 0;
    m_mcs = ( isSignalMC ) ? 1 : 0;

    //reweight(m_mc_gam_isr_boost_e);
    //cout << "m_mc = " << m_mc << endl;
    //cout << "m_msc = " << m_mcs << endl;

    /*
    cout << "VS-D-number of charged tracks for B tagging (pion) = " << pic_list.size() << endl;
    cout << "VS-D-number of charged tracks for B tagging (muon) = " << muc_list.size() << endl;
    cout << "VS-D-number of high-quality positively charged pions = " << pip_list.size() << endl;
    cout << "VS-D-number of high-quality negatively charged pions = " << pim_list.size() << endl;
    cout << "VS-D-number of high-quality positively charged muons = " << mup_list.size() << endl;
    cout << "VS-D-number of high-quality negatively charged muons = " << mum_list.size() << endl;
    cout << "VS-D-number of photon candidates = " << gam_list.size() << endl;
    cout << "VS-D-number of rho->pi+pi- candidates = " << rho_list.size() << endl;
    cout << "VS-D-number of YnS->mu+mu- candidates = " << YnS_list.size() << endl;
    cout << "VS-D-number of Wbj->YnSrho candidates = " << Wbj_list.size() << endl;
    cout << "VS-D-number of Y5S candidates = " << Y5S_list.size() << endl;
    cout << "VS-D-number of signal candidates = " << candidatesRec.size() << endl;
    */

    //photonSelectionEfficiencies();
    //chargedSelectionEfficiencies();
    //overallEfficiencies();
    
    prepareTreeRec(); //Copy candidates to STL vectors for writing out to ROOT file and prepare MC truth TTree

    diagnoseTreeRec();

    if ( isSignalMC || isISRMC || isNonISRMC ||  candidatesRec.size() > 0 ) writeTreesOut();

    /*
    cout << "count_Y5S = " << count_Y5S << endl;  // Use this to make sure your MC truth turns out okay! 
    cout << "count_Wbj = " << count_Wbj << endl;    
    cout << "count_gam = " << count_gam << endl; 
    cout << "count_YnS = " << count_YnS << endl;     
    cout << "count_rho = " << count_rho << endl;    
    cout << "count_mup = " << count_mup << endl;    
    cout << "count_mum = " << count_mum << endl;     
    cout << "count_pip = " << count_pip << endl;    
    cout << "count_pim = " << count_pim << endl;     
    */
    
    // vs  skimming: the event will be copied to output stream if there is at least one signal candidate
    // non-zero status will trigger skimming which occurs only for events with at least one signal candidate
    *status = candidatesRec.size();    

    // 11/29/2017:  for data:  skim only if the event has no candidates in the signal region (blinding)
    
    if ( !isMC ) {
      *status = 0;
      if ( m_ncand_sr == 0 ) *status = candidatesRec.size();    
    }

    // total number of events that pass skimming selection
    if ( *status ) m_event_seq_skim++;
    
  }




  
  void Analysis::analyzeHepEvtRecordSignal() {
    //cout << "analyzeHepEvtRecordSignal()" << endl;
    // @TODO make this significantly less redundant by turning it into a tree
    isSignalMC = false;
    isISRMC = false;
    isNonISRMC = false;
    
    m_mc_fsr_YnS = 0; // keeping track of Upsilon(nS) final state radiation
    m_mc_fsr_rho = 0; // and rho0 final state radiation
    m_mc_fsr_YnS_e = 0;
    m_mc_fsr_YnS_boost_e = 0;
    m_mc_fsr_rho_e = 0;
    m_mc_fsr_rho_boost_e = 0;
    nt_mc_fsr_YnS_ids.clear();
    cout << "testing1... " << m_mc_fsr_YnS_e << endl;

    if ( isMC ) {
      //cout << "if isMC" << endl;

      Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
      Gen_hepevt_Index gen_i = GenHepevtMgr.index( "mother" );
      gen_i.update(); // sort by pointer "mother" 	  
      
      int mc_Y5S_id = 0;
      int mc_Wbj_id = 0;
      int mc_gam_id = 0;
      int mc_YnS_id = 0;
      int mc_rho_id = 0;
      int mc_pip_id = 0;
      int mc_pim_id = 0;
      int mc_mup_id = 0;
      int mc_mum_id = 0;
      
      int mc_vgam_id = 0;
      
      int mc_YNS_isr_id = 0;
      int mc_gam_isr_id = 0;
      int mc_YnS_isr_id = 0;
      int mc_pip_isr_id = 0;
      int mc_pim_isr_id = 0;
      int mc_mup_isr_id = 0;
      int mc_mum_isr_id = 0;      
      
      for(std::vector<Gen_hepevt>::iterator it = GenHepevtMgr.begin(); it != GenHepevtMgr.end(); it++){

	Gen_hepevt* mc_Y5S   = 0;
	Gen_hepevt* mc_Wbj   = 0; 
	Gen_hepevt* mc_gam   = 0; 
	Gen_hepevt* mc_YnS   = 0;
	Gen_hepevt* mc_rho   = 0;  
	Gen_hepevt* mc_pip   = 0;
	Gen_hepevt* mc_pim   = 0;
	Gen_hepevt* mc_mup   = 0;
	Gen_hepevt* mc_mum   = 0;

	Gen_hepevt* mc_vgam  = 0;

	Gen_hepevt* mc_YNS_isr   = 0;
	Gen_hepevt* mc_gam_isr   = 0; 
	Gen_hepevt* mc_YnS_isr   = 0;
	Gen_hepevt* mc_pip_isr   = 0;
	Gen_hepevt* mc_pim_isr   = 0;
	Gen_hepevt* mc_mup_isr   = 0;
	Gen_hepevt* mc_mum_isr   = 0;


	Gen_hepevt& mc_particle = *it;
	int id_hep = mc_particle.idhep();

	//---------- MC TRUTH ----------//
	// Check for ISR 
	// ISR is simulated as e+e- -> virtual photon -> Y(5S)isr
	if ( id_hep == 10022) {
        
	  std::vector<Gen_hepevt> mc_vgam_daughters = point_from( mc_particle.get_ID(), gen_i ); // fetch virtual photon daughters

	  mc_vgam = &mc_particle;
	  mc_vgam_id = mc_vgam->get_ID();
	  
	  for(std::vector<Gen_hepevt>::iterator it_vgam = mc_vgam_daughters.begin(); it_vgam != mc_vgam_daughters.end(); it_vgam++){
	    Gen_hepevt &mc_daughter = *it_vgam;
	    
	    int mc_daughter_id_hep = mc_daughter.idhep();
	    int mc_daughter_id = mc_daughter.get_ID();
	    
	    if (mc_daughter_id_hep == 22) {
	      mc_gam_isr_id = mc_daughter_id;
	      mc_gam_isr = &mc_daughter;
	      m_isr = 1;
	    }
	    
	    // MC truth for ISR and RR to Y(2, 3S) 
	    else if (mc_daughter_id_hep == 9000553 || mc_daughter_id_hep == 200553 || mc_daughter_id_hep == 100553) {
	      mc_YNS_isr_id = mc_daughter_id;
	      mc_YNS_isr = &mc_daughter;
	      
	      std::vector<Gen_hepevt> mc_YNS_isr_daughters = point_from( mc_YNS_isr->get_ID(), gen_i );
	      for(std::vector<Gen_hepevt>::iterator it_YNS_isr = mc_YNS_isr_daughters.begin(); it_YNS_isr != mc_YNS_isr_daughters.end(); it_YNS_isr++){
		Gen_hepevt &mc_daughter = *it_YNS_isr;
		
		int mc_daughter_id_hep = mc_daughter.idhep();
		int mc_daughter_id = mc_daughter.get_ID();
		
		if (mc_daughter_id_hep == 211) {
		    mc_pip_isr_id = mc_daughter_id;
		    mc_pip_isr = &mc_daughter;
		}
		else if (mc_daughter_id_hep == -211) {
		  mc_pim_isr_id = mc_daughter_id;
		  mc_pim_isr = &mc_daughter;
		}
		
		else if (mc_daughter_id_hep == 553 || mc_daughter_id_hep == 100553) {
		  mc_YnS_isr_id = mc_daughter_id;
		  mc_YnS_isr = &mc_daughter;
		  
		  std::vector<Gen_hepevt> mc_YnS_isr_daughters = point_from( mc_YnS_isr->get_ID(), gen_i );
		  for(std::vector<Gen_hepevt>::iterator it_Y1S_isr = mc_YnS_isr_daughters.begin(); it_Y1S_isr != mc_YnS_isr_daughters.end(); it_Y1S_isr++){
		    Gen_hepevt &mc_daughter = *it_Y1S_isr;
		    
		    int mc_daughter_id_hep = mc_daughter.idhep();
		    int mc_daughter_id = mc_daughter.get_ID();
		    
		    if (mc_daughter_id_hep == -13) {
		      mc_mup_isr_id = mc_daughter_id;
		      mc_mup_isr = &mc_daughter;
		    }
		    else if (mc_daughter_id_hep == 13) {
		      mc_mum_isr_id = mc_daughter_id;
		      mc_mum_isr = &mc_daughter;
		    }
		  } //for YnS_isr
		} //if Y(1S)
	      } //for YNS_isr
	    } // if Y(5S)
	  }
	}

	// Non-ISR MC, i.e. Y5StoY1Smumu2pic
	// Signal MC truth
	if ( id_hep == 9000553 ) {
	  std::vector<Gen_hepevt> mc_Y5S_daughters = point_from( mc_particle.get_ID(), gen_i ); // fetch Y(5S) daughters

	  mc_Y5S = &mc_particle;
	  mc_Y5S_id = mc_Y5S->get_ID();

	  for(std::vector<Gen_hepevt>::iterator it_Y5S = mc_Y5S_daughters.begin(); it_Y5S != mc_Y5S_daughters.end(); it_Y5S++){
	    Gen_hepevt &mc_daughter = *it_Y5S;
	    
	    int mc_daughter_id_hep = mc_daughter.idhep();
	    int mc_daughter_id = mc_daughter.get_ID();
	    
	    if ( mc_daughter_id_hep == 553 ) {
	      mc_YnS_id = mc_daughter_id;
	      mc_YnS = &mc_daughter;
	    }

	    else if ( mc_daughter_id_hep == 211 ){
	      mc_pip_id = mc_daughter_id;
	      mc_pip = &mc_daughter;
	    }

	    else if ( mc_daughter_id_hep == -211 ){
	      mc_pim_id = mc_daughter_id;
	      mc_pim = &mc_daughter;
	    }

	    if ( mc_YnS_id != 0 && mc_pip_id != 0 && mc_pim_id != 0 ) {
	      std::vector<Gen_hepevt> mc_YnS_daughters = point_from( mc_YnS->get_ID(), gen_i );
	      
	      for ( std::vector<Gen_hepevt>::iterator it_YnS = mc_YnS_daughters.begin(); it_YnS != mc_YnS_daughters.end(); it_YnS++ ) {
		Gen_hepevt &mc_daughter = *it_YnS;
		
		int mc_daughter_id_hep = mc_daughter.idhep();
		int mc_daughter_id = mc_daughter.get_ID();
		
		if ( mc_daughter_id_hep == 13) { 
		  mc_mum_id = mc_daughter_id;
		  mc_mum = &mc_daughter;
		  count_mum++;
		}
		
		else if ( mc_daughter_id_hep == -13) {
		  mc_mup_id = mc_daughter_id;
		  mc_mup = &mc_daughter;
		  count_mup++;
		}
	      }
	    }
	  }
	}

	// Signal MC truth
	if ( id_hep == 9000553 ) {
	  std::vector<Gen_hepevt> mc_Y5S_daughters = point_from( mc_particle.get_ID(), gen_i ); // fetch Y(5S) daughters

	  if (mc_Y5S_daughters.size() == 2) { 
	    mc_Y5S = &mc_particle;
	    mc_Y5S_id = mc_Y5S->get_ID();

	    for(std::vector<Gen_hepevt>::iterator it_Y5S = mc_Y5S_daughters.begin(); it_Y5S != mc_Y5S_daughters.end(); it_Y5S++){
	      Gen_hepevt &mc_daughter = *it_Y5S;
	      
	      int mc_daughter_id_hep = mc_daughter.idhep();
	      int mc_daughter_id = mc_daughter.get_ID();
	      
	      if ( mc_daughter_id_hep == 50557 || mc_daughter_id_hep == 50558 || mc_daughter_id_hep == 50559 || mc_daughter_id_hep == 50560 ) {
		mc_Wbj_id = mc_daughter_id;
		mc_Wbj = &mc_daughter;
		count_Wbj++;
	      }

	      else if ( mc_daughter_id_hep == 22 ){
		mc_gam_id = mc_daughter_id;
		mc_gam = &mc_daughter;
		count_gam++;

	      }
	    } //for

	    if ( mc_Wbj_id != 0 && mc_gam_id != 0 ) {
	      std::vector<Gen_hepevt> mc_Wbj_daughters = point_from( mc_Wbj->get_ID(), gen_i );

	      if (mc_Wbj_daughters.size() == 2) {
	 	count_Y5S++;

		for ( std::vector<Gen_hepevt>::iterator it_Wbj = mc_Wbj_daughters.begin(); it_Wbj != mc_Wbj_daughters.end(); it_Wbj++ ) {
		  Gen_hepevt &mc_daughter = *it_Wbj;
	      
		  int mc_daughter_id_hep = mc_daughter.idhep();
		  int mc_daughter_id = mc_daughter.get_ID();
		 
		  if ( mc_daughter_id_hep == 553 || mc_daughter_id_hep == 100553 || mc_daughter_id_hep == 200553) { // include Upsilon(1, 2, 3S)
		    mc_YnS_id = mc_daughter_id;
		    mc_YnS = &mc_daughter;
		    count_YnS++;
		  }

		  else if ( mc_daughter_id_hep == 113 ) {
		    mc_rho_id = mc_daughter_id;
		    mc_rho = &mc_daughter;
		    count_rho++;
		  }
		} //for
	      
		if ( mc_YnS_id != 0 ) {
		  std::vector<Gen_hepevt> mc_YnS_daughters = point_from( mc_YnS->get_ID(), gen_i );
		  int n_other_daughters_YnS = 0;

		  if (mc_YnS_daughters.size() >= 2) {
		    for ( std::vector<Gen_hepevt>::iterator it_YnS = mc_YnS_daughters.begin(); it_YnS != mc_YnS_daughters.end(); it_YnS++ ) {
		      Gen_hepevt &mc_daughter = *it_YnS;
	      
		      int mc_daughter_id_hep = mc_daughter.idhep();
		      int mc_daughter_id = mc_daughter.get_ID();
		
		      if ( mc_daughter_id_hep == 13) { 
			mc_mum_id = mc_daughter_id;
			mc_mum = &mc_daughter;
			count_mum++;
		      }

		      else if ( mc_daughter_id_hep == -13) {
			mc_mup_id = mc_daughter_id;
			mc_mup = &mc_daughter;
			count_mup++;
			}

		      else if ( mc_daughter_id_hep == 22 ) {
			m_mc_fsr_YnS++;
			nt_mc_fsr_YnS_ids.insert(mc_daughter_id);
		      }

		      else {
			n_other_daughters_YnS++;
		      }
		    } //for
		  }
	      
		  if (mc_rho_id != 0 && n_other_daughters_YnS == 0 ) {
		    std::vector<Gen_hepevt> mc_rho_daughters = point_from( mc_rho->get_ID(), gen_i );
		    int n_other_daughters_rho = 0;

		    if (mc_rho_daughters.size() >= 2 ) {
		      for ( std::vector<Gen_hepevt>::iterator it_rho = mc_rho_daughters.begin(); it_rho != mc_rho_daughters.end(); it_rho++ ) {
			  
			Gen_hepevt &mc_daughter = *it_rho;
	      
		        int mc_daughter_id_hep = mc_daughter.idhep();
		        int mc_daughter_id = mc_daughter.get_ID();
	       
			if ( mc_daughter_id_hep == 211) { 
			  mc_pip_id = mc_daughter_id;
			  mc_pip = &mc_daughter;
			  count_pip++;
			}

			else if ( mc_daughter_id_hep == -211) {
			  mc_pim_id = mc_daughter_id;
			  mc_pim = &mc_daughter;
			  count_pim++;
			}

			else if ( mc_daughter_id_hep == 22 ) {
			  m_mc_fsr_rho++;
			  nt_mc_fsr_rho_ids.insert(mc_daughter_id);
			}

			else {
			  n_other_daughters_rho++;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	//cout<< "test ids : " << mc_Y5S_isr_id << " " << mc_gam_isr_id << endl;
	/*
    	if (m_mc_fsr_YnS > 0 && mc_YnS_id != 0) { //redundant check for YnS
	  
	  HepLorentzVector p4_mc_fsr_YnS_total(0,0,0,0);

	  Hep3Vector p_boost;
	  std::set<int>::iterator it;
	  cout << "SIZE = " << nt_mc_fsr_YnS_ids.size() << endl;
	  for (it=nt_mc_fsr_YnS_ids.begin(); it!=nt_mc_fsr_YnS_ids.end(); ++it) {
	    int id = *it;
	    cout << "###id = " << id << endl;
	    Gen_hepevt &mc_fsr_YnS = GenHepevtMgr[id-1];
	    Gen_hepevt &mc_YnS = GenHepevtMgr[m_mc_YnS_id-1];
	    HepLorentzVector p4_mc_fsr_YnS( mc_fsr_YnS.PX(), mc_fsr_YnS.PY(), mc_fsr_YnS.PZ(), mc_fsr_YnS.E() );
	    HepLorentzVector p4_mc_YnS( mc_YnS.PX(), mc_YnS.PY(), mc_YnS.PZ(), mc_YnS.E() );

	    p_boost = p4_mc_YnS.boostVector(); // calcuated several times, but oh well.
	    p4_mc_fsr_YnS_total += p4_mc_fsr_YnS;
	  }

	  HepLorentzVector p4_mc_fsr_YnS_total_boost = p4_mc_fsr_YnS_total;
	  p4_mc_fsr_YnS_total_boost.boost(-1.0*p_boost);
	  m_mc_fsr_YnS_e = p4_mc_fsr_YnS_total.e();
	  m_mc_fsr_YnS_boost_e = p4_mc_fsr_YnS_total_boost.e();
	  
	  nt_mc_fsr_YnS_ids.clear();
	}
	*/
	/*
	if (m_mc_fsr_rho > 0 && mc_rho_id !=0) { //redundant check for rho

	  HepLorentzVector p4_mc_fsr_rho_total;
	  std::set<int>::iterator it;
	  for (it=nt_mc_fsr_rho_ids.begin(); it!=nt_mc_fsr_rho_ids.end(); ++it) {
	    int id = *it;
	    Gen_hepevt &mc_fsr_rho = GenHepevtMgr[id-1];
	    Gen_hepevt &mc_rho = GenHepevtMgr[m_mc_rho_id-1];
	    HepLorentzVector p4_mc_fsr_rho( mc_fsr_rho.PX(), mc_fsr_rho.PY(), mc_fsr_rho.PZ(), mc_fsr_rho.E() );
	    HepLorentzVector p4_mc_rho( mc_rho.PX(), mc_rho.PY(), mc_rho.PZ(), mc_rho.E() );
	    p4_mc_fsr_rho_total += p4_mc_fsr_rho;
	  }

	  m_mc_fsr_rho_e = p4_mc_fsr_rho_total.e();
	  nt_mc_fsr_rho_ids.clear();
	}
	*/
	// ISR MC
	if (//mc_vgam_id    != 0 &&
	    mc_YNS_isr_id != 0 &&
	    //mc_gam_isr_id != 0 &&
	    mc_YnS_isr_id != 0 &&
	    mc_pip_isr_id != 0 &&
	    mc_pim_isr_id != 0 &&
	    mc_mup_isr_id != 0 &&
	    mc_mum_isr_id != 0    ) {

	  isISRMC = true;
	  isNonISRMC = false;
	  isSignalMC = false;
	  

	  m_mc_YNS_isr_id  = mc_YNS_isr_id;
	  m_mc_YnS_isr_id  = mc_YnS_isr_id;
	  m_mc_pip_isr_id  = mc_pip_isr_id;
	  m_mc_pim_isr_id  = mc_pim_isr_id;
	  m_mc_mup_isr_id  = mc_mup_isr_id;
	  m_mc_mum_isr_id  = mc_mum_isr_id;


	  Gen_hepevt &mc_YNS_isr = GenHepevtMgr[m_mc_YNS_isr_id-1];
	  Gen_hepevt &mc_YnS_isr = GenHepevtMgr[m_mc_YnS_isr_id-1];
	  Gen_hepevt &mc_pip_isr = GenHepevtMgr[m_mc_pip_isr_id-1];
	  Gen_hepevt &mc_pim_isr = GenHepevtMgr[m_mc_pim_isr_id-1];
	  Gen_hepevt &mc_mup_isr = GenHepevtMgr[m_mc_mup_isr_id-1];
	  Gen_hepevt &mc_mum_isr = GenHepevtMgr[m_mc_mum_isr_id-1];

	  HepLorentzVector p4_mc_YNS_isr( mc_YNS_isr.PX(), mc_YNS_isr.PY(), mc_YNS_isr.PZ(), mc_YNS_isr.E() );           
	  HepLorentzVector p4_mc_YnS_isr( mc_YnS_isr.PX(), mc_YnS_isr.PY(), mc_YnS_isr.PZ(), mc_YnS_isr.E() );
	  HepLorentzVector p4_mc_pip_isr( mc_pip_isr.PX(), mc_pip_isr.PY(), mc_pip_isr.PZ(), mc_pip_isr.E() );
	  HepLorentzVector p4_mc_pim_isr( mc_pim_isr.PX(), mc_pim_isr.PY(), mc_pim_isr.PZ(), mc_pim_isr.E() );
	  HepLorentzVector p4_mc_mup_isr( mc_mup_isr.PX(), mc_mup_isr.PY(), mc_mup_isr.PZ(), mc_mup_isr.E() );
	  HepLorentzVector p4_mc_mum_isr( mc_mum_isr.PX(), mc_mum_isr.PY(), mc_mum_isr.PZ(), mc_mum_isr.E() );

	  HepLorentzVector p4_mc_pip_pim_isr = p4_mc_pip_isr + p4_mc_pim_isr;

	  if (mc_vgam_id != 0 && mc_gam_isr_id != 0) {

	    m_mc_vgam_id    = mc_vgam_id;
	    m_mc_gam_isr_id = mc_gam_isr_id;
	    Gen_hepevt &mc_vgam    = GenHepevtMgr[m_mc_vgam_id-1];
	    Gen_hepevt &mc_gam_isr = GenHepevtMgr[m_mc_gam_isr_id-1];
	    HepLorentzVector p4_mc_vgam( mc_YNS_isr.PX(), mc_YNS_isr.PY(), mc_YNS_isr.PZ(), mc_YNS_isr.E() );
	    HepLorentzVector p4_mc_gam_isr( mc_gam_isr.PX(), mc_gam_isr.PY(), mc_gam_isr.PZ(), mc_gam_isr.E() ); 

	    Hep3Vector p_boost = p4_mc_vgam.boostVector();
	    HepLorentzVector p4_mc_gam_isr_boost = p4_mc_gam_isr;
	    p4_mc_gam_isr_boost.boost(-1.0*p_boost);

	    m_mc_gam_isr_e     = p4_mc_gam_isr.e();
	    m_mc_gam_isr_p     = p4_mc_gam_isr.rho();
	    m_mc_gam_isr_m     = p4_mc_gam_isr.m();
	    m_mc_gam_isr_pt    = p4_mc_gam_isr.perp();
	    m_mc_gam_isr_phi   = p4_mc_gam_isr.phi();
	    m_mc_gam_isr_costh = p4_mc_gam_isr.cosTheta();
	    
	    m_mc_gam_isr_boost_e     = p4_mc_gam_isr_boost.e();
	    m_mc_gam_isr_boost_p     = p4_mc_gam_isr_boost.rho();
	    m_mc_gam_isr_boost_m     = p4_mc_gam_isr_boost.m();
	    m_mc_gam_isr_boost_pt    = p4_mc_gam_isr_boost.perp();
	    m_mc_gam_isr_boost_phi   = p4_mc_gam_isr_boost.phi();
	    m_mc_gam_isr_boost_costh = p4_mc_gam_isr_boost.cosTheta();

	    m_mc_angle_bt_mup_gam_isr = p4_mc_mup_isr.angle( p4_mc_gam_isr );
	    m_mc_angle_bt_mum_gam_isr = p4_mc_mum_isr.angle( p4_mc_gam_isr );
	    
	    m_mc_angle_bt_pip_gam_isr = p4_mc_pip_isr.angle( p4_mc_gam_isr );
	    m_mc_angle_bt_pim_gam_isr = p4_mc_pim_isr.angle( p4_mc_gam_isr );

	  }

	  m_mc_YNS_isr_e     = p4_mc_YNS_isr.e();
	  m_mc_YNS_isr_p     = p4_mc_YNS_isr.rho();
	  m_mc_YNS_isr_m     = p4_mc_YNS_isr.m();
	  m_mc_YNS_isr_pt    = p4_mc_YNS_isr.perp();
	  m_mc_YNS_isr_phi   = p4_mc_YNS_isr.phi();
	  m_mc_YNS_isr_costh = p4_mc_YNS_isr.cosTheta();

	  m_mc_YnS_isr_e     = p4_mc_YnS_isr.e();
	  m_mc_YnS_isr_p     = p4_mc_YnS_isr.rho();
	  m_mc_YnS_isr_m     = p4_mc_YnS_isr.m();
	  m_mc_YnS_isr_pt    = p4_mc_YnS_isr.perp();
	  m_mc_YnS_isr_phi   = p4_mc_YnS_isr.phi();
	  m_mc_YnS_isr_costh = p4_mc_YnS_isr.cosTheta();

	  m_mc_pip_pim_isr_e     = p4_mc_pip_pim_isr.e();
	  m_mc_pip_pim_isr_p     = p4_mc_pip_pim_isr.rho();
	  m_mc_pip_pim_isr_m     = p4_mc_pip_pim_isr.m();
	  m_mc_pip_pim_isr_pt    = p4_mc_pip_pim_isr.perp();
	  m_mc_pip_pim_isr_phi   = p4_mc_pip_pim_isr.phi();
	  m_mc_pip_pim_isr_costh = p4_mc_pip_pim_isr.cosTheta();
	  
	  m_mc_pip_isr_e     = p4_mc_pip_isr.e();
	  m_mc_pip_isr_p     = p4_mc_pip_isr.rho();
	  m_mc_pip_isr_m     = p4_mc_pip_isr.m();
	  m_mc_pip_isr_pt    = p4_mc_pip_isr.perp();
	  m_mc_pip_isr_phi   = p4_mc_pip_isr.phi();
	  m_mc_pip_isr_costh = p4_mc_pip_isr.cosTheta();

	  m_mc_pim_isr_e     = p4_mc_pim_isr.e();
	  m_mc_pim_isr_p     = p4_mc_pim_isr.rho();
	  m_mc_pim_isr_m     = p4_mc_pim_isr.m();
	  m_mc_pim_isr_pt    = p4_mc_pim_isr.perp();
	  m_mc_pim_isr_phi   = p4_mc_pim_isr.phi();
	  m_mc_pim_isr_costh = p4_mc_pim_isr.cosTheta();

	  m_mc_mup_isr_e     = p4_mc_mup_isr.e();
	  m_mc_mup_isr_p     = p4_mc_mup_isr.rho();
	  m_mc_mup_isr_m     = p4_mc_mup_isr.m();
	  m_mc_mup_isr_pt    = p4_mc_mup_isr.perp();
	  m_mc_mup_isr_phi   = p4_mc_mup_isr.phi();
	  m_mc_mup_isr_costh = p4_mc_mup_isr.cosTheta();
	  
	  m_mc_mum_isr_e     = p4_mc_mum_isr.e();
	  m_mc_mum_isr_p     = p4_mc_mum_isr.rho();
	  m_mc_mum_isr_m     = p4_mc_mum_isr.m();
	  m_mc_mum_isr_pt    = p4_mc_mum_isr.perp();
	  m_mc_mum_isr_phi   = p4_mc_mum_isr.phi();
	  m_mc_mum_isr_costh = p4_mc_mum_isr.cosTheta();

	
	}

	// Non-ISR MC
	// Note that we require that Wbj, gam, and rho are NOT in the decay!!!
	// This could be seen as redundant because non-ISR MC is in a sense a "subset" of signal, because it just excludes some of the final state particles,
	// but keeping the signal/non-ISR analyses separate will help with clarity and debugging.
	if( mc_Y5S_id != 0 && 
	    mc_Wbj_id == 0 && 
	    mc_gam_id == 0 && 
	    mc_YnS_id != 0 && 
	    mc_rho_id == 0 && 
	    mc_pip_id != 0 && 
	    mc_pim_id != 0 && 
	    mc_mup_id != 0 && 
	    mc_mum_id != 0     ) {

	  isSignalMC = false;;
	  isISRMC = false;
	  isNonISRMC = true;

	  m_mc_Y5S_id  = mc_Y5S_id;
	  //m_mc_Wbj_id  = mc_Wbj_id;
	  //m_mc_gam_id  = mc_gam_id;
	  m_mc_YnS_id  = mc_YnS_id;
	  //m_mc_rho_id  = mc_rho_id;
	  m_mc_pip_id  = mc_pip_id;
	  m_mc_pim_id  = mc_pim_id;
	  m_mc_mup_id  = mc_mup_id;
	  m_mc_mum_id  = mc_mum_id;

	  Gen_hepevt &mc_Y5S = GenHepevtMgr[m_mc_Y5S_id-1];
	  //Gen_hepevt &mc_Wbj = GenHepevtMgr[m_mc_Wbj_id-1];
	  //Gen_hepevt &mc_gam = GenHepevtMgr[m_mc_gam_id-1];
	  Gen_hepevt &mc_YnS = GenHepevtMgr[m_mc_YnS_id-1];
	  //Gen_hepevt &mc_rho = GenHepevtMgr[m_mc_rho_id-1];
	  Gen_hepevt &mc_pip = GenHepevtMgr[m_mc_pip_id-1];
	  Gen_hepevt &mc_pim = GenHepevtMgr[m_mc_pim_id-1];
	  Gen_hepevt &mc_mup = GenHepevtMgr[m_mc_mup_id-1];
	  Gen_hepevt &mc_mum = GenHepevtMgr[m_mc_mum_id-1];

	  HepLorentzVector p4_mc_Y5S( mc_Y5S.PX(), mc_Y5S.PY(), mc_Y5S.PZ(), mc_Y5S.E() );
	  //HepLorentzVector p4_mc_Wbj( mc_Wbj.PX(), mc_Wbj.PY(), mc_Wbj.PZ(), mc_Wbj.E() );
	  //HepLorentzVector p4_mc_gam( mc_gam.PX(), mc_gam.PY(), mc_gam.PZ(), mc_gam.E() );
	  HepLorentzVector p4_mc_YnS( mc_YnS.PX(), mc_YnS.PY(), mc_YnS.PZ(), mc_YnS.E() );
	  //HepLorentzVector p4_mc_rho( mc_rho.PX(), mc_rho.PY(), mc_rho.PZ(), mc_rho.E() );            
	  HepLorentzVector p4_mc_pip( mc_pip.PX(), mc_pip.PY(), mc_pip.PZ(), mc_pip.E() );
	  HepLorentzVector p4_mc_pim( mc_pim.PX(), mc_pim.PY(), mc_pim.PZ(), mc_pim.E() );
	  HepLorentzVector p4_mc_mup( mc_mup.PX(), mc_mup.PY(), mc_mup.PZ(), mc_mup.E() );
	  HepLorentzVector p4_mc_mum( mc_mum.PX(), mc_mum.PY(), mc_mum.PZ(), mc_mum.E() );

	  HepLorentzVector p4_mc_pip_pim = p4_mc_pip + p4_mc_pim;

	  m_mc_Y5S_e     = p4_mc_Y5S.e();
	  m_mc_Y5S_p     = p4_mc_Y5S.rho();
	  m_mc_Y5S_m     = p4_mc_Y5S.m();
	  m_mc_Y5S_pt    = p4_mc_Y5S.perp();
	  m_mc_Y5S_phi   = p4_mc_Y5S.phi();
	  m_mc_Y5S_costh = p4_mc_Y5S.cosTheta();
	  
	  m_mc_YnS_e     = p4_mc_YnS.e();
	  m_mc_YnS_p     = p4_mc_YnS.rho();
	  m_mc_YnS_m     = p4_mc_YnS.m();
	  m_mc_YnS_pt    = p4_mc_YnS.perp();
	  m_mc_YnS_phi   = p4_mc_YnS.phi();
	  m_mc_YnS_costh = p4_mc_YnS.cosTheta();
	  
	  m_mc_pip_e     = p4_mc_pip.e();
	  m_mc_pip_p     = p4_mc_pip.rho();
	  m_mc_pip_m     = p4_mc_pip.m();
	  m_mc_pip_pt    = p4_mc_pip.perp();
	  m_mc_pip_phi   = p4_mc_pip.phi();
	  m_mc_pip_costh = p4_mc_pip.cosTheta();
	  
	  m_mc_pim_e     = p4_mc_pim.e();
	  m_mc_pim_p     = p4_mc_pim.rho();
	  m_mc_pim_m     = p4_mc_pim.m();
	  m_mc_pim_pt    = p4_mc_pim.perp();
	  m_mc_pim_phi   = p4_mc_pim.phi();
	  m_mc_pim_costh = p4_mc_pim.cosTheta();
	  
	  m_mc_mup_e     = p4_mc_mup.e();
	  m_mc_mup_p     = p4_mc_mup.rho();
	  m_mc_mup_m     = p4_mc_mup.m();
	  m_mc_mup_pt    = p4_mc_mup.perp();
	  m_mc_mup_phi   = p4_mc_mup.phi();
	  m_mc_mup_costh = p4_mc_mup.cosTheta();
	  
	  m_mc_mum_e     = p4_mc_mum.e();
	  m_mc_mum_p     = p4_mc_mum.rho();
	  m_mc_mum_m     = p4_mc_mum.m();
	  m_mc_mum_pt    = p4_mc_mum.perp();
	  m_mc_mum_phi   = p4_mc_mum.phi();
	  m_mc_mum_costh = p4_mc_mum.cosTheta();

	  m_mc_pip_pim_e     = p4_mc_pip_pim.e();
	  m_mc_pip_pim_p     = p4_mc_pip_pim.rho();
	  m_mc_pip_pim_m     = p4_mc_pip_pim.m();
	  m_mc_pip_pim_pt    = p4_mc_pip_pim.perp();
	  m_mc_pip_pim_phi   = p4_mc_pip_pim.phi();
	  m_mc_pip_pim_costh = p4_mc_pip_pim.cosTheta();

	}


	if( mc_Y5S_id != 0 && 
	    mc_Wbj_id != 0 && 
	    mc_gam_id != 0 && 
	    mc_YnS_id != 0 && 
	    mc_rho_id != 0 && 
	    mc_pip_id != 0 && 
	    mc_pim_id != 0 && 
	    mc_mup_id != 0 && 
	    mc_mum_id != 0     ) {

	  if ( isSignalMC ) {
	    //cout << "VS-W-Something is fishy about this event (more than one signal candidate found?) - will report its entire MC history now:" << endl;
	    //printCompleteMCHistory();
	  }

	  isSignalMC = true;
	  isISRMC = false;
	  isNonISRMC = false;
	  
	  m_mc_Y5S_id  = mc_Y5S_id;
	  m_mc_Wbj_id  = mc_Wbj_id;
	  m_mc_gam_id  = mc_gam_id;
	  m_mc_YnS_id  = mc_YnS_id;
	  m_mc_rho_id  = mc_rho_id;
	  m_mc_pip_id  = mc_pip_id;
	  m_mc_pim_id  = mc_pim_id;
	  m_mc_mup_id  = mc_mup_id;
	  m_mc_mum_id  = mc_mum_id;

	  Gen_hepevt &mc_Y5S = GenHepevtMgr[m_mc_Y5S_id-1];
	  Gen_hepevt &mc_Wbj = GenHepevtMgr[m_mc_Wbj_id-1];
	  Gen_hepevt &mc_gam = GenHepevtMgr[m_mc_gam_id-1];
	  Gen_hepevt &mc_YnS = GenHepevtMgr[m_mc_YnS_id-1];
	  Gen_hepevt &mc_rho = GenHepevtMgr[m_mc_rho_id-1];
	  Gen_hepevt &mc_pip = GenHepevtMgr[m_mc_pip_id-1];
	  Gen_hepevt &mc_pim = GenHepevtMgr[m_mc_pim_id-1];
	  Gen_hepevt &mc_mup = GenHepevtMgr[m_mc_mup_id-1];
	  Gen_hepevt &mc_mum = GenHepevtMgr[m_mc_mum_id-1];

	  HepLorentzVector p4_mc_Y5S( mc_Y5S.PX(), mc_Y5S.PY(), mc_Y5S.PZ(), mc_Y5S.E() );
	  HepLorentzVector p4_mc_Wbj( mc_Wbj.PX(), mc_Wbj.PY(), mc_Wbj.PZ(), mc_Wbj.E() );
	  HepLorentzVector p4_mc_gam( mc_gam.PX(), mc_gam.PY(), mc_gam.PZ(), mc_gam.E() );
	  HepLorentzVector p4_mc_rho( mc_rho.PX(), mc_rho.PY(), mc_rho.PZ(), mc_rho.E() );            
	  HepLorentzVector p4_mc_YnS( mc_YnS.PX(), mc_YnS.PY(), mc_YnS.PZ(), mc_YnS.E() );
	  HepLorentzVector p4_mc_pip( mc_pip.PX(), mc_pip.PY(), mc_pip.PZ(), mc_pip.E() );
	  HepLorentzVector p4_mc_pim( mc_pim.PX(), mc_pim.PY(), mc_pim.PZ(), mc_pim.E() );
	  HepLorentzVector p4_mc_mup( mc_mup.PX(), mc_mup.PY(), mc_mup.PZ(), mc_mup.E() );
	  HepLorentzVector p4_mc_mum( mc_mum.PX(), mc_mum.PY(), mc_mum.PZ(), mc_mum.E() );
	  
	  // Boost gamma to Y(5S) CM reference frame
	  Hep3Vector p_boost = p4_mc_Y5S.boostVector();

	  m_p4_mc_Y5S = p4_mc_Y5S;
	  m_p4_mc_gam = p4_mc_gam;

	  HepLorentzVector p4_mc_gam_boost = p4_mc_gam;
	  p4_mc_gam_boost.boost(-1.0*p_boost);
	  
	  // Calculate missing p4
	  //HepLorentzVector p4_mc_gam_recoil_4 = m_pExp - (p4_mc_gam);
	  //HepLorentzVector p4_mc_gam_recoil_1 = p4_mc_Y5S - (p4_mc_gam);
	  HepLorentzVector p4_mc_mup_mum_gam_recoil = m_pExp - (p4_mc_gam + p4_mc_mup + p4_mc_mum);
	  HepLorentzVector p4_mc_pip_pim_gam_recoil = m_pExp - (p4_mc_gam + p4_mc_pip + p4_mc_pim);
	  HepLorentzVector p4_mc_pip_pim_mup_mum_recoil = m_pExp - (p4_mc_mup + p4_mc_mum + p4_mc_pip + p4_mc_pim);

	  HepLorentzVector p4_mc_YnS_boost = p4_mc_YnS;
	  HepLorentzVector p4_mc_Wbj_boost = p4_mc_Wbj;
	  HepLorentzVector p4_mc_pip_pim_gam_recoil_boost = p4_mc_pip_pim_gam_recoil;
	  HepLorentzVector p4_mc_mup_mum_gam_recoil_boost = p4_mc_mup_mum_gam_recoil;	  
	  HepLorentzVector p4_mc_pip_pim_mup_mum_recoil_boost = p4_mc_pip_pim_mup_mum_recoil;

	  p4_mc_YnS_boost.boost(-1.0*p_boost);
	  p4_mc_Wbj_boost.boost(-1.0*p_boost);
	  p4_mc_pip_pim_gam_recoil_boost.boost(-1.0*p_boost);
	  p4_mc_mup_mum_gam_recoil_boost.boost(-1.0*p_boost);
	  p4_mc_pip_pim_mup_mum_recoil_boost.boost(-1.0*p_boost);	  

	  HepLorentzVector p4_mup_mum_gam = p4_mc_mup + p4_mc_mum + p4_mc_gam;
	  HepLorentzVector p4_pip_pim_gam = p4_mc_pip + p4_mc_pim + p4_mc_gam;

	  double mc_pip_pim_gam_recoil_m2 = p4_mc_pip_pim_gam_recoil.m2();
	  double mc_mup_mum_gam_recoil_m2 = p4_mc_mup_mum_gam_recoil.m2();
	  double mc_pip_pim_gam_recoil_boost_m2 = p4_mc_pip_pim_gam_recoil.m2();
	  double mc_mup_mum_gam_recoil_boost_m2 = p4_mc_mup_mum_gam_recoil.m2();

	  m_mc_angle_bt_mup_gam = p4_mc_mup.angle( p4_mc_gam );
	  m_mc_angle_bt_mum_gam = p4_mc_mum.angle( p4_mc_gam );

	  m_mc_angle_bt_pip_gam = p4_mc_pip.angle( p4_mc_gam );
	  m_mc_angle_bt_pim_gam = p4_mc_pim.angle( p4_mc_gam );

	  m_mc_Y5S_e     = p4_mc_Y5S.e();
	  m_mc_Y5S_p     = p4_mc_Y5S.rho();
	  m_mc_Y5S_m     = p4_mc_Y5S.m();
	  m_mc_Y5S_pt    = p4_mc_Y5S.perp();
	  m_mc_Y5S_phi   = p4_mc_Y5S.phi();
	  m_mc_Y5S_costh = p4_mc_Y5S.cosTheta();
	  
	  m_mc_Wbj_e     = p4_mc_Wbj.e();
	  m_mc_Wbj_p     = p4_mc_Wbj.rho();
	  m_mc_Wbj_m     = p4_mc_Wbj.m();
	  m_mc_Wbj_pt    = p4_mc_Wbj.perp();
	  m_mc_Wbj_phi   = p4_mc_Wbj.phi();
	  m_mc_Wbj_costh = p4_mc_Wbj.cosTheta();

	  m_mc_Wbj_boost_e     = p4_mc_Wbj_boost.e();
	  m_mc_Wbj_boost_p     = p4_mc_Wbj_boost.rho();
	  m_mc_Wbj_boost_m     = p4_mc_Wbj_boost.m();
	  m_mc_Wbj_boost_pt    = p4_mc_Wbj_boost.perp();
	  m_mc_Wbj_boost_phi   = p4_mc_Wbj_boost.phi();
	  m_mc_Wbj_boost_costh = p4_mc_Wbj_boost.cosTheta();
	 
	  m_mc_gam_e     = p4_mc_gam.e();
	  m_mc_gam_p     = p4_mc_gam.rho();
	  m_mc_gam_m     = p4_mc_gam.m();
	  m_mc_gam_pt    = p4_mc_gam.perp();
	  m_mc_gam_phi   = p4_mc_gam.phi();
	  m_mc_gam_costh = p4_mc_gam.cosTheta();
	  
	  m_mc_YnS_e     = p4_mc_YnS.e();
	  m_mc_YnS_p     = p4_mc_YnS.rho();
	  m_mc_YnS_m     = p4_mc_YnS.m();
	  m_mc_YnS_pt    = p4_mc_YnS.perp();
	  m_mc_YnS_phi   = p4_mc_YnS.phi();
	  m_mc_YnS_costh = p4_mc_YnS.cosTheta();
	  
	  m_mc_YnS_boost_e     = p4_mc_YnS_boost.e();
	  m_mc_YnS_boost_p     = p4_mc_YnS_boost.rho();
	  m_mc_YnS_boost_m     = p4_mc_YnS_boost.m();
	  m_mc_YnS_boost_pt    = p4_mc_YnS_boost.perp();
	  m_mc_YnS_boost_phi   = p4_mc_YnS_boost.phi();
	  m_mc_YnS_boost_costh = p4_mc_YnS_boost.cosTheta();

	  m_mc_rho_e     = p4_mc_rho.e();
	  m_mc_rho_p     = p4_mc_rho.rho();
	  m_mc_rho_m     = p4_mc_rho.m();
	  m_mc_rho_pt    = p4_mc_rho.perp();
	  m_mc_rho_phi   = p4_mc_rho.phi();
	  m_mc_rho_costh = p4_mc_rho.cosTheta();
	  
	  m_mc_pip_e     = p4_mc_pip.e();
	  m_mc_pip_p     = p4_mc_pip.rho();
	  m_mc_pip_m     = p4_mc_pip.m();
	  m_mc_pip_pt    = p4_mc_pip.perp();
	  m_mc_pip_phi   = p4_mc_pip.phi();
	  m_mc_pip_costh = p4_mc_pip.cosTheta();
	  
	  m_mc_pim_e     = p4_mc_pim.e();
	  m_mc_pim_p     = p4_mc_pim.rho();
	  m_mc_pim_m     = p4_mc_pim.m();
	  m_mc_pim_pt    = p4_mc_pim.perp();
	  m_mc_pim_phi   = p4_mc_pim.phi();
	  m_mc_pim_costh = p4_mc_pim.cosTheta();
	  
	  m_mc_mup_e     = p4_mc_mup.e();
	  m_mc_mup_p     = p4_mc_mup.rho();
	  m_mc_mup_m     = p4_mc_mup.m();
	  m_mc_mup_pt    = p4_mc_mup.perp();
	  m_mc_mup_phi   = p4_mc_mup.phi();
	  m_mc_mup_costh = p4_mc_mup.cosTheta();
	  
	  m_mc_mum_e     = p4_mc_mum.e();
	  m_mc_mum_p     = p4_mc_mum.rho();
	  m_mc_mum_m     = p4_mc_mum.m();
	  m_mc_mum_pt    = p4_mc_mum.perp();
	  m_mc_mum_phi   = p4_mc_mum.phi();
	  m_mc_mum_costh = p4_mc_mum.cosTheta();

	  m_mc_pip_pim_gam_recoil_e     = p4_mc_pip_pim_gam_recoil.e();
	  m_mc_pip_pim_gam_recoil_p     = p4_mc_pip_pim_gam_recoil.rho();
	  m_mc_pip_pim_gam_recoil_m     = (mc_pip_pim_gam_recoil_m2 > 0.) ? sqrt(mc_pip_pim_gam_recoil_m2) : -sqrt(-mc_pip_pim_gam_recoil_m2);
	  m_mc_pip_pim_gam_recoil_pt    = p4_mc_pip_pim_gam_recoil.perp();
	  m_mc_pip_pim_gam_recoil_phi   = p4_mc_pip_pim_gam_recoil.phi();
	  m_mc_pip_pim_gam_recoil_costh = p4_mc_pip_pim_gam_recoil.cosTheta();
	  
	  m_mc_pip_pim_gam_recoil_boost_e     = p4_mc_pip_pim_gam_recoil_boost.e();
	  m_mc_pip_pim_gam_recoil_boost_p     = p4_mc_pip_pim_gam_recoil_boost.rho();
	  m_mc_pip_pim_gam_recoil_boost_m     = (mc_pip_pim_gam_recoil_boost_m2 > 0.) ? sqrt(mc_pip_pim_gam_recoil_boost_m2) : -sqrt(-mc_pip_pim_gam_recoil_boost_m2);
	  m_mc_pip_pim_gam_recoil_boost_pt    = p4_mc_pip_pim_gam_recoil_boost.perp();
	  m_mc_pip_pim_gam_recoil_boost_phi   = p4_mc_pip_pim_gam_recoil_boost.phi();
	  m_mc_pip_pim_gam_recoil_boost_costh = p4_mc_pip_pim_gam_recoil_boost.cosTheta();

	  
	  m_mc_mup_mum_gam_recoil_e     = p4_mc_mup_mum_gam_recoil.e();
	  m_mc_mup_mum_gam_recoil_p     = p4_mc_mup_mum_gam_recoil.rho();
	  m_mc_mup_mum_gam_recoil_m     = (mc_mup_mum_gam_recoil_m2 > 0.) ? sqrt(mc_mup_mum_gam_recoil_m2) : -sqrt(-mc_mup_mum_gam_recoil_m2);
	  m_mc_mup_mum_gam_recoil_pt    = p4_mc_mup_mum_gam_recoil.perp();
	  m_mc_mup_mum_gam_recoil_phi   = p4_mc_mup_mum_gam_recoil.phi();
	  m_mc_mup_mum_gam_recoil_costh = p4_mc_mup_mum_gam_recoil.cosTheta();

	  m_mc_mup_mum_gam_recoil_boost_e     = p4_mc_mup_mum_gam_recoil_boost.e();
	  m_mc_mup_mum_gam_recoil_boost_p     = p4_mc_mup_mum_gam_recoil_boost.rho();
	  m_mc_mup_mum_gam_recoil_boost_m     = (mc_mup_mum_gam_recoil_boost_m2 > 0.) ? sqrt(mc_mup_mum_gam_recoil_boost_m2) : -sqrt(-mc_mup_mum_gam_recoil_boost_m2);
	  m_mc_mup_mum_gam_recoil_boost_pt    = p4_mc_mup_mum_gam_recoil_boost.perp();
	  m_mc_mup_mum_gam_recoil_boost_phi   = p4_mc_mup_mum_gam_recoil_boost.phi();
	  m_mc_mup_mum_gam_recoil_boost_costh = p4_mc_mup_mum_gam_recoil_boost.cosTheta();
	  
	  m_mc_pip_pim_mup_mum_recoil_e     = p4_mc_pip_pim_mup_mum_recoil.e();
	  m_mc_pip_pim_mup_mum_recoil_p     = p4_mc_pip_pim_mup_mum_recoil.rho();
	  m_mc_pip_pim_mup_mum_recoil_m     = p4_mc_pip_pim_mup_mum_recoil.m();
	  m_mc_pip_pim_mup_mum_recoil_pt    = p4_mc_pip_pim_mup_mum_recoil.perp();
	  m_mc_pip_pim_mup_mum_recoil_phi   = p4_mc_pip_pim_mup_mum_recoil.phi();
	  m_mc_pip_pim_mup_mum_recoil_costh = p4_mc_pip_pim_mup_mum_recoil.cosTheta();

	  m_mc_pip_pim_mup_mum_recoil_m2    = p4_mc_pip_pim_mup_mum_recoil.m2();

	  m_mc_pip_pim_mup_mum_recoil_boost_e     = p4_mc_pip_pim_mup_mum_recoil_boost.e();
	  m_mc_pip_pim_mup_mum_recoil_boost_p     = p4_mc_pip_pim_mup_mum_recoil_boost.rho();
	  m_mc_pip_pim_mup_mum_recoil_boost_m     = p4_mc_pip_pim_mup_mum_recoil_boost.m();
	  m_mc_pip_pim_mup_mum_recoil_boost_pt    = p4_mc_pip_pim_mup_mum_recoil_boost.perp();
	  m_mc_pip_pim_mup_mum_recoil_boost_phi   = p4_mc_pip_pim_mup_mum_recoil_boost.phi();
	  m_mc_pip_pim_mup_mum_recoil_boost_costh = p4_mc_pip_pim_mup_mum_recoil_boost.cosTheta();

	  m_mc_gam_boost_e     = p4_mc_gam_boost.e();
	  m_mc_gam_boost_p     = p4_mc_gam_boost.rho();
	  m_mc_gam_boost_m     = p4_mc_gam_boost.m();
	  m_mc_gam_boost_pt    = p4_mc_gam_boost.perp();
	  m_mc_gam_boost_phi   = p4_mc_gam_boost.phi();
	  m_mc_gam_boost_costh = p4_mc_gam_boost.cosTheta(); 
	  
	  m_mc_mup_mum_gam_m = p4_mup_mum_gam.m();
	  m_mc_pip_pim_gam_m = p4_pip_pim_gam.m();

	  //cout << "mc_vgam_id = " << mc_vgam_id << endl;
	  //cout << "mc_gam_isr_id = " << mc_gam_isr_id << endl;
	  if ( mc_vgam_id != 0 && mc_gam_isr_id != 0 ) {

	    m_mc_gam_isr_id  = mc_gam_isr_id;

	    Gen_hepevt &mc_gam_isr = GenHepevtMgr[m_mc_gam_isr_id-1];

	    HepLorentzVector p4_mc_gam_isr( mc_gam_isr.PX(), mc_gam_isr.PY(), mc_gam_isr.PZ(), mc_gam_isr.E() );

	    HepLorentzVector p4_mc_gam_isr_boost = p4_mc_gam_isr;
	    p4_mc_gam_isr_boost.boost(-1.0*p_boost);

	    m_mc_gam_isr_e     = p4_mc_gam_isr.e();
	    m_mc_gam_isr_p     = p4_mc_gam_isr.rho();
	    m_mc_gam_isr_m     = p4_mc_gam_isr.m();
	    m_mc_gam_isr_pt    = p4_mc_gam_isr.perp();
	    m_mc_gam_isr_phi   = p4_mc_gam_isr.phi();
	    m_mc_gam_isr_costh = p4_mc_gam_isr.cosTheta();

	    m_mc_gam_isr_boost_e     = p4_mc_gam_isr_boost.e();
	    m_mc_gam_isr_boost_p     = p4_mc_gam_isr_boost.rho();
	    m_mc_gam_isr_boost_m     = p4_mc_gam_isr_boost.m();
	    m_mc_gam_isr_boost_pt    = p4_mc_gam_isr_boost.perp();
	    m_mc_gam_isr_boost_phi   = p4_mc_gam_isr_boost.phi();
	    m_mc_gam_isr_boost_costh = p4_mc_gam_isr_boost.cosTheta();
	  }

	  if (nt_mc_fsr_YnS_ids.size() != 0 && mc_YnS_id != 0) { //redundant check for YnS
	    
	    HepLorentzVector p4_mc_fsr_YnS_total(0,0,0,0);
	    
	    Hep3Vector p_boost;
	    std::set<int>::iterator it;
	    cout << "SIZE = " << nt_mc_fsr_YnS_ids.size() << endl;
	    for (it=nt_mc_fsr_YnS_ids.begin(); it!=nt_mc_fsr_YnS_ids.end(); ++it) {
	      int id = *it;
	      cout << "###id = " << id << endl;
	      Gen_hepevt &mc_fsr_YnS = GenHepevtMgr[id-1];
	      Gen_hepevt &mc_YnS = GenHepevtMgr[m_mc_YnS_id-1];
	      HepLorentzVector p4_mc_fsr_YnS( mc_fsr_YnS.PX(), mc_fsr_YnS.PY(), mc_fsr_YnS.PZ(), mc_fsr_YnS.E() );
	      HepLorentzVector p4_mc_YnS( mc_YnS.PX(), mc_YnS.PY(), mc_YnS.PZ(), mc_YnS.E() );
	      
	      p_boost = p4_mc_YnS.boostVector(); // calcuated several times, but oh well.
	      p4_mc_fsr_YnS_total += p4_mc_fsr_YnS;
	    }
	    nt_mc_fsr_YnS_ids.clear();
	    HepLorentzVector p4_mc_fsr_YnS_total_boost = p4_mc_fsr_YnS_total;
	    p4_mc_fsr_YnS_total_boost.boost(-1.0*p_boost);
	    m_mc_fsr_YnS_e = p4_mc_fsr_YnS_total.e();
	    cout << "energY? = " << m_mc_fsr_YnS_e << endl;
	    m_mc_fsr_YnS_boost_e = p4_mc_fsr_YnS_total_boost.e();
	  }

	  if (nt_mc_fsr_rho_ids.size() != 0 && mc_rho_id != 0) { //redundant check for rho
	    
	    HepLorentzVector p4_mc_fsr_rho_total(0,0,0,0);
	    
	    Hep3Vector p_boost;
	    std::set<int>::iterator it;
	    cout << "SIZE = " << nt_mc_fsr_rho_ids.size() << endl;
	    for (it=nt_mc_fsr_rho_ids.begin(); it!=nt_mc_fsr_rho_ids.end(); ++it) {
	      int id = *it;
	      cout << "###id = " << id << endl;
	      Gen_hepevt &mc_fsr_rho = GenHepevtMgr[id-1];
	      Gen_hepevt &mc_rho = GenHepevtMgr[m_mc_rho_id-1];
	      HepLorentzVector p4_mc_fsr_rho( mc_fsr_rho.PX(), mc_fsr_rho.PY(), mc_fsr_rho.PZ(), mc_fsr_rho.E() );
	      HepLorentzVector p4_mc_rho( mc_rho.PX(), mc_rho.PY(), mc_rho.PZ(), mc_rho.E() );
	      
	      p_boost = p4_mc_rho.boostVector(); // calcuated several times, but oh well.
	      p4_mc_fsr_rho_total += p4_mc_fsr_rho;
	    }
	    nt_mc_fsr_rho_ids.clear();
	    HepLorentzVector p4_mc_fsr_rho_total_boost = p4_mc_fsr_rho_total;
	    p4_mc_fsr_rho_total_boost.boost(-1.0*p_boost);
	    m_mc_fsr_rho_e = p4_mc_fsr_rho_total.e();
	    cout << "energY? = " << m_mc_fsr_rho_e << endl;
	    m_mc_fsr_rho_boost_e = p4_mc_fsr_rho_total_boost.e();
	  }

	}
      }//for
    } //if
  } //func


  void Analysis::selectChargedPions() {
    //cout << "selectChargedPions" << endl;
    
    pic_list.clear();
    pip_list.clear();
    pim_list.clear();

    Mdst_charged_Manager &ChargedMgr = Mdst_charged_Manager::get_manager();
    for ( std::vector<Mdst_charged>::iterator it = ChargedMgr.begin(); it != ChargedMgr.end(); it++ ) {
      
      Mdst_charged &chrg = *it;
      
      Particle pion( chrg, Ptype( chrg.charge() > 0 ? "PI+" : "PI-" ) );
      
      if ( goodTrackQuality( chrg, "pi" ) ) {
	
	pic_list.push_back( pion );
	set_mdst_chrg_id.insert( chrg.get_ID() );

	if ( goodChargedPion( chrg ) ) {
	  
	  if ( chrg.charge() > 0 ) {
	    pip_list.push_back( pion );
	  }
	  else {
	    pim_list.push_back( pion );
	  }
	}
      }
    }
    if ( isMC ) {
      setGenHepInfoF( pic_list );
      setGenHepInfoF( pip_list );
      setGenHepInfoF( pim_list );
    }
  }



  

  void Analysis::selectChargedMuons() {
    //cout << "selectChargedMuons" << endl;

    muc_list.clear();
    mup_list.clear();
    mum_list.clear();

    Mdst_charged_Manager &ChargedMgr = Mdst_charged_Manager::get_manager();

    for ( std::vector<Mdst_charged>::iterator it = ChargedMgr.begin(); it != ChargedMgr.end(); it++ ) {

      Mdst_charged &chrg = *it;

      Particle muon( chrg, Ptype( chrg.charge() > 0 ? "MU+" : "MU-"));

      if ( goodTrackQuality( chrg, "mu" ) ) {
	
	set_mdst_chrg_id.insert( chrg.get_ID() );

	if ( goodChargedMuon( chrg ) ) {
	  
	  if ( chrg.charge() > 0 ) {
	    mup_list.push_back( muon );
	  }
	  else {
	    mum_list.push_back( muon );
	  }
	}
      }
    }
    if ( isMC ) {
      setGenHepInfoF( mup_list );
      setGenHepInfoF( mum_list );
    } 
  }


      
  //I use "gamma" as my variable rather than gam because these functions should work for any photon in any process. So, it's a little misleading to use gam since thats a specific photon.
  void Analysis::selectPhotons() {
    //cout << "selectPhotons" << endl;

    gam_list.clear();

    Mdst_gamma_Manager &GammaMgr = Mdst_gamma_Manager::get_manager();

    for( std::vector<Mdst_gamma>::iterator it = GammaMgr.begin(); it != GammaMgr.end(); it++ ) {

      Mdst_gamma &mdst_gamma = *it;

      if( goodGamma( mdst_gamma ) ) {

	// Create photon particle 
	Particle photon( mdst_gamma );

	if ( photon.e() >= 0.1 && photon.e() <= 0.6 ) {
	  
	  gam_list.push_back( photon );
	  m_ncand_gam = gam_list.size(); //for signal background analysis
	}
      }
    }
    if ( isMC ) {
      setGenHepInfoG( gam_list );
    }
  }






  void Analysis::selectRhos() {
    //cout << "selectRhos" << endl;

    rho_list.clear();

    for ( std::vector<Particle>::iterator it_pip = pip_list.begin(); it_pip !=  pip_list.end(); it_pip++ ) {
      for ( std::vector<Particle>::iterator it_pim = pim_list.begin(); it_pim !=  pim_list.end(); it_pim++ ) {

	  Particle &pip = (*it_pip);
	  Particle &pim = (*it_pim);

	  HepLorentzVector p4_rho = pip.p()+pim.p();
	  if ( abs(p4_rho.m() - 0.770 ) < 0.35) {
	  //if ( p4_rho.m() >= 0.280 || p4_rho.m() <= 1.0 ) {
	    
	    // create rho
	    Particle rho( p4_rho, Ptype("RHO0"));
	    rho.relation().append(pip);
	    rho.relation().append(pim);
	    
	    rho_list.push_back(rho);
	    m_total_rho++;
	  }
      }
    }
  }


  void Analysis::countChargedTracks() {

    m_NChargedTracks = 0;

    std::set<int> trk_ids;

    for ( std::vector<Particle>::iterator it_trk = pip_list.begin(); it_trk !=  pip_list.end(); it_trk++ ) {
      const Mdst_charged &mdst_trk = (*it_trk).mdstCharged();
      int id = mdst_trk.get_ID();
      //      cout << "-D-VS-pip trk id = " << id << endl;
      trk_ids.insert(id);
    }

    for ( std::vector<Particle>::iterator it_trk = pim_list.begin(); it_trk !=  pim_list.end(); it_trk++ ) {
      const Mdst_charged &mdst_trk = (*it_trk).mdstCharged();
      int id = mdst_trk.get_ID();
      //      cout << "-D-VS-pim trk id = " << id << endl;
      trk_ids.insert(id);
    }

    for ( std::vector<Particle>::iterator it_trk = mup_list.begin(); it_trk !=  mup_list.end(); it_trk++ ) {
      const Mdst_charged &mdst_trk = (*it_trk).mdstCharged();
      int id = mdst_trk.get_ID();
      //      cout << "-D-VS-mup trk id = " << id << endl;
      trk_ids.insert(id);
    }

    for ( std::vector<Particle>::iterator it_trk = mum_list.begin(); it_trk !=  mum_list.end(); it_trk++ ) {
      const Mdst_charged &mdst_trk = (*it_trk).mdstCharged();
      int id = mdst_trk.get_ID();
      //      cout << "-D-VS-mum trk id = " << id << endl;
      trk_ids.insert(id);
    }

    //    cout << "-D-VS-size of set = " << trk_ids.size() << endl;

    m_NChargedTracks = trk_ids.size();

  }


  void Analysis::selectYnS() {
    //cout << "selectYnS" << endl;

    YnS_list.clear();

    for ( std::vector<Particle>::iterator it_mup = mup_list.begin(); it_mup !=  mup_list.end(); it_mup++ ) {
      for ( std::vector<Particle>::iterator it_mum = mum_list.begin(); it_mum !=  mum_list.end(); it_mum++ ) {

	Particle &mup = (*it_mup);
	Particle &mum = (*it_mum);
	  
	HepLorentzVector p4_YnS = mup.p()+mum.p();

	if ( p4_YnS.m() > 9.3 && p4_YnS.m() < 19.6 )  {
	  // create YnS
	  Particle YnS( p4_YnS, Ptype("UPS1")); // It doesn't matter if we use UPS1, UPS2, or UPS3 here
	  YnS.relation().append(mup);
	  YnS.relation().append(mum);
	  
	  YnS_list.push_back(YnS);
	  m_total_YnS++;
	}
      }
    }
  }

  void Analysis::selectWbj() {
    //cout << "selectWbj" << endl;

    Wbj_list.clear();

    for ( std::vector<Particle>::iterator it_rho = rho_list.begin(); it_rho !=  rho_list.end(); it_rho++ ) {
      for ( std::vector<Particle>::iterator it_YnS = YnS_list.begin(); it_YnS !=  YnS_list.end(); it_YnS++ ) {

	Particle &rho = (*it_rho);
	Particle &YnS = (*it_YnS);

	HepLorentzVector p4_Wbj = YnS.p()+rho.p();

	if ( p4_Wbj.m() > 9.2 ) {
	  Particle Wbj( p4_Wbj, Ptype() ); // This particle doesn't exist in the .pdl file, so we must use the default Ptype constructor--no concerns here.

	  Wbj.relation().append(YnS);
	  Wbj.relation().append(rho);
	    
	  Wbj_list.push_back(Wbj);
	  m_total_Wbj++;
	}
      }
    }
  }


  void Analysis::selectY5S() {
    //cout << "selectY5S" << endl;
    
    Y5S_list.clear();

    for ( std::vector<Particle>::iterator it_Wbj = Wbj_list.begin(); it_Wbj !=  Wbj_list.end(); it_Wbj++ ) {
      for ( std::vector<Particle>::iterator it_gam = gam_list.begin(); it_gam !=  gam_list.end(); it_gam++ ) {

	Particle &Wbj = (*it_Wbj);
     	Particle &gam = (*it_gam);

	// Retrieve daughters
	Mdst_charged const &mdst_mup = Wbj.relation().child(0).relation().child(0).mdstCharged();
	Mdst_charged const &mdst_mum = Wbj.relation().child(0).relation().child(1).mdstCharged();
	Mdst_charged const &mdst_pip = Wbj.relation().child(1).relation().child(0).mdstCharged();
	Mdst_charged const &mdst_pim = Wbj.relation().child(1).relation().child(1).mdstCharged();
	
	// Check if all charged tracks are unique
	if ( mdst_mup != mdst_pip &&
	     mdst_mum != mdst_pim    ) {
	  
	  HepLorentzVector p4_Y5S = Wbj.p()+gam.p();
	  HepLorentzVector p4_Y5S_boost = Wbj.p()+gam.p();
	  Hep3Vector boost_CM = p4_Y5S.boostVector();

	  p4_Y5S_boost.boost(-1.0*boost_CM);

	  m_Ntrk_good = set_mdst_chrg_id.size();

	  // Cut
	  // DON'T DO THIS: it's a bad idea!!!	  
	  //if ( p4_Y5S.m() > 10.2 && p4_Y5S.m() < 11.5 && m_evt_Ntrk == 4) {

	  if ( p4_Y5S.m() > 10.2 && p4_Y5S.m() < 11.5 && m_NChargedTracks == 4) {

	  // This is a better (but NOT perfect!!!) "exactly four-tracks cut":
	  //	  if ( p4_Y5S.m() > 10.2 && p4_Y5S.m() < 11.5 && pic_list.size() == 4 ) {

	  // IMPORTANT NOTE: There is an additional cut requiring deltaE < 50 MeV for Y5S candidate in COM frame with
	  // muons fit to Y1S. This occurs in prepareCandidates()
	  // DELTA E CUT IS 
	  //if ( m_evt_Ntrk == 4) { // additionally constrain the number of allowed reconstructed charged tracks
	  
	    
	    // create Y5S
	    Particle Y5S(p4_Y5S, Ptype(9000553));
	    Y5S.relation().append(Wbj);
	    Y5S.relation().append(gam);
	    
	    Y5S_list.push_back(Y5S);
	  }
	}
      }
    }
  }



  void Analysis::prepareCandidates() {
    //cout << "prepareCandidates" << endl;

    candidatesRec.clear();
    
    int n_mctag = 0;
    int n_cand = 0;

    //    int index_candidate = 0;
    //    int index_candidate_mct = 0;

    for ( std::vector<Particle>::iterator it_Y5S = Y5S_list.begin(); it_Y5S !=  Y5S_list.end(); it_Y5S++ ) {
      
      //      index_candidate++;

      Particle &Y5S = (*it_Y5S);
      
      //Retrieve all daughter particles
      Particle Wbj(Y5S.relation().child(0));
      Particle gam(Y5S.relation().child(1));
      
      Particle YnS(Wbj.relation().child(0));
      Particle rho(Wbj.relation().child(1));
      Particle mup(YnS.relation().child(0));
      Particle mum(YnS.relation().child(1));
      Particle pip(rho.relation().child(0));
      Particle pim(rho.relation().child(1));
      
      candidateRec mycand;

      HepLorentzVector p4_gam = gam.p();
      HepLorentzVector p4_mup = mup.p();
      HepLorentzVector p4_mum = mum.p();
      HepLorentzVector p4_pip = pip.p();
      HepLorentzVector p4_pim = pim.p();
      
      HepLorentzVector p4_YnS = YnS.p();
      HepLorentzVector p4_rho = rho.p();
      HepLorentzVector p4_Wbj = Wbj.p();
      HepLorentzVector p4_Y5S = Y5S.p();

      HepLorentzVector p4_YnS_tmp = p4_YnS; // temp used in mass fitting

      double YnS_m = p4_YnS.m();

      if ( YnS_m <= 9.75 ) {
	m_ptype_param = "UPS1";
	m_resonance_fit = 1;
      }
      else if ( YnS_m > 9.75 && YnS_m <= 10.2 ) {
	m_ptype_param = "UPS2";
	m_resonance_fit = 2;
      }
      else {
	m_ptype_param = "UPS3";
	m_resonance_fit = 3;
      }
      
      // Mass fitting
      Particle YnS_fit( p4_YnS_tmp, Ptype(m_ptype_param)) ;        
      YnS_fit.relation().append(mup);
      YnS_fit.relation().append(mum);

      mycand.chi2_YnS = massFit2(YnS_fit, 2);
      HepLorentzVector p4_YnS_fit = YnS_fit.p();
      
      // Boost gamma to Y(5S) CM reference frame
      Hep3Vector p_boost = p4_Y5S.boostVector();
      //      cout << "-D-NC-p_boost = " << p_boost << endl;
      
      HepLorentzVector p4_gam_boost = p4_gam;
      p4_gam_boost.boost(-1.0*p_boost);     

      // Calculate missing p4
      HepLorentzVector p4_gam_recoil_1 = m_p4_mc_Y5S - (m_p4_mc_gam);
      HepLorentzVector p4_gam_recoil_2 = m_p4_mc_Y5S - (p4_gam);
      HepLorentzVector p4_gam_recoil_3 = m_pExp - (m_p4_mc_gam);
      HepLorentzVector p4_gam_recoil = m_pExp - (p4_gam);
      HepLorentzVector p4_pip_pim_recoil = m_pExp - (p4_pip + p4_pim);
      HepLorentzVector p4_mup_mum_recoil = m_pExp - (p4_mup + p4_mum);
      HepLorentzVector p4_pip_pim_gam_recoil = m_pExp - (p4_gam + p4_pip + p4_pim);
      HepLorentzVector p4_mup_mum_gam_recoil = m_pExp - (p4_gam + p4_mup + p4_mum);
      HepLorentzVector p4_mup_mum_gam_fit_recoil = m_pExp - (p4_gam + p4_YnS_fit);
      HepLorentzVector p4_pip_pim_mup_mum_recoil = m_pExp - (p4_mup + p4_mum + p4_pip + p4_pim);
      HepLorentzVector p4_pip_pim_mup_mum_fit_recoil = m_pExp - (p4_YnS_fit + p4_pip + p4_pim);

      HepLorentzVector p4_Wbj_fit = p4_Wbj - p4_YnS + p4_YnS_fit;
      HepLorentzVector p4_Y5S_fit = p4_Y5S - p4_YnS + p4_YnS_fit;

      HepLorentzVector p4_YnS_boost = p4_YnS;
      HepLorentzVector p4_Y5S_boost = p4_Y5S;
      HepLorentzVector p4_YnS_fit_boost = p4_YnS_fit;
      HepLorentzVector p4_Wbj_fit_boost = p4_Wbj_fit;
      HepLorentzVector p4_Y5S_fit_boost = p4_Y5S_fit;
      //HepLorentzVector p4_gam_recoil_boost = p4_gam_recoil;
      HepLorentzVector p4_pip_pim_gam_recoil_boost = p4_pip_pim_gam_recoil;
      HepLorentzVector p4_mup_mum_gam_recoil_boost = p4_mup_mum_gam_recoil;	  
      HepLorentzVector p4_pip_pim_mup_mum_recoil_boost = p4_pip_pim_mup_mum_recoil;
      
      p4_YnS_boost.boost(-1.0*p_boost);
      p4_Y5S_boost.boost(-1.0*p_boost);
      p4_YnS_fit_boost.boost(-1.0*p_boost);
      p4_Wbj_fit_boost.boost(-1.0*p_boost);
      p4_Y5S_fit_boost.boost(-1.0*p_boost);
      //p4_gam_recoil_boost.boost(-1.0*p_boost);
      p4_pip_pim_gam_recoil_boost.boost(-1.0*p_boost);
      p4_mup_mum_gam_recoil_boost.boost(-1.0*p_boost);
      p4_pip_pim_mup_mum_recoil_boost.boost(-1.0*p_boost);	  
      
      HepLorentzVector p4_mup_mum_gam = p4_mup + p4_mum + p4_gam;
      HepLorentzVector p4_pip_pim_gam = p4_pip + p4_pim + p4_gam;

      Hep3Vector p_boost_Y1D = p4_mup_mum_gam.boostVector();
      
      double mup_mum_recoil_m2 = p4_mup_mum_recoil.m2();
      double pip_pim_gam_recoil_m2 = p4_pip_pim_gam_recoil.m2();
      double mup_mum_gam_recoil_m2 = p4_mup_mum_gam_recoil.m2();
      double mup_mum_gam_fit_recoil_m2 = p4_mup_mum_gam_fit_recoil.m2();
      double pip_pim_gam_recoil_boost_m2 = p4_pip_pim_gam_recoil.m2();
      double mup_mum_gam_recoil_boost_m2 = p4_mup_mum_gam_recoil.m2();
      double pip_pim_mup_mum_recoil_m2 = p4_pip_pim_mup_mum_recoil.m2();
      double pip_pim_mup_mum_fit_recoil_m2 = p4_pip_pim_mup_mum_fit_recoil.m2();

      // Adjust the energy of the boosted photon, leaving momentum unchanged
      HepLorentzVector p4_gam_boost_adjust = p4_gam; 
      p4_gam_boost_adjust.setE( p4_pip_pim_mup_mum_fit_recoil.e() );
      HepLorentzVector p4_gam_boost_adjust_recoil = m_pExp - (p4_gam_boost_adjust);

 
      double deltaE = p4_Y5S_fit_boost.e() - m_eBeamCM;

      HepLorentzVector deltaP4_Y1D = p4_Y5S_fit - m_pExp;
      deltaP4_Y1D.boost(-1.0*p_boost_Y1D);
      double deltaE_Y1D = deltaP4_Y1D.e();

      // angle between charged tracks and photon candidate
      // ---
      /*
      Hep3Vector p_boost_mupmum = p4_YnS.boostVector(); 
      HepLorentzVector p4_gam_boost_mupmum = p4_gam;
      p4_gam_boost_mupmum.boost(-1.0*p_boost_mupmum);
      HepLorentzVector p4_mup_boost_mupmum = p4_mup;
      p4_mup_boost_mupmum.boost(-1.0*p_boost_mupmum);
      HepLorentzVector p4_mum_boost_mupmum = p4_mum;
      p4_mum_boost_mupmum.boost(-1.0*p_boost_mupmum);
      */
      mycand.angle_bt_mup_gam = p4_mup.angle( p4_gam );
      mycand.angle_bt_mum_gam = p4_mum.angle( p4_gam );
      /*
      Hep3Vector p_boost_pippim = p4_rho.boostVector(); 
      HepLorentzVector p4_gam_boost_pippim = p4_gam;
      p4_gam_boost_pippim.boost(-1.0*p_boost_pippim);
      HepLorentzVector p4_pip_boost_pippim = p4_pip;
      p4_pip_boost_pippim.boost(-1.0*p_boost_pippim);
      HepLorentzVector p4_pim_boost_pippim = p4_pim;
      p4_pim_boost_pippim.boost(-1.0*p_boost_pippim);
      */
      mycand.angle_bt_pip_gam = p4_pip.angle( p4_gam );
      mycand.angle_bt_pim_gam = p4_pim.angle( p4_gam );

      // ---

      mycand.angle_bt_Wbj_gam = p4_Wbj_fit_boost.angle( p4_gam_boost.vect() );

      mycand.gam_boost_adjust_e = p4_gam_boost_adjust.e();
      mycand.gam_boost_adjust_recoil_m = p4_gam_boost_adjust_recoil.m();

      mycand.Y5S_e = p4_Y5S.e();
      mycand.Y5S_p = p4_Y5S.rho();
      mycand.Y5S_m = p4_Y5S.m();
      mycand.Y5S_pt = p4_Y5S.perp();
      mycand.Y5S_phi = p4_Y5S.phi();
      mycand.Y5S_costh = p4_Y5S.cosTheta();

      mycand.Y5S_fit_e = p4_Y5S_fit.e();
      mycand.Y5S_fit_p = p4_Y5S_fit.rho();
      mycand.Y5S_fit_m = p4_Y5S_fit.m();
      mycand.Y5S_fit_pt = p4_Y5S_fit.perp();
      mycand.Y5S_fit_phi = p4_Y5S_fit.phi();
      mycand.Y5S_fit_costh = p4_Y5S_fit.cosTheta();

      mycand.Y5S_boost_e = p4_Y5S_boost.e();
      mycand.Y5S_boost_p = p4_Y5S_boost.rho();
      mycand.Y5S_boost_m = p4_Y5S_boost.m();
      mycand.Y5S_boost_pt = p4_Y5S_boost.perp();
      mycand.Y5S_boost_phi = p4_Y5S_boost.phi();
      mycand.Y5S_boost_costh = p4_Y5S_boost.cosTheta();

      mycand.Y5S_fit_boost_e = p4_Y5S_fit_boost.e();
      mycand.Y5S_fit_boost_p = p4_Y5S_fit_boost.rho();
      mycand.Y5S_fit_boost_m = p4_Y5S_fit_boost.m();
      mycand.Y5S_fit_boost_pt = p4_Y5S_fit_boost.perp();
      mycand.Y5S_fit_boost_phi = p4_Y5S_fit_boost.phi();
      mycand.Y5S_fit_boost_costh = p4_Y5S_fit_boost.cosTheta();
      
      mycand.Wbj_e = p4_Wbj.e();
      mycand.Wbj_p = p4_Wbj.rho();
      mycand.Wbj_m = p4_Wbj.m();
      mycand.Wbj_pt = p4_Wbj.perp();
      mycand.Wbj_phi = p4_Wbj.phi();
      mycand.Wbj_costh = p4_Wbj.cosTheta();

      mycand.Wbj_fit_e = p4_Wbj_fit.e();
      mycand.Wbj_fit_p = p4_Wbj_fit.rho();
      mycand.Wbj_fit_m = p4_Wbj_fit.m();
      mycand.Wbj_fit_pt = p4_Wbj_fit.perp();
      mycand.Wbj_fit_phi = p4_Wbj_fit.phi();
      mycand.Wbj_fit_costh = p4_Wbj_fit.cosTheta();

      mycand.Wbj_fit_boost_e = p4_Wbj_fit_boost.e();
      mycand.Wbj_fit_boost_p = p4_Wbj_fit_boost.rho();
      mycand.Wbj_fit_boost_m = p4_Wbj_fit_boost.m();
      mycand.Wbj_fit_boost_pt = p4_Wbj_fit_boost.perp();
      mycand.Wbj_fit_boost_phi = p4_Wbj_fit_boost.phi();
      mycand.Wbj_fit_boost_costh = p4_Wbj_fit_boost.cosTheta();

      mycand.gam_recoil_e_1 = p4_gam_recoil_1.e();
      mycand.gam_recoil_p_1 = p4_gam_recoil_1.rho();
      mycand.gam_recoil_m_1 = p4_gam_recoil_1.m();
      mycand.gam_recoil_pt_1 = p4_gam_recoil_1.perp();
      mycand.gam_recoil_phi_1 = p4_gam_recoil_1.phi();
      mycand.gam_recoil_costh_1 = p4_gam_recoil_1.cosTheta();

      mycand.gam_recoil_e_2 = p4_gam_recoil_2.e();
      mycand.gam_recoil_p_2 = p4_gam_recoil_2.rho();
      mycand.gam_recoil_m_2 = p4_gam_recoil_2.m();
      mycand.gam_recoil_pt_2 = p4_gam_recoil_2.perp();
      mycand.gam_recoil_phi_2 = p4_gam_recoil_2.phi();
      mycand.gam_recoil_costh_2 = p4_gam_recoil_2.cosTheta();

      mycand.gam_recoil_e_3 = p4_gam_recoil_3.e();
      mycand.gam_recoil_p_3 = p4_gam_recoil_3.rho();
      mycand.gam_recoil_m_3 = p4_gam_recoil_3.m();
      mycand.gam_recoil_pt_3 = p4_gam_recoil_3.perp();
      mycand.gam_recoil_phi_3 = p4_gam_recoil_3.phi();
      mycand.gam_recoil_costh_3 = p4_gam_recoil_3.cosTheta();

      mycand.gam_recoil_e = p4_gam_recoil.e();
      mycand.gam_recoil_p = p4_gam_recoil.rho();
      mycand.gam_recoil_m = p4_gam_recoil.m();
      mycand.gam_recoil_pt = p4_gam_recoil.perp();
      mycand.gam_recoil_phi = p4_gam_recoil.phi();
      mycand.gam_recoil_costh = p4_gam_recoil.cosTheta();

      //mycand.gam_recoil_boost_m = p4_gam_recoil_boost.m();

      mycand.YnS_e = p4_YnS.e();
      mycand.YnS_p = p4_YnS.rho();
      mycand.YnS_m = p4_YnS.m();
      mycand.YnS_pt = p4_YnS.perp();
      mycand.YnS_phi = p4_YnS.phi();
      mycand.YnS_costh = p4_YnS.cosTheta();

      mycand.YnS_boost_e = p4_YnS_boost.e();
      mycand.YnS_boost_p = p4_YnS_boost.rho();
      mycand.YnS_boost_m = p4_YnS_boost.m();
      mycand.YnS_boost_pt = p4_YnS_boost.perp();
      mycand.YnS_boost_phi = p4_YnS_boost.phi();
      mycand.YnS_boost_costh = p4_YnS_boost.cosTheta();

      mycand.YnS_fit_e = p4_YnS_fit.e();
      mycand.YnS_fit_p = p4_YnS_fit.rho();
      mycand.YnS_fit_m = p4_YnS_fit.m();
      mycand.YnS_fit_pt = p4_YnS_fit.perp();
      mycand.YnS_fit_phi = p4_YnS_fit.phi();
      mycand.YnS_fit_costh = p4_YnS_fit.cosTheta();
     
      mycand.rho_e = p4_rho.e();
      mycand.rho_p = p4_rho.rho();
      mycand.rho_m = p4_rho.m();
      mycand.rho_pt = p4_rho.perp();
      mycand.rho_phi = p4_rho.phi();
      mycand.rho_costh = p4_rho.cosTheta();

      mycand.gam_e = p4_gam.e();
      mycand.gam_p = p4_gam.rho();
      mycand.gam_m = p4_gam.m();
      mycand.gam_pt = p4_gam.perp();
      mycand.gam_phi = p4_gam.phi();
      mycand.gam_costh = p4_gam.cosTheta();

      Mdst_gamma const& mdst_gam = gam.mdstGamma();
      Mdst_ecl const& mdst_shower = mdst_gam.ecl();
      
      Mdst_ecl_aux_Manager &auxShowerMgr = Mdst_ecl_aux_Manager::get_manager();
      Mdst_ecl_aux const& auxShowerInfo = auxShowerMgr(mdst_shower.get_ID());
	
      int gam_nhit = auxShowerInfo.nhits();
      double gam_width = auxShowerInfo.width();
      double gam_e9oe25 = auxShowerInfo.e9oe25();

      //????????????????????????

      mycand.gam_nhit = gam_nhit;
      mycand.gam_width = gam_width;
      mycand.gam_e9oe25 = gam_e9oe25;

      mycand.mup_e = p4_mup.e();
      mycand.mup_p = p4_mup.rho();
      mycand.mup_m = p4_mup.m();
      mycand.mup_pt = p4_mup.perp();
      mycand.mup_phi = p4_mup.phi();
      mycand.mup_costh = p4_mup.cosTheta();
      
      mycand.mum_e = p4_mum.e();
      mycand.mum_p = p4_mum.rho();
      mycand.mum_m = p4_mum.m();
      mycand.mum_pt = p4_mum.perp();
      mycand.mum_phi = p4_mum.phi();
      mycand.mum_costh = p4_mum.cosTheta();
      
      mycand.pip_e = p4_pip.e();
      mycand.pip_p = p4_pip.rho();
      mycand.pip_m = p4_pip.m();
      mycand.pip_pt = p4_pip.perp();
      mycand.pip_phi = p4_pip.phi();
      mycand.pip_costh = p4_pip.cosTheta();

      mycand.pim_e = p4_pim.e();
      mycand.pim_p = p4_pim.rho();
      mycand.pim_m = p4_pim.m();
      mycand.pim_pt = p4_pim.perp();
      mycand.pim_phi = p4_pim.phi();
      mycand.pim_costh = p4_pim.cosTheta();
      
      mycand.pip_pim_recoil_m         = p4_pip_pim_recoil.m();
      mycand.pip_pim_recoil_e         = p4_pip_pim_recoil.e();

      mycand.pip_pim_gam_recoil_e     = p4_pip_pim_gam_recoil.e();
      mycand.pip_pim_gam_recoil_p     = p4_pip_pim_gam_recoil.rho();
      mycand.pip_pim_gam_recoil_m     = (pip_pim_gam_recoil_m2 > 0.) ? sqrt(pip_pim_gam_recoil_m2) : -sqrt(-pip_pim_gam_recoil_m2); 
      mycand.pip_pim_gam_recoil_pt    = p4_pip_pim_gam_recoil.perp();
      mycand.pip_pim_gam_recoil_phi   = p4_pip_pim_gam_recoil.phi();
      mycand.pip_pim_gam_recoil_costh = p4_pip_pim_gam_recoil.cosTheta();
      
      mycand.pip_pim_gam_recoil_boost_e     = p4_pip_pim_gam_recoil_boost.e();
      mycand.pip_pim_gam_recoil_boost_p     = p4_pip_pim_gam_recoil_boost.rho();
      mycand.pip_pim_gam_recoil_boost_m     = (pip_pim_gam_recoil_boost_m2 > 0.) ? sqrt(pip_pim_gam_recoil_boost_m2) : -sqrt(-pip_pim_gam_recoil_boost_m2); 
      mycand.pip_pim_gam_recoil_boost_pt    = p4_pip_pim_gam_recoil_boost.perp();
      mycand.pip_pim_gam_recoil_boost_phi   = p4_pip_pim_gam_recoil_boost.phi();
      mycand.pip_pim_gam_recoil_boost_costh = p4_pip_pim_gam_recoil_boost.cosTheta();

      //mycand.mup_mum_recoil_m         = p4_mup_mum_recoil.m();
      mycand.mup_mum_recoil_m         = (p4_mup_mum_recoil.m2() > 0.) ? sqrt(mup_mum_recoil_m2) : -sqrt(-mup_mum_recoil_m2);
      mycand.mup_mum_recoil_e         = p4_mup_mum_recoil.e();

      mycand.mup_mum_gam_recoil_e     = p4_mup_mum_gam_recoil.e();
      mycand.mup_mum_gam_recoil_p     = p4_mup_mum_gam_recoil.rho();
      mycand.mup_mum_gam_recoil_m     = (mup_mum_gam_recoil_m2 > 0.) ? sqrt(mup_mum_gam_recoil_m2) : -sqrt(-mup_mum_gam_recoil_m2); 
      mycand.mup_mum_gam_recoil_pt    = p4_mup_mum_gam_recoil.perp();
      mycand.mup_mum_gam_recoil_phi   = p4_mup_mum_gam_recoil.phi();
      mycand.mup_mum_gam_recoil_costh = p4_mup_mum_gam_recoil.cosTheta();

      mycand.mup_mum_gam_fit_recoil_m     = (mup_mum_gam_fit_recoil_m2 > 0.) ? sqrt(mup_mum_gam_fit_recoil_m2) : -sqrt(-mup_mum_gam_fit_recoil_m2);

      mycand.mup_mum_gam_recoil_boost_e     = p4_mup_mum_gam_recoil_boost.e();
      mycand.mup_mum_gam_recoil_boost_p     = p4_mup_mum_gam_recoil_boost.rho();
      mycand.mup_mum_gam_recoil_boost_m     = (mup_mum_gam_recoil_boost_m2 > 0.) ? sqrt(mup_mum_gam_recoil_boost_m2) : -sqrt(-mup_mum_gam_recoil_boost_m2); 
      mycand.mup_mum_gam_recoil_boost_pt    = p4_mup_mum_gam_recoil_boost.perp();
      mycand.mup_mum_gam_recoil_boost_phi   = p4_mup_mum_gam_recoil_boost.phi();
      mycand.mup_mum_gam_recoil_boost_costh = p4_mup_mum_gam_recoil_boost.cosTheta();

      mycand.pip_pim_mup_mum_recoil_e     = p4_pip_pim_mup_mum_recoil.e();
      mycand.pip_pim_mup_mum_recoil_p     = p4_pip_pim_mup_mum_recoil.rho();
      mycand.pip_pim_mup_mum_recoil_m     = (pip_pim_mup_mum_recoil_m2 > 0.) ? sqrt(pip_pim_mup_mum_recoil_m2) : -sqrt(-pip_pim_mup_mum_recoil_m2);
      mycand.pip_pim_mup_mum_recoil_pt    = p4_pip_pim_mup_mum_recoil.perp();
      mycand.pip_pim_mup_mum_recoil_phi   = p4_pip_pim_mup_mum_recoil.phi();
      mycand.pip_pim_mup_mum_recoil_costh = p4_pip_pim_mup_mum_recoil.cosTheta();

      mycand.pip_pim_mup_mum_recoil_m2    = p4_pip_pim_mup_mum_recoil.m2();

      mycand.pip_pim_mup_mum_fit_recoil_m  = (pip_pim_mup_mum_fit_recoil_m2 > 0.) ? sqrt(pip_pim_mup_mum_fit_recoil_m2) : -sqrt(-pip_pim_mup_mum_fit_recoil_m2);
      mycand.pip_pim_mup_mum_fit_recoil_m2 = pip_pim_mup_mum_fit_recoil_m2;
      mycand.pip_pim_mup_mum_fit_recoil_e  = p4_pip_pim_mup_mum_fit_recoil.e();
      mycand.pip_pim_mup_mum_fit_recoil_costh  = p4_pip_pim_mup_mum_fit_recoil.cosTheta();

      mycand.pip_pim_mup_mum_recoil_boost_e     = p4_pip_pim_mup_mum_recoil_boost.e();
      mycand.pip_pim_mup_mum_recoil_boost_p     = p4_pip_pim_mup_mum_recoil_boost.rho();
      mycand.pip_pim_mup_mum_recoil_boost_m     = p4_pip_pim_mup_mum_recoil_boost.m();
      mycand.pip_pim_mup_mum_recoil_boost_pt    = p4_pip_pim_mup_mum_recoil_boost.perp();
      mycand.pip_pim_mup_mum_recoil_boost_phi   = p4_pip_pim_mup_mum_recoil_boost.phi();
      mycand.pip_pim_mup_mum_recoil_boost_costh = p4_pip_pim_mup_mum_recoil_boost.cosTheta();
      
      mycand.gam_boost_e = p4_gam_boost.e();
      mycand.gam_boost_p = p4_gam_boost.rho();
      mycand.gam_boost_m = p4_gam_boost.m();
      mycand.gam_boost_pt = p4_gam_boost.perp();
      mycand.gam_boost_phi = p4_gam_boost.phi();
      mycand.gam_boost_costh = p4_gam_boost.cosTheta();
     
      mycand.deltaE = deltaE;
      mycand.deltaE_Y1D = deltaE_Y1D;

      mycand.mup_mum_gam_m = p4_mup_mum_gam.m();
      mycand.mup_mum_gam_e = p4_mup_mum_gam.e();     
      mycand.pip_pim_gam_m = p4_pip_pim_gam.m(); 

      // @TODO: gam_id and gam_parent_id
      /*


      */

      if ( isMC ) {
	Mdst_gamma const &mdst_gam = Y5S.relation().child(1).mdstGamma();

	Mdst_charged const &mdst_mup  = YnS.relation().child(0).mdstCharged();
	Mdst_charged const &mdst_mum  = YnS.relation().child(1).mdstCharged();
	
	Mdst_charged const &mdst_pip = rho.relation().child(0).mdstCharged();
	Mdst_charged const &mdst_pim = rho.relation().child(1).mdstCharged();

	const Gen_hepevt& hep_gam(gen_level(get_hepevt(mdst_gam)));
	const Gen_hepevt& hep_mup(gen_level(get_hepevt(mdst_mup)));
	const Gen_hepevt& hep_mum(gen_level(get_hepevt(mdst_mum)));
	const Gen_hepevt& hep_pip(gen_level(get_hepevt(mdst_pip)));
	const Gen_hepevt& hep_pim(gen_level(get_hepevt(mdst_pim)));

	mycand.gam_id = mdst_gam.get_ID();    
	mycand.pip_id = mdst_pip.get_ID(); 
	mycand.pim_id = mdst_pim.get_ID(); 
	mycand.mup_id = mdst_mup.get_ID(); 
	mycand.mum_id = mdst_mum.get_ID();

	mycand.gam_from_pi0 = vetoPi0(mycand.gam_id);

	leftover_charged_list.clear();
	leftover_neutral_list.clear();

	Mdst_charged_Manager &ChargedMgr = Mdst_charged_Manager::get_manager();
	for ( std::vector<Mdst_charged>::iterator it = ChargedMgr.begin(); it != ChargedMgr.end(); it++ ) {

	  Mdst_charged &chrg = *it;
	  int chrg_id = chrg.get_ID();

	  if ( chrg_id != mycand.pip_id &&
	       chrg_id != mycand.pim_id &&
	       chrg_id != mycand.mup_id &&
	       chrg_id != mycand.mum_id    ) {

	    // the program crashes if I use the default constructor, so I just call all of them pions--this means nothing, though.
	    Particle leftover_chrg( chrg, Ptype( chrg.charge() > 0 ? "PI+" : "PI-" ) ); 
	    leftover_charged_list.push_back(leftover_chrg);
	  }
	}

	Mdst_gamma_Manager &GammaMgr = Mdst_gamma_Manager::get_manager();
	for( std::vector<Mdst_gamma>::iterator it = GammaMgr.begin(); it != GammaMgr.end(); it++ ) {

	  Mdst_gamma &mdst_gamma = *it;
	  int gamma_id = mdst_gamma.get_ID();

	  
	  if ( gamma_id != mycand.gam_id) {

	    Particle leftover_gam( mdst_gamma );
	    leftover_neutral_list.push_back(leftover_gam);
	  }
	}

	mycand.neutral_e = getNeutralEnergy();
	mycand.charged_e = getChargedEnergy();

	//cout<<"SIZES: " << leftover_charged_list.size() << "," << leftover_neutral_list.size() << endl;

	/*
	  for ( int i = 0; i < pi0.nChildren(); i++) {
	  Particle pi0_daughter(pi0.relation().child(i));
	  }
	*/
      }
      
      // candidate is not MC tagged (yet)
      mycand.mct = 0;

      //MC TAGGING
      if ( isSignalMC ) {
	//cout << "MC Tagging" << endl;

	Mdst_gamma const &mdst_gam = Y5S.relation().child(1).mdstGamma();
	
	Mdst_charged const &mdst_mup  = YnS.relation().child(0).mdstCharged();
	Mdst_charged const &mdst_mum  = YnS.relation().child(1).mdstCharged();
	
	Mdst_charged const &mdst_pip = rho.relation().child(0).mdstCharged();
	Mdst_charged const &mdst_pim = rho.relation().child(1).mdstCharged();

	const Gen_hepevt& hep_gam(gen_level(get_hepevt(mdst_gam)));
	const Gen_hepevt& hep_mup(gen_level(get_hepevt(mdst_mup)));
	const Gen_hepevt& hep_mum(gen_level(get_hepevt(mdst_mum)));
	const Gen_hepevt& hep_pip(gen_level(get_hepevt(mdst_pip)));
	const Gen_hepevt& hep_pim(gen_level(get_hepevt(mdst_pim)));
	
	mycand.gam_mct = 0;

	if ( hep_gam ) {

	  if( hep_gam.get_ID() == m_mc_gam_id ) {
	    if ( hep_gam.mother() ) {
	      if ( hep_gam.mother().get_ID() == m_mc_Y5S_id  ) {
	    		mycand.gam_mct = 1;
	      }
	    }
	  }
	}

	mycand.mup_mct = 0;
	
	if(  hep_mup.get_ID() == m_mc_mup_id && 
	     hep_mup.mother().get_ID() == m_mc_YnS_id &&
	     hep_mup.mother().mother().get_ID() == m_mc_Wbj_id &&
	     hep_mup.mother().mother().mother().get_ID() == m_mc_Y5S_id  ) {
	  
	  mycand.mup_mct = 1;
	}
	
	mycand.mum_mct = 0;
	
	if(  hep_mum.get_ID() == m_mc_mum_id && 
	     hep_mum.mother().get_ID() == m_mc_YnS_id &&
	     hep_mum.mother().mother().get_ID() == m_mc_Wbj_id &&
	     hep_mum.mother().mother().mother().get_ID() == m_mc_Y5S_id  ) {

	  mycand.mum_mct = 1;
	}	     
	
	mycand.pip_mct = 0;
	
	if(  hep_pip.get_ID() == m_mc_pip_id && 
	     hep_pip.mother().get_ID() == m_mc_rho_id &&
	     hep_pip.mother().mother().get_ID() == m_mc_Wbj_id &&
	     hep_pip.mother().mother().mother().get_ID() == m_mc_Y5S_id  ) {
	  
	  mycand.pip_mct = 1;
	}	  
	
	mycand.pim_mct = 0;
	
	if(  hep_pim.get_ID() == m_mc_pim_id && 
	     hep_pim.mother().get_ID() == m_mc_rho_id &&
	     hep_pim.mother().mother().get_ID() == m_mc_Wbj_id &&
	     hep_pim.mother().mother().mother().get_ID() == m_mc_Y5S_id  ) {
	  
	  mycand.pim_mct = 1;
	}	   	
	
	if ( mycand.gam_mct   == 1 &&
	     mycand.mup_mct   == 1 &&
	     mycand.mum_mct   == 1 && 
	     mycand.pip_mct   == 1 && 
	     mycand.pim_mct   == 1     ) {

	  //	  index_candidate_mct++;

	  mycand.mct = 1;
	  //	      cout <<"TAGGED"<< endl;
	}

	/*	      
	  cout << "VS-D-MC tagging: eta_g1_mct = " << mycand.eta_g1_mct << endl;
	  cout << "VS-D-MC tagging: eta_g2_mct = " << mycand.eta_g2_mct << endl;
	  cout << "VS-D-MC tagging: ks_pip_mct = " << mycand.ks_pip_mct << endl;
	  cout << "VS-D-MC tagging: ks_pim_mct = " << mycand.ks_pim_mct << endl;
	  cout << "VS-D-MC tagging: ks_mct     = " << mycand.pip_mct    << endl;
	  cout << "VS-D-MC tagging: ks_mct     = " << mycand.pim_mct    << endl;
	  cout << "VS-D-MC tagging: total mct = "  << mycand.mct        << endl;
	*/	      
	
      }
      
      // Nov. 29, 2017:  blinding the signal region (see e-mails exchanged with Nicholas on Nov. 29, 2017)

      bool pass_cut_1 = mycand.Wbj_fit_m >= 10.49 && mycand.Wbj_fit_m <= 10.72;
      bool pass_cut_2 = mycand.Wbj_fit_m <= mycand.gam_recoil_m + 0.04;
      bool pass_cut_3 = mycand.gam_recoil_m <= 10.72;

      mycand.in_signal_region = 0;

      if ( pass_cut_1 && pass_cut_2 && pass_cut_3 ) mycand.in_signal_region = 1;

      n_cand++;
      //	  cout << "n_cand = " << n_cand  << endl;
      // place signal event candidates into the list of candidates
      candidatesRec.push_back(mycand);
           
    }
  }
  
 
  void Analysis::prepareTreeRec() {
    //cout << "prepareTreeRec" << endl;


    nt_angle_bt_pip_gam.clear();
    nt_angle_bt_pim_gam.clear();
    nt_angle_bt_mup_gam.clear();
    nt_angle_bt_mum_gam.clear();




    m_ncand = candidatesRec.size(); 

    m_ncand_sr = 0; 

    set_mdst_chrg_id.clear();

    nt_gam_nhit.clear();
    nt_gam_width.clear();
    nt_gam_e9oe25.clear();

    nt_gam_id.clear();
    nt_gam_parent_id.clear();

    nt_index.clear();
    nt_gam_from_pi0.clear();

    nt_neutral_e.clear();
    nt_charged_e.clear();

    nt_gam_boost_adjust_e.clear();
    nt_gam_boost_adjust_recoil_m.clear();

    // Rec values
    nt_gam_mct.clear();
    nt_mup_mct.clear();
    nt_mum_mct.clear();
    nt_pip_mct.clear();
    nt_pim_mct.clear();
    nt_mct.clear();

    nt_sr.clear();

    nt_angle_bt_Wbj_gam.clear();

    nt_Y5S_e.clear();
    nt_Y5S_p.clear();
    nt_Y5S_m.clear();
    nt_Y5S_pt.clear();
    nt_Y5S_phi.clear();
    nt_Y5S_costh.clear();

    nt_Y5S_fit_e.clear();
    nt_Y5S_fit_p.clear();
    nt_Y5S_fit_m.clear();
    nt_Y5S_fit_pt.clear();
    nt_Y5S_fit_phi.clear();
    nt_Y5S_fit_costh.clear();

    nt_Y5S_boost_e.clear();
    nt_Y5S_boost_p.clear();
    nt_Y5S_boost_m.clear();
    nt_Y5S_boost_pt.clear();
    nt_Y5S_boost_phi.clear();
    nt_Y5S_boost_costh.clear();

    nt_Y5S_fit_boost_e.clear();
    nt_Y5S_fit_boost_p.clear();
    nt_Y5S_fit_boost_m.clear();
    nt_Y5S_fit_boost_pt.clear();
    nt_Y5S_fit_boost_phi.clear();
    nt_Y5S_fit_boost_costh.clear();
    
    nt_Wbj_e.clear();
    nt_Wbj_p.clear();
    nt_Wbj_m.clear();
    nt_Wbj_pt.clear();
    nt_Wbj_phi.clear();
    nt_Wbj_costh.clear();

    nt_Wbj_fit_e.clear();
    nt_Wbj_fit_p.clear();
    nt_Wbj_fit_m.clear();
    nt_Wbj_fit_pt.clear();
    nt_Wbj_fit_phi.clear();
    nt_Wbj_fit_costh.clear();

    nt_Wbj_fit_boost_e.clear();
    nt_Wbj_fit_boost_p.clear();
    nt_Wbj_fit_boost_m.clear();
    nt_Wbj_fit_boost_pt.clear();
    nt_Wbj_fit_boost_phi.clear();
    nt_Wbj_fit_boost_costh.clear();

    nt_gam_recoil_e_1.clear();
    nt_gam_recoil_p_1.clear();
    nt_gam_recoil_m_1.clear();
    nt_gam_recoil_pt_1.clear();
    nt_gam_recoil_phi_1.clear();
    nt_gam_recoil_costh_1.clear();

    nt_gam_recoil_e_2.clear();
    nt_gam_recoil_p_2.clear();
    nt_gam_recoil_m_2.clear();
    nt_gam_recoil_pt_2.clear();
    nt_gam_recoil_phi_2.clear();
    nt_gam_recoil_costh_2.clear();

    nt_gam_recoil_e_3.clear();
    nt_gam_recoil_p_3.clear();
    nt_gam_recoil_m_3.clear();
    nt_gam_recoil_pt_3.clear();
    nt_gam_recoil_phi_3.clear();
    nt_gam_recoil_costh_3.clear();

    nt_gam_recoil_e.clear();
    nt_gam_recoil_p.clear();
    nt_gam_recoil_m.clear();
    nt_gam_recoil_pt.clear();
    nt_gam_recoil_phi.clear();
    nt_gam_recoil_costh.clear();
    
    nt_gam_recoil_boost_m.clear();

    nt_YnS_e.clear();
    nt_YnS_p.clear();
    nt_YnS_m.clear();
    nt_YnS_pt.clear();
    nt_YnS_phi.clear();
    nt_YnS_costh.clear();

    nt_YnS_boost_e.clear();
    nt_YnS_boost_p.clear();
    nt_YnS_boost_m.clear();
    nt_YnS_boost_pt.clear();
    nt_YnS_boost_phi.clear();
    nt_YnS_boost_costh.clear();
    
    nt_YnS_fit_e.clear();
    nt_YnS_fit_p.clear();
    nt_YnS_fit_m.clear();
    nt_YnS_fit_pt.clear();
    nt_YnS_fit_phi.clear();
    nt_YnS_fit_costh.clear();

    nt_pip_pim_recoil_m.clear();
    nt_pip_pim_recoil_e.clear();
    
    nt_pip_pim_gam_recoil_e.clear();
    nt_pip_pim_gam_recoil_p.clear();
    nt_pip_pim_gam_recoil_m.clear();
    nt_pip_pim_gam_recoil_pt.clear();
    nt_pip_pim_gam_recoil_phi.clear();
    nt_pip_pim_gam_recoil_costh.clear();

    nt_pip_pim_gam_recoil_boost_e.clear();
    nt_pip_pim_gam_recoil_boost_p.clear();
    nt_pip_pim_gam_recoil_boost_m.clear();
    nt_pip_pim_gam_recoil_boost_pt.clear();
    nt_pip_pim_gam_recoil_boost_phi.clear();
    nt_pip_pim_gam_recoil_boost_costh.clear();

    nt_mup_mum_recoil_m.clear();
    nt_mup_mum_recoil_e.clear();

    nt_mup_mum_gam_recoil_e.clear();
    nt_mup_mum_gam_recoil_p.clear();
    nt_mup_mum_gam_recoil_m.clear();
    nt_mup_mum_gam_recoil_pt.clear();
    nt_mup_mum_gam_recoil_phi.clear();
    nt_mup_mum_gam_recoil_costh.clear();

    nt_mup_mum_gam_fit_recoil_m.clear();

    nt_mup_mum_gam_recoil_boost_e.clear();
    nt_mup_mum_gam_recoil_boost_p.clear();
    nt_mup_mum_gam_recoil_boost_m.clear();
    nt_mup_mum_gam_recoil_boost_pt.clear();
    nt_mup_mum_gam_recoil_boost_phi.clear();
    nt_mup_mum_gam_recoil_boost_costh.clear();

    nt_rho_e.clear();
    nt_rho_p.clear();
    nt_rho_m.clear();
    nt_rho_pt.clear();
    nt_rho_phi.clear();
    nt_rho_costh.clear();

    nt_mup_e.clear();
    nt_mup_p.clear();
    nt_mup_m.clear();
    nt_mup_pt.clear();
    nt_mup_phi.clear();
    nt_mup_costh.clear();

    nt_mum_e.clear();
    nt_mum_p.clear();
    nt_mum_m.clear();
    nt_mum_pt.clear();
    nt_mum_phi.clear();
    nt_mum_costh.clear();

    nt_pip_e.clear();
    nt_pip_p.clear();
    nt_pip_m.clear();
    nt_pip_pt.clear();
    nt_pip_phi.clear();
    nt_pip_costh.clear();

    nt_pim_e.clear();
    nt_pim_p.clear();
    nt_pim_m.clear();
    nt_pim_pt.clear();
    nt_pim_phi.clear();
    nt_pim_costh.clear();

    nt_gam_e.clear();
    nt_gam_p.clear();
    nt_gam_m.clear();
    nt_gam_pt.clear();
    nt_gam_phi.clear();
    nt_gam_costh.clear();

    nt_gam_boost_e.clear();
    nt_gam_boost_p.clear();
    nt_gam_boost_m.clear();
    nt_gam_boost_pt.clear();
    nt_gam_boost_phi.clear();
    nt_gam_boost_costh.clear();

    nt_pip_pim_mup_mum_recoil_e.clear();
    nt_pip_pim_mup_mum_recoil_p.clear();
    nt_pip_pim_mup_mum_recoil_m.clear();
    nt_pip_pim_mup_mum_recoil_pt.clear();
    nt_pip_pim_mup_mum_recoil_phi.clear();
    nt_pip_pim_mup_mum_recoil_costh.clear();

    nt_pip_pim_mup_mum_recoil_m2.clear();

    nt_pip_pim_mup_mum_recoil_boost_e.clear();
    nt_pip_pim_mup_mum_recoil_boost_p.clear();
    nt_pip_pim_mup_mum_recoil_boost_m.clear();
    nt_pip_pim_mup_mum_recoil_boost_pt.clear();
    nt_pip_pim_mup_mum_recoil_boost_phi.clear();
    nt_pip_pim_mup_mum_recoil_boost_costh.clear();

    nt_pip_pim_mup_mum_fit_recoil_m.clear();
    nt_pip_pim_mup_mum_fit_recoil_m2.clear();
    nt_pip_pim_mup_mum_fit_recoil_e.clear();
    nt_pip_pim_mup_mum_fit_recoil_costh.clear();

    // Stand-alone variables
    nt_mup_mum_gam_m.clear();
    nt_mup_mum_gam_e.clear();
    nt_pip_pim_gam_m.clear();
    nt_pip_pim_gam_e.clear();
    nt_deltaE.clear();
    nt_deltaE_Y1D.clear();
    nt_chi2_YnS.clear();
    
    m_ncand_mct = 0;

    // count the number of MC-tagged candidates (normally there should be no more than 1 but stuff happens)
    //@TODO: can't this be combined with the other for loop below? Do they loop over the same thing?
    for ( std::vector<candidateRec>::iterator it_cand = candidatesRec.begin(); it_cand !=  candidatesRec.end(); it_cand++ ) {
      //      cout << "loop over candidiates" << endl;

      candidateRec cand = *it_cand;

      if ( cand.mct == 1 ) m_ncand_mct++;

      if ( cand.in_signal_region == 1 ) m_ncand_sr++;

    }

    /*
    if ( m_ncand_mct > 1 ) {

      cout << "-D-VS-MC tagging===========================================================" << endl;

      int index = 0;

      for ( std::vector<candidateRec>::iterator it_cand = candidatesRec.begin(); it_cand !=  candidatesRec.end(); it_cand++ ) {
	//      cout << "loop over candidiates" << endl;

	index++;
	
	candidateRec cand = *it_cand;

	cout << "-D-VS-MC-tagged candidate # " << index << endl;
	cout << "-D-VS-MC-gam_id = " << cand.gam_id << endl;

	Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
	Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

	Mdst_gamma &mdst_gamma = GammaMgr[cand.gam_id-1]; 

	const Gen_hepevt& mc_particle(gen_level(get_hepevt(mdst_gamma)));

	cout << "-D-VS-MC-&mc_particle = " << &mc_particle << endl;
	
	if ( mc_particle ) {
	  cout << "-D-VS-MC tagging:  YES" << endl;
	  cout << "-D-VS-MC-mc_particle.get_ID() = " << mc_particle.get_ID() << endl;
	  cout << "-D-VS-MC-mc_particle.idhep() = " << mc_particle.idhep() << endl;
	  
	  if ( mc_particle.mother() ) {
	    cout << "-D-VS-MC tagging:  MOTHER: YES" << endl;
	    cout << "-D-VS-MC-mc_particle.mother().get_ID() = " << mc_particle.mother().get_ID() << endl;
	    cout << "-D-VS-MC-mc_particle.mother().idhep() = " << mc_particle.mother().idhep() << endl;
	  }
	  else {
	    cout << "-D-VS-MC tagging:  NO MOTHER--------------" << endl;
	  }

	}
	else {
	  cout << "-D-VS-MC tagging:  NO--------------------------------" << endl;
	}
	
      }
    }
    */

    /*
    cout << "VS-D-runseq, eventseq = " << m_run_seq << ", " << m_event_seq << endl;
    cout << "VS-D-number of charged tracks of any quality = " << pic_list.size() << endl;
    cout << "VS-D-number of high-quality positively charged pions = " << pip_list.size() << endl;
    cout << "VS-D-number of high-quality negatively charged pions = " << pim_list.size() << endl;
    cout << "VS-D-number of high-quality positively charged muons = " << mup_list.size() << endl;
    cout << "VS-D-number of high-quality negatively charged muons = " << mum_list.size() << endl;
    cout << "VS-D-number of photon candidates = " << gam_list.size() << endl;
    cout << "VS-D-number of rho->pi+pi- candidates = " << rho_list.size() << endl;
    cout << "VS-D-number of YnS->mu+mu- candidates = " << YnS_list.size() << endl;
    cout << "VS-D-number of Wbj->YnSrho candidates = " << Wbj_list.size() << endl;
    cout << "VS-D-number of Y5S candidates = " << Y5S_list.size() << endl;
    cout << "VS-D-number of signal candidates = " << candidatesRec.size() << endl;
    */

    int index = 0;
    int icand_mct = 0;
    
    for ( std::vector<candidateRec>::iterator it_cand = candidatesRec.begin(); it_cand !=  candidatesRec.end(); it_cand++ ) {
      //      cout << "loop over candidates again" << endl;
      index++;
      
      candidateRec cand = *it_cand;

      nt_angle_bt_pip_gam.push_back(cand.angle_bt_pip_gam);
      nt_angle_bt_pim_gam.push_back(cand.angle_bt_pim_gam);
      nt_angle_bt_mup_gam.push_back(cand.angle_bt_mup_gam);
      nt_angle_bt_mum_gam.push_back(cand.angle_bt_mum_gam);

      nt_gam_nhit.push_back(cand.gam_nhit);
      nt_gam_width.push_back(cand.gam_width);
      nt_gam_e9oe25.push_back(cand.gam_e9oe25);

      nt_index.push_back(index);
      nt_gam_from_pi0.push_back(cand.gam_from_pi0);

      nt_neutral_e.push_back(cand.neutral_e);
      nt_charged_e.push_back(cand.charged_e);

      nt_gam_boost_adjust_e.push_back(cand.gam_boost_adjust_e);
      nt_gam_boost_adjust_recoil_m.push_back(cand.gam_boost_adjust_recoil_m);

      // Rec values
      nt_gam_mct.push_back(cand.gam_mct);
      nt_mup_mct.push_back(cand.mup_mct);
      nt_mum_mct.push_back(cand.mum_mct);
      nt_pip_mct.push_back(cand.pip_mct);
      nt_pim_mct.push_back(cand.pim_mct);
      nt_mct.push_back(cand.mct);

      nt_sr.push_back(cand.in_signal_region);

      nt_angle_bt_Wbj_gam.push_back(cand.angle_bt_Wbj_gam);

      nt_Y5S_e.push_back(cand.Y5S_e);
      nt_Y5S_p.push_back(cand.Y5S_p);
      nt_Y5S_m.push_back(cand.Y5S_m);
      nt_Y5S_pt.push_back(cand.Y5S_pt);
      nt_Y5S_phi.push_back(cand.Y5S_phi);
      nt_Y5S_costh.push_back(cand.Y5S_costh);

      nt_Y5S_fit_e.push_back(cand.Y5S_fit_e);
      nt_Y5S_fit_p.push_back(cand.Y5S_fit_p);
      nt_Y5S_fit_m.push_back(cand.Y5S_fit_m);
      nt_Y5S_fit_pt.push_back(cand.Y5S_fit_pt);
      nt_Y5S_fit_phi.push_back(cand.Y5S_fit_phi);
      nt_Y5S_fit_costh.push_back(cand.Y5S_fit_costh);

      nt_Y5S_boost_e.push_back(cand.Y5S_boost_e);
      nt_Y5S_boost_p.push_back(cand.Y5S_boost_p);
      nt_Y5S_boost_m.push_back(cand.Y5S_boost_m);
      nt_Y5S_boost_pt.push_back(cand.Y5S_boost_pt);
      nt_Y5S_boost_phi.push_back(cand.Y5S_boost_phi);
      nt_Y5S_boost_costh.push_back(cand.Y5S_boost_costh);

      nt_Y5S_fit_boost_e.push_back(cand.Y5S_fit_boost_e);
      nt_Y5S_fit_boost_p.push_back(cand.Y5S_fit_boost_p);
      nt_Y5S_fit_boost_m.push_back(cand.Y5S_fit_boost_m);
      nt_Y5S_fit_boost_pt.push_back(cand.Y5S_fit_boost_pt);
      nt_Y5S_fit_boost_phi.push_back(cand.Y5S_fit_boost_phi);
      nt_Y5S_fit_boost_costh.push_back(cand.Y5S_fit_boost_costh);

      nt_Wbj_e.push_back(cand.Wbj_e);
      nt_Wbj_p.push_back(cand.Wbj_p);
      nt_Wbj_m.push_back(cand.Wbj_m);
      nt_Wbj_pt.push_back(cand.Wbj_pt);
      nt_Wbj_phi.push_back(cand.Wbj_phi);
      nt_Wbj_costh.push_back(cand.Wbj_costh);

      nt_Wbj_fit_e.push_back(cand.Wbj_fit_e);
      nt_Wbj_fit_p.push_back(cand.Wbj_fit_p);
      nt_Wbj_fit_m.push_back(cand.Wbj_fit_m);
      nt_Wbj_fit_pt.push_back(cand.Wbj_fit_pt);
      nt_Wbj_fit_phi.push_back(cand.Wbj_fit_phi);
      nt_Wbj_fit_costh.push_back(cand.Wbj_fit_costh);

      nt_Wbj_fit_boost_e.push_back(cand.Wbj_fit_boost_e);
      nt_Wbj_fit_boost_p.push_back(cand.Wbj_fit_boost_p);
      nt_Wbj_fit_boost_m.push_back(cand.Wbj_fit_boost_m);
      nt_Wbj_fit_boost_pt.push_back(cand.Wbj_fit_boost_pt);
      nt_Wbj_fit_boost_phi.push_back(cand.Wbj_fit_boost_phi);
      nt_Wbj_fit_boost_costh.push_back(cand.Wbj_fit_boost_costh);
      
      nt_gam_recoil_e_1.push_back(cand.gam_recoil_e_1);
      nt_gam_recoil_p_1.push_back(cand.gam_recoil_p_1);
      nt_gam_recoil_m_1.push_back(cand.gam_recoil_m_1);
      nt_gam_recoil_pt_1.push_back(cand.gam_recoil_pt_1);
      nt_gam_recoil_phi_1.push_back(cand.gam_recoil_phi_1);
      nt_gam_recoil_costh_1.push_back(cand.gam_recoil_costh_1);

      nt_gam_recoil_e_2.push_back(cand.gam_recoil_e_2);
      nt_gam_recoil_p_2.push_back(cand.gam_recoil_p_2);
      nt_gam_recoil_m_2.push_back(cand.gam_recoil_m_2);
      nt_gam_recoil_pt_2.push_back(cand.gam_recoil_pt_2);
      nt_gam_recoil_phi_2.push_back(cand.gam_recoil_phi_2);
      nt_gam_recoil_costh_2.push_back(cand.gam_recoil_costh_2);

      nt_gam_recoil_e_3.push_back(cand.gam_recoil_e_3);
      nt_gam_recoil_p_3.push_back(cand.gam_recoil_p_3);
      nt_gam_recoil_m_3.push_back(cand.gam_recoil_m_3);
      nt_gam_recoil_pt_3.push_back(cand.gam_recoil_pt_3);
      nt_gam_recoil_phi_3.push_back(cand.gam_recoil_phi_3);
      nt_gam_recoil_costh_3.push_back(cand.gam_recoil_costh_3);

      nt_gam_recoil_e.push_back(cand.gam_recoil_e);
      nt_gam_recoil_p.push_back(cand.gam_recoil_p);
      nt_gam_recoil_m.push_back(cand.gam_recoil_m);
      nt_gam_recoil_pt.push_back(cand.gam_recoil_pt);
      nt_gam_recoil_phi.push_back(cand.gam_recoil_phi);
      nt_gam_recoil_costh.push_back(cand.gam_recoil_costh);
      
      //nt_gam_recoil_boost_m.push_back(cand.gam_recoil_boost_m);

      nt_YnS_e.push_back(cand.YnS_e);
      nt_YnS_p.push_back(cand.YnS_p);
      nt_YnS_m.push_back(cand.YnS_m);
      nt_YnS_pt.push_back(cand.YnS_pt);
      nt_YnS_phi.push_back(cand.YnS_phi);
      nt_YnS_costh.push_back(cand.YnS_costh);

      nt_YnS_boost_e.push_back(cand.YnS_boost_e);
      nt_YnS_boost_p.push_back(cand.YnS_boost_p);
      nt_YnS_boost_m.push_back(cand.YnS_boost_m);
      nt_YnS_boost_pt.push_back(cand.YnS_boost_pt);
      nt_YnS_boost_phi.push_back(cand.YnS_boost_phi);
      nt_YnS_boost_costh.push_back(cand.YnS_boost_costh);

      nt_YnS_fit_e.push_back(cand.YnS_fit_e);
      nt_YnS_fit_p.push_back(cand.YnS_fit_p);
      nt_YnS_fit_m.push_back(cand.YnS_fit_m);
      nt_YnS_fit_pt.push_back(cand.YnS_fit_pt);
      nt_YnS_fit_phi.push_back(cand.YnS_fit_phi);
      nt_YnS_fit_costh.push_back(cand.YnS_fit_costh);
      
      nt_rho_e.push_back(cand.rho_e);
      nt_rho_p.push_back(cand.rho_p);
      nt_rho_m.push_back(cand.rho_m);
      nt_rho_pt.push_back(cand.rho_pt);
      nt_rho_phi.push_back(cand.rho_phi);
      nt_rho_costh.push_back(cand.rho_costh);

      nt_mup_e.push_back(cand.mup_e);
      nt_mup_p.push_back(cand.mup_p);
      nt_mup_m.push_back(cand.mup_m);
      nt_mup_pt.push_back(cand.mup_pt);
      nt_mup_phi.push_back(cand.mup_phi);
      nt_mup_costh.push_back(cand.mup_costh);
      
      nt_mum_e.push_back(cand.mum_e);
      nt_mum_p.push_back(cand.mum_p);
      nt_mum_m.push_back(cand.mum_m);
      nt_mum_pt.push_back(cand.mum_pt);
      nt_mum_phi.push_back(cand.mum_phi);
      nt_mum_costh.push_back(cand.mum_costh);
      
      nt_pip_e.push_back(cand.pip_e);
      nt_pip_p.push_back(cand.pip_p);
      nt_pip_m.push_back(cand.pip_m);
      nt_pip_pt.push_back(cand.pip_pt);
      nt_pip_phi.push_back(cand.pip_phi);
      nt_pip_costh.push_back(cand.pip_costh);
      
      nt_pim_e.push_back(cand.pim_e);
      nt_pim_p.push_back(cand.pim_p);
      nt_pim_m.push_back(cand.pim_m);
      nt_pim_pt.push_back(cand.pim_pt);
      nt_pim_phi.push_back(cand.pim_phi);
      nt_pim_costh.push_back(cand.pim_costh);
      
      nt_gam_e.push_back(cand.gam_e);
      nt_gam_p.push_back(cand.gam_p);
      nt_gam_m.push_back(cand.gam_m);
      nt_gam_pt.push_back(cand.gam_pt);
      nt_gam_phi.push_back(cand.gam_phi);
      nt_gam_costh.push_back(cand.gam_costh);

      // Missing rec values
      nt_pip_pim_recoil_m.push_back(cand.pip_pim_recoil_m);
      nt_pip_pim_recoil_e.push_back(cand.pip_pim_recoil_e);

      nt_pip_pim_gam_recoil_e.push_back(cand.pip_pim_gam_recoil_e);
      nt_pip_pim_gam_recoil_p.push_back(cand.pip_pim_gam_recoil_p);
      nt_pip_pim_gam_recoil_m.push_back(cand.pip_pim_gam_recoil_m);
      nt_pip_pim_gam_recoil_pt.push_back(cand.pip_pim_gam_recoil_pt);
      nt_pip_pim_gam_recoil_phi.push_back(cand.pip_pim_gam_recoil_phi);
      nt_pip_pim_gam_recoil_costh.push_back(cand.pip_pim_gam_recoil_costh);

      nt_pip_pim_gam_recoil_boost_e.push_back(cand.pip_pim_gam_recoil_boost_e);
      nt_pip_pim_gam_recoil_boost_p.push_back(cand.pip_pim_gam_recoil_boost_p);
      nt_pip_pim_gam_recoil_boost_m.push_back(cand.pip_pim_gam_recoil_boost_m);
      nt_pip_pim_gam_recoil_boost_pt.push_back(cand.pip_pim_gam_recoil_boost_pt);
      nt_pip_pim_gam_recoil_boost_phi.push_back(cand.pip_pim_gam_recoil_boost_phi);
      nt_pip_pim_gam_recoil_boost_costh.push_back(cand.pip_pim_gam_recoil_boost_costh);

      nt_mup_mum_recoil_m.push_back(cand.mup_mum_recoil_m);
      nt_mup_mum_recoil_e.push_back(cand.mup_mum_recoil_e);

      nt_mup_mum_gam_recoil_e.push_back(cand.mup_mum_gam_recoil_e);
      nt_mup_mum_gam_recoil_p.push_back(cand.mup_mum_gam_recoil_p);
      nt_mup_mum_gam_recoil_m.push_back(cand.mup_mum_gam_recoil_m);
      nt_mup_mum_gam_recoil_pt.push_back(cand.mup_mum_gam_recoil_pt);
      nt_mup_mum_gam_recoil_phi.push_back(cand.mup_mum_gam_recoil_phi);
      nt_mup_mum_gam_recoil_costh.push_back(cand.mup_mum_gam_recoil_costh);

      nt_mup_mum_gam_fit_recoil_m.push_back(cand.mup_mum_gam_fit_recoil_m);

      nt_mup_mum_gam_recoil_boost_e.push_back(cand.mup_mum_gam_recoil_boost_e);
      nt_mup_mum_gam_recoil_boost_p.push_back(cand.mup_mum_gam_recoil_boost_p);
      nt_mup_mum_gam_recoil_boost_m.push_back(cand.mup_mum_gam_recoil_boost_m);
      nt_mup_mum_gam_recoil_boost_pt.push_back(cand.mup_mum_gam_recoil_boost_pt);
      nt_mup_mum_gam_recoil_boost_phi.push_back(cand.mup_mum_gam_recoil_boost_phi);
      nt_mup_mum_gam_recoil_boost_costh.push_back(cand.mup_mum_gam_recoil_boost_costh);

      nt_pip_pim_mup_mum_recoil_e.push_back(cand.pip_pim_mup_mum_recoil_e);
      nt_pip_pim_mup_mum_recoil_p.push_back(cand.pip_pim_mup_mum_recoil_p);
      nt_pip_pim_mup_mum_recoil_m.push_back(cand.pip_pim_mup_mum_recoil_m);
      nt_pip_pim_mup_mum_recoil_pt.push_back(cand.pip_pim_mup_mum_recoil_pt);
      nt_pip_pim_mup_mum_recoil_phi.push_back(cand.pip_pim_mup_mum_recoil_phi);
      nt_pip_pim_mup_mum_recoil_costh.push_back(cand.pip_pim_mup_mum_recoil_costh);

      nt_pip_pim_mup_mum_recoil_m2.push_back(cand.pip_pim_mup_mum_recoil_m2);
      nt_pip_pim_mup_mum_fit_recoil_m.push_back(cand.pip_pim_mup_mum_fit_recoil_m);
      nt_pip_pim_mup_mum_fit_recoil_m2.push_back(cand.pip_pim_mup_mum_fit_recoil_m2);
      nt_pip_pim_mup_mum_fit_recoil_e.push_back(cand.pip_pim_mup_mum_fit_recoil_e);
      nt_pip_pim_mup_mum_fit_recoil_costh.push_back(cand.pip_pim_mup_mum_fit_recoil_costh);

      nt_pip_pim_mup_mum_recoil_boost_e.push_back(cand.pip_pim_mup_mum_recoil_boost_e);
      nt_pip_pim_mup_mum_recoil_boost_p.push_back(cand.pip_pim_mup_mum_recoil_boost_p);
      nt_pip_pim_mup_mum_recoil_boost_m.push_back(cand.pip_pim_mup_mum_recoil_boost_m);
      nt_pip_pim_mup_mum_recoil_boost_pt.push_back(cand.pip_pim_mup_mum_recoil_boost_pt);
      nt_pip_pim_mup_mum_recoil_boost_phi.push_back(cand.pip_pim_mup_mum_recoil_boost_phi);
      nt_pip_pim_mup_mum_recoil_boost_costh.push_back(cand.pip_pim_mup_mum_recoil_boost_costh);

      nt_gam_boost_e.push_back(cand.gam_boost_e);
      nt_gam_boost_p.push_back(cand.gam_boost_p);
      nt_gam_boost_m.push_back(cand.gam_boost_m);
      nt_gam_boost_pt.push_back(cand.gam_boost_pt);
      nt_gam_boost_phi.push_back(cand.gam_boost_phi);
      nt_gam_boost_costh.push_back(cand.gam_boost_costh);
    
      nt_deltaE.push_back(cand.deltaE);
      nt_deltaE_Y1D.push_back(cand.deltaE_Y1D);

      nt_mup_mum_gam_m.push_back(cand.mup_mum_gam_m);
      nt_mup_mum_gam_e.push_back(cand.mup_mum_gam_e);
      nt_pip_pim_gam_m.push_back(cand.pip_pim_gam_m);      
      nt_pip_pim_gam_e.push_back(cand.pip_pim_gam_e); 
      nt_chi2_YnS.push_back(cand.chi2_YnS);

    }

    //    cout << "end prepareTreeRec" << endl;
  }


  bool Analysis::goodGamma(const Mdst_gamma &gamma) {
    //    cout << "goodGamma" << endl;

    bool passed = false;

    // Baseline requirements: 150>theta>17, E>20 MeV, e9/e25>0.75, width<6, no match
    if( good_gamma( gamma ) ) {

      Mdst_ecl &shower = gamma.ecl();
      double energy = shower.energy();
      double theta = shower.theta() * 180.0/M_PI;
      
      // require at least 50MeV/100MeV in barrel/endcap
      if( ( theta >= 32 && theta <= 129 && energy >= 0.050 ) || energy >= 0.100 ) {
	passed = true;
      }
    }
    return passed;
  }



  bool Analysis::goodChargedPion( const Mdst_charged &chrg ) {
    //    cout << "goodChargedPion" << endl;

    // initialize
    eid eid( chrg );
    atc_pid kid(3, 1, 5, 3, 2);
    
    // retrieve probabilities
    double prob_K = kid.prob( chrg );
    double prob_e = eid.prob( 3, -1, 5 );

    // We'll say it's a pion unless were are fairly certain it looks like a kaon or electron/positron.
    if( prob_K < 0.90 && prob_e < 0.90 ) {

      return true;
    }
    return false;
  }


  
  // Obtain helix wrt the IP
  Helix Analysis::getChargedTrackHelixAtIP( const Mdst_charged &chrg, int mass_hypothesis) {
    //    cout << "getPionHelixAtIP" << endl;
   
    const Mdst_trk_fit& trkFit = chrg.trk().mhyp(mass_hypothesis);

    const HepPoint3D pivot( trkFit.pivot(0), trkFit.pivot(1), trkFit.pivot(2) );
    
    HepVector a(5);
    a[0] = trkFit.helix(0);
    a[1] = trkFit.helix(1);
    a[2] = trkFit.helix(2);
    a[3] = trkFit.helix(3);
    a[4] = trkFit.helix(4);

    /*
      HepSymMatrix Ea(5,0);
      Ea[0][0] = trkFit.error(0);
      Ea[1][0] = trkFit.error(1);
      Ea[1][1] = trkFit.error(2);
      Ea[2][0] = trkFit.error(3);
      Ea[2][1] = trkFit.error(4);
      Ea[2][2] = trkFit.error(5);
      Ea[3][0] = trkFit.error(6);
      Ea[3][1] = trkFit.error(7);
      Ea[3][2] = trkFit.error(8);
      Ea[3][3] = trkFit.error(9);
      Ea[4][0] = trkFit.error(10);
      Ea[4][1] = trkFit.error(11);
      Ea[4][2] = trkFit.error(12);
      Ea[4][3] = trkFit.error(13);
      Ea[4][4] = trkFit.error(14);
    */

    //Helix helix(pivot, a, Ea );
    Helix helix(pivot, a);
    
    if( m_ip.mag() ) helix.pivot( m_ip );
    
    return helix;
  }

	
      
  bool Analysis::goodTrackQuality( const Mdst_charged &chrg, string particle_type ) {

    // Assume e hypothesis initially, but change it if it's wrong
    int mass_hypothesis = 0;

    if      (particle_type == "mu" ) mass_hypothesis = 1; 
    else if (particle_type == "pi" ) mass_hypothesis = 2;
    else if (particle_type == "K"  ) mass_hypothesis = 3; 
    else if (particle_type == "rho") mass_hypothesis = 4; 

    // baseline selection
    if( good_charged( chrg, 1.0e-25, 10., 5. ) ) {
    
      // require trk_fit_prob > 0,  dz < 10.cm,  dr < 5.0cm wrt origin
      Helix helix = getChargedTrackHelixAtIP( chrg, mass_hypothesis);
      double dr = helix.dr();
      double dz = helix.dz();
      
      // dz and dr wrt IP
      //if( fabs( dz ) <= 4.0 && fabs( dr ) <= 0.2 ) {
      if( fabs( dz ) <= 2.0 && fabs( dr ) <= 0.3 ) {
    
	double px = chrg.px();
	double py = chrg.py();
	double pt = sqrt( px*px + py*py );

	// require transverse momentum >= 0.1 GeV/c
	if( pt >= 0.100 ) {
	  return true;
	}
      }
    }
    return false;
  }


  bool Analysis::goodChargedMuon( const Mdst_charged &chrg) {
    //    cout << "goodChargedMuon" << endl;
    
    // initialize
    Muid_mdst muID(chrg);
    
    // retrieve probability
    double prob_mu = ( muID.Chi_2() <= 0 ) ? -2 : muID.Muon_likelihood();

    if ( prob_mu >= 0.10 ) {

      return true;   
    }
    return false;
  }



  void Analysis::writeTreesOut() {
    //cout << "writeTreesOut" << endl;
    
    tree_rec->Fill();
    
    if ( isSignalMC ) {
      tree_mc->Fill();
      //    tree_mc->Print();
    }
    
  }







// Kinematic Mass Fitter Attempt #2
  double Analysis::massFit2(Particle &p, int numDaughters)
  {
    kmassfitter kmf;
    //  kmf.invariantMass(PI0_MASS);
    kmf.invariantMass(p.pType().mass());
    
    // add tracks to vertex fit
    for ( int i = 0; i < numDaughters; i++) {
      addTrack2fit(kmf, p.relation().child(i));
    }
    //addTrack2fit(kmf, p.relation().child(0));
    //addTrack2fit(kmf, p.relation().child(1));
    
    HepPoint3D origin(0.,0.,0.);
    kmf.vertex(origin); // set position where "PI0" is reconstructed.
    
    unsigned kmf_error = kmf.fit();
    
    if ( kmf_error != 0 ) {
      // cout << "Non-zero RC from the fit! kmf_error = " << kmf_error << endl; 
    }
    
    setUserInfo(p);
    
    if ( kmf_error == 0 ) {
      dynamic_cast<UserInfo&>(p.userInfo()).cl(kmf.cl());
      dynamic_cast<UserInfo&>(p.userInfo()).chisq(kmf.chisq());
      dynamic_cast<UserInfo&>(p.userInfo()).ndf(kmf.dgf());
    }
    else {
      dynamic_cast<UserInfo&>(p.userInfo()).cl(-999.);
      dynamic_cast<UserInfo&>(p.userInfo()).chisq(-999.);
      dynamic_cast<UserInfo&>(p.userInfo()).ndf(-999);
    }
    
    double chi2 = kmf.chisq();
    /*
      p.userInfo(&double chi2);
      // cout << "chi2 = " << chi2 << endl;
      */  
    
    kmakemother kmm;
    
    makeMother(kmm,kmf,p,1);
    
    kmm.vertex(origin);
    
    unsigned kmm_error = kmm.make();
    
    if ( kmm_error != 0 ) {
      // cout << "Non-zero RC from kmm,make()! kmm_error = " << kmm_error << endl; 
    }
    /*
    // cout << "----------------Before Fit---------------------" << endl;
    // cout << "pi: px = " << p.p().x() << endl;
    // cout << "pi: py = " << p.p().y() << endl;
    // cout << "pi: pz = " << p.p().z() << endl;
    */
    
    p.momentum().momentumPosition(kmm.momentum(), // set "PI0" information.
				  kmm.position(), // 4-momentum, position and these error.
				  kmm.error()); 
    /*
    // cout << "----------------After Fit---------------------" << endl;
    // cout << "pi: px = " << p.p().x() << endl;
    // cout << "pi: py = " << p.p().y() << endl;
    // cout << "pi: pz = " << p.p().z() << endl;
    */
    HepSymMatrix errVtx(3,0);
    p.momentum().decayVertex(kmf.vertex(),errVtx); // set decay point. (error matrix is meaningless.)
    
    //  Particle mother;
    //  mother.pType(Ptype("PI0"));
    
    //  mother.relation().append(p.relation().child(0));
    //  mother.relation().append(p.relation().child(1));
    //  makeMother(km, mother);
    
    // MORE WORK IS NEEDED HERE
    
    
    
    //  double mother_mass = mother.p().mag();
    //  double mother_energy = mother.p().t();
    //
    //  // cout << "mother_mass = " << mother_mass << endl;
    //  // cout << "mother_energy = " << mother_energy << endl;
    
    // VS:   DANGEROUS - NEED TO DISCUSS!!!  p = mother;
    
  return chi2;

}



 



  void Analysis::printCompleteMCHistory() {

    if ( isMC ) {

      cout << "VS-I-Reporting the entire MC record for the current event:" << endl;
      
      Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
      
      Gen_hepevt_Index gen_i = GenHepevtMgr.index( "mother" );
      gen_i.update(); // sort by pointer "mother" 	  
      
      for ( std::vector<Gen_hepevt>::iterator it = GenHepevtMgr.begin(); it != GenHepevtMgr.end(); it++ ) {
	
	Gen_hepevt &mc_particle = *it;
	
	int id = mc_particle.get_ID();
	int idhep = 0;
	if ( id != 0 ) idhep = mc_particle.idhep();
	
	cout << "    VS-I-ID, idhep = " << id << ", " << idhep << endl;
	
	if ( idhep == 911 ) {
	  cout << "Noise photon energy (as merged to data) = " << mc_particle.E() << endl;
	}

	std::vector<Gen_hepevt> mc_daughters = point_from( mc_particle.get_ID(), gen_i );
	cout << "        The number of daughter particles = "<< mc_daughters.size() << endl;
	
	if ( mc_daughters.size() != 0 ) {
	  std::vector<Gen_hepevt>::iterator it = mc_daughters.begin();
	  Gen_hepevt &mc_daughter = *it;
	  double vtx_radius = sqrt( mc_daughter.VX()*mc_daughter.VX() + mc_daughter.VY()*mc_daughter.VY() );
	  cout << "Radius (in cm) where the first daughter of this particle originated from = " << vtx_radius << endl;
	}

	for( std::vector<Gen_hepevt>::iterator it = mc_daughters.begin(); it != mc_daughters.end(); it++ ) {
	  
	  Gen_hepevt &mc_daughter = *it;
	  
	  int id = mc_daughter.get_ID();
	  int idhep = mc_daughter.idhep();
	  
	  cout << "            VS-I-daughter: ID, idhep = " << id << ", " << idhep << endl;

	}      
      }	
    }
  }

 



  void Analysis::photonSelectionEfficiencies() {
    //cout << "photonEfficiencyRec" << endl;

    m_count_gam_rec = 0;

    Mdst_gamma_Manager &GammaMgr = Mdst_gamma_Manager::get_manager();

    for( std::vector<Mdst_gamma>::iterator it = GammaMgr.begin(); it != GammaMgr.end(); it++ ) {
      
      Mdst_gamma &mdst_gamma = *it;
      const Gen_hepevt& hep_gamma(gen_level(get_hepevt(mdst_gamma)));

      if ( hep_gamma.get_ID() == m_mc_gam_id) {

	m_count_gam_rec++;

      }
    }
  }


  void Analysis::chargedSelectionEfficiencies() {
    //    cout << "chargedPionAndMuonEfficiencyRec" << endl;

    m_count_mup_rec = 0;
    m_count_mum_rec = 0;
    m_count_pip_rec = 0;
    m_count_pim_rec = 0;

    m_count_mup_good_trk = 0;
    m_count_mum_good_trk = 0;
    m_count_pip_good_trk = 0;
    m_count_pim_good_trk = 0;

    m_count_mup_good_dz = 0;
    m_count_mum_good_dz = 0;
    m_count_pip_good_dz = 0;
    m_count_pim_good_dz = 0;  

    m_count_mup_good_dr = 0;
    m_count_mum_good_dr = 0;
    m_count_pip_good_dr = 0;
    m_count_pim_good_dr = 0;  

    m_count_mup_good_pt = 0;
    m_count_mum_good_pt = 0;
    m_count_pip_good_pt = 0;
    m_count_pim_good_pt = 0;  

    Mdst_charged_Manager &ChargedMgr = Mdst_charged_Manager::get_manager();
    
    for ( std::vector<Mdst_charged>::iterator it = ChargedMgr.begin(); it != ChargedMgr.end(); it++ ) {

      Mdst_charged &mdst_chrg = *it;
      const Gen_hepevt& hep_chrg(gen_level(get_hepevt(mdst_chrg)));

      Helix helix_pi = getChargedTrackHelixAtIP( mdst_chrg, 2);
      double dr_pi = helix_pi.dr();
      double dz_pi = helix_pi.dz();

      Helix helix_mu = getChargedTrackHelixAtIP( mdst_chrg, 1);
      double dr_mu = helix_mu.dr();
      double dz_mu = helix_mu.dz();

      double px = mdst_chrg.px();
      double py = mdst_chrg.py();
      double pt = sqrt( px*px + py*py );

      double dr_cut = 0.3;
      double dz_cut = 2.0;
      double pt_cut = 0.1;

      if ( mdst_chrg.charge() > 0 ) {

	if ( hep_chrg.get_ID() == m_mc_pip_id ) {
	  m_count_pip_rec++;
	  
	  if ( fabs(dr_pi) <= dr_cut ) {
	    m_count_pip_good_dr++;
	  }
	  
	  if ( fabs(dz_pi) <= dz_cut ) {
	    m_count_pip_good_dz++;
	  }
	  
	  if ( pt > pt_cut ) {
	    m_count_pip_good_pt++;
	  }
	  
	  if ( goodTrackQuality( mdst_chrg, "pi" ) ) {
	    m_count_pip_good_trk++;
	  }
	}
 
	else if ( hep_chrg.get_ID() == m_mc_mup_id ) {
	  m_count_mup_rec++;
	  
	  if ( fabs(dr_mu) <= dr_cut ) {
	    m_count_mup_good_dr++;
	  }
	  
	  if ( fabs(dz_mu) <= dz_cut ) {
	    m_count_mup_good_dz++;
	  }
	  
	  if ( pt > pt_cut ) {
	    m_count_mup_good_pt++;
	  }
	  
	  if ( goodTrackQuality( mdst_chrg, "mu" ) ) {
	    m_count_mup_good_trk++;
	  }
	}
      }

      else if ( mdst_chrg.charge() < 0 ) {

	if ( hep_chrg.get_ID() == m_mc_pim_id ) {
	  m_count_pim_rec++;
	  
	  if ( fabs(dr_pi) <= dr_cut ) {
	    m_count_pim_good_dr++;
	  }
	  
	  if ( fabs(dz_pi) <= dz_cut ) {
	    m_count_pim_good_dz++;
	  }
	  
	  if ( pt > pt_cut ) {
	    m_count_pim_good_pt++;
	  }
	  
	  if ( goodTrackQuality( mdst_chrg, "pi" ) ) {
	    m_count_pim_good_trk++;
	  }
	}
 
	else if ( hep_chrg.get_ID() == m_mc_mum_id ) {
	  m_count_mum_rec++;
	  
	  if ( fabs(dr_mu) <= dr_cut ) {
	    m_count_mum_good_dr++;
	  }
	  
	  if ( fabs(dz_mu) <= dz_cut ) {
	    m_count_mum_good_dz++;
	  }
	  
	  if ( pt > pt_cut ) {
	    m_count_mum_good_pt++;
	  }
	  
	  if ( goodTrackQuality( mdst_chrg, "mu" ) ) {
	    m_count_mum_good_trk++;
	  }
	}
      }
    }
  }

  void Analysis::overallEfficiencies() {
    //    cout << "chargedPionAndMuonEfficiencyPID" << endl;

    m_count_gam_eff = 0;
    m_count_mup_eff = 0;
    m_count_mum_eff = 0;
    m_count_pip_eff = 0;
    m_count_pim_eff = 0;
      
    for ( std::vector<Particle>::iterator it_gam = gam_list.begin(); it_gam != gam_list.end(); it_gam++ ) {

      const Mdst_gamma &mdst_gamma = (*it_gam).mdstGamma();
      const Gen_hepevt& hep_gam(gen_level(get_hepevt(mdst_gamma)));

      if ( hep_gam.get_ID() == m_mc_gam_id) {
	m_count_gam_eff++;
      }
    }  
    
    for ( std::vector<Particle>::iterator it_pip = pip_list.begin(); it_pip !=  pip_list.end(); it_pip++ ) {

      const Mdst_charged &mdst_pip = (*it_pip).mdstCharged();
      const Gen_hepevt& hep_pip(gen_level(get_hepevt(mdst_pip)));

      if ( hep_pip.get_ID() == m_mc_pip_id ) {
	m_count_pip_eff++;
      }
    }

    for ( std::vector<Particle>::iterator it_pim = pim_list.begin(); it_pim !=  pim_list.end(); it_pim++ ) {

      const Mdst_charged &mdst_pim = (*it_pim).mdstCharged();
      const Gen_hepevt& hep_pim(gen_level(get_hepevt(mdst_pim)));

      if ( hep_pim.get_ID() == m_mc_pim_id ) {
	m_count_pim_eff++;
      }
    }

    for ( std::vector<Particle>::iterator it_mup = mup_list.begin(); it_mup !=  mup_list.end(); it_mup++ ) {
 
      const Mdst_charged &mdst_mup = (*it_mup).mdstCharged();
      const Gen_hepevt& hep_mup(gen_level(get_hepevt(mdst_mup)));

      if ( hep_mup.get_ID() == m_mc_mup_id ) {
	m_count_mup_eff++;
      }
    }

    for ( std::vector<Particle>::iterator it_mum = mum_list.begin(); it_mum !=  mum_list.end(); it_mum++ ) {

      const Mdst_charged &mdst_mum = (*it_mum).mdstCharged();
      const Gen_hepevt& hep_mum(gen_level(get_hepevt(mdst_mum)));

      if ( hep_mum.get_ID() == m_mc_mum_id ) {
	m_count_mum_eff++;
      }
    }
  }


  void Analysis::getRelevantPhotonIDHEPs(int id) {

    int parent_idhep = 0;

    Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
    Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

    Mdst_gamma &mdst_gamma = GammaMgr[id-1]; 

    const Gen_hepevt& mc_particle(gen_level(get_hepevt(mdst_gamma)));

    if ( mc_particle ) {

      nt_gam_id.push_back(mc_particle.idhep());

      cout << "Signal photon candidate is tagged as MC particle with ID and idhep = " << mc_particle.get_ID() << ", " << mc_particle.idhep() << endl;

      Gen_hepevt_Index gen_i = GenHepevtMgr.index( "mother" );
      gen_i.update(); // sort by pointer "mother"
      std::vector<Gen_hepevt> mc_daughters = point_from( mc_particle.get_ID(), gen_i );
      cout << "The number of daughters: " << mc_daughters.size() << endl;

      if ( mc_daughters.size() != 0 ) {

	std::vector<Gen_hepevt>::iterator it = mc_daughters.begin();
	Gen_hepevt &mc_daughter = *it;
	double vtx_radius = sqrt( mc_daughter.VX()*mc_daughter.VX() + mc_daughter.VY()*mc_daughter.VY() );
	cout << "Radius (in cm) where the first daughter of this particle originated from = " << vtx_radius << endl;

      }

	Gen_hepevt &mc_mother = mc_particle.mother();
	if (mc_mother) {
	  nt_gam_parent_id.push_back(mc_mother.idhep());
	}
	else {
	  nt_gam_parent_id.push_back(-2);
	  //	  cout << "-D-VS- Signal photon candidate was MC tagged but has no MC parent" << endl;
	  //printCompleteMCHistory();
	}
    }
    else {       
      nt_gam_id.push_back(-1);
      nt_gam_parent_id.push_back(-1);

      //      cout << "-D-VS- Signal photon candidate was not MC tagged" << endl;
      //printCompleteMCHistory();
    }
  }


  void Analysis::getRelevantChargedIDHEPs(int id) {

    int parent_idhep = 0;

    Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
    Mdst_charged_Manager &ChargedMgr = Mdst_charged_Manager::get_manager();
    Gen_hepevt_Index gen_i = GenHepevtMgr.index( "mother" );
    gen_i.update(); // sort by pointer "mother"

    Mdst_charged &mdst_chrg = ChargedMgr[id-1]; 

    const Gen_hepevt& mc_particle(gen_level(get_hepevt(mdst_chrg)));

    if ( mc_particle ) {
      int chrg_idhep = mc_particle.idhep();
      Gen_hepevt &mc_mother = mc_particle.mother();

      if ( mc_mother ) {
	int mc_mother_id = mc_mother.idhep();
	
	if      ( chrg_idhep == -13 ) m_mup_parent_id = mc_mother_id;
        else if ( chrg_idhep ==  13 ) m_mum_parent_id = mc_mother_id;

	else if ( fabs(chrg_idhep) == 211 ) {
	  int pip_count = 0;
	  int pim_count = 0;
	  std::vector<Gen_hepevt> mc_daughters = point_from( mc_mother.get_ID(), gen_i );
	  for( std::vector<Gen_hepevt>::iterator it = mc_daughters.begin(); it != mc_daughters.end(); it++ ) {
       
	    Gen_hepevt &mc_daughter = *it;
	    int mc_recoil_idhep = mc_daughter.idhep();

	    if ( chrg_idhep ==  211) {

	      pip_count++;
	      if ( fabs(mc_recoil_idhep) != 211) m_pip_recoil_id = mc_recoil_idhep;
	      m_pip_parent_id = mc_mother_id;              
	    }
	    else if( chrg_idhep == -211) {

	      pim_count++;
	      if ( fabs(mc_recoil_idhep) != 211) m_pim_recoil_id = mc_recoil_idhep;
	      m_pim_parent_id = mc_mother_id;
	    }
	  }
	  /*
	  if (pip_count != 3 && pip_count != 0) {
	    cout <<"pip_count daughters: " << pip_count << endl;
	    cout << "pip_recoil_id: " << m_pip_recoil_id << endl;
	    printCompleteMCHistory();
	  }

	  if (pim_count != 3 && pim_count != 0) {
	    cout <<"pim_count daughters: " << pim_count << endl;
	    cout << "pim_recoil_id: " << m_pim_recoil_id << endl;
	    printCompleteMCHistory();
	  }
	  */
	}
      }
    }
  }
  
  void Analysis::diagnoseTreeRec() {
    //cout << "diagnoseTreeRec" << endl;

    m_bkg_A = 0;
    m_bkg_B = 0;
    m_bkg_C = 0;
    m_bkg_D = 0;
    m_bkg_E = 0;
    m_bkg_F = 0;
    m_bkg_G = 0;
    m_bkg_H = 0;
    m_bkg_I = 0;
    m_bkg_J = 0;
    m_bkg_K = 0;

    if ( !isMC ) return;

    // photon_parent_id* lists are NOT cleared here. The idea is that we want them to contain ALL decays detected. So, we clear them in begin_run only when m_run_seq == 0
    candidatesBkg_01.clear();
    candidatesBkg_02.clear();
    candidatesBkg_03.clear();
    candidatesBkg_04.clear();
    candidatesBkg_05.clear();
    candidatesBkg_06.clear();
    candidatesBkg_07.clear();
    candidatesBkg_08.clear();
    candidatesBkg_09.clear();

    int index = 0;

    for ( std::vector<candidateRec>::iterator it_cand = candidatesRec.begin(); it_cand !=  candidatesRec.end(); it_cand++ ) {
      //      cout << "loop over candidates again" << endl;

      index++;

      candidateRec cand = *it_cand;

      //if ( not(isSignalMC) ) {

      m_pip_parent_id = 0;
      m_pip_recoil_id = 0;
      
      m_pim_parent_id = 0;
      m_pim_recoil_id = 0;
      
      m_mup_parent_id = 0;
      m_mum_parent_id = 0;
      
      //      m_gam_parent_id = 0;
      
      getRelevantChargedIDHEPs(cand.pip_id);
      getRelevantChargedIDHEPs(cand.pim_id);
      getRelevantChargedIDHEPs(cand.mup_id);
      getRelevantChargedIDHEPs(cand.mum_id);
      
      //getRelevantPhotonIDHEPs(cand.gam_id);
      if (index == 1) {

	if (cand.Wbj_fit_m > 10.72 && cand.Wbj_fit_m < 10.8 ) {
	  cout << "R1:" << endl;
	  printCompleteMCHistory();
	}
	
	if ( cand.Wbj_fit_m > 10.38 && cand.Wbj_fit_m < 10.49 ) {
	  cout << "R2:" << endl;
	  printCompleteMCHistory();
	}
	
	if (cand.Wbj_fit_m > 10.49 && cand.gam_recoil_m < 10.72 && cand.deltaE > 0.03 ) {
	  cout << "R3:" << endl;
	  printCompleteMCHistory();
	}
	
	if (cand.Wbj_fit_m > 10.8 && cand.Wbj_fit_m < 11) {
	  cout << "Y1Spipi:" << endl;
	  printCompleteMCHistory();
	}

	if (cand.mup_mum_recoil_m < 0 ) {
	  cout << "Mrec(mumu) < 0:" << endl;
	  printCompleteMCHistory();
	}
      }

      if (index == 1) {

	if ( cand.deltaE < 0.03 && cand.deltaE > -0.05 && cand.Wbj_fit_m > 10.42 && cand.Wbj_fit_m < 10.72) {
	    cout << "deltaE cut" << endl;
	    //printCompleteMCHistory();
	  }

	if ( cand.Wbj_fit_m > 10.4 && cand.Wbj_fit_m < 10.8 && m_Y1S) {
	    cout << "blahblah" << endl;
	    //getRelevantPhotonIDHEPs(cand.gam_id);
	    //printCompleteMCHistory();
	    cout << "Photon energy (as reconstructed) = " << cand.gam_e << endl;
	  }
      }


      /*
      if ( index == 1 ) {
	
	double x = cand.pip_pim_recoil_m;
	double y = cand.pip_pim_gam_recoil_m-cand.gam_recoil_m+cand.Wbj_fit_m;

	double x2 = cand.pip_pim_gam_recoil_m;
	double y2 = cand.Wbj_fit_m;
	
	if ( x >= 9.44 && x <= 9.47 ) {
	  if ( y >= 9.44 && y <= 9.47) {
	    cout << "bkgA" << endl;
	    printCompleteMCHistory();
	    m_bkg_A = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	if ( x >= 9.95 && x <= 10.00 ) {
	  if ( y >= 9.45 && y <= 9.47) {
	    cout << "bkgB" << endl;
	    printCompleteMCHistory();
	    m_bkg_B = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.01 && x <= 10.04 ) {
	  if ( y >= 9.45 && y <= 9.47) {
	    cout << "bkgC" << endl;
	    printCompleteMCHistory();
	    m_bkg_C = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.28 && x <= 10.31 ) {
	  if ( y >= 9.44 && y <= 9.47) {
	    cout << "bkgD" << endl;
	    printCompleteMCHistory();
	    m_bkg_D = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.35 && x <= 10.36 ) {
	  if ( y >= 9.45 && y <= 9.47) {
	    cout << "bkgE" << endl;
	    printCompleteMCHistory();
	    m_bkg_E = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	// Bkg_F is very obvious and also extremely small...
	if ( x >= 9.70 && x <= 10.30 ) {
	  if ( y >= 10.015 && y <= 10.03) {
	    cout << "bkgF" << endl;
	    printCompleteMCHistory();
	    m_bkg_F = 1;
	  }
	}
	
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.03 && x <= 10.40 ) {
	  if ( y >= 10.01 && y <= 10.03) {
	    cout << "bkgG" << endl;
	    printCompleteMCHistory();
	    m_bkg_G = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.02 && x <= 10.03 ) {
	  if ( y >= 10.35 && y <= 10.37) {
	    cout << "bkgH" << endl;
	    printCompleteMCHistory();
	    m_bkg_H = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	// Bkg_I is also very obvious...
	if ( x >= 9.80 && x <= 10.40 ) {
	  if ( y >= 10.5 && y <= 10.60) {
	    cout << "bkgI" << endl;
	    printCompleteMCHistory();
	    m_bkg_I = 1;
	  }
	}
       
	//-------------------------------------------------------------------------------------------
	if ( x >= 10.34 && x <= 10.36 ) {
	  if ( y2 >= 10.34 && y2 <= 10.37) {
	    cout << "bkgJ" << endl;
	    printCompleteMCHistory();
	    m_bkg_J = 1;
	  }
	}
	//-------------------------------------------------------------------------------------------
	// The events in the bottom right corner of the signal region...
	if ( x2 >= 9.2 && x2 <= 9.75 ) {
	  if ( y2 >= 10.40 && y2 <= 10.75) {
	    cout << "bkgK" << endl;
	    printCompleteMCHistory();
	    m_bkg_K = 1;
	  }
	}

	if ( cand.pip_pim_gam_recoil_m >= 9.3 && cand.pip_pim_gam_recoil_m <= 9.5 ) {
	  if ( cand.Wbj_fit_m > 10 && cand.Wbj_fit_m < 10.4 ) {
	    cout << "Weird bkg...Reporting MC history" << endl;
	    printCompleteMCHistory();
	  }
	}

	if ( cand.pip_pim_gam_recoil_m >= 9.45 && cand.pip_pim_gam_recoil_m <= 9.55 ) {
	  if ( cand.Wbj_fit_m > 10.6 && cand.Wbj_fit_m < 10.62 ) {
	    cout << "Background in signal region after extending vertical range:" << endl;
	    printCompleteMCHistory();
	  }
	}
	  
      }
      */
    }
  }

  double Analysis::getNeutralEnergy() {

    HepLorentzVector p4_neutral(0, 0, 0, 0); 

    for ( std::vector<Particle>::iterator it_neutral = leftover_neutral_list.begin(); it_neutral !=  leftover_neutral_list.end(); it_neutral++ ) {
      
	Particle &neutral = (*it_neutral);
	p4_neutral += neutral.p();
    }
    return p4_neutral.e();
  }

  double Analysis::getChargedEnergy() {
 
    HepLorentzVector p4_charged(0, 0, 0, 0); 

    for ( std::vector<Particle>::iterator it_charged = leftover_charged_list.begin(); it_charged !=  leftover_charged_list.end(); it_charged++ ) {
	Particle &charged = (*it_charged);      
	p4_charged += charged.p();
    }
    return p4_charged.e();
  }

  int Analysis::vetoPi0(int gam_id) {

    int i = 0;
    Mdst_pi0_Manager& Pi0Mgr = Mdst_pi0_Manager::get_manager();
    for ( std::vector<Mdst_pi0>::iterator it_pi0 = Pi0Mgr.begin(); it_pi0 != Pi0Mgr.end(); it_pi0++ ) { 

      Mdst_pi0 &mdst_pi0 = (*it_pi0); 
      Particle pi0(mdst_pi0);
      
      Mdst_gamma const &mdst_gam_1 = pi0.relation().child(0).mdstGamma();
      Mdst_gamma const &mdst_gam_2 = pi0.relation().child(1).mdstGamma();

      if ( gam_id == mdst_gam_1.get_ID() || gam_id == mdst_gam_2.get_ID()     ) {
	m_count_gam_from_pi0++;
	i++;
	//	cout << "-D-NC-selected gam: " << gam_id << endl;
	//	cout << "-D-NC-mdst_gam_1, mdst_gam_2: " << mdst_gam_1.get_ID() << ", " << mdst_gam_2.get_ID() << endl;
	//cout<<"pi0 count: " <<m_count_gam_from_pi0 << endl;
      }
    }
    return i;
  }
  

  vector<int> Analysis::getPhotonChain( int id ) {
    cout << "getPhotonChain" << endl;

    photon_parent_id_list.clear();
    
    Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
    Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

    Mdst_gamma &mdst_gamma = GammaMgr[id-1]; 

    const Gen_hepevt& mc_particle(gen_level(get_hepevt(mdst_gamma)));

    if ( mc_particle ) {
      int gamma_tag_idhep = mc_particle.idhep();

      if ( gamma_tag_idhep == 22 ) {
	Gen_hepevt &mc_mother = mc_particle.mother();
	int mc_mother_id = mc_mother.idhep();

	while (mc_mother_id && mc_mother_id != 9000553) {
	  photon_parent_id_list.push_back(mc_mother_id);

	  mc_mother = mc_mother.mother();
	  mc_mother_id = mc_mother.idhep();
	}	

	std::reverse( photon_parent_id_list.begin(), photon_parent_id_list.end() );
	/*
	for ( int j=0; j < photon_parent_id_list.size(); j++) {
	  cout << "Particle " << j+1 << ": " << photon_parent_id_list[j]  << endl;    
	}
	       
	if (photon_parent_id_bkg_06.size() == 0) {
	  photon_parent_id_bkg_06.push_back(photon_parent_id_list);
	}
	
	else if ( find(photon_parent_id_bkg_01.begin(), photon_parent_id_bkg_01.end(), photon_parent_id_list) == photon_parent_id_bkg_01.end() ) {
	  photon_parent_id_bkg_06.push_back(photon_parent_id_list); 
	}
	
	for ( int i=0; i < photon_parent_id_bkg_06.size(); i++) {
	  cout << "BKG_06 DECAY CHAIN " << i+1 << endl;
	  for ( int j=0; j < photon_parent_id_bkg_06.size(); j++) {
	    cout << "Particle " << j+1 << ": " << photon_parent_id_bkg_06[i][j]  << endl;    
	  }
	}
	*/
      }
    }
  }

  void Analysis::analyzeBackgroundSignal() {

    if ( isMC == true) {

      Gen_hepevt_Manager &GenHepevtMgr = Gen_hepevt_Manager::get_manager();
      Gen_hepevt_Index gen_i = GenHepevtMgr.index( "mother" );
      gen_i.update(); // sort by pointer "mother" 

      m_Y1S = 0;
      
      m_Y2S        = 0;   
      m_Y2S_to_Y1S = 0;
      
      m_Y3S               = 0; 
      m_Y3S_to_Y1S        = 0;
      m_Y3S_to_Y2S        = 0;
      m_Y3S_to_Y2S_to_Y1S = 0;	  

      m_Y2S_to_pi0pi0 = 0;
      m_Y3S_to_pi0pi0 = 0;   
      m_Y5S_to_pi0pi0 = 0; 

      m_Y2S_to_pippim = 0;
      m_Y3S_to_pippim = 0;   
      m_Y5S_to_pippim = 0;

      m_Y3S_to_X2P = 0;

      int count_Y5S_to_pi0 = 0;
      int count_Y3S_to_pi0 = 0;
      int count_Y2S_to_pi0 = 0;

      int Y5S_to_pip = 0;
      int Y5S_to_pim = 0;
      int Y3S_to_pip = 0;
      int Y3S_to_pim = 0;
      int Y2S_to_pip = 0;
      int Y2S_to_pim = 0;

      int Y5S_to_K0 = 0;
      int Y5S_to_K0bar = 0;

      int count_Y3S_to_gam = 0;

      m_bkg_A1 = 0;
      
      m_bkg_B1 = 0;
      m_bkg_B2 = 0;
      
      m_bkg_C1 = 0;
      m_bkg_C2 = 0;
      m_bkg_C3 = 0;
      
      m_bkg_D1 = 0;
      
      m_bkg_E1 = 0;
      
      m_bkg_F1 = 0;
      
      m_bkg_G1 = 0;
      m_bkg_G2 = 0;
      m_bkg_G3 = 0;
      
      m_bkg_H1 = 0;
      m_bkg_H2 = 0;
      
      m_bkg_I1 = 0;
      
      m_bkg_J1 = 0;

      m_bkg_K1 = 0;
      m_Y5S_ndaughters = 0;
      
      for(std::vector<Gen_hepevt>::iterator it = GenHepevtMgr.begin(); it != GenHepevtMgr.end(); it++){
	
	Gen_hepevt* bkg_Y3S = 0;
	Gen_hepevt* bkg_Y2S = 0; 
	Gen_hepevt* bkg_Y1S = 0; 
	
	Gen_hepevt& bkg_particle = *it;
	int idhep = bkg_particle.idhep();
	
	//---------- MC TRUTH ----------//
	// @TODO: Add MC tagging for bkg?
	if ( idhep == 9000553 ) {
	  
	  std::vector<Gen_hepevt> bkg_Y5S_daughters = point_from( bkg_particle.get_ID(), gen_i );
	  if ( bkg_Y5S_daughters.size() != 0 ) m_Y5S_ndaughters = bkg_Y5S_daughters.size();
	  for(std::vector<Gen_hepevt>::iterator it_Y5S = bkg_Y5S_daughters.begin(); it_Y5S != bkg_Y5S_daughters.end(); it_Y5S++){
	    
	    Gen_hepevt &bkg_daughter = *it_Y5S;
	    int bkg_daughter_idhep = bkg_daughter.idhep();
	    
	    if ( bkg_daughter_idhep == 111 ) count_Y5S_to_pi0++;
	    else if ( bkg_daughter_idhep ==  211) Y5S_to_pip = 1;
	    else if ( bkg_daughter_idhep == -211) Y5S_to_pim = 1; 
	    else if ( bkg_daughter_idhep ==  311) Y5S_to_K0 = 1;
	    else if ( bkg_daughter_idhep ==  311) Y5S_to_K0bar = 1;

	    else if ( bkg_daughter_idhep == 200553 ) { // check Y(3S)
	      bkg_Y3S = &bkg_daughter;
	      m_Y3S = 1; 

	      std::vector<Gen_hepevt> bkg_Y3S_daughters = point_from( bkg_Y3S->get_ID(), gen_i );
	      for(std::vector<Gen_hepevt>::iterator it_Y3S = bkg_Y3S_daughters.begin(); it_Y3S != bkg_Y3S_daughters.end(); it_Y3S++){
		
		Gen_hepevt &bkg_daughter = *it_Y3S;
		int bkg_daughter_idhep = bkg_daughter.idhep();

		if ( bkg_daughter_idhep == 111 ) count_Y3S_to_pi0++;
		else if ( bkg_daughter_idhep ==  211 ) Y3S_to_pip = 1;
		else if ( bkg_daughter_idhep == -211 ) Y3S_to_pim = 1;
		else if ( bkg_daughter_idhep ==   22 ) count_Y3S_to_gam++;

		else if ( bkg_daughter_idhep == 100553 ) { // check Y(3S) -> Y(2S)
		  bkg_Y2S = &bkg_daughter;
		  m_Y3S_to_Y2S = 1;
		  m_Y3S = 0;
  		
		  std::vector<Gen_hepevt> bkg_Y2S_daughters = point_from( bkg_Y2S->get_ID(), gen_i );
		  for(std::vector<Gen_hepevt>::iterator it_Y2S = bkg_Y2S_daughters.begin(); it_Y2S != bkg_Y2S_daughters.end(); it_Y2S++){
		    
		    Gen_hepevt &bkg_daughter = *it_Y2S;
		    int bkg_daughter_idhep = bkg_daughter.idhep();

		    if ( bkg_daughter_idhep == 111 ) count_Y2S_to_pi0++; 
		    else if ( bkg_daughter_idhep ==  211) Y2S_to_pip = 1;
		    else if ( bkg_daughter_idhep == -211) Y2S_to_pim = 1; 

		    else if ( bkg_daughter_idhep == 553 ) { // check Y(3S) -> Y(2S) -> Y(1S)
		      m_Y3S_to_Y2S_to_Y1S = 1;
		      m_Y3S_to_Y2S = 0;
		    }
		  } // loop Y(2S)
		} // if Y(2S)

		else if ( bkg_daughter_idhep == 110551 || bkg_daughter_idhep == 120553 || bkg_daughter_idhep == 100555 ) {
		  m_Y3S_to_X2P = 1;
		  m_Y3S = 0;
		}
		else if ( bkg_daughter_idhep == 553 ) { // check Y(3S) -> Y(1S)
		  bkg_Y1S = &bkg_daughter;
		  m_Y3S_to_Y1S = 1;
		  m_Y3S = 0;
		}
	      } // loop Y(3S)
	    } // if Y(3S)

	    else if ( bkg_daughter_idhep == 100553 ) { // check Y(2S)
	      bkg_Y2S = &bkg_daughter;
	      m_Y2S = 1; 

	      std::vector<Gen_hepevt> bkg_Y2S_daughters = point_from( bkg_Y2S->get_ID(), gen_i );
	      for(std::vector<Gen_hepevt>::iterator it_Y2S = bkg_Y2S_daughters.begin(); it_Y2S != bkg_Y2S_daughters.end(); it_Y2S++){

		Gen_hepevt &bkg_daughter = *it_Y2S;
		int bkg_daughter_idhep = bkg_daughter.idhep();

		if ( bkg_daughter_idhep == 111 ) count_Y2S_to_pi0++; 
		else if ( bkg_daughter_idhep ==  211) Y2S_to_pip = 1;
		else if ( bkg_daughter_idhep == -211) Y2S_to_pim = 1; 

		else if ( bkg_daughter_idhep == 553 ) { // check Y(2S) -> Y(1S)
		  bkg_Y1S = &bkg_daughter;
		  m_Y2S_to_Y1S = 1;
		  m_Y2S = 0;
		}
	      } // loop Y(2S)
	    } // if Y(2S)
	    
	    else if ( bkg_daughter_idhep == 553 ) { // check Y(1S)
	      bkg_Y1S = &bkg_daughter;
	      m_Y1S = 1; 
	    }
	  } // loop Y(5S)

	  if ( count_Y5S_to_pi0 == 2 ) m_Y5S_to_pi0pi0 = 1;
	  if ( count_Y3S_to_pi0 == 2 ) m_Y3S_to_pi0pi0 = 1;
	  if ( count_Y2S_to_pi0 == 2 ) m_Y2S_to_pi0pi0 = 1;
	  
	  if ( Y5S_to_pip && Y5S_to_pim ) m_Y5S_to_pippim = 1;
	  if ( Y3S_to_pip && Y3S_to_pim ) m_Y3S_to_pippim = 1;
	  if ( Y2S_to_pip && Y2S_to_pim ) m_Y2S_to_pippim = 1;
   
    
	  if ( (m_Y3S + m_Y3S_to_Y2S + m_Y3S_to_Y1S + m_Y3S_to_Y2S_to_Y1S + m_Y2S + m_Y2S_to_Y1S + m_Y1S + m_Y3S_to_X2P == 1) &&
	       !(m_Y5S_to_pi0pi0 == 1 && m_Y5S_to_pippim == 1) &&
	       !(m_Y3S_to_pi0pi0 == 1 && m_Y3S_to_pippim == 1) &&
	       !(m_Y2S_to_pi0pi0 == 1 && m_Y2S_to_pippim == 1)    ) {
	    
	    // Assigning labels
	    // NOTE: You *must* specificy cand == 1 when plotting these; there can be multiple values for m_Y*S_pippim and m_Y*S_pi0pi0!!!
	    // There will be ~22,000 background events if you don't specify cand == 1 and ~14,000 if you do.
	    
	    m_bkg_A1 = (m_Y1S && m_Y5S_to_pippim) && !(m_Y3S_to_pippim || m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0); //good
	    
	    m_bkg_B1 = (m_Y3S_to_Y1S && m_Y5S_to_pippim && m_Y3S_to_pippim) && !(m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0); //not good
	    m_bkg_B2 = (m_Y3S_to_Y1S && m_Y5S_to_pi0pi0 && m_Y3S_to_pippim) && !(m_Y5S_to_pippim || m_Y2S_to_pippim || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0); //not good
	    
	    m_bkg_C1 = (m_Y2S_to_Y1S && m_Y5S_to_pippim && m_Y2S_to_pippim) && !(m_Y3S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0); //good
	    m_bkg_C2 = (m_Y2S_to_Y1S && m_Y5S_to_pippim && m_Y2S_to_pi0pi0) && !(m_Y3S_to_pippim || m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0); //good
	    //m_bkg_C3 = m_Y2S && m_resonance_fit == 1; 
	    
	    m_bkg_D1 = (m_Y2S_to_Y1S && m_Y5S_to_pi0pi0 && m_Y2S_to_pippim) && !(m_Y5S_to_pippim || m_Y3S_to_pippim || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0); 
	    
	    m_bkg_E1 = (m_Y3S_to_Y1S && m_Y5S_to_pippim && m_Y3S_to_pi0pi0) && !(m_Y3S_to_pippim || m_Y3S_to_pippim || m_Y5S_to_pi0pi0 || m_Y2S_to_pi0pi0); 
	    
	    //m_bkg_F1 = ( m_Y1S +  (m_resonance_fit == 2) == 2 ) ? 1 : 0 ;
	    
	    m_bkg_G1 = m_Y2S && m_Y5S_to_pippim && !(m_Y3S_to_pippim || m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0);
	    m_bkg_G2 = m_Y3S_to_Y2S && m_Y5S_to_pippim && count_Y3S_to_gam == 2 && !(m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0);
	    m_bkg_G3 = m_Y3S_to_X2P && m_Y5S_to_pippim && count_Y3S_to_gam == 1 && !(m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0);
	    
	    //m_bkg_H1 = 0; 
	    m_bkg_H2 = m_Y3S_to_Y2S && m_Y5S_to_pippim && m_Y3S_to_pippim && !(m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0);
	    
	    //m_bkg_I1 = m_Y2S && m_resonance_fit == 3;
	    
	    m_bkg_J1 = m_Y3S && m_Y5S_to_pippim && !(m_Y3S_to_pippim || m_Y2S_to_pippim || m_Y5S_to_pi0pi0 || m_Y3S_to_pi0pi0 || m_Y2S_to_pi0pi0);
	    /*
	    cout << "count gam = " << count_Y3S_to_gam << endl;
	    if (m_Y3S_to_X2P == 1) {
	      cout <<"X2P" <<endl;
	      printCompleteMCHistory();
	    }
	    */


	  }
	  else {
	    //	    cout << "New bkg: "<<endl;
	    //printCompleteMCHistory();
	    // @TODO: Keep track of the decay that you domn't account for in MC Truth (yet)
	  }	  
	}
      } // if Y(5S)
    } // loop MC
  } // isMC 
 
  void Analysis::binary2char(unsigned n) {

    char buf[32];
    
    int index = 0;
    for (unsigned i = 1 << 31; i > 0; i >>= 1) {
      sprintf(&buf[index],"%u",!!(n & i));
      index++;
    }

    std::cout << buf << std::endl;

  }

  void Analysis::reweight(double mc_gam_isr_boost_e) {

    double weights[1200] = {0.49712, 0.0312597, 0.019187, 0.0140328, 0.0111305, 0.00925674,	
			    0.00794209, 0.00696643, 0.00621229, 0.00561113, 0.00512018, 
			    0.00471134, 0.00436537, 0.00406863, 0.00381119, 0.00358562, 
			    0.00338629, 0.0032088, 0.00304972, 0.00290627, 0.00277624, 0.0026578, 
			    0.00254945, 0.00244994, 0.0023582, 0.00227336, 0.00219465, 
			    0.00212142, 0.00205311, 0.00198924, 0.00192937, 0.00187315, 
			    0.00182023, 0.00177034, 0.00172322, 0.00167864, 0.0016364, 
			    0.00159631, 0.00155821, 0.00152196, 0.00148743, 0.00145448, 
			    0.00142302, 0.00139295, 0.00136416, 0.00133659, 0.00131015, 
			    0.00128478, 0.0012604, 0.00123697, 0.00121443, 0.00119272, 
			    0.00117181, 0.00115164, 0.00113218, 0.00111339, 0.00109524, 
			    0.00107769, 0.00106071, 0.00104428, 0.00102837, 0.00101295, 
			    0.00099801, 0.000983515, 0.000969451, 0.000955797, 0.000942536, 
			    0.00092965, 0.000917124, 0.000904943, 0.000893092, 0.000881558, 
			    0.000870328, 0.00085939, 0.000848733, 0.000838346, 0.000828219, 
			    0.000818341, 0.000808704, 0.000799298, 0.000790116, 0.00078115,
			    0.000772391, 0.000763832, 0.000755468, 0.00074729, 0.000739293, 
			    0.00073147, 0.000723817, 0.000716327, 0.000708995, 0.000701817, 
			    0.000694787, 0.0006879, 0.000681153, 0.000674541, 0.000668059, 
			    0.000661705, 0.000655474, 0.000649363, 0.000643368, 0.000637486, 
			    0.000631714, 0.000626048, 0.000620486, 0.000615024, 0.000609661, 
			    0.000604393, 0.000599217, 0.000594132, 0.000589135, 0.000584224, 
			    0.000579396, 0.000574649, 0.000569982, 0.000565392, 0.000560877, 
			    0.000556435, 0.000552065, 0.000547765, 0.000543533, 0.000539368, 
			    0.000535268, 0.000531231, 0.000527256, 0.000523342, 0.000519486, 
			    0.000515689, 0.000511948, 0.000508262, 0.00050463, 0.000501051, 
			    0.000497524, 0.000494047, 0.000490619, 0.00048724, 0.000483908, 
			    0.000480622, 0.000477382, 0.000474186, 0.000471033, 0.000467923, 
			    0.000464855, 0.000461828, 0.00045884, 0.000455892, 0.000452982, 
			    0.00045011, 0.000447275, 0.000444477, 0.000441713, 0.000438985, 
			    0.00043629, 0.00043363, 0.000431002, 0.000428406, 0.000425842, 
			    0.00042331, 0.000420807, 0.000418335, 0.000415892, 0.000413478, 
			    0.000411093, 0.000408735, 0.000406404, 0.000404101, 0.000401824, 
			    0.000399573, 0.000397347, 0.000395147, 0.000392971, 0.00039082, 
			    0.000388692, 0.000386587, 0.000384506, 0.000382448, 0.000380411, 
			    0.000378397, 0.000376404, 0.000374432, 0.000372481, 0.000370551, 
			    0.000368641, 0.000366751, 0.00036488, 0.000363029, 0.000361196, 
			    0.000359383, 0.000357587, 0.00035581, 0.000354051, 0.000352309, 
			    0.000350584, 0.000348876, 0.000347186, 0.000345511, 0.000343853, 
			    0.000342211, 0.000340585, 0.000338974, 0.000337379, 0.000335799, 
			    0.000334234, 0.000332683, 0.000331147, 0.000329625, 0.000328117, 
			    0.000326623, 0.000325143, 0.000323677, 0.000322223, 0.000320783, 
			    0.000319356, 0.000317941, 0.000316539, 0.00031515, 0.000313772, 
			    0.000312407, 0.000311054, 0.000309712, 0.000308383, 0.000307064, 
			    0.000305757, 0.000304461, 0.000303176, 0.000301902, 0.000300639, 
			    0.000299386, 0.000298144, 0.000296912, 0.000295691, 0.000294479, 
			    0.000293278, 0.000292086, 0.000290904, 0.000289731, 0.000288568, 
			    0.000287414, 0.00028627, 0.000285135, 0.000284008, 0.000282891, 
			    0.000281782, 0.000280683, 0.000279591, 0.000278508, 0.000277434, 
			    0.000276368, 0.00027531, 0.00027426, 0.000273218, 0.000272184, 
			    0.000271158, 0.00027014, 0.000269129, 0.000268126, 0.00026713, 
			    0.000266142, 0.000265161, 0.000264187, 0.00026322, 0.000262261, 
			    0.000261308, 0.000260362, 0.000259423, 0.000258491, 0.000257566, 
			    0.000256647, 0.000255734, 0.000254829, 0.000253929, 0.000253036, 
			    0.000252149, 0.000251268, 0.000250394, 0.000249525, 0.000248663, 
			    0.000247806, 0.000246955, 0.000246111, 0.000245271, 0.000244438, 
			    0.00024361, 0.000242788, 0.000241971, 0.00024116, 0.000240354, 
			    0.000239554, 0.000238758, 0.000237968, 0.000237184, 0.000236404, 
			    0.00023563, 0.00023486, 0.000234096, 0.000233336, 0.000232581, 
			    0.000231831, 0.000231086, 0.000230346, 0.000229611, 0.00022888, 
			    0.000228153, 0.000227432, 0.000226714, 0.000226002, 0.000225293, 
			    0.00022459, 0.00022389, 0.000223195, 0.000222504, 0.000221817, 
			    0.000221135, 0.000220456, 0.000219782, 0.000219112, 0.000218446, 
			    0.000217784, 0.000217126, 0.000216471, 0.000215821, 0.000215175, 
			    0.000214532, 0.000213893, 0.000213258, 0.000212627, 0.000211999, 
			    0.000211375, 0.000210755, 0.000210138, 0.000209525, 0.000208915, 
			    0.000208309, 0.000207706, 0.000207107, 0.000206511, 0.000205919, 
			    0.00020533, 0.000204744, 0.000204162, 0.000203582, 0.000203006, 
			    0.000202434, 0.000201864, 0.000201298, 0.000200734, 0.000200174, 
			    0.000199617, 0.000199063, 0.000198512, 0.000197964, 0.000197419, 
			    0.000196877, 0.000196337, 0.000195801, 0.000195268, 0.000194737, 
			    0.00019421, 0.000193685, 0.000193163, 0.000192643, 0.000192127, 
			    0.000191613, 0.000191102, 0.000190593, 0.000190087, 0.000189584, 
			    0.000189084, 0.000188586, 0.00018809, 0.000187597, 0.000187107, 
			    0.000186619, 0.000186134, 0.000185651, 0.000185171, 0.000184693, 
			    0.000184217, 0.000183744, 0.000183273, 0.000182805, 0.000182339, 
			    0.000181875, 0.000181414, 0.000180955, 0.000180498, 0.000180043, 
			    0.000179591, 0.000179141, 0.000178693, 0.000178247, 0.000177803, 
			    0.000177362, 0.000176923, 0.000176486, 0.000176051, 0.000175618, 
			    0.000175187, 0.000174758, 0.000174331, 0.000173907, 0.000173484, 
			    0.000173063, 0.000172645, 0.000172228, 0.000171813, 0.000171401, 
			    0.00017099, 0.000170581, 0.000170174, 0.000169769, 0.000169365,
			    0.000168964, 0.000168565, 0.000168167, 0.000167771, 0.000167377, 
			    0.000166985, 0.000166594, 0.000166206, 0.000165819, 0.000165434, 
			    0.00016505, 0.000164668, 0.000164289, 0.00016391, 0.000163534, 
			    0.000163159, 0.000162786, 0.000162414, 0.000162044, 0.000161676, 
			    0.000161309, 0.000160944, 0.000160581, 0.000160219, 0.000159859, 
			    0.0001595, 0.000159143, 0.000158787, 0.000158433, 0.000158081, 
			    0.00015773, 0.000157381, 0.000157033, 0.000156686, 0.000156341, 
			    0.000155998, 0.000155656, 0.000155315, 0.000154976, 0.000154638, 
			    0.000154302, 0.000153967, 0.000153634, 0.000153302, 0.000152971, 
			    0.000152642, 0.000152314, 0.000151987, 0.000151662, 0.000151338, 
			    0.000151016, 0.000150694, 0.000150375, 0.000150056, 0.000149739, 
			    0.000149423, 0.000149108, 0.000148795, 0.000148483, 0.000148172, 
			    0.000147862, 0.000147554, 0.000147247, 0.000146941, 0.000146636, 
			    0.000146333, 0.000146031, 0.00014573, 0.00014543, 0.000145131, 
			    0.000144834, 0.000144538, 0.000144243, 0.000143949, 0.000143656, 
			    0.000143364, 0.000143074, 0.000142785, 0.000142496, 0.000142209, 
			    0.000141923, 0.000141638, 0.000141355, 0.000141072, 0.00014079, 
			    0.00014051, 0.00014023, 0.000139952, 0.000139675, 0.000139398, 
			    0.000139123, 0.000138849, 0.000138576, 0.000138304, 0.000138033, 
			    0.000137763, 0.000137494, 0.000137226, 0.000136959, 0.000136693, 
			    0.000136428, 0.000136164, 0.000135901, 0.000135639, 0.000135378, 
			    0.000135117, 0.000134858, 0.0001346, 0.000134343, 0.000134086, 
			    0.000133831, 0.000133576, 0.000133323, 0.00013307, 0.000132819, 
			    0.000132568, 0.000132318, 0.000132069, 0.000131821, 0.000131574, 
			    0.000131327, 0.000131082, 0.000130837, 0.000130594, 0.000130351, 
			    0.000130109, 0.000129868, 0.000129628, 0.000129388, 0.00012915, 
			    0.000128912, 0.000128675, 0.000128439, 0.000128204, 0.00012797, 
			    0.000127736, 0.000127503, 0.000127271, 0.00012704, 0.00012681,
			    0.00012658, 0.000126352, 0.000126124, 0.000125897, 0.00012567,
			    0.000125445, 0.00012522, 0.000124996, 0.000124772, 0.00012455, 
			    0.000124328, 0.000124107, 0.000123887, 0.000123667, 0.000123448, 
			    0.00012323, 0.000123013, 0.000122796, 0.000122581, 0.000122365, 
			    0.000122151, 0.000121937, 0.000121724, 0.000121512, 0.0001213, 
			    0.000121089, 0.000120879, 0.00012067, 0.000120461, 0.000120253, 
			    0.000120045, 0.000119839, 0.000119633, 0.000119427, 0.000119222, 
			    0.000119018, 0.000118815, 0.000118612, 0.00011841, 0.000118209, 
			    0.000118008, 0.000117808, 0.000117608, 0.000117409, 0.000117211, 
			    0.000117013, 0.000116817, 0.00011662, 0.000116425, 0.000116229, 
			    0.000116035, 0.000115841, 0.000115648, 0.000115455, 0.000115263, 
			    0.000115072, 0.000114881, 0.000114691, 0.000114501, 0.000114312, 
			    0.000114124, 0.000113936, 0.000113749, 0.000113562, 0.000113376, 
			    0.00011319, 0.000113005, 0.000112821, 0.000112637, 0.000112454, 
			    0.000112271, 0.000112089, 0.000111908, 0.000111727, 0.000111546, 
			    0.000111366, 0.000111187, 0.000111008, 0.00011083, 0.000110652,
			    0.000110475, 0.000110298, 0.000110122, 0.000109947, 0.000109771,
			    0.000109597, 0.000109423, 0.000109249, 0.000109076, 0.000108904, 
			    0.000108732, 0.000108561, 0.00010839, 0.000108219, 0.000108049, 
			    0.00010788, 0.000107711, 0.000107543, 0.000107375, 0.000107207, 
			    0.00010704, 0.000106874, 0.000106708, 0.000106542, 0.000106377, 
			    0.000106213, 0.000106049, 0.000105885, 0.000105722, 0.00010556, 
			    0.000105398, 0.000105236, 0.000105075, 0.000104914, 0.000104754, 
			    0.000104594, 0.000104434, 0.000104276, 0.000104117, 0.000103959, 
			    0.000103801, 0.000103644, 0.000103488, 0.000103331, 0.000103176, 
			    0.00010302, 0.000102865, 0.000102711, 0.000102557, 0.000102403, 
			    0.00010225, 0.000102097, 0.000101945, 0.000101793, 0.000101642, 
			    0.000101491, 0.00010134, 0.00010119, 0.00010104, 0.000100891, 
			    0.000100742, 0.000100593, 0.000100445, 0.000100297, 0.00010015, 
			    0.000100003, 0.0000998567, 0.0000997106, 0.0000995649, 0.0000994196, 
			    0.0000992747, 0.0000991302, 0.000098986, 0.0000988423, 0.0000986989, 
			    0.000098556, 0.0000984134, 0.0000982712, 0.0000981294, 0.000097988, 
			    0.0000978469, 0.0000977063, 0.000097566, 0.000097426, 0.0000972865, 
			    0.0000971473, 0.0000970085, 0.0000968701, 0.0000967321, 0.0000965944, 
			    0.000096457, 0.0000963201, 0.0000961835, 0.0000960472, 0.0000959113, 
			    0.0000957758, 0.0000956406, 0.0000955058, 0.0000953714, 0.0000952373, 
			    0.0000951035, 0.0000949701, 0.000094837, 0.0000947043, 0.000094572, 
			    0.0000944399, 0.0000943083, 0.0000941769, 0.0000940459, 0.0000939153, 
			    0.0000937849, 0.0000936549, 0.0000935253, 0.000093396, 0.000093267, 
			    0.0000931383, 0.00009301, 0.000092882, 0.0000927543, 0.000092627, 
			    0.0000924999, 0.0000923732, 0.0000922469, 0.0000921208, 0.0000919951, 
			    0.0000918696, 0.0000917445, 0.0000916197, 0.0000914953, 0.0000913711, 
			    0.0000912473, 0.0000911237, 0.0000910005, 0.0000908776, 0.000090755, 
			    0.0000906326, 0.0000905106, 0.0000903889, 0.0000902676, 0.0000901465, 
			    0.0000900257, 0.0000899052, 0.000089785, 0.0000896651, 0.0000895455, 
			    0.0000894262, 0.0000893072, 0.0000891884, 0.00008907, 0.0000889519, 
			    0.000088834, 0.0000887165, 0.0000885992, 0.0000884822, 0.0000883655,
			    0.0000882491, 0.0000881329, 0.0000880171, 0.0000879015, 0.0000877862, 
			    0.0000876712, 0.0000875565, 0.000087442, 0.0000873278, 0.0000872139, 
			    0.0000871003, 0.0000869869, 0.0000868739, 0.000086761, 0.0000866485, 
			    0.0000865362, 0.0000864242, 0.0000863125, 0.000086201, 0.0000860898, 
			    0.0000859788, 0.0000858681, 0.0000857577, 0.0000856475, 0.0000855376, 
			    0.000085428, 0.0000853186, 0.0000852095, 0.0000851006, 0.000084992, 
			    0.0000848837, 0.0000847756, 0.0000846677, 0.0000845601, 0.0000844528, 
			    0.0000843457, 0.0000842388, 0.0000841322, 0.0000840259, 0.0000839198, 
			    0.0000838139, 0.0000837083, 0.0000836029, 0.0000834978, 0.0000833929, 
			    0.0000832883, 0.0000831839, 0.0000830797, 0.0000829758, 0.0000828721, 
			    0.0000827687, 0.0000826654, 0.0000825625, 0.0000824597, 0.0000823572, 
			    0.000082255, 0.0000821529, 0.0000820511, 0.0000819496, 0.0000818482, 
			    0.0000817471, 0.0000816462, 0.0000815456, 0.0000814451, 0.0000813449,
			    0.0000812449, 0.0000811452, 0.0000810457, 0.0000809464, 0.0000808473,
			    0.0000807484, 0.0000806498, 0.0000805513, 0.0000804531, 0.0000803552,
			    0.0000802574, 0.0000801599, 0.0000800625, 0.0000799654, 0.0000798685, 
			    0.0000797718, 0.0000796754, 0.0000795791, 0.0000794831, 0.0000793872, 
			    0.0000792916, 0.0000791962, 0.000079101, 0.000079006, 0.0000789112, 
			    0.0000788166, 0.0000787223, 0.0000786281, 0.0000785341, 0.0000784404, 
			    0.0000783468, 0.0000782535, 0.0000781603, 0.0000780674, 0.0000779747, 
			    0.0000778821, 0.0000777898, 0.0000776976, 0.0000776057, 0.0000775139, 
			    0.0000774224, 0.000077331, 0.0000772399, 0.0000771489, 0.0000770582, 
			    0.0000769676, 0.0000768772, 0.000076787, 0.0000766971, 0.0000766073, 
			    0.0000765177, 0.0000764282, 0.000076339, 0.00007625, 0.0000761611, 
			    0.0000760725, 0.000075984, 0.0000758957, 0.0000758076, 0.0000757197, 
			    0.000075632, 0.0000755444, 0.0000754571, 0.0000753699, 0.0000752829, 
			    0.0000751961, 0.0000751094, 0.000075023, 0.0000749367, 0.0000748506, 
			    0.0000747647, 0.000074679, 0.0000745935, 0.0000745081, 0.0000744229, 
			    0.0000743379, 0.000074253, 0.0000741683, 0.0000740839, 0.0000739995, 
			    0.0000739154, 0.0000738314, 0.0000737476, 0.000073664, 0.0000735805, 
			    0.0000734972, 0.0000734141, 0.0000733312, 0.0000732484, 0.0000731658, 
			    0.0000730834, 0.0000730011, 0.000072919, 0.0000728371, 0.0000727553, 
			    0.0000726737, 0.0000725923, 0.000072511, 0.0000724299, 0.0000723489, 
			    0.0000722682, 0.0000721875, 0.0000721071, 0.0000720268, 0.0000719467, 
			    0.0000718667, 0.0000717869, 0.0000717072, 0.0000716278, 0.0000715484, 
			    0.0000714692, 0.0000713902, 0.0000713114, 0.0000712327, 0.0000711541, 
			    0.0000710758, 0.0000709975, 0.0000709194, 0.0000708415, 0.0000707638, 
			    0.0000706861, 0.0000706087, 0.0000705314, 0.0000704542, 0.0000703772, 
			    0.0000703004, 0.0000702237, 0.0000701471, 0.0000700707, 0.0000699945, 
			    0.0000699184, 0.0000698424, 0.0000697666, 0.000069691, 0.0000696155, 
			    0.0000695401, 0.0000694649, 0.0000693898, 0.0000693149, 0.0000692402,
			    0.0000691655, 0.000069091, 0.0000690167, 0.0000689425, 0.0000688684, 
			    0.0000687945, 0.0000687208, 0.0000686471, 0.0000685737, 0.0000685003, 
			    0.0000684271, 0.0000683541, 0.0000682811, 0.0000682084, 0.0000681357, 
			    0.0000680632, 0.0000679908, 0.0000679186, 0.0000678465, 0.0000677746, 
			    0.0000677028, 0.0000676311, 0.0000675595, 0.0000674881, 0.0000674169, 
			    0.0000673457, 0.0000672747, 0.0000672039, 0.0000671331, 0.0000670625, 
			    0.0000669921, 0.0000669217, 0.0000668515, 0.0000667815, 0.0000667115, 
			    0.0000666417, 0.000066572, 0.0000665025, 0.0000664331, 0.0000663638, 
			    0.0000662946, 0.0000662256, 0.0000661567, 0.0000660879, 0.0000660193, 
			    0.0000659508, 0.0000658824, 0.0000658141, 0.000065746, 0.000065678, 
			    0.0000656101, 0.0000655423, 0.0000654747, 0.0000654072, 0.0000653398, 
			    0.0000652725, 0.0000652054, 0.0000651384, 0.0000650715, 0.0000650047, 
			    0.0000649381, 0.0000648716, 0.0000648052, 0.0000647389, 0.0000646727,
			    0.0000646067, 0.0000645408, 0.000064475, 0.0000644093, 0.0000643437, 
			    0.0000642783, 0.000064213, 0.0000641478, 0.0000640827, 0.0000640177, 
			    0.0000639528, 0.0000638881, 0.0000638235, 0.000063759, 0.0000636946, 
			    0.0000636303, 0.0000635662, 0.0000635022, 0.0000634382, 0.0000633744,
			    0.0000633107, 0.0000632471, 0.0000631837, 0.0000631203, 0.0000630571, 
			    0.000062994, 0.0000629309, 0.000062868, 0.0000628052, 0.0000627426, 
			    0.00006268, 0.0000626175, 0.0000625552, 0.000062493, 0.0000624308, 
			    0.0000623688, 0.0000623069, 0.0000622451, 0.0000621834, 0.0000621218, 
			    0.0000620604, 0.000061999, 0.0000619377, 0.0000618766, 0.0000618156,
			    0.0000617546, 0.0000616938, 0.0000616331, 0.0000615725, 0.0000615119, 
			    0.0000614515, 0.0000613912, 0.000061331, 0.0000612709, 0.000061211, 
			    0.0000611511, 0.0000610913, 0.0000610316, 0.0000609721, 0.0000609126, 
			    0.0000608532, 0.000060794, 0.0000607348, 0.0000606757, 0.0000606168, 
			    0.0000605579, 0.0000604992, 0.0000604405, 0.000060382, 0.0000603235, 
			    0.0000602652, 0.0000602069, 0.0000601488, 0.0000600907, 0.0000600328, 
			    0.0000599749, 0.0000599171, 0.0000598595, 0.0000598019, 0.0000597445, 
			    0.0000596871, 0.0000596298, 0.0000595727, 0.0000595156, 0.0000594586, 
			    0.0000594017, 0.000059345, 0.0000592883, 0.0000592317, 0.0000591752, 
			    0.0000591188, 0.0000590625, 0.0000590063, 0.0000589501, 0.0000588941, 
			    0.0000588382, 0.0000587823, 0.0000587266, 0.000058671, 0.0000586154, 
			    0.0000585599, 0.0000585046, 0.0000584493, 0.0000583941, 0.000058339, 
			    0.000058284, 0.0000582291, 0.0000581743, 0.0000581195, 0.0000580649, 
			    0.0000580104, 0.0000579559, 0.0000579015, 0.0000578473, 0.0000577931, 
			    0.000057739, 0.000057685, 0.000057631, 0.0000575772, 0.0000575235, 
			    0.0000574698, 0.0000574162, 0.0000573628, 0.0000573094, 0.0000572561, 
			    0.0000572028, 0.0000571497, 0.0000570967, 0.0000570437, 0.0000569908, 
			    0.0000569381, 0.0000568854, 0.0000568328, 0.0000567802, 0.0000567278, 
			    0.0000566754, 0.0000566232, 0.000056571, 0.0000565189, 0.0000564669, 
			    0.0000564149, 0.0000563631, 0.0000563113};

    
    double r = rand() % 1000000000;
    r /= 1000000000;

    // Determine which bin the event falls into in the ISR spectrum. Each bin is 1 MeV wide.
    int bin = (int) floor(mc_gam_isr_boost_e/(0.2*m_eBeam/2) * 1200) + 1; //if energy = 1, make sure you use the correct bin
    double prob = weights[bin-1]/weights[0];
    /*
    cout << "#$#$#" << endl;
    cout << "e =    " << mc_gam_isr_boost_e << endl;
    cout << "bin =  " << bin << endl;
    cout << "r =    " << r << endl;
    cout << "prob = " << prob << endl;
    */

    if (r < prob) {
      //cout << "included" << endl;
      m_mc_reweight = 1;
    }
    else {
      m_mc_reweight = 0;
    }
  }








  void Analysis::clearValuesMC() {


      m_mc_fsr_YnS_e = 0;
      m_mc_fsr_YnS_boost_e = 0;
      m_mc_fsr_rho_e = 0;
      m_mc_fsr_rho_boost_e = 0;

      //ISR MC
      m_isr = 0;
      m_mc_reweight = 0;

      m_mc_YNS_isr_e     = 0;
      m_mc_YNS_isr_p     = 0;
      m_mc_YNS_isr_m     = 0;
      m_mc_YNS_isr_pt    = 0;
      m_mc_YNS_isr_phi   = 0;
      m_mc_YNS_isr_costh = 0;

      m_mc_gam_isr_recoil_e     = 0;
      m_mc_gam_isr_recoil_p     = 0;
      m_mc_gam_isr_recoil_m     = 0;
      m_mc_gam_isr_recoil_pt    = 0;
      m_mc_gam_isr_recoil_phi   = 0;
      m_mc_gam_isr_recoil_costh = 0;
      
      m_mc_gam_isr_e     = 0;
      m_mc_gam_isr_p     = 0;
      m_mc_gam_isr_m     = 0;
      m_mc_gam_isr_pt    = 0;
      m_mc_gam_isr_phi   = 0;
      m_mc_gam_isr_costh = 0;

      m_mc_gam_isr_boost_e     = 0;
      m_mc_gam_isr_boost_p     = 0;
      m_mc_gam_isr_boost_m     = 0;
      m_mc_gam_isr_boost_pt    = 0;
      m_mc_gam_isr_boost_phi   = 0;
      m_mc_gam_isr_boost_costh = 0;
      
      m_mc_YnS_isr_e     = 0;
      m_mc_YnS_isr_p     = 0;
      m_mc_YnS_isr_m     = 0;
      m_mc_YnS_isr_pt    = 0;
      m_mc_YnS_isr_phi   = 0;
      m_mc_YnS_isr_costh = 0;
      
      m_mc_pip_pim_isr_e     = 0;
      m_mc_pip_pim_isr_p     = 0;
      m_mc_pip_pim_isr_m     = 0;
      m_mc_pip_pim_isr_pt    = 0;
      m_mc_pip_pim_isr_phi   = 0;
      m_mc_pip_pim_isr_costh = 0;
      
      m_mc_pip_isr_e     = 0;
      m_mc_pip_isr_p     = 0;
      m_mc_pip_isr_m     = 0;
      m_mc_pip_isr_pt    = 0;
      m_mc_pip_isr_phi   = 0;
      m_mc_pip_isr_costh = 0;
      
      m_mc_pim_isr_e     = 0;
      m_mc_pim_isr_p     = 0;
      m_mc_pim_isr_m     = 0;
      m_mc_pim_isr_pt    = 0;
      m_mc_pim_isr_phi   = 0;
      m_mc_pim_isr_costh = 0;
      
      m_mc_mup_isr_e     = 0;
      m_mc_mup_isr_p     = 0;
      m_mc_mup_isr_m     = 0;
      m_mc_mup_isr_pt    = 0;
      m_mc_mup_isr_phi   = 0;
      m_mc_mup_isr_costh = 0;
      
      m_mc_mum_isr_e     = 0;
      m_mc_mum_isr_p     = 0;
      m_mc_mum_isr_m     = 0;
      m_mc_mum_isr_pt    = 0;
      m_mc_mum_isr_phi   = 0;
      m_mc_mum_isr_costh = 0;

      m_mc_pip_pim_gam_isr_recoil_e     = 0;
      m_mc_pip_pim_gam_isr_recoil_p     = 0;
      m_mc_pip_pim_gam_isr_recoil_m     = 0;
      m_mc_pip_pim_gam_isr_recoil_pt    = 0;
      m_mc_pip_pim_gam_isr_recoil_phi   = 0;
      m_mc_pip_pim_gam_isr_recoil_costh = 0;

      m_mc_mup_mum_gam_isr_recoil_e     = 0;
      m_mc_mup_mum_gam_isr_recoil_p     = 0;
      m_mc_mup_mum_gam_isr_recoil_m     = 0;
      m_mc_mup_mum_gam_isr_recoil_pt    = 0;
      m_mc_mup_mum_gam_isr_recoil_phi   = 0;
      m_mc_mup_mum_gam_isr_recoil_costh = 0;

      m_mc_pip_pim_mup_mum_isr_recoil_e     = 0;
      m_mc_pip_pim_mup_mum_isr_recoil_p     = 0;
      m_mc_pip_pim_mup_mum_isr_recoil_m     = 0;
      m_mc_pip_pim_mup_mum_isr_recoil_pt    = 0;
      m_mc_pip_pim_mup_mum_isr_recoil_phi   = 0;
      m_mc_pip_pim_mup_mum_isr_recoil_costh = 0;

      m_mc_pip_pim_mup_mum_isr_recoil_m2 = 0;

      m_mc_pip_pim_e     = 0;
      m_mc_pip_pim_p     = 0;
      m_mc_pip_pim_m     = 0;
      m_mc_pip_pim_pt    = 0;
      m_mc_pip_pim_phi   = 0;
      m_mc_pip_pim_costh = 0;


      // Non-ISR MC
      // Y5StoY1Smumu2pic
      m_mc_pip_pim_e     = 0;
      m_mc_pip_pim_p     = 0;
      m_mc_pip_pim_m     = 0;
      m_mc_pip_pim_pt    = 0;
      m_mc_pip_pim_phi   = 0;
      m_mc_pip_pim_costh = 0;

      // Signal MC
      m_mc_Y5S_e     = 0;
      m_mc_Y5S_p     = 0;
      m_mc_Y5S_m     = 0;
      m_mc_Y5S_pt    = 0;
      m_mc_Y5S_phi   = 0;
      m_mc_Y5S_costh = 0;

      m_mc_Wbj_e     = 0;
      m_mc_Wbj_p     = 0;
      m_mc_Wbj_m     = 0;
      m_mc_Wbj_pt    = 0;
      m_mc_Wbj_phi   = 0;
      m_mc_Wbj_costh = 0;

      m_mc_Wbj_boost_e     = 0;
      m_mc_Wbj_boost_p     = 0;
      m_mc_Wbj_boost_m     = 0;
      m_mc_Wbj_boost_pt    = 0;
      m_mc_Wbj_boost_phi   = 0;
      m_mc_Wbj_boost_costh = 0;

      m_mc_gam_recoil_e     = 0;
      m_mc_gam_recoil_p     = 0;
      m_mc_gam_recoil_m     = 0;
      m_mc_gam_recoil_pt    = 0;
      m_mc_gam_recoil_phi   = 0;
      m_mc_gam_recoil_costh = 0;
      
      m_mc_gam_e     = 0;
      m_mc_gam_p     = 0;
      m_mc_gam_m     = 0;
      m_mc_gam_pt    = 0;
      m_mc_gam_phi   = 0;
      m_mc_gam_costh = 0;
      
      m_mc_YnS_e     = 0;
      m_mc_YnS_p     = 0;
      m_mc_YnS_m     = 0;
      m_mc_YnS_pt    = 0;
      m_mc_YnS_phi   = 0;
      m_mc_YnS_costh = 0;

      m_mc_YnS_boost_e     = 0;
      m_mc_YnS_boost_p     = 0;
      m_mc_YnS_boost_m     = 0;
      m_mc_YnS_boost_pt    = 0;
      m_mc_YnS_boost_phi   = 0;
      m_mc_YnS_boost_costh = 0;
      
      m_mc_rho_e     = 0;
      m_mc_rho_p     = 0;
      m_mc_rho_m     = 0;
      m_mc_rho_pt    = 0;
      m_mc_rho_phi   = 0;
      m_mc_rho_costh = 0;
      
      m_mc_pip_e     = 0;
      m_mc_pip_p     = 0;
      m_mc_pip_m     = 0;
      m_mc_pip_pt    = 0;
      m_mc_pip_phi   = 0;
      m_mc_pip_costh = 0;
      
      m_mc_pim_e     = 0;
      m_mc_pim_p     = 0;
      m_mc_pim_m     = 0;
      m_mc_pim_pt    = 0;
      m_mc_pim_phi   = 0;
      m_mc_pim_costh = 0;
      
      m_mc_mup_e     = 0;
      m_mc_mup_p     = 0;
      m_mc_mup_m     = 0;
      m_mc_mup_pt    = 0;
      m_mc_mup_phi   = 0;
      m_mc_mup_costh = 0;
      
      m_mc_mum_e     = 0;
      m_mc_mum_p     = 0;
      m_mc_mum_m     = 0;
      m_mc_mum_pt    = 0;
      m_mc_mum_phi   = 0;
      m_mc_mum_costh = 0;

      m_mc_pip_pim_gam_recoil_e     = 0;
      m_mc_pip_pim_gam_recoil_p     = 0;
      m_mc_pip_pim_gam_recoil_m     = 0;
      m_mc_pip_pim_gam_recoil_pt    = 0;
      m_mc_pip_pim_gam_recoil_phi   = 0;
      m_mc_pip_pim_gam_recoil_costh = 0;

      m_mc_pip_pim_gam_recoil_boost_e     = 0;
      m_mc_pip_pim_gam_recoil_boost_p     = 0;
      m_mc_pip_pim_gam_recoil_boost_m     = 0;
      m_mc_pip_pim_gam_recoil_boost_pt    = 0;
      m_mc_pip_pim_gam_recoil_boost_phi   = 0;
      m_mc_pip_pim_gam_recoil_boost_costh = 0;

      m_mc_mup_mum_gam_recoil_e     = 0;
      m_mc_mup_mum_gam_recoil_p     = 0;
      m_mc_mup_mum_gam_recoil_m     = 0;
      m_mc_mup_mum_gam_recoil_pt    = 0;
      m_mc_mup_mum_gam_recoil_phi   = 0;
      m_mc_mup_mum_gam_recoil_costh = 0;

      m_mc_mup_mum_gam_fit_recoil_m     = 0;

      m_mc_mup_mum_gam_recoil_boost_e     = 0;
      m_mc_mup_mum_gam_recoil_boost_p     = 0;
      m_mc_mup_mum_gam_recoil_boost_m     = 0;
      m_mc_mup_mum_gam_recoil_boost_pt    = 0;
      m_mc_mup_mum_gam_recoil_boost_phi   = 0;
      m_mc_mup_mum_gam_recoil_boost_costh = 0;

      m_mc_pip_pim_mup_mum_recoil_e     = 0;
      m_mc_pip_pim_mup_mum_recoil_p     = 0;
      m_mc_pip_pim_mup_mum_recoil_m     = 0;
      m_mc_pip_pim_mup_mum_recoil_pt    = 0;
      m_mc_pip_pim_mup_mum_recoil_phi   = 0;
      m_mc_pip_pim_mup_mum_recoil_costh = 0;

      m_mc_pip_pim_mup_mum_recoil_m2 = 0;

      m_mc_pip_pim_mup_mum_recoil_boost_e     = 0;
      m_mc_pip_pim_mup_mum_recoil_boost_p     = 0;
      m_mc_pip_pim_mup_mum_recoil_boost_m     = 0;
      m_mc_pip_pim_mup_mum_recoil_boost_pt    = 0;
      m_mc_pip_pim_mup_mum_recoil_boost_phi   = 0;
      m_mc_pip_pim_mup_mum_recoil_boost_costh = 0;

      m_mc_gam_boost_e     = 0;
      m_mc_gam_boost_p     = 0;
      m_mc_gam_boost_m     = 0;
      m_mc_gam_boost_pt    = 0;
      m_mc_gam_boost_phi   = 0;
      m_mc_gam_boost_costh = 0;

      //m_mc_gam_recoil_m     = 0;

      m_mc_mup_mum_gam_m = 0;
      m_mc_pip_pim_gam_m = 0;

      m_mc_fsr_YnS= 0;
      m_mc_fsr_rho= 0;

  }


#if defined(BELLE_NAMESPACE)
}
#endif
