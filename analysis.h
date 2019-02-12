#include <iostream>                           

#include "belle.h"                           

#include "event/BelleEvent.h"                
#include "tuple/BelleTupleManager.h"          
#include "basf/module.h"                      
#include "basf/module_descr.h"                

#ifdef HAVE_BOOST_FOREACH                     
#undef HAVE_BOOST_FOREACH
#endif
#include "panther/panther.h"                  

#include MDST_H                               
#include BELLETDF_H                           
#include HEPEVT_H                             

#include "particle/utility.h"                 
#include "particle/combination.h"             
#include "kid/atc_pid.h"                     
#include "mdst/mdst.h"                        
#include "mdst/Muid_mdst.h"                   
#include "eid/eid.h"                          
#include "ip/IpProfile.h"                     
#include "benergy/BeamEnergy.h"              

#include "hamlet/Hamlet.h"                    
#include "hamlet/Fbtag_MultDimLikelihood0.h"  

#include <toolbox/FoxWolfr.h>                 
#include "toolbox/Thrust.h"                   
#include "toolbox/FuncPtr.h"                  
#include "particle/gammac.h"                  

#include <TFile.h>                            
#include <TTree.h>                            
#include <TH1D.h>                             


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  
  using namespace std;
  
  class Analysis : public Module {
  public:
    Analysis () : m_event_seq(0), m_run_seq(0), count_Y5S(0), count_Wbj(0), count_YnS(0), count_rho(0), count_pip(0), count_pim(0), count_mup(0), count_mum(0), m_count_gam_from_pi0(0), m_total_rho(0), m_total_YnS(0), m_total_Wbj(0) {
      cout << "Analysis constructor" << endl;
    }
      ~Analysis (void) {};
      
      void init(int *status){
	
	Ptype ptype_dummy("VPHO");
	
	cout << "input file name: " << inpath << endl;
	cout << "writing ROOT output into a file with name " << outpath << endl;

	// ROOT File 
	m_outputRootFile = new TFile(outpath, "RECREATE", "Analysis ntuple");
	
	// Create TTree for reconstructed and MC data
	tree_rec = new TTree("rec",  "Y5S candidates" );
	tree_mc  = new TTree("mc",   "Monte Carlo Truth");
	tree_bkg  = new TTree("bkg",   "Background Stuff");
	
	// see, e.g., https://root.cern.ch/root/html/tutorials/math/mathcoreVectorCollection.C.html
	
	//---------------------------------------------------------------------------------------------------------------

	tree_rec->Branch("angle_bt_pip_gam"         , "std::vector<double>"              , &nt_angle_bt_pip_gam   );
	tree_rec->Branch("angle_bt_pim_gam"         , "std::vector<double>"              , &nt_angle_bt_pim_gam   );
	tree_rec->Branch("angle_bt_mup_gam"         , "std::vector<double>"              , &nt_angle_bt_mup_gam   );
	tree_rec->Branch("angle_bt_mum_gam"         , "std::vector<double>"              , &nt_angle_bt_mum_gam   );

	//---------------------------------------------------------------------------------------------------------------
	// Relevant run information
 
	tree_rec->Branch("istrig"         , &m_istrig              , "istrig/I"             );
	tree_rec->Branch("trigbits"       , &m_trigbits            , "trigbits/I"           );

	tree_rec->Branch("exp"            , &m_exp                  , "exp/I"                );
	tree_rec->Branch("run"            , &m_event_run            , "m_event_run/I"        );
	tree_rec->Branch("runseq"         , &m_run_seq              , "runseq/I"             );
	tree_rec->Branch("event"          , &m_event                , "event/I"              );
	tree_rec->Branch("eventseq"       , &m_event_seq            , "eventseq/I"           );

	tree_rec->Branch("ebeam"          , &m_eBeam                , "ebeam/D"              );
	tree_rec->Branch("ebeamCM"        , &m_eBeamCM              , "ebeamCM/D"            );
	tree_rec->Branch("pbeam"          , &m_pBeamCM              , "pbeam/D"              );

	tree_rec->Branch("NtrkGood"       , &m_NChargedTracks       , "NtrkGood/I"           );

	tree_rec->Branch("Ntrk"           , &m_evt_Ntrk             , "Ntrk/I"               );
	tree_rec->Branch("Ncls"           , &m_evt_Ncls             , "Ncls/I"               );
	tree_rec->Branch("Psum"           , &m_evt_Psum             , "Psum/D"               );
	tree_rec->Branch("Esum"           , &m_evt_Esum             , "Esum/D"               );
	tree_rec->Branch("Evis"           , &m_evt_Evis             , "Evis/D"               );
	tree_rec->Branch("Pz"             , &m_evt_Pz               , "Pz/D"                 );
	tree_rec->Branch("HeavyJetMass"   , &m_evt_HeavyJetMass     , "HeavyJetMass/D"       );
	tree_rec->Branch("Thrust"         , &m_evt_Thrust           , "Thrust/D"             );
	tree_rec->Branch("R2"             , &m_evt_R2               , "R2/D"                 );
	tree_rec->Branch("HadronBflag"    , &m_HadronB_flag         , "HadronBflag/I"        );
	tree_rec->Branch("HadronAflag"    , &m_HadronA_flag         , "HadronAflag/I"        );
	tree_rec->Branch("Tauflag"        , &m_Tau_flag             , "Tauflag/I"            );

	tree_rec->Branch("mc"             , &m_mc                   , "mc/I"                 );
	tree_rec->Branch("mcs"            , &m_mcs                  , "mcs/I"                );

	//---------------------------------------------------------------------------------------------------------------
	// MC counting statistics

	tree_rec->Branch( "cand"           , "std::vector<int>"     , &nt_index              );
	tree_rec->Branch( "ncand"          , &m_ncand               , "ncand/I"              );
	tree_rec->Branch( "ncand_sr"       , &m_ncand_sr            , "ncand_sr/I"           );
	tree_rec->Branch( "ncand_mct"      , &m_ncand_mct           , "ncand_mct/I"          );
	tree_rec->Branch( "ncand_gam"      , &m_ncand_gam           , "ncand_gam/I"          );
	tree_rec->Branch( "ncand_pi0"      , &m_ncand_pi0           , "ncand_pi0/I"          );

	tree_rec->Branch( "total_rho"          , &m_total_rho               , "total_rho/I"              );
	tree_rec->Branch( "total_YnS"          , &m_total_YnS               , "total_YnS/I"              );
	tree_rec->Branch( "total_Wbj"          , &m_total_Wbj               , "total_Wbj/I"              );
	
	tree_rec->Branch( "gam_from_pi0"   , "std::vector<int>"     , &nt_gam_from_pi0            );

	tree_rec->Branch( "sr"             , "std::vector<int>"     , &nt_sr             );
  
	tree_rec->Branch( "gam_rec"      , &m_count_gam_rec     , "gam_rec/I"      );
	tree_rec->Branch( "mup_rec"      , &m_count_mup_rec     , "mup_rec/I"    );
	tree_rec->Branch( "mum_rec"      , &m_count_mum_rec     , "mum_rec/I"    );
	tree_rec->Branch( "pip_rec"      , &m_count_pip_rec     , "pip_rec/I"    );
	tree_rec->Branch( "pim_rec"      , &m_count_pim_rec     , "pim_rec/I"    );

	tree_rec->Branch( "gam_good_trk"      , &m_count_gam_good_trk     , "gam_good_trk/I"      );
	tree_rec->Branch( "mup_good_trk"      , &m_count_mup_good_trk     , "mup_good_trk/I"    );
	tree_rec->Branch( "mum_good_trk"      , &m_count_mum_good_trk     , "mum_good_trk/I"    );
	tree_rec->Branch( "pip_good_trk"      , &m_count_pip_good_trk     , "pip_good_trk/I"    );
	tree_rec->Branch( "pim_good_trk"      , &m_count_pim_good_trk     , "pim_good_trk/I"    );

	tree_rec->Branch( "gam_good_dr"      , &m_count_gam_good_dr     , "gam_good_dr/I"      );
	tree_rec->Branch( "mup_good_dr"      , &m_count_mup_good_dr     , "mup_good_dr/I"    );
	tree_rec->Branch( "mum_good_dr"      , &m_count_mum_good_dr     , "mum_good_dr/I"    );
	tree_rec->Branch( "pip_good_dr"      , &m_count_pip_good_dr     , "pip_good_dr/I"    );
	tree_rec->Branch( "pim_good_dr"      , &m_count_pim_good_dr     , "pim_good_dr/I"    );

	tree_rec->Branch( "gam_good_dz"      , &m_count_gam_good_dz     , "gam_good_dz/I"      );
	tree_rec->Branch( "mup_good_dz"      , &m_count_mup_good_dz     , "mup_good_dz/I"    );
	tree_rec->Branch( "mum_good_dz"      , &m_count_mum_good_dz     , "mum_good_dz/I"    );
	tree_rec->Branch( "pip_good_dz"      , &m_count_pip_good_dz     , "pip_good_dz/I"    );
	tree_rec->Branch( "pim_good_dz"      , &m_count_pim_good_dz     , "pim_good_dz/I"    );

	tree_rec->Branch( "gam_good_pt"      , &m_count_gam_good_pt     , "gam_good_pt/I"      );
	tree_rec->Branch( "mup_good_pt"      , &m_count_mup_good_pt     , "mup_good_pt/I"    );
	tree_rec->Branch( "mum_good_pt"      , &m_count_mum_good_pt     , "mum_good_pt/I"    );
	tree_rec->Branch( "pip_good_pt"      , &m_count_pip_good_pt     , "pip_good_pt/I"    );
	tree_rec->Branch( "pim_good_pt"      , &m_count_pim_good_pt     , "pim_good_pt/I"    );

	tree_rec->Branch( "gam_eff"      , &m_count_gam_eff     , "gam_eff/I"      );
	tree_rec->Branch( "mup_eff"      , &m_count_mup_eff     , "mup_eff/I"    );
	tree_rec->Branch( "mum_eff"      , &m_count_mum_eff     , "mum_eff/I"    );
	tree_rec->Branch( "pip_eff"      , &m_count_pip_eff     , "pip_eff/I"    );
	tree_rec->Branch( "pim_eff"      , &m_count_pim_eff     , "pim_eff/I"    );

	tree_rec->Branch( "mc_fsr_YnS"  , &m_mc_fsr_YnS        , "mc_fsr_YnS/I"      );
	tree_rec->Branch( "mc_fsr_rho"  , &m_mc_fsr_rho        , "mc_fsr_rho/I"      );

	//---------------------------------------------------------------------------------------------------------------
	// MC tagging

	tree_rec->Branch( "mct"            , "std::vector<int>"     , &nt_mct            );
	tree_rec->Branch( "gam_mct"        , "std::vector<int>"     , &nt_gam_mct        );
	tree_rec->Branch( "mup_mct"        , "std::vector<int>"     , &nt_mup_mct        );
	tree_rec->Branch( "mum_mct"        , "std::vector<int>"     , &nt_mum_mct        );
	tree_rec->Branch( "pip_mct"        , "std::vector<int>"     , &nt_pip_mct        );
	tree_rec->Branch( "pim_mct"        , "std::vector<int>"     , &nt_pim_mct        );

	//---------------------------------------------------------------------------------------------------------------
	// Quantities for background MC samples with initial state radiation / radiative returns to lower mass Upsilon

	tree_rec->Branch( "mc_YNS_isr_e"        , &m_mc_YNS_isr_e              , "mc_YNS_isr_e/D"            );
	tree_rec->Branch( "mc_YNS_isr_p"        , &m_mc_YNS_isr_p              , "mc_YNS_isr_p/D"            );
	tree_rec->Branch( "mc_YNS_isr_m"        , &m_mc_YNS_isr_m              , "mc_YNS_isr_m/D"            );
	tree_rec->Branch( "mc_YNS_isr_pt"       , &m_mc_YNS_isr_pt             , "mc_YNS_isr_pt/D"           );
	tree_rec->Branch( "mc_YNS_isr_phi"      , &m_mc_YNS_isr_phi            , "mc_YNS_isr_phi/D"          );
	tree_rec->Branch( "mc_YNS_isr_costh"    , &m_mc_YNS_isr_costh          , "mc_YNS_isr_costh/D"        );

	tree_rec->Branch( "mc_YnS_isr_e"        , &m_mc_YnS_isr_e              , "mc_YnS_isr_e/D"            );
	tree_rec->Branch( "mc_YnS_isr_p"        , &m_mc_YnS_isr_p              , "mc_YnS_isr_p/D"            );
	tree_rec->Branch( "mc_YnS_isr_m"        , &m_mc_YnS_isr_m              , "mc_YnS_isr_m/D"            );
	tree_rec->Branch( "mc_YnS_isr_pt"       , &m_mc_YnS_isr_pt             , "mc_YnS_isr_pt/D"           );
	tree_rec->Branch( "mc_YnS_isr_phi"      , &m_mc_YnS_isr_phi            , "mc_YnS_isr_phi/D"          );
	tree_rec->Branch( "mc_YnS_isr_costh"    , &m_mc_YnS_isr_costh          , "mc_YnS_isr_costh/D"        );

	tree_rec->Branch( "mc_pip_isr_e"        , &m_mc_pip_isr_e              , "mc_pip_isr_e/D"            );
	tree_rec->Branch( "mc_pip_isr_p"        , &m_mc_pip_isr_p              , "mc_pip_isr_p/D"            );
	tree_rec->Branch( "mc_pip_isr_m"        , &m_mc_pip_isr_m              , "mc_pip_isr_m/D"            );
	tree_rec->Branch( "mc_pip_isr_pt"       , &m_mc_pip_isr_pt             , "mc_pip_isr_pt/D"           );
	tree_rec->Branch( "mc_pip_isr_phi"      , &m_mc_pip_isr_phi            , "mc_pip_isr_phi/D"          );
	tree_rec->Branch( "mc_pip_isr_costh"    , &m_mc_pip_isr_costh          , "mc_pip_isr_costh/D"        );

	tree_rec->Branch( "mc_pim_isr_e"        , &m_mc_pim_isr_e              , "mc_pim_isr_e/D"            );
	tree_rec->Branch( "mc_pim_isr_p"        , &m_mc_pim_isr_p              , "mc_pim_isr_p/D"            );
	tree_rec->Branch( "mc_pim_isr_m"        , &m_mc_pim_isr_m              , "mc_pim_isr_m/D"            );
	tree_rec->Branch( "mc_pim_isr_pt"       , &m_mc_pim_isr_pt             , "mc_pim_isr_pt/D"           );
	tree_rec->Branch( "mc_pim_isr_phi"      , &m_mc_pim_isr_phi            , "mc_pim_isr_phi/D"          );
	tree_rec->Branch( "mc_pim_isr_costh"    , &m_mc_pim_isr_costh          , "mc_pim_isr_costh/D"        );

	tree_rec->Branch( "mc_mup_isr_e"        , &m_mc_mup_isr_e              , "mc_mup_isr_e/D"            );
	tree_rec->Branch( "mc_mup_isr_p"        , &m_mc_mup_isr_p              , "mc_mup_isr_p/D"            );
	tree_rec->Branch( "mc_mup_isr_m"        , &m_mc_mup_isr_m              , "mc_mup_isr_m/D"            );
	tree_rec->Branch( "mc_mup_isr_pt"       , &m_mc_mup_isr_pt             , "mc_mup_isr_pt/D"           );
	tree_rec->Branch( "mc_mup_isr_phi"      , &m_mc_mup_isr_phi            , "mc_mup_isr_phi/D"          );
	tree_rec->Branch( "mc_mup_isr_costh"    , &m_mc_mup_isr_costh          , "mc_mup_isr_costh/D"        );

	tree_rec->Branch( "mc_mum_isr_e"        , &m_mc_mum_isr_e              , "mc_mum_isr_e/D"            );
	tree_rec->Branch( "mc_mum_isr_p"        , &m_mc_mum_isr_p              , "mc_mum_isr_p/D"            );
	tree_rec->Branch( "mc_mum_isr_m"        , &m_mc_mum_isr_m              , "mc_mum_isr_m/D"            );
	tree_rec->Branch( "mc_mum_isr_pt"       , &m_mc_mum_isr_pt             , "mc_mum_isr_pt/D"           );
	tree_rec->Branch( "mc_mum_isr_phi"      , &m_mc_mum_isr_phi            , "mc_mum_isr_phi/D"          );
	tree_rec->Branch( "mc_mum_isr_costh"    , &m_mc_mum_isr_costh          , "mc_mum_isr_costh/D"        );

	tree_rec->Branch( "mc_pip_pim_isr_e"        , &m_mc_pip_pim_isr_e              , "mc_pip_pim_isr_e/D"            );
	tree_rec->Branch( "mc_pip_pim_isr_p"        , &m_mc_pip_pim_isr_p              , "mc_pip_pim_isr_p/D"            );
	tree_rec->Branch( "mc_pip_pim_isr_m"        , &m_mc_pip_pim_isr_m              , "mc_pip_pim_isr_m/D"            );
	tree_rec->Branch( "mc_pip_pim_isr_pt"       , &m_mc_pip_pim_isr_pt             , "mc_pip_pim_isr_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_isr_phi"      , &m_mc_pip_pim_isr_phi            , "mc_pip_pim_isr_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_isr_costh"    , &m_mc_pip_pim_isr_costh          , "mc_pip_pim_isr_costh/D"        );

	tree_rec->Branch( "mc_gam_isr_e"        , &m_mc_gam_isr_e              , "mc_gam_isr_e/D"            );
	tree_rec->Branch( "mc_gam_isr_p"        , &m_mc_gam_isr_p              , "mc_gam_isr_p/D"            );
	tree_rec->Branch( "mc_gam_isr_m"        , &m_mc_gam_isr_m              , "mc_gam_isr_m/D"            );
	tree_rec->Branch( "mc_gam_isr_pt"       , &m_mc_gam_isr_pt             , "mc_gam_isr_pt/D"           );
	tree_rec->Branch( "mc_gam_isr_phi"      , &m_mc_gam_isr_phi            , "mc_gam_isr_phi/D"          );
	tree_rec->Branch( "mc_gam_isr_costh"    , &m_mc_gam_isr_costh          , "mc_gam_isr_costh/D"        );

	tree_rec->Branch( "mc_gam_isr_boost_e"        , &m_mc_gam_isr_boost_e              , "mc_gam_isr_boost_e/D"            );
	tree_rec->Branch( "mc_gam_isr_boost_p"        , &m_mc_gam_isr_boost_p              , "mc_gam_isr_boost_p/D"            );
	tree_rec->Branch( "mc_gam_isr_boost_m"        , &m_mc_gam_isr_boost_m              , "mc_gam_isr_boost_m/D"            );
	tree_rec->Branch( "mc_gam_isr_boost_pt"       , &m_mc_gam_isr_boost_pt             , "mc_gam_isr_boost_pt/D"           );
	tree_rec->Branch( "mc_gam_isr_boost_phi"      , &m_mc_gam_isr_boost_phi            , "mc_gam_isr_boost_phi/D"          );
	tree_rec->Branch( "mc_gam_isr_boost_costh"    , &m_mc_gam_isr_boost_costh          , "mc_gam_isr_boost_costh/D"        );

	tree_rec->Branch( "mc_reweight"    , &m_mc_reweight        , "mc_reweight/I"        );
	tree_rec->Branch( "isr"            , &m_isr                , "isr/I"                );

	//---------------------------------------------------------------------------------------------------------------
	// Quantities for background MC of dipion transitions to Upsilon(nS) (no ISR, no radiative return to lower mass Upsilon -- e.g sample Y5StoY1Smumu2pic)
	// We reuse the variables for signal MC for these decays, so we only introduce one new variable (which is really just a renaming of "rho" to pip_pim since we don't have rho is these decays)

	tree_rec->Branch( "mc_pip_pim_e"        , &m_mc_pip_pim_e              , "mc_pip_pim_e/D"            );
	tree_rec->Branch( "mc_pip_pim_p"        , &m_mc_pip_pim_p              , "mc_pip_pim_p/D"            );
	tree_rec->Branch( "mc_pip_pim_m"        , &m_mc_pip_pim_m              , "mc_pip_pim_m/D"            );
	tree_rec->Branch( "mc_pip_pim_pt"       , &m_mc_pip_pim_pt             , "mc_pip_pim_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_phi"      , &m_mc_pip_pim_phi            , "mc_pip_pim_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_costh"    , &m_mc_pip_pim_costh          , "mc_pip_pim_costh/D"        );

	//---------------------------------------------------------------------------------------------------------------
	// Quantities for signal MC and generic MC

	tree_rec->Branch( "mc_Y5S_e"        , &m_mc_Y5S_e              , "mc_Y5S_e/D"            );
	tree_rec->Branch( "mc_Y5S_p"        , &m_mc_Y5S_p              , "mc_Y5S_p/D"            );
	tree_rec->Branch( "mc_Y5S_m"        , &m_mc_Y5S_m              , "mc_Y5S_m/D"            );
	tree_rec->Branch( "mc_Y5S_pt"       , &m_mc_Y5S_pt             , "mc_Y5S_pt/D"           );
	tree_rec->Branch( "mc_Y5S_phi"      , &m_mc_Y5S_phi            , "mc_Y5S_phi/D"          );
	tree_rec->Branch( "mc_Y5S_costh"    , &m_mc_Y5S_costh          , "mc_Y5S_costh/D"        );

	tree_rec->Branch( "mc_Wbj_e"      , &m_mc_Wbj_e            , "mc_Wbj_e/D"          );
	tree_rec->Branch( "mc_Wbj_p"      , &m_mc_Wbj_p            , "mc_Wbj_p/D"          );
	tree_rec->Branch( "mc_Wbj_m"      , &m_mc_Wbj_m            , "mc_Wbj_m/D"          );
	tree_rec->Branch( "mc_Wbj_pt"     , &m_mc_Wbj_pt           , "mc_Wbj_pt/D"         );
	tree_rec->Branch( "mc_Wbj_phi"    , &m_mc_Wbj_phi          , "mc_Wbj_phi/D"        );
	tree_rec->Branch( "mc_Wbj_costh"  , &m_mc_Wbj_costh        , "mc_Wbj_costh/D"      );

	tree_rec->Branch( "mc_Wbj_boost_e"      , &m_mc_Wbj_boost_e            , "mc_Wbj_boost_e/D"          );
	tree_rec->Branch( "mc_Wbj_boost_p"      , &m_mc_Wbj_boost_p            , "mc_Wbj_boost_p/D"          );
	tree_rec->Branch( "mc_Wbj_boost_m"      , &m_mc_Wbj_boost_m            , "mc_Wbj_boost_m/D"          );
	tree_rec->Branch( "mc_Wbj_boost_pt"     , &m_mc_Wbj_boost_pt           , "mc_Wbj_boost_pt/D"         );
	tree_rec->Branch( "mc_Wbj_boost_phi"    , &m_mc_Wbj_boost_phi          , "mc_Wbj_boost_phi/D"        );
	tree_rec->Branch( "mc_Wbj_boost_costh"  , &m_mc_Wbj_boost_costh        , "mc_Wbj_boost_costh/D"      );
	
	tree_rec->Branch( "mc_gam_e"       , &m_mc_gam_e             , "mc_gam_e/D"           );
	tree_rec->Branch( "mc_gam_p"       , &m_mc_gam_p             , "mc_gam_p/D"           );
	tree_rec->Branch( "mc_gam_m"       , &m_mc_gam_m             , "mc_gam_m/D"           );
	tree_rec->Branch( "mc_gam_pt"      , &m_mc_gam_pt            , "mc_gam_pt/D"          );
	tree_rec->Branch( "mc_gam_phi"     , &m_mc_gam_phi           , "mc_gam_phi/D"         );
	tree_rec->Branch( "mc_gam_costh"   , &m_mc_gam_costh         , "mc_gam_costh/D"       );
	
	tree_rec->Branch( "mc_YnS_e"    , &m_mc_YnS_e          , "mc_YnS_e/D"        );
	tree_rec->Branch( "mc_YnS_p"    , &m_mc_YnS_p          , "mc_YnS_p/D"        );
	tree_rec->Branch( "mc_YnS_m"    , &m_mc_YnS_m          , "mc_YnS_m/D"        );
	tree_rec->Branch( "mc_YnS_pt"   , &m_mc_YnS_pt         , "mc_YnS_pt/D"       );
	tree_rec->Branch( "mc_YnS_phi"  , &m_mc_YnS_phi        , "mc_YnS_phi/D"      );
	tree_rec->Branch( "mc_YnS_costh", &m_mc_YnS_costh      , "mc_YnS_costh/D"    );

	tree_rec->Branch( "mc_YnS_boost_e"    , &m_mc_YnS_boost_e          , "mc_YnS_boost_e/D"        );
	tree_rec->Branch( "mc_YnS_boost_p"    , &m_mc_YnS_boost_p          , "mc_YnS_boost_p/D"        );
	tree_rec->Branch( "mc_YnS_boost_m"    , &m_mc_YnS_boost_m          , "mc_YnS_boost_m/D"        );
	tree_rec->Branch( "mc_YnS_boost_pt"   , &m_mc_YnS_boost_pt         , "mc_YnS_boost_pt/D"       );
	tree_rec->Branch( "mc_YnS_boost_phi"  , &m_mc_YnS_boost_phi        , "mc_YnS_boost_phi/D"      );
	tree_rec->Branch( "mc_YnS_boost_costh", &m_mc_YnS_boost_costh      , "mc_YnS_boost_costh/D"    );
	
	tree_rec->Branch( "mc_rho_e"    , &m_mc_rho_e          , "mc_rho_e/D"        );
	tree_rec->Branch( "mc_rho_p"    , &m_mc_rho_p          , "mc_rho_p/D"        );
	tree_rec->Branch( "mc_rho_m"    , &m_mc_rho_m          , "mc_rho_m/D"        );
	tree_rec->Branch( "mc_rho_pt"   , &m_mc_rho_pt         , "mc_rho_pt/D"       );
	tree_rec->Branch( "mc_rho_phi"  , &m_mc_rho_phi        , "mc_rho_phi/D"      );
	tree_rec->Branch( "mc_rho_costh", &m_mc_rho_costh      , "mc_rho_costh/D"    );
	
	tree_rec->Branch( "mc_pip_e"       , &m_mc_pip_e             , "mc_pip_e/D"           );
	tree_rec->Branch( "mc_pip_p"       , &m_mc_pip_p             , "mc_pip_p/D"           );
	tree_rec->Branch( "mc_pip_m"       , &m_mc_pip_m             , "mc_pip_m/D"           );
	tree_rec->Branch( "mc_pip_pt"      , &m_mc_pip_pt            , "mc_pip_pt/D"          );
	tree_rec->Branch( "mc_pip_phi"     , &m_mc_pip_phi           , "mc_pip_phi/D"         );
	tree_rec->Branch( "mc_pip_costh"   , &m_mc_pip_costh         , "mc_pip_costh/D"       );
	
	tree_rec->Branch( "mc_pim_e"       , &m_mc_pim_e             , "mc_pim_e/D"           );
	tree_rec->Branch( "mc_pim_p"       , &m_mc_pim_p             , "mc_pim_p/D"           );
	tree_rec->Branch( "mc_pim_m"       , &m_mc_pim_m             , "mc_pim_m/D"           );
	tree_rec->Branch( "mc_pim_pt"      , &m_mc_pim_pt            , "mc_pim_pt/D"          );
	tree_rec->Branch( "mc_pim_phi"     , &m_mc_pim_phi           , "mc_pim_phi/D"         );
	tree_rec->Branch( "mc_pim_costh"   , &m_mc_pim_costh         , "mc_pim_costh/D"       );
	
	tree_rec->Branch( "mc_mup_e"        , &m_mc_mup_e              , "mc_mup_e/D"            );
	tree_rec->Branch( "mc_mup_p"        , &m_mc_mup_p              , "mc_mup_p/D"            );
	tree_rec->Branch( "mc_mup_m"        , &m_mc_mup_m              , "mc_mup_m/D"            );
	tree_rec->Branch( "mc_mup_pt"       , &m_mc_mup_pt             , "mc_mup_pt/D"           );
	tree_rec->Branch( "mc_mup_phi"      , &m_mc_mup_phi            , "mc_mup_phi/D"          );
	tree_rec->Branch( "mc_mup_costh"    , &m_mc_mup_costh          , "mc_mup_costh/D"        );
	
	tree_rec->Branch( "mc_mum_e"    , &m_mc_mum_e          , "mc_mum_e/D"        );
	tree_rec->Branch( "mc_mum_p"    , &m_mc_mum_p          , "mc_mum_p/D"        );
	tree_rec->Branch( "mc_mum_m"    , &m_mc_mum_m          , "mc_mum_m/D"        );
	tree_rec->Branch( "mc_mum_pt"   , &m_mc_mum_pt         , "mc_mum_pt/D"       );
	tree_rec->Branch( "mc_mum_phi"  , &m_mc_mum_phi        , "mc_mum_phi/D"      );
	tree_rec->Branch( "mc_mum_costh", &m_mc_mum_costh      , "mc_mum_costh/D"    );

	tree_rec->Branch( "mc_pip_pim_gam_recoil_e"        , &m_mc_pip_pim_gam_recoil_e              , "mc_pip_pim_gam_recoil_e/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_p"        , &m_mc_pip_pim_gam_recoil_p              , "mc_pip_pim_gam_recoil_p/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_m"        , &m_mc_pip_pim_gam_recoil_m              , "mc_pip_pim_gam_recoil_m/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_pt"       , &m_mc_pip_pim_gam_recoil_pt             , "mc_pip_pim_gam_recoil_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_phi"      , &m_mc_pip_pim_gam_recoil_phi            , "mc_pip_pim_gam_recoil_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_costh"    , &m_mc_pip_pim_gam_recoil_costh          , "mc_pip_pim_gam_recoil_costh/D"        );

	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_e"        , &m_mc_pip_pim_gam_recoil_boost_e              , "mc_pip_pim_gam_recoil_boost_e/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_p"        , &m_mc_pip_pim_gam_recoil_boost_p              , "mc_pip_pim_gam_recoil_boost_p/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_m"        , &m_mc_pip_pim_gam_recoil_boost_m              , "mc_pip_pim_gam_recoil_boost_m/D"            );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_pt"       , &m_mc_pip_pim_gam_recoil_boost_pt             , "mc_pip_pim_gam_recoil_boost_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_phi"      , &m_mc_pip_pim_gam_recoil_boost_phi            , "mc_pip_pim_gam_recoil_boost_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_gam_recoil_boost_costh"    , &m_mc_pip_pim_gam_recoil_boost_costh          , "mc_pip_pim_gam_recoil_boost_costh/D"        );

	tree_rec->Branch( "mc_mup_mum_gam_recoil_e"        , &m_mc_mup_mum_gam_recoil_e              , "mc_mup_mum_gam_recoil_e/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_p"        , &m_mc_mup_mum_gam_recoil_p              , "mc_mup_mum_gam_recoil_p/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_m"        , &m_mc_mup_mum_gam_recoil_m              , "mc_mup_mum_gam_recoil_m/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_pt"       , &m_mc_mup_mum_gam_recoil_pt             , "mc_mup_mum_gam_recoil_pt/D"           );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_phi"      , &m_mc_mup_mum_gam_recoil_phi            , "mc_mup_mum_gam_recoil_phi/D"          );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_costh"    , &m_mc_mup_mum_gam_recoil_costh          , "mc_mup_mum_gam_recoil_costh/D"        );

	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_e"        , &m_mc_mup_mum_gam_recoil_boost_e              , "mc_mup_mum_gam_recoil_boost_e/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_p"        , &m_mc_mup_mum_gam_recoil_boost_p              , "mc_mup_mum_gam_recoil_boost_p/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_m"        , &m_mc_mup_mum_gam_recoil_boost_m              , "mc_mup_mum_gam_recoil_boost_m/D"            );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_pt"       , &m_mc_mup_mum_gam_recoil_boost_pt             , "mc_mup_mum_gam_recoil_boost_pt/D"           );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_phi"      , &m_mc_mup_mum_gam_recoil_boost_phi            , "mc_mup_mum_gam_recoil_boost_phi/D"          );
	tree_rec->Branch( "mc_mup_mum_gam_recoil_boost_costh"    , &m_mc_mup_mum_gam_recoil_boost_costh          , "mc_mup_mum_gam_recoil_boost_costh/D"        );

	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_e"        , &m_mc_pip_pim_mup_mum_recoil_e              , "mc_pip_pim_mup_mum_recoil_e/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_p"        , &m_mc_pip_pim_mup_mum_recoil_p              , "mc_pip_pim_mup_mum_recoil_p/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_m"        , &m_mc_pip_pim_mup_mum_recoil_m              , "mc_pip_pim_mup_mum_recoil_m/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_pt"       , &m_mc_pip_pim_mup_mum_recoil_pt             , "mc_pip_pim_mup_mum_recoil_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_phi"      , &m_mc_pip_pim_mup_mum_recoil_phi            , "mc_pip_pim_mup_mum_recoil_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_costh"    , &m_mc_pip_pim_mup_mum_recoil_costh          , "mc_pip_pim_mup_mum_recoil_costh/D"        );

	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_e"        , &m_mc_pip_pim_mup_mum_recoil_boost_e              , "mc_pip_pim_mup_mum_recoil_boost_e/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_p"        , &m_mc_pip_pim_mup_mum_recoil_boost_p              , "mc_pip_pim_mup_mum_recoil_boost_p/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_m"        , &m_mc_pip_pim_mup_mum_recoil_boost_m              , "mc_pip_pim_mup_mum_recoil_boost_m/D"            );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_pt"       , &m_mc_pip_pim_mup_mum_recoil_boost_pt             , "mc_pip_pim_mup_mum_recoil_boost_pt/D"           );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_phi"      , &m_mc_pip_pim_mup_mum_recoil_boost_phi            , "mc_pip_pim_mup_mum_recoil_boost_phi/D"          );
	tree_rec->Branch( "mc_pip_pim_mup_mum_recoil_boost_costh"    , &m_mc_pip_pim_mup_mum_recoil_boost_costh          , "mc_pip_pim_mup_mum_recoil_boost_costh/D"        );

	tree_rec->Branch( "mc_gam_recoil_m"       , &m_mc_gam_recoil_m             , "mc_gam_recoil_m/D"           );

	tree_rec->Branch( "mc_gam_boost_e"       , &m_mc_gam_boost_e             , "mc_gam_boost_e/D"           );
	tree_rec->Branch( "mc_gam_boost_p"       , &m_mc_gam_boost_p             , "mc_gam_boost_p/D"           );
	tree_rec->Branch( "mc_gam_boost_m"       , &m_mc_gam_boost_m             , "mc_gam_boost_m/D"           );
	tree_rec->Branch( "mc_gam_boost_pt"      , &m_mc_gam_boost_pt            , "mc_gam_boost_pt/D"          );
	tree_rec->Branch( "mc_gam_boost_phi"     , &m_mc_gam_boost_phi           , "mc_gam_boost_phi/D"         );
	tree_rec->Branch( "mc_gam_boost_costh"   , &m_mc_gam_boost_costh         , "mc_gam_boost_costh/D"       );

	tree_rec->Branch( "mc_mup_mum_gam_m"   , &m_mc_mup_mum_gam_m         , "mc_mup_mum_gam_m/D"       );
	tree_rec->Branch( "mc_pip_pim_gam_m"   , &m_mc_pip_pim_gam_m         , "mc_pip_pim_gam_m/D"       );

	tree_rec->Branch("mc_angle_bt_pip_gam_isr"         , &m_mc_angle_bt_pip_gam_isr              , "mc_angle_bt_pip_gam_isr/D"   );
	tree_rec->Branch("mc_angle_bt_pim_gam_isr"         , &m_mc_angle_bt_pim_gam_isr              , "mc_angle_bt_pim_gam_isr/D"   );
	tree_rec->Branch("mc_angle_bt_mup_gam_isr"         , &m_mc_angle_bt_mup_gam_isr              , "mc_angle_bt_mup_gam_isr/D"   );
	tree_rec->Branch("mc_angle_bt_mum_gam_isr"         , &m_mc_angle_bt_mum_gam_isr              , "mc_angle_bt_mum_gam_isr/D"   );

	tree_rec->Branch("mc_angle_bt_pip_gam"         , &m_mc_angle_bt_pip_gam             , "mc_angle_bt_pip_gam/D"   );
	tree_rec->Branch("mc_angle_bt_pim_gam"         , &m_mc_angle_bt_pim_gam             , "mc_angle_bt_pim_gam/D"   );
	tree_rec->Branch("mc_angle_bt_mup_gam"         , &m_mc_angle_bt_mup_gam             , "mc_angle_bt_mup_gam/D"   );
	tree_rec->Branch("mc_angle_bt_mum_gam"         , &m_mc_angle_bt_mum_gam             , "mc_angle_bt_mum_gam/D"   );

	tree_rec->Branch( "mc_fsr_YnS_e"              , &m_mc_fsr_YnS_e              , "mc_fsr_YnS_e/D"                  );
	tree_rec->Branch( "mc_fsr_YnS_boost_e"        , &m_mc_fsr_YnS_boost_e        , "mc_fsr_YnS_boost_e/D"            );
	tree_rec->Branch( "mc_fsr_rho_e"              , &m_mc_fsr_rho_e              , "mc_fsr_rho_e/D"                  );
	tree_rec->Branch( "mc_fsr_rho_boost_e"        , &m_mc_fsr_rho_boost_e        , "mc_fsr_rho_boost_e/D"            );

	//---------------------------------------------------------------------------------------------------------------
	// Reconstructed quantities

	tree_rec->Branch( "Y5S_e"        , "std::vector<double>"      , &nt_Y5S_e            );
	tree_rec->Branch( "Y5S_p"        , "std::vector<double>"      , &nt_Y5S_p            );
	tree_rec->Branch( "Y5S_m"        , "std::vector<double>"      , &nt_Y5S_m            );
	tree_rec->Branch( "Y5S_pt"       , "std::vector<double>"      , &nt_Y5S_pt           );
	tree_rec->Branch( "Y5S_phi"      , "std::vector<double>"      , &nt_Y5S_phi          );
	tree_rec->Branch( "Y5S_costh"    , "std::vector<double>"      , &nt_Y5S_costh        );

	tree_rec->Branch( "Y5S_fit_e"        , "std::vector<double>"      , &nt_Y5S_fit_e            );
	tree_rec->Branch( "Y5S_fit_p"        , "std::vector<double>"      , &nt_Y5S_fit_p            );
	tree_rec->Branch( "Y5S_fit_m"        , "std::vector<double>"      , &nt_Y5S_fit_m            );
	tree_rec->Branch( "Y5S_fit_pt"       , "std::vector<double>"      , &nt_Y5S_fit_pt           );
	tree_rec->Branch( "Y5S_fit_phi"      , "std::vector<double>"      , &nt_Y5S_fit_phi          );
	tree_rec->Branch( "Y5S_fit_costh"    , "std::vector<double>"      , &nt_Y5S_fit_costh        );

	tree_rec->Branch( "Y5S_boost_e"        , "std::vector<double>"      , &nt_Y5S_boost_e            );
	tree_rec->Branch( "Y5S_boost_p"        , "std::vector<double>"      , &nt_Y5S_boost_p            );
	tree_rec->Branch( "Y5S_boost_m"        , "std::vector<double>"      , &nt_Y5S_boost_m            );
	tree_rec->Branch( "Y5S_boost_pt"       , "std::vector<double>"      , &nt_Y5S_boost_pt           );
	tree_rec->Branch( "Y5S_boost_phi"      , "std::vector<double>"      , &nt_Y5S_boost_phi          );
	tree_rec->Branch( "Y5S_boost_costh"    , "std::vector<double>"      , &nt_Y5S_boost_costh        );

	tree_rec->Branch( "Y5S_fit_boost_e"        , "std::vector<double>"      , &nt_Y5S_fit_boost_e            );
	tree_rec->Branch( "Y5S_fit_boost_p"        , "std::vector<double>"      , &nt_Y5S_fit_boost_p            );
	tree_rec->Branch( "Y5S_fit_boost_m"        , "std::vector<double>"      , &nt_Y5S_fit_boost_m            );
	tree_rec->Branch( "Y5S_fit_boost_pt"       , "std::vector<double>"      , &nt_Y5S_fit_boost_pt           );
	tree_rec->Branch( "Y5S_fit_boost_phi"      , "std::vector<double>"      , &nt_Y5S_fit_boost_phi          );
	tree_rec->Branch( "Y5S_fit_boost_costh"    , "std::vector<double>"      , &nt_Y5S_fit_boost_costh        );
	
	tree_rec->Branch( "Wbj_e"      , "std::vector<double>"        , &nt_Wbj_e          );
	tree_rec->Branch( "Wbj_p"      , "std::vector<double>"        , &nt_Wbj_p          );
	tree_rec->Branch( "Wbj_m"      , "std::vector<double>"        , &nt_Wbj_m          );
	tree_rec->Branch( "Wbj_pt"     , "std::vector<double>"        , &nt_Wbj_pt         );
	tree_rec->Branch( "Wbj_phi"    , "std::vector<double>"        , &nt_Wbj_phi        );
	tree_rec->Branch( "Wbj_costh"  , "std::vector<double>"        , &nt_Wbj_costh      );

	tree_rec->Branch( "Wbj_fit_e"      , "std::vector<double>"        , &nt_Wbj_fit_e          );
	tree_rec->Branch( "Wbj_fit_p"      , "std::vector<double>"        , &nt_Wbj_fit_p          );
	tree_rec->Branch( "Wbj_fit_m"      , "std::vector<double>"        , &nt_Wbj_fit_m          );
	tree_rec->Branch( "Wbj_fit_pt"     , "std::vector<double>"        , &nt_Wbj_fit_pt         );
	tree_rec->Branch( "Wbj_fit_phi"    , "std::vector<double>"        , &nt_Wbj_fit_phi        );
	tree_rec->Branch( "Wbj_fit_costh"  , "std::vector<double>"        , &nt_Wbj_fit_costh      );

	tree_rec->Branch( "Wbj_fit_boost_e"      , "std::vector<double>"        , &nt_Wbj_fit_boost_e          );
	tree_rec->Branch( "Wbj_fit_boost_p"      , "std::vector<double>"        , &nt_Wbj_fit_boost_p          );
	tree_rec->Branch( "Wbj_fit_boost_m"      , "std::vector<double>"        , &nt_Wbj_fit_boost_m          );
	tree_rec->Branch( "Wbj_fit_boost_pt"     , "std::vector<double>"        , &nt_Wbj_fit_boost_pt         );
	tree_rec->Branch( "Wbj_fit_boost_phi"    , "std::vector<double>"        , &nt_Wbj_fit_boost_phi        );
	tree_rec->Branch( "Wbj_fit_boost_costh"  , "std::vector<double>"        , &nt_Wbj_fit_boost_costh      );

	tree_rec->Branch( "gam_recoil_e_1"      , "std::vector<double>"        , &nt_gam_recoil_e_1          );
	tree_rec->Branch( "gam_recoil_p_1"      , "std::vector<double>"        , &nt_gam_recoil_p_1          );
	tree_rec->Branch( "gam_recoil_m_1"      , "std::vector<double>"        , &nt_gam_recoil_m_1          );
	tree_rec->Branch( "gam_recoil_pt_1"     , "std::vector<double>"        , &nt_gam_recoil_pt_1         );
	tree_rec->Branch( "gam_recoil_phi_1"    , "std::vector<double>"        , &nt_gam_recoil_phi_1        );
	tree_rec->Branch( "gam_recoil_costh_1"  , "std::vector<double>"        , &nt_gam_recoil_costh_1      );

	tree_rec->Branch( "gam_recoil_e_2"      , "std::vector<double>"        , &nt_gam_recoil_e_2          );
	tree_rec->Branch( "gam_recoil_p_2"      , "std::vector<double>"        , &nt_gam_recoil_p_2          );
	tree_rec->Branch( "gam_recoil_m_2"      , "std::vector<double>"        , &nt_gam_recoil_m_2          );
	tree_rec->Branch( "gam_recoil_pt_2"     , "std::vector<double>"        , &nt_gam_recoil_pt_2         );
	tree_rec->Branch( "gam_recoil_phi_2"    , "std::vector<double>"        , &nt_gam_recoil_phi_2        );
	tree_rec->Branch( "gam_recoil_costh_2"  , "std::vector<double>"        , &nt_gam_recoil_costh_2      );

	tree_rec->Branch( "gam_recoil_e_3"      , "std::vector<double>"        , &nt_gam_recoil_e_3          );
	tree_rec->Branch( "gam_recoil_p_3"      , "std::vector<double>"        , &nt_gam_recoil_p_3          );
	tree_rec->Branch( "gam_recoil_m_3"      , "std::vector<double>"        , &nt_gam_recoil_m_3          );
	tree_rec->Branch( "gam_recoil_pt_3"     , "std::vector<double>"        , &nt_gam_recoil_pt_3         );
	tree_rec->Branch( "gam_recoil_phi_3"    , "std::vector<double>"        , &nt_gam_recoil_phi_3        );
	tree_rec->Branch( "gam_recoil_costh_3"  , "std::vector<double>"        , &nt_gam_recoil_costh_3      );

	tree_rec->Branch( "gam_recoil_e"      , "std::vector<double>"        , &nt_gam_recoil_e          );
	tree_rec->Branch( "gam_recoil_p"      , "std::vector<double>"        , &nt_gam_recoil_p          );
	tree_rec->Branch( "gam_recoil_m"      , "std::vector<double>"        , &nt_gam_recoil_m          );
	tree_rec->Branch( "gam_recoil_pt"     , "std::vector<double>"        , &nt_gam_recoil_pt         );
	tree_rec->Branch( "gam_recoil_phi"    , "std::vector<double>"        , &nt_gam_recoil_phi        );
	tree_rec->Branch( "gam_recoil_costh"  , "std::vector<double>"        , &nt_gam_recoil_costh      );

	tree_rec->Branch( "gam_recoil_boost_m"      , "std::vector<double>"        , &nt_gam_recoil_boost_m          );

	tree_rec->Branch( "YnS_e"    , "std::vector<double>"     , &nt_YnS_e        );
	tree_rec->Branch( "YnS_p"    , "std::vector<double>"     , &nt_YnS_p        );
	tree_rec->Branch( "YnS_m"    , "std::vector<double>"     , &nt_YnS_m        );
	tree_rec->Branch( "YnS_pt"   , "std::vector<double>"     , &nt_YnS_pt       );
	tree_rec->Branch( "YnS_phi"  , "std::vector<double>"     , &nt_YnS_phi      );
	tree_rec->Branch( "YnS_costh", "std::vector<double>"     , &nt_YnS_costh    );

	tree_rec->Branch( "YnS_boost_e"    , "std::vector<double>"     , &nt_YnS_boost_e        );
	tree_rec->Branch( "YnS_boost_p"    , "std::vector<double>"     , &nt_YnS_boost_p        );
	tree_rec->Branch( "YnS_boost_m"    , "std::vector<double>"     , &nt_YnS_boost_m        );
	tree_rec->Branch( "YnS_boost_pt"   , "std::vector<double>"     , &nt_YnS_boost_pt       );
	tree_rec->Branch( "YnS_boost_phi"  , "std::vector<double>"     , &nt_YnS_boost_phi      );
	tree_rec->Branch( "YnS_boost_costh", "std::vector<double>"     , &nt_YnS_boost_costh    );

	tree_rec->Branch( "YnS_fit_e"    , "std::vector<double>"     , &nt_YnS_fit_e        );
	tree_rec->Branch( "YnS_fit_p"    , "std::vector<double>"     , &nt_YnS_fit_p        );
	tree_rec->Branch( "YnS_fit_m"    , "std::vector<double>"     , &nt_YnS_fit_m        );
	tree_rec->Branch( "YnS_fit_pt"   , "std::vector<double>"     , &nt_YnS_fit_pt       );
	tree_rec->Branch( "YnS_fit_phi"  , "std::vector<double>"     , &nt_YnS_fit_phi      );
	tree_rec->Branch( "YnS_fit_costh", "std::vector<double>"     , &nt_YnS_fit_costh    );
	
	tree_rec->Branch( "rho_e"    , "std::vector<double>"       , &nt_rho_e        );
	tree_rec->Branch( "rho_p"    , "std::vector<double>"       , &nt_rho_p        );
	tree_rec->Branch( "rho_m"    , "std::vector<double>"       , &nt_rho_m        );
	tree_rec->Branch( "rho_pt"   , "std::vector<double>"       , &nt_rho_pt       );
	tree_rec->Branch( "rho_phi"  , "std::vector<double>"       , &nt_rho_phi      );
	tree_rec->Branch( "rho_costh", "std::vector<double>"       , &nt_rho_costh    );

	tree_rec->Branch( "mup_e"      , "std::vector<double>"        , &nt_mup_e          );
	tree_rec->Branch( "mup_p"      , "std::vector<double>"        , &nt_mup_p          );
	tree_rec->Branch( "mup_m"      , "std::vector<double>"        , &nt_mup_m          );
	tree_rec->Branch( "mup_pt"     , "std::vector<double>"        , &nt_mup_pt         );
	tree_rec->Branch( "mup_phi"    , "std::vector<double>"        , &nt_mup_phi        );
	tree_rec->Branch( "mup_costh"  , "std::vector<double>"        , &nt_mup_costh      );

	tree_rec->Branch( "mum_e"      , "std::vector<double>"        , &nt_mum_e          );
	tree_rec->Branch( "mum_p"      , "std::vector<double>"        , &nt_mum_p          );
	tree_rec->Branch( "mum_m"      , "std::vector<double>"        , &nt_mum_m          );
	tree_rec->Branch( "mum_pt"     , "std::vector<double>"        , &nt_mum_pt         );
	tree_rec->Branch( "mum_phi"    , "std::vector<double>"        , &nt_mum_phi        );
	tree_rec->Branch( "mum_costh"  , "std::vector<double>"        , &nt_mum_costh      );

	tree_rec->Branch( "pip_e"      , "std::vector<double>"        , &nt_pip_e          );
	tree_rec->Branch( "pip_p"      , "std::vector<double>"        , &nt_pip_p          );
	tree_rec->Branch( "pip_m"      , "std::vector<double>"        , &nt_pip_m          );
	tree_rec->Branch( "pip_pt"     , "std::vector<double>"        , &nt_pip_pt         );
	tree_rec->Branch( "pip_phi"    , "std::vector<double>"        , &nt_pip_phi        );
	tree_rec->Branch( "pip_costh"  , "std::vector<double>"        , &nt_pip_costh      );

	tree_rec->Branch( "pim_e"      , "std::vector<double>"        , &nt_pim_e          );
	tree_rec->Branch( "pim_p"      , "std::vector<double>"        , &nt_pim_p          );
	tree_rec->Branch( "pim_m"      , "std::vector<double>"        , &nt_pim_m          );
	tree_rec->Branch( "pim_pt"     , "std::vector<double>"        , &nt_pim_pt         );
	tree_rec->Branch( "pim_phi"    , "std::vector<double>"        , &nt_pim_phi        );
	tree_rec->Branch( "pim_costh"  , "std::vector<double>"        , &nt_pim_costh      );

	tree_rec->Branch( "gam_e"      , "std::vector<double>"        , &nt_gam_e          );
	tree_rec->Branch( "gam_p"      , "std::vector<double>"        , &nt_gam_p          );
	tree_rec->Branch( "gam_m"      , "std::vector<double>"        , &nt_gam_m          );
	tree_rec->Branch( "gam_pt"     , "std::vector<double>"        , &nt_gam_pt         );
	tree_rec->Branch( "gam_phi"    , "std::vector<double>"        , &nt_gam_phi        );
	tree_rec->Branch( "gam_costh"  , "std::vector<double>"        , &nt_gam_costh      );

	tree_rec->Branch( "pip_pim_recoil_m"      , "std::vector<double>"        , &nt_pip_pim_recoil_m          );
	tree_rec->Branch( "pip_pim_recoil_e"      , "std::vector<double>"        , &nt_pip_pim_recoil_e          );

	tree_rec->Branch( "pip_pim_gam_recoil_e"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_e          );
	tree_rec->Branch( "pip_pim_gam_recoil_p"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_p          );
	tree_rec->Branch( "pip_pim_gam_recoil_m"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_m          );
	tree_rec->Branch( "pip_pim_gam_recoil_pt"     , "std::vector<double>"        , &nt_pip_pim_gam_recoil_pt         );
	tree_rec->Branch( "pip_pim_gam_recoil_phi"    , "std::vector<double>"        , &nt_pip_pim_gam_recoil_phi        );
	tree_rec->Branch( "pip_pim_gam_recoil_costh"  , "std::vector<double>"        , &nt_pip_pim_gam_recoil_costh      );

	tree_rec->Branch( "pip_pim_gam_recoil_boost_e"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_e          );
	tree_rec->Branch( "pip_pim_gam_recoil_boost_p"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_p          );
	tree_rec->Branch( "pip_pim_gam_recoil_boost_m"      , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_m          );
	tree_rec->Branch( "pip_pim_gam_recoil_boost_pt"     , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_pt         );
	tree_rec->Branch( "pip_pim_gam_recoil_boost_phi"    , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_phi        );
	tree_rec->Branch( "pip_pim_gam_recoil_boost_costh"  , "std::vector<double>"        , &nt_pip_pim_gam_recoil_boost_costh      );

	tree_rec->Branch( "mup_mum_recoil_m"      , "std::vector<double>"        , &nt_mup_mum_recoil_m          );
	tree_rec->Branch( "mup_mum_recoil_e"      , "std::vector<double>"        , &nt_mup_mum_recoil_e          );

	tree_rec->Branch( "mup_mum_gam_recoil_e"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_e          );
	tree_rec->Branch( "mup_mum_gam_recoil_p"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_p          );
	tree_rec->Branch( "mup_mum_gam_recoil_m"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_m          );
	tree_rec->Branch( "mup_mum_gam_recoil_pt"     , "std::vector<double>"        , &nt_mup_mum_gam_recoil_pt         );
	tree_rec->Branch( "mup_mum_gam_recoil_phi"    , "std::vector<double>"        , &nt_mup_mum_gam_recoil_phi        );
	tree_rec->Branch( "mup_mum_gam_recoil_costh"  , "std::vector<double>"        , &nt_mup_mum_gam_recoil_costh      );

	tree_rec->Branch( "mup_mum_gam_fit_recoil_m"      , "std::vector<double>"        , &nt_mup_mum_gam_fit_recoil_m          );

	tree_rec->Branch( "mup_mum_gam_recoil_boost_e"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_e          );
	tree_rec->Branch( "mup_mum_gam_recoil_boost_p"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_p          );
	tree_rec->Branch( "mup_mum_gam_recoil_boost_m"      , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_m          );
	tree_rec->Branch( "mup_mum_gam_recoil_boost_pt"     , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_pt         );
	tree_rec->Branch( "mup_mum_gam_recoil_boost_phi"    , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_phi        );
	tree_rec->Branch( "mup_mum_gam_recoil_boost_costh"  , "std::vector<double>"        , &nt_mup_mum_gam_recoil_boost_costh      );

	tree_rec->Branch( "pip_pim_mup_mum_recoil_e"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_e          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_p"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_p          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_m"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_m          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_pt"     , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_pt         );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_phi"    , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_phi        );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_costh"  , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_costh      );

	tree_rec->Branch( "pip_pim_mup_mum_recoil_m2"     , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_m2         );
	tree_rec->Branch( "pip_pim_mup_mum_fit_recoil_m"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_fit_recoil_m          );
	tree_rec->Branch( "pip_pim_mup_mum_fit_recoil_m2"     , "std::vector<double>"        , &nt_pip_pim_mup_mum_fit_recoil_m2         );
	tree_rec->Branch( "pip_pim_mup_mum_fit_recoil_e"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_fit_recoil_e          );
	tree_rec->Branch( "pip_pim_mup_mum_fit_recoil_costh"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_fit_recoil_costh          );

	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_e"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_e          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_p"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_p          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_m"      , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_m          );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_pt"     , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_pt         );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_phi"    , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_phi        );
	tree_rec->Branch( "pip_pim_mup_mum_recoil_boost_costh"  , "std::vector<double>"        , &nt_pip_pim_mup_mum_recoil_boost_costh      );

	tree_rec->Branch( "gam_boost_e"      , "std::vector<double>"        , &nt_gam_boost_e          );
	tree_rec->Branch( "gam_boost_p"      , "std::vector<double>"        , &nt_gam_boost_p          );
	tree_rec->Branch( "gam_boost_m"      , "std::vector<double>"        , &nt_gam_boost_m          );
	tree_rec->Branch( "gam_boost_pt"     , "std::vector<double>"        , &nt_gam_boost_pt         );
	tree_rec->Branch( "gam_boost_phi"    , "std::vector<double>"        , &nt_gam_boost_phi        );
	tree_rec->Branch( "gam_boost_costh"  , "std::vector<double>"        , &nt_gam_boost_costh      );

	tree_rec->Branch( "mup_mum_gam_m"  , "std::vector<double>"   , &nt_mup_mum_gam_m );
	tree_rec->Branch( "mup_mum_gam_e"  , "std::vector<double>"   , &nt_mup_mum_gam_e );
	tree_rec->Branch( "pip_pim_gam_m"  , "std::vector<double>"   , &nt_pip_pim_gam_m );
	tree_rec->Branch( "pip_pim_gam_e"  , "std::vector<double>"   , &nt_pip_pim_gam_e );
	tree_rec->Branch( "deltaE"         , "std::vector<double>"   , &nt_deltaE        );
	tree_rec->Branch( "deltaE_Y1D"         , "std::vector<double>"   , &nt_deltaE_Y1D        );
	tree_rec->Branch( "chi2_YnS"       , "std::vector<double>"   , &nt_chi2_YnS      );
    
	tree_rec->Branch( "neutral_e"       , "std::vector<double>"   , &nt_neutral_e      );
	tree_rec->Branch( "charged_e"       , "std::vector<double>"   , &nt_charged_e      );

	tree_rec->Branch( "Ntrk_good",       &m_Ntrk_good      , "Ntrk_good/I"    );

	tree_rec->Branch( "angle_bt_Wbj_gam"                 , "std::vector<double>"      , &nt_angle_bt_Wbj_gam            );

	tree_rec->Branch( "gam_boost_adjust_e"               , "std::vector<double>"      , &nt_gam_boost_adjust_e            );
	tree_rec->Branch( "gam_boost_adjust_recoil_m"        , "std::vector<double>"      , &nt_gam_boost_adjust_recoil_m            );

	tree_rec->Branch( "gam_nhit"    , "std::vector<double>"   , &nt_gam_nhit );
	tree_rec->Branch( "gam_width"   , "std::vector<double>"   , &nt_gam_width );
	tree_rec->Branch( "gam_e9oe25"  , "std::vector<double>"   , &nt_gam_e9oe25 );

	tree_rec->Branch( "resonance_fit"  , &m_resonance_fit       , "m_resonance_fit/I" );

	//---------------------------------------------------------------------------------------------------------------
	// Parent MC truth information for reconstructed quantities

	tree_rec->Branch( "gam_id"         , "std::vector<int>"  , &nt_gam_id        );
	tree_rec->Branch( "gam_parent_id"  , "std::vector<int>"  , &nt_gam_parent_id );

	tree_rec->Branch( "pip_parent_id"  , &m_pip_parent_id       , "pip_parent_id/I" );
	tree_rec->Branch( "pip_recoil_id"  , &m_pip_recoil_id       , "pip_recoil_id/I" );
	tree_rec->Branch( "pim_parent_id"  , &m_pim_parent_id       , "pim_parent_id/I" );
	tree_rec->Branch( "pim_recoil_id"  , &m_pim_recoil_id       , "pim_recoil_id/I" );
	tree_rec->Branch( "mup_parent_id"  , &m_mup_parent_id       , "mup_parent_id/I" );
	tree_rec->Branch( "mum_parent_id"  , &m_mum_parent_id       , "mum_parent_id/I" );

	//---------------------------------------------------------------------------------------------------------------
	// Misc. quantities, primarily used for investigation of background sources

	tree_rec->Branch( "bkg_A", &m_bkg_A, "bkg_A/I" );
	tree_rec->Branch( "bkg_B", &m_bkg_B, "bkg_B/I" );
	tree_rec->Branch( "bkg_C", &m_bkg_C, "bkg_C/I" );
	tree_rec->Branch( "bkg_D", &m_bkg_D, "bkg_D/I" );
	tree_rec->Branch( "bkg_E", &m_bkg_E, "bkg_E/I" );
	tree_rec->Branch( "bkg_F", &m_bkg_F, "bkg_F/I" );
	tree_rec->Branch( "bkg_G", &m_bkg_G, "bkg_G/I" );
	tree_rec->Branch( "bkg_H", &m_bkg_H, "bkg_H/I" );
	tree_rec->Branch( "bkg_I", &m_bkg_I, "bkg_I/I" );
	tree_rec->Branch( "bkg_J", &m_bkg_J, "bkg_J/I" );
	tree_rec->Branch( "bkg_K", &m_bkg_K, "bkg_K/I" );

	tree_rec->Branch( "Y3S"              , &m_Y3S              , "Y3S/I" );
	tree_rec->Branch( "Y3S_to_Y1S"       , &m_Y3S_to_Y1S       , "Y3S_to_Y1S/I" );
	tree_rec->Branch( "Y3S_to_Y2S"       , &m_Y3S_to_Y2S       , "Y3S_to_Y2S/I" );
	tree_rec->Branch( "Y3S_to_Y2S_to_Y1S", &m_Y3S_to_Y2S_to_Y1S, "Y3S_to_Y2S_to_Y1S/I" );

	tree_rec->Branch( "Y2S"       , &m_Y2S, "Y2S/I" );
	tree_rec->Branch( "Y2S_to_Y1S", &m_Y2S_to_Y1S, "Y2S_to_Y1S/I" );

	tree_rec->Branch( "Y1S", &m_Y1S, "Y1S/I" );
	tree_rec->Branch( "Y3S_to_X2P", &m_Y3S_to_X2P, "Y3S_to_X2P/I" );
	/*
	tree_rec->Branch( "Y1S_pi0pi0", &m_Y1S_pi0pi0, "Y1S_pi0pi0/I" );
	tree_rec->Branch( "Y2S_pi0pi0", &m_Y2S_pi0pi0, "Y2S_pi0pi0/I" );
	tree_rec->Branch( "Y3S_pi0pi0", &m_Y3S_pi0pi0, "Y3S_pi0pi0/I" );

	tree_rec->Branch( "Y1S_pippim", &m_Y1S_pippim, "Y1S_pippim/I" );
	tree_rec->Branch( "Y2S_pippim", &m_Y2S_pippim, "Y2S_pippim/I" );
	tree_rec->Branch( "Y3S_pippim", &m_Y3S_pippim, "Y3S_pippim/I" );
	*/
	tree_rec->Branch( "bkg_A1", &m_bkg_A1, "bkg_A1/I" );
	tree_rec->Branch( "bkg_B1", &m_bkg_B1, "bkg_B1/I" );
	tree_rec->Branch( "bkg_B2", &m_bkg_B2, "bkg_B2/I" );
	tree_rec->Branch( "bkg_C1", &m_bkg_C1, "bkg_C1/I" );
	tree_rec->Branch( "bkg_C2", &m_bkg_C2, "bkg_C2/I" );
	tree_rec->Branch( "bkg_C3", &m_bkg_C3, "bkg_C3/I" );
	tree_rec->Branch( "bkg_D1", &m_bkg_D1, "bkg_D1/I" );
	tree_rec->Branch( "bkg_E1", &m_bkg_E1, "bkg_E1/I" );
	tree_rec->Branch( "bkg_F1", &m_bkg_F1, "bkg_F1/I" );
	tree_rec->Branch( "bkg_G1", &m_bkg_G1, "bkg_G1/I" );
	tree_rec->Branch( "bkg_G2", &m_bkg_G2, "bkg_G2/I" );
	tree_rec->Branch( "bkg_G3", &m_bkg_G3, "bkg_G3/I" );
	tree_rec->Branch( "bkg_H1", &m_bkg_H1, "bkg_H1/I" );
	tree_rec->Branch( "bkg_H2", &m_bkg_H2, "bkg_H2/I" );
	tree_rec->Branch( "bkg_I1", &m_bkg_I1, "bkg_I1/I" );
	tree_rec->Branch( "bkg_J1", &m_bkg_J1, "bkg_J1/I" );
	tree_rec->Branch( "bkg_K1", &m_bkg_K1, "bkg_K1/I" );

	tree_rec->Branch("Y5S_ndaughters", &m_Y5S_ndaughters, "Y5S_ndaughters/I");
	
	//---------------------------------------------------------------------------------------------------------------

      }
      
      void term(void) {
	m_outputRootFile->Write();
	m_outputRootFile->Close();
	cout << "-I-term(): " << m_event_seq << " events processed (total)" << endl;
      };
      
      void disp_stat(const char*) {};
      void hist_def(void);
      void event(BelleEvent*, int*);
      void begin_run(BelleEvent*, int*);
      void end_run(BelleEvent*, int*) {
	cout << "-I-end_run(): " << eventcount << " events processed in this run" << endl;
      };
      
      void other(int*, BelleEvent*, int*) {};
      
      // arrays for the names of input and output files
      char inpath[4096];
      char outpath[4096];
      
      // member functions  
      
  private:
      

      //---------------------------------------------------------------------------------------------------------------
      // pointers to ROOT objects

      TFile* m_outputRootFile;

      TTree* tree_rec;
      TTree* tree_mc;
      TTree* tree_bkg;

      //---------------------------------------------------------------------------------------------------------------
      // declaring ntuples for reconstructed candidates

      std::vector<double> nt_angle_bt_pip_gam;
      std::vector<double> nt_angle_bt_pim_gam;
      std::vector<double> nt_angle_bt_mup_gam;
      std::vector<double> nt_angle_bt_mum_gam;

      std::set<int> set_mdst_chrg_id ; //unique mdst chrg ids of high quality tracks

      std::vector<int> nt_gam_nhit;
      std::vector<double> nt_gam_width;
      std::vector<double> nt_gam_e9oe25;

      std::vector<int>    nt_gam_from_pi0       ;
      std::vector<int>    nt_index       ;

      std::vector<int>    nt_gam_mct  ;
      std::vector<int>    nt_mup_mct  ;
      std::vector<int>    nt_mum_mct  ;
      std::vector<int>    nt_pip_mct     ;
      std::vector<int>    nt_pim_mct     ;
      std::vector<int>    nt_mct         ;

      std::vector<int>    nt_sr         ;

      std::vector<int> nt_gam_id;
      std::vector<int> nt_gam_parent_id;

      std::vector<double> nt_angle_bt_Wbj_gam;

      std::vector<double> nt_gam_boost_adjust_e;
      std::vector<double> nt_gam_boost_adjust_recoil_m;

      std::vector<double> nt_Y5S_e       ;
      std::vector<double> nt_Y5S_p       ;
      std::vector<double> nt_Y5S_m       ;
      std::vector<double> nt_Y5S_pt  ;
      std::vector<double> nt_Y5S_phi  ;
      std::vector<double> nt_Y5S_costh  ;

      std::vector<double> nt_Y5S_fit_e       ;
      std::vector<double> nt_Y5S_fit_p       ;
      std::vector<double> nt_Y5S_fit_m       ;
      std::vector<double> nt_Y5S_fit_pt  ;
      std::vector<double> nt_Y5S_fit_phi  ;
      std::vector<double> nt_Y5S_fit_costh  ;

      std::vector<double> nt_Y5S_boost_e       ;
      std::vector<double> nt_Y5S_boost_p       ;
      std::vector<double> nt_Y5S_boost_m       ;
      std::vector<double> nt_Y5S_boost_pt  ;
      std::vector<double> nt_Y5S_boost_phi  ;
      std::vector<double> nt_Y5S_boost_costh  ;

      std::vector<double> nt_Y5S_fit_boost_e       ;
      std::vector<double> nt_Y5S_fit_boost_p       ;
      std::vector<double> nt_Y5S_fit_boost_m       ;
      std::vector<double> nt_Y5S_fit_boost_pt  ;
      std::vector<double> nt_Y5S_fit_boost_phi  ;
      std::vector<double> nt_Y5S_fit_boost_costh  ;

      std::vector<double> nt_Wbj_e       ;
      std::vector<double> nt_Wbj_p       ;
      std::vector<double> nt_Wbj_m       ;
      std::vector<double> nt_Wbj_pt  ;
      std::vector<double> nt_Wbj_phi  ;
      std::vector<double> nt_Wbj_costh  ;

      std::vector<double> nt_Wbj_fit_e       ;
      std::vector<double> nt_Wbj_fit_p       ;
      std::vector<double> nt_Wbj_fit_m       ;
      std::vector<double> nt_Wbj_fit_pt  ;
      std::vector<double> nt_Wbj_fit_phi  ;
      std::vector<double> nt_Wbj_fit_costh  ;

      std::vector<double> nt_Wbj_fit_boost_e       ;
      std::vector<double> nt_Wbj_fit_boost_p       ;
      std::vector<double> nt_Wbj_fit_boost_m       ;
      std::vector<double> nt_Wbj_fit_boost_pt  ;
      std::vector<double> nt_Wbj_fit_boost_phi  ;
      std::vector<double> nt_Wbj_fit_boost_costh  ;

      std::vector<double> nt_gam_recoil_e_1       ;
      std::vector<double> nt_gam_recoil_p_1       ;
      std::vector<double> nt_gam_recoil_m_1       ;
      std::vector<double> nt_gam_recoil_pt_1      ;
      std::vector<double> nt_gam_recoil_phi_1     ;
      std::vector<double> nt_gam_recoil_costh_1   ;

      std::vector<double> nt_gam_recoil_e_2       ;
      std::vector<double> nt_gam_recoil_p_2       ;
      std::vector<double> nt_gam_recoil_m_2       ;
      std::vector<double> nt_gam_recoil_pt_2      ;
      std::vector<double> nt_gam_recoil_phi_2     ;
      std::vector<double> nt_gam_recoil_costh_2   ;

      std::vector<double> nt_gam_recoil_e_3       ;
      std::vector<double> nt_gam_recoil_p_3       ;
      std::vector<double> nt_gam_recoil_m_3       ;
      std::vector<double> nt_gam_recoil_pt_3      ;
      std::vector<double> nt_gam_recoil_phi_3     ;
      std::vector<double> nt_gam_recoil_costh_3   ;

      std::vector<double> nt_gam_recoil_e       ;
      std::vector<double> nt_gam_recoil_p       ;
      std::vector<double> nt_gam_recoil_m       ;
      std::vector<double> nt_gam_recoil_pt      ;
      std::vector<double> nt_gam_recoil_phi     ;
      std::vector<double> nt_gam_recoil_costh   ;

      std::vector<double> nt_gam_recoil_boost_m       ;

      std::vector<double> nt_YnS_e  ;
      std::vector<double> nt_YnS_p       ;
      std::vector<double> nt_YnS_m       ;
      std::vector<double> nt_YnS_pt  ;
      std::vector<double> nt_YnS_phi  ;
      std::vector<double> nt_YnS_costh  ;

      std::vector<double> nt_YnS_boost_e  ;
      std::vector<double> nt_YnS_boost_p       ;
      std::vector<double> nt_YnS_boost_m       ;
      std::vector<double> nt_YnS_boost_pt  ;
      std::vector<double> nt_YnS_boost_phi  ;
      std::vector<double> nt_YnS_boost_costh  ;

      std::vector<double> nt_YnS_fit_e  ;
      std::vector<double> nt_YnS_fit_p       ;
      std::vector<double> nt_YnS_fit_m       ;
      std::vector<double> nt_YnS_fit_pt  ;
      std::vector<double> nt_YnS_fit_phi  ;
      std::vector<double> nt_YnS_fit_costh  ;

      std::vector<double> nt_rho_e  ;
      std::vector<double> nt_rho_p       ;
      std::vector<double> nt_rho_m       ;
      std::vector<double> nt_rho_pt  ;
      std::vector<double> nt_rho_phi   ;
      std::vector<double> nt_rho_costh  ;
      
      std::vector<double> nt_mup_e  ;
      std::vector<double> nt_mup_p       ;
      std::vector<double> nt_mup_m       ;
      std::vector<double> nt_mup_pt  ;
      std::vector<double> nt_mup_phi   ;
      std::vector<double> nt_mup_costh  ;

      std::vector<double> nt_mum_e  ;
      std::vector<double> nt_mum_p       ;
      std::vector<double> nt_mum_m       ;
      std::vector<double> nt_mum_pt  ;
      std::vector<double> nt_mum_phi   ;
      std::vector<double> nt_mum_costh  ;

      std::vector<double> nt_pip_e  ;
      std::vector<double> nt_pip_p       ;
      std::vector<double> nt_pip_m       ;
      std::vector<double> nt_pip_pt  ;
      std::vector<double> nt_pip_phi   ;
      std::vector<double> nt_pip_costh  ;

      std::vector<double> nt_pim_e  ;
      std::vector<double> nt_pim_p       ;
      std::vector<double> nt_pim_m       ;
      std::vector<double> nt_pim_pt  ;
      std::vector<double> nt_pim_phi   ;
      std::vector<double> nt_pim_costh  ;

      std::vector<double> nt_gam_e  ;
      std::vector<double> nt_gam_p       ;
      std::vector<double> nt_gam_m       ;
      std::vector<double> nt_gam_pt  ;
      std::vector<double> nt_gam_phi   ;
      std::vector<double> nt_gam_costh  ;

      std::vector<double> nt_pip_pim_recoil_m       ;
      std::vector<double> nt_pip_pim_recoil_e       ;

      std::vector<double> nt_pip_pim_gam_recoil_e  ;
      std::vector<double> nt_pip_pim_gam_recoil_p       ;
      std::vector<double> nt_pip_pim_gam_recoil_m       ;
      std::vector<double> nt_pip_pim_gam_recoil_pt  ;
      std::vector<double> nt_pip_pim_gam_recoil_phi   ;
      std::vector<double> nt_pip_pim_gam_recoil_costh  ;

      std::vector<double> nt_pip_pim_gam_recoil_boost_e  ;
      std::vector<double> nt_pip_pim_gam_recoil_boost_p       ;
      std::vector<double> nt_pip_pim_gam_recoil_boost_m       ;
      std::vector<double> nt_pip_pim_gam_recoil_boost_pt  ;
      std::vector<double> nt_pip_pim_gam_recoil_boost_phi   ;
      std::vector<double> nt_pip_pim_gam_recoil_boost_costh  ;

      std::vector<double> nt_mup_mum_recoil_m       ;
      std::vector<double> nt_mup_mum_recoil_e       ;

      std::vector<double> nt_mup_mum_gam_recoil_e  ;
      std::vector<double> nt_mup_mum_gam_recoil_p       ;
      std::vector<double> nt_mup_mum_gam_recoil_m       ;
      std::vector<double> nt_mup_mum_gam_recoil_pt  ;
      std::vector<double> nt_mup_mum_gam_recoil_phi   ;
      std::vector<double> nt_mup_mum_gam_recoil_costh  ;

      std::vector<double> nt_mup_mum_gam_fit_recoil_m       ;

      std::vector<double> nt_mup_mum_gam_recoil_boost_e  ;
      std::vector<double> nt_mup_mum_gam_recoil_boost_p       ;
      std::vector<double> nt_mup_mum_gam_recoil_boost_m       ;
      std::vector<double> nt_mup_mum_gam_recoil_boost_pt  ;
      std::vector<double> nt_mup_mum_gam_recoil_boost_phi   ;
      std::vector<double> nt_mup_mum_gam_recoil_boost_costh  ;

      std::vector<double> nt_pip_pim_mup_mum_recoil_e  ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_p       ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_m       ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_pt  ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_phi   ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_costh  ;

      std::vector<double> nt_pip_pim_mup_mum_recoil_m2       ;
      std::vector<double> nt_pip_pim_mup_mum_fit_recoil_m       ;
      std::vector<double> nt_pip_pim_mup_mum_fit_recoil_m2      ;
      std::vector<double> nt_pip_pim_mup_mum_fit_recoil_e       ;
      std::vector<double> nt_pip_pim_mup_mum_fit_recoil_costh   ;

      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_e  ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_p       ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_m       ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_pt  ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_phi   ;
      std::vector<double> nt_pip_pim_mup_mum_recoil_boost_costh  ;

      std::vector<double> nt_gam_boost_e  ;
      std::vector<double> nt_gam_boost_p       ;
      std::vector<double> nt_gam_boost_m       ;
      std::vector<double> nt_gam_boost_pt  ;
      std::vector<double> nt_gam_boost_phi   ;
      std::vector<double> nt_gam_boost_costh  ;
  
      std::vector<double> nt_mup_mum_gam_m     ;       
      std::vector<double> nt_mup_mum_gam_e     ;    
      std::vector<double> nt_pip_pim_gam_m     ;
      std::vector<double> nt_pip_pim_gam_e     ;
      std::vector<double> nt_deltaE     ;
      std::vector<double> nt_deltaE_Y1D     ;

      std::vector<double> nt_chi2_YnS     ;

      std::vector<double> nt_neutral_e     ;
      std::vector<double> nt_charged_e     ;

      std::set<int> nt_mc_fsr_YnS_ids     ;
      std::set<int> nt_mc_fsr_rho_ids     ;

      // We will use vectors of particles for preselected candidates according to some criteria
      std::vector<Particle> Y5S_list;
      std::vector<Particle> Wbj_list;
      std::vector<Particle> YnS_list;
      std::vector<Particle> rho_list;
      std::vector<Particle> pip_list;
      std::vector<Particle> pim_list;
      std::vector<Particle> pic_list;
      std::vector<Particle> mup_list;
      std::vector<Particle> mum_list;
      std::vector<Particle> muc_list;
      std::vector<Particle> gam_list;

      std::vector<Particle> leftover_charged_list;
      std::vector<Particle> leftover_neutral_list;
      std::vector<Particle> signal_part_list;

      std::vector< int > photon_parent_id_bkg_01;
      std::vector< int > photon_parent_id_bkg_02;
      std::vector< int > photon_parent_id_bkg_03;
      std::vector< int > photon_parent_id_bkg_04;
      std::vector< int > photon_parent_id_bkg_05;
      std::vector< int > photon_parent_id_bkg_06;
      std::vector< int > photon_parent_id_bkg_07;
      std::vector< int > photon_parent_id_bkg_08;
      std::vector< int > photon_parent_id_bkg_09;

      std::vector< int > pip_parent_id_bkg_01;
      std::vector< int > pip_parent_id_bkg_02;
      std::vector< int > pip_parent_id_bkg_03;
      std::vector< int > pip_parent_id_bkg_04;
      std::vector< int > pip_parent_id_bkg_05;
      std::vector< int > pip_parent_id_bkg_06;
      std::vector< int > pip_parent_id_bkg_07;
      std::vector< int > pip_parent_id_bkg_08;
      std::vector< int > pip_parent_id_bkg_09;

      std::vector< int > pip_recoil_id_bkg_01;
      std::vector< int > pip_recoil_id_bkg_02;
      std::vector< int > pip_recoil_id_bkg_03;
      std::vector< int > pip_recoil_id_bkg_04;
      std::vector< int > pip_recoil_id_bkg_05;
      std::vector< int > pip_recoil_id_bkg_06;
      std::vector< int > pip_recoil_id_bkg_07;
      std::vector< int > pip_recoil_id_bkg_08;
      std::vector< int > pip_recoil_id_bkg_09;

      std::vector< int > mup_parent_id_bkg_01;
      std::vector< int > mup_parent_id_bkg_02;
      std::vector< int > mup_parent_id_bkg_03;
      std::vector< int > mup_parent_id_bkg_04;
      std::vector< int > mup_parent_id_bkg_05;
      std::vector< int > mup_parent_id_bkg_06;
      std::vector< int > mup_parent_id_bkg_07;
      std::vector< int > mup_parent_id_bkg_08;
      std::vector< int > mup_parent_id_bkg_09;

      std::vector<int> temp;
      std::vector<int> photon_parent_id_list ;

      std::vector<candidateRec> candidatesRec; 

      std::vector<candidateRec> candidatesBkg_01;
      std::vector<candidateRec> candidatesBkg_02;
      std::vector<candidateRec> candidatesBkg_03;
      std::vector<candidateRec> candidatesBkg_04;
      std::vector<candidateRec> candidatesBkg_05;
      std::vector<candidateRec> candidatesBkg_06;
      std::vector<candidateRec> candidatesBkg_07;
      std::vector<candidateRec> candidatesBkg_08;
      std::vector<candidateRec> candidatesBkg_09;
      
      //---------------------------------------------------------------------------------------------------------------
      // reconstructed candidate information

      struct candidateRec{

	double angle_bt_pip_gam;
	double angle_bt_pim_gam;
	double angle_bt_mup_gam;
	double angle_bt_mum_gam;

	int gam_nhit;
	double gam_width;
	double gam_e9oe25;

	double angle_bt_Wbj_gam;
	double gam_boost_adjust_e;
	double gam_boost_adjust_recoil_m;

	double Y5S_e     ;
	double Y5S_p     ;
	double Y5S_m     ;
	double Y5S_pt    ;
	double Y5S_phi   ;
	double Y5S_costh ;

	double Y5S_fit_e     ;
	double Y5S_fit_p     ;
	double Y5S_fit_m     ;
	double Y5S_fit_pt    ;
	double Y5S_fit_phi   ;
	double Y5S_fit_costh ;

	double Y5S_boost_e     ;
	double Y5S_boost_p     ;
	double Y5S_boost_m     ;
	double Y5S_boost_pt    ;
	double Y5S_boost_phi   ;
	double Y5S_boost_costh ;

	double Y5S_fit_boost_e     ;
	double Y5S_fit_boost_p     ;
	double Y5S_fit_boost_m     ;
	double Y5S_fit_boost_pt    ;
	double Y5S_fit_boost_phi   ;
	double Y5S_fit_boost_costh ;
	
	double Wbj_e     ;
	double Wbj_p     ;
	double Wbj_m     ;
	double Wbj_pt    ;
	double Wbj_phi   ;
	double Wbj_costh ;

	double Wbj_fit_e     ;
	double Wbj_fit_p ;
	double Wbj_fit_m ;
	double Wbj_fit_pt    ;
	double Wbj_fit_phi   ;
	double Wbj_fit_costh ;

	double Wbj_fit_boost_e     ;
	double Wbj_fit_boost_p ;
	double Wbj_fit_boost_m ;
	double Wbj_fit_boost_pt    ;
	double Wbj_fit_boost_phi   ;
	double Wbj_fit_boost_costh ;
	
	double gam_recoil_e_1     ;
	double gam_recoil_p_1     ;
	double gam_recoil_m_1     ;
	double gam_recoil_pt_1    ;
	double gam_recoil_phi_1   ;
	double gam_recoil_costh_1 ;

	double gam_recoil_e_2     ;
	double gam_recoil_p_2     ;
	double gam_recoil_m_2     ;
	double gam_recoil_pt_2    ;
	double gam_recoil_phi_2   ;
	double gam_recoil_costh_2 ;

	double gam_recoil_e_3     ;
	double gam_recoil_p_3     ;
	double gam_recoil_m_3     ;
	double gam_recoil_pt_3    ;
	double gam_recoil_phi_3   ;
	double gam_recoil_costh_3 ;

	double gam_recoil_e     ;
	double gam_recoil_p     ;
	double gam_recoil_m     ;
	double gam_recoil_pt    ;
	double gam_recoil_phi   ;
	double gam_recoil_costh ;

	double gam_recoil_boost_m;
      
	double gam_e     ;
	double gam_p     ;
	double gam_m     ;
	double gam_pt    ;
	double gam_phi   ;
	double gam_costh ;
	
	double YnS_e     ;
	double YnS_p     ;
	double YnS_m     ;
	double YnS_pt    ;
	double YnS_phi   ;
	double YnS_costh ;

	double YnS_boost_e     ;
	double YnS_boost_p     ;
	double YnS_boost_m     ;
	double YnS_boost_pt    ;
	double YnS_boost_phi   ;
	double YnS_boost_costh ;
	
	double YnS_fit_e     ;
	double YnS_fit_p     ;
	double YnS_fit_m     ;
	double YnS_fit_pt    ;
	double YnS_fit_phi   ;
	double YnS_fit_costh ;

	double rho_e     ;
	double rho_p     ;
	double rho_m     ;
	double rho_pt    ;
	double rho_phi   ;
	double rho_costh ;
	
	double pip_e     ;
	double pip_p     ;
	double pip_m     ;
	double pip_pt    ;
	double pip_phi   ;
	double pip_costh ;
	
	double pim_e     ;
	double pim_p     ;
	double pim_m     ;
	double pim_pt    ;
	double pim_phi   ;
	double pim_costh ;
      
	double mup_e     ;
	double mup_p     ;
	double mup_m     ;
	double mup_pt    ;
	double mup_phi   ;
	double mup_costh ;
      
	double mum_e     ;
	double mum_p     ;
	double mum_m     ;
	double mum_pt    ;
	double mum_phi   ;
	double mum_costh ;

	double pip_pim_recoil_m     ;
	double pip_pim_recoil_e     ;

	double pip_pim_gam_recoil_e     ;
	double pip_pim_gam_recoil_p     ;
	double pip_pim_gam_recoil_m     ;
	double pip_pim_gam_recoil_pt    ;
	double pip_pim_gam_recoil_phi   ;
	double pip_pim_gam_recoil_costh ;

	double pip_pim_gam_recoil_boost_e     ;
	double pip_pim_gam_recoil_boost_p     ;
	double pip_pim_gam_recoil_boost_m     ;
	double pip_pim_gam_recoil_boost_pt    ;
	double pip_pim_gam_recoil_boost_phi   ;
	double pip_pim_gam_recoil_boost_costh ;

	double mup_mum_recoil_m     ;
	double mup_mum_recoil_e     ;

	double mup_mum_gam_recoil_e     ;
	double mup_mum_gam_recoil_p     ;
	double mup_mum_gam_recoil_m     ;
	double mup_mum_gam_recoil_pt    ;
	double mup_mum_gam_recoil_phi   ;
	double mup_mum_gam_recoil_costh ;

	double mup_mum_gam_fit_recoil_m     ;
	double mup_mum_gam_fit_recoil_m2     ;

	double mup_mum_gam_recoil_boost_e     ;
	double mup_mum_gam_recoil_boost_p     ;
	double mup_mum_gam_recoil_boost_m     ;
	double mup_mum_gam_recoil_boost_pt    ;
	double mup_mum_gam_recoil_boost_phi   ;
	double mup_mum_gam_recoil_boost_costh ;

	double pip_pim_mup_mum_recoil_e     ;
	double pip_pim_mup_mum_recoil_p     ;
	double pip_pim_mup_mum_recoil_m     ;
	double pip_pim_mup_mum_recoil_pt    ;
	double pip_pim_mup_mum_recoil_phi   ;
	double pip_pim_mup_mum_recoil_costh ;

	double pip_pim_mup_mum_recoil_m2    ;
	double pip_pim_mup_mum_fit_recoil_m     ;
	double pip_pim_mup_mum_fit_recoil_m2    ;
	double pip_pim_mup_mum_fit_recoil_e     ;
	double pip_pim_mup_mum_fit_recoil_costh ;

	double pip_pim_mup_mum_recoil_boost_e     ;
	double pip_pim_mup_mum_recoil_boost_p     ;
	double pip_pim_mup_mum_recoil_boost_m     ;
	double pip_pim_mup_mum_recoil_boost_pt    ;
	double pip_pim_mup_mum_recoil_boost_phi   ;
	double pip_pim_mup_mum_recoil_boost_costh ;

	double gam_boost_e     ;
	double gam_boost_p     ;
	double gam_boost_m     ;
	double gam_boost_pt    ;
	double gam_boost_phi   ;
	double gam_boost_costh ;

	double mup_mum_gam_m ;
	double mup_mum_gam_e ;
	double pip_pim_gam_m ;
	double pip_pim_gam_e ;
	double deltaE;
	double deltaE_Y1D;
	double chi2_YnS;

	int gam_parent_id;

	int mct;

	int gam_mct;
	int mup_mct;
	int mum_mct;
	int pip_mct;
	int pim_mct;

	int gam_id;
	int pip_id;
	int pim_id;
	int mup_id;
	int mum_id;

	int gam_from_pi0;

	double neutral_e;
	double charged_e;

	int in_signal_region;

	
      };


      //---------------------------------------------------------------------------------------------------------------
      // Declaring MC simulation info and MC truth quantities

      HepLorentzVector m_p4_mc_Y5S;
      HepLorentzVector m_p4_mc_gam;
      
      int m_exp;  
      int m_run;  
      int m_run_seq;  
      int m_event;
      int m_event_seq;

      int m_event_run;
      int m_event_seq_skim;

      int m_Upsilon_n;

      int m_ncand;
      int m_ncand_sr;
      int m_ncand_mct;
      int m_ncand_gam;
      int m_ncand_pi0;
      
      int    m_evt_Ntrk;
      int    m_evt_Ncls;
      double m_evt_Psum;
      double m_evt_Esum;
      double m_evt_Evis;
      double m_evt_Pz;
      double m_evt_HeavyJetMass;
      double m_evt_Thrust;
      double m_evt_R2;
      
      int m_HadronB_flag;
      int m_HadronA_flag;
      
      int m_Tau_flag;

      int m_mc;
      int m_mcs;
      
      double m_vr2;
      double m_cosbt;
      double m_vthrust;

      int count_Y5S;
      int count_Wbj;
      int count_gam;
      int count_YnS;
      int count_rho;
      int count_mup;
      int count_mum;
      int count_pip;
      int count_pim;

      int m_count_gam_rec;
      int m_count_mup_rec;
      int m_count_mum_rec;
      int m_count_pip_rec;
      int m_count_pim_rec;

      int m_count_gam_good_trk;
      int m_count_mup_good_trk;
      int m_count_mum_good_trk;
      int m_count_pip_good_trk;
      int m_count_pim_good_trk;

      int m_count_gam_good_pt;
      int m_count_mup_good_pt;
      int m_count_mum_good_pt;
      int m_count_pip_good_pt;
      int m_count_pim_good_pt;   

      int m_count_gam_good_dr;
      int m_count_mup_good_dr;
      int m_count_mum_good_dr;
      int m_count_pip_good_dr;
      int m_count_pim_good_dr; 

      int m_count_gam_good_dz;
      int m_count_mup_good_dz;
      int m_count_mum_good_dz;
      int m_count_pip_good_dz;
      int m_count_pim_good_dz;    

      int m_count_gam_eff;
      int m_count_mup_eff;
      int m_count_mum_eff;
      int m_count_pip_eff;
      int m_count_pim_eff;

      int m_count_gam_from_pi0;
      
      int m_Ntrk_good;
      int m_Y5S_ndaughters;

      double m_mc_angle_bt_pip_gam_isr;
      double m_mc_angle_bt_pim_gam_isr;
      double m_mc_angle_bt_mup_gam_isr;
      double m_mc_angle_bt_mum_gam_isr;

      double m_mc_angle_bt_pip_gam;
      double m_mc_angle_bt_pim_gam;
      double m_mc_angle_bt_mup_gam;
      double m_mc_angle_bt_mum_gam;

      const char* m_ptype_param;

      int m_total_rho;
      int m_total_YnS;
      int m_total_Wbj;

      //      int m_gam_id;
      //      int m_gam_parent_id;

      int m_pip_parent_id;
      int m_pip_recoil_id;
      int m_pim_parent_id;
      int m_pim_recoil_id;
      int m_mup_parent_id;
      int m_mum_parent_id;

      int m_resonance_fit;
	
      int eventcount;
      
      bool isMC;
      bool isSignalMC;
      bool isISRMC;
      bool isNonISRMC;

      int m_mc_Y5S_id   ;
      int m_mc_Wbj_id   ;
      int m_mc_gam_id   ;
      int m_mc_YnS_id   ;
      int m_mc_rho_id   ;
      int m_mc_mup_id   ;
      int m_mc_mum_id   ;
      int m_mc_pip_id   ;
      int m_mc_pim_id   ;

      int m_mc_vgam_id ;

      int m_mc_YNS_isr_id   ;
      int m_mc_gam_isr_id   ;
      int m_mc_YnS_isr_id   ;
      int m_mc_mup_isr_id   ;
      int m_mc_mum_isr_id   ;
      int m_mc_pip_isr_id   ;
      int m_mc_pim_isr_id   ;

      double m_mc_fsr_YnS_e;
      double m_mc_fsr_YnS_boost_e;
      double m_mc_fsr_rho_e;
      double m_mc_fsr_rho_boost_e;

      // ISR MC

      int m_mc_reweight;

      double m_mc_YNS_isr_e     ;
      double m_mc_YNS_isr_p     ;
      double m_mc_YNS_isr_m     ;
      double m_mc_YNS_isr_pt    ;
      double m_mc_YNS_isr_phi   ;
      double m_mc_YNS_isr_costh ;

      double m_mc_gam_isr_recoil_e     ;
      double m_mc_gam_isr_recoil_p     ;
      double m_mc_gam_isr_recoil_m     ;
      double m_mc_gam_isr_recoil_pt    ;
      double m_mc_gam_isr_recoil_phi   ;
      double m_mc_gam_isr_recoil_costh ;
      
      double m_mc_gam_isr_e     ;
      double m_mc_gam_isr_p     ;
      double m_mc_gam_isr_m     ;
      double m_mc_gam_isr_pt    ;
      double m_mc_gam_isr_phi   ;
      double m_mc_gam_isr_costh ;

      double m_mc_gam_isr_boost_e     ;
      double m_mc_gam_isr_boost_p     ;
      double m_mc_gam_isr_boost_m     ;
      double m_mc_gam_isr_boost_pt    ;
      double m_mc_gam_isr_boost_phi   ;
      double m_mc_gam_isr_boost_costh ;
      
      double m_mc_YnS_isr_e     ;
      double m_mc_YnS_isr_p     ;
      double m_mc_YnS_isr_m     ;
      double m_mc_YnS_isr_pt    ;
      double m_mc_YnS_isr_phi   ;
      double m_mc_YnS_isr_costh ;
      
      double m_mc_pip_pim_isr_e     ;
      double m_mc_pip_pim_isr_p     ;
      double m_mc_pip_pim_isr_m     ;
      double m_mc_pip_pim_isr_pt    ;
      double m_mc_pip_pim_isr_phi   ;
      double m_mc_pip_pim_isr_costh ;
      
      double m_mc_pip_isr_e     ;
      double m_mc_pip_isr_p     ;
      double m_mc_pip_isr_m     ;
      double m_mc_pip_isr_pt    ;
      double m_mc_pip_isr_phi   ;
      double m_mc_pip_isr_costh ;
      
      double m_mc_pim_isr_e     ;
      double m_mc_pim_isr_p     ;
      double m_mc_pim_isr_m     ;
      double m_mc_pim_isr_pt    ;
      double m_mc_pim_isr_phi   ;
      double m_mc_pim_isr_costh ;
      
      double m_mc_mup_isr_e     ;
      double m_mc_mup_isr_p     ;
      double m_mc_mup_isr_m     ;
      double m_mc_mup_isr_pt    ;
      double m_mc_mup_isr_phi   ;
      double m_mc_mup_isr_costh ;
      
      double m_mc_mum_isr_e     ;
      double m_mc_mum_isr_p     ;
      double m_mc_mum_isr_m     ;
      double m_mc_mum_isr_pt    ;
      double m_mc_mum_isr_phi   ;
      double m_mc_mum_isr_costh ;

      double m_mc_pip_pim_gam_isr_recoil_e     ;
      double m_mc_pip_pim_gam_isr_recoil_p     ;
      double m_mc_pip_pim_gam_isr_recoil_m     ;
      double m_mc_pip_pim_gam_isr_recoil_pt    ;
      double m_mc_pip_pim_gam_isr_recoil_phi   ;
      double m_mc_pip_pim_gam_isr_recoil_costh ;

      double m_mc_mup_mum_gam_isr_recoil_e     ;
      double m_mc_mup_mum_gam_isr_recoil_p     ;
      double m_mc_mup_mum_gam_isr_recoil_m     ;
      double m_mc_mup_mum_gam_isr_recoil_pt    ;
      double m_mc_mup_mum_gam_isr_recoil_phi   ;
      double m_mc_mup_mum_gam_isr_recoil_costh ;

      double m_mc_pip_pim_mup_mum_isr_recoil_e     ;
      double m_mc_pip_pim_mup_mum_isr_recoil_p     ;
      double m_mc_pip_pim_mup_mum_isr_recoil_m     ;
      double m_mc_pip_pim_mup_mum_isr_recoil_pt    ;
      double m_mc_pip_pim_mup_mum_isr_recoil_phi   ;
      double m_mc_pip_pim_mup_mum_isr_recoil_costh ;

      double m_mc_pip_pim_mup_mum_isr_recoil_m2 ;

      // Non-ISR MC
      // Y5StoY1Smumu2pic
      double m_mc_pip_pim_e     ;
      double m_mc_pip_pim_p     ;
      double m_mc_pip_pim_m     ;
      double m_mc_pip_pim_pt    ;
      double m_mc_pip_pim_phi   ;
      double m_mc_pip_pim_costh ;

      // Signal MC
      double m_mc_Y5S_e     ;
      double m_mc_Y5S_p     ;
      double m_mc_Y5S_m     ;
      double m_mc_Y5S_pt    ;
      double m_mc_Y5S_phi   ;
      double m_mc_Y5S_costh ;

      double m_mc_Wbj_e     ;
      double m_mc_Wbj_p     ;
      double m_mc_Wbj_m     ;
      double m_mc_Wbj_pt    ;
      double m_mc_Wbj_phi   ;
      double m_mc_Wbj_costh ;

      double m_mc_Wbj_boost_e     ;
      double m_mc_Wbj_boost_p     ;
      double m_mc_Wbj_boost_m     ;
      double m_mc_Wbj_boost_pt    ;
      double m_mc_Wbj_boost_phi   ;
      double m_mc_Wbj_boost_costh ;

      double m_mc_gam_recoil_e     ;
      double m_mc_gam_recoil_p     ;
      double m_mc_gam_recoil_m     ;
      double m_mc_gam_recoil_pt    ;
      double m_mc_gam_recoil_phi   ;
      double m_mc_gam_recoil_costh ;
      
      double m_mc_gam_e     ;
      double m_mc_gam_p     ;
      double m_mc_gam_m     ;
      double m_mc_gam_pt    ;
      double m_mc_gam_phi   ;
      double m_mc_gam_costh ;
      
      double m_mc_YnS_e     ;
      double m_mc_YnS_p     ;
      double m_mc_YnS_m     ;
      double m_mc_YnS_pt    ;
      double m_mc_YnS_phi   ;
      double m_mc_YnS_costh ;

      double m_mc_YnS_boost_e     ;
      double m_mc_YnS_boost_p     ;
      double m_mc_YnS_boost_m     ;
      double m_mc_YnS_boost_pt    ;
      double m_mc_YnS_boost_phi   ;
      double m_mc_YnS_boost_costh ;
      
      double m_mc_rho_e     ;
      double m_mc_rho_p     ;
      double m_mc_rho_m     ;
      double m_mc_rho_pt    ;
      double m_mc_rho_phi   ;
      double m_mc_rho_costh ;
      
      double m_mc_pip_e     ;
      double m_mc_pip_p     ;
      double m_mc_pip_m     ;
      double m_mc_pip_pt    ;
      double m_mc_pip_phi   ;
      double m_mc_pip_costh ;
      
      double m_mc_pim_e     ;
      double m_mc_pim_p     ;
      double m_mc_pim_m     ;
      double m_mc_pim_pt    ;
      double m_mc_pim_phi   ;
      double m_mc_pim_costh ;
      
      double m_mc_mup_e     ;
      double m_mc_mup_p     ;
      double m_mc_mup_m     ;
      double m_mc_mup_pt    ;
      double m_mc_mup_phi   ;
      double m_mc_mup_costh ;
      
      double m_mc_mum_e     ;
      double m_mc_mum_p     ;
      double m_mc_mum_m     ;
      double m_mc_mum_pt    ;
      double m_mc_mum_phi   ;
      double m_mc_mum_costh ;

      double m_mc_pip_pim_gam_recoil_e     ;
      double m_mc_pip_pim_gam_recoil_p     ;
      double m_mc_pip_pim_gam_recoil_m     ;
      double m_mc_pip_pim_gam_recoil_pt    ;
      double m_mc_pip_pim_gam_recoil_phi   ;
      double m_mc_pip_pim_gam_recoil_costh ;

      double m_mc_pip_pim_gam_recoil_boost_e     ;
      double m_mc_pip_pim_gam_recoil_boost_p     ;
      double m_mc_pip_pim_gam_recoil_boost_m     ;
      double m_mc_pip_pim_gam_recoil_boost_pt    ;
      double m_mc_pip_pim_gam_recoil_boost_phi   ;
      double m_mc_pip_pim_gam_recoil_boost_costh ;

      double m_mc_mup_mum_gam_recoil_e     ;
      double m_mc_mup_mum_gam_recoil_p     ;
      double m_mc_mup_mum_gam_recoil_m     ;
      double m_mc_mup_mum_gam_recoil_pt    ;
      double m_mc_mup_mum_gam_recoil_phi   ;
      double m_mc_mup_mum_gam_recoil_costh ;

      double m_mc_mup_mum_gam_fit_recoil_m     ;

      double m_mc_mup_mum_gam_recoil_boost_e     ;
      double m_mc_mup_mum_gam_recoil_boost_p     ;
      double m_mc_mup_mum_gam_recoil_boost_m     ;
      double m_mc_mup_mum_gam_recoil_boost_pt    ;
      double m_mc_mup_mum_gam_recoil_boost_phi   ;
      double m_mc_mup_mum_gam_recoil_boost_costh ;

      double m_mc_pip_pim_mup_mum_recoil_e     ;
      double m_mc_pip_pim_mup_mum_recoil_p     ;
      double m_mc_pip_pim_mup_mum_recoil_m     ;
      double m_mc_pip_pim_mup_mum_recoil_pt    ;
      double m_mc_pip_pim_mup_mum_recoil_phi   ;
      double m_mc_pip_pim_mup_mum_recoil_costh ;

      double m_mc_pip_pim_mup_mum_recoil_m2 ;

      double m_mc_pip_pim_mup_mum_recoil_boost_e     ;
      double m_mc_pip_pim_mup_mum_recoil_boost_p     ;
      double m_mc_pip_pim_mup_mum_recoil_boost_m     ;
      double m_mc_pip_pim_mup_mum_recoil_boost_pt    ;
      double m_mc_pip_pim_mup_mum_recoil_boost_phi   ;
      double m_mc_pip_pim_mup_mum_recoil_boost_costh ;

      double m_mc_gam_boost_e     ;
      double m_mc_gam_boost_p     ;
      double m_mc_gam_boost_m     ;
      double m_mc_gam_boost_pt    ;
      double m_mc_gam_boost_phi   ;
      double m_mc_gam_boost_costh ;

      //double m_mc_gam_recoil_m     ;

      double m_mc_mup_mum_gam_m ;
      double m_mc_pip_pim_gam_m ;

      int m_mc_fsr_YnS;
      int m_mc_fsr_rho;

      int m_isr;

      // background variables

      int m_bkg_A ;
      int m_bkg_B ;
      int m_bkg_C ;
      int m_bkg_D ;
      int m_bkg_E ;
      int m_bkg_F ;
      int m_bkg_G ;
      int m_bkg_H ;
      int m_bkg_I ;
      int m_bkg_J ;
      int m_bkg_K ;


      int m_bkg_A1;
      
      int m_bkg_B1;
      int m_bkg_B2;
      
      int m_bkg_C1;
      int m_bkg_C2;
      int m_bkg_C3;
      
      int m_bkg_D1;
      
      int m_bkg_E1;
      
      int m_bkg_F1;
      
      int m_bkg_G1;
      int m_bkg_G2;
      int m_bkg_G3;
      
      int m_bkg_H1;
      int m_bkg_H2;
      
      int m_bkg_I1;
      
      int m_bkg_J1;

      int m_bkg_K1;

      double m_bkg_Y1S_1_Wbj_fit_m;
      double m_bkg_Y1S_2_Wbj_fit_m;
      double m_bkg_Y1S_3_Wbj_fit_m;
      double m_bkg_Y1S_4_Wbj_fit_m;
      double m_bkg_Y1S_5_Wbj_fit_m;
      double m_bkg_Y1S_6_Wbj_fit_m;
      double m_bkg_Y1S_7_Wbj_fit_m;

      double m_bkg_Y2S_1_Wbj_fit_m;
      double m_bkg_Y2S_2_Wbj_fit_m;

      double m_bkg_Y1S_1_pip_pim_recoil_m;
      double m_bkg_Y1S_2_pip_pim_recoil_m;
      double m_bkg_Y1S_3_pip_pim_recoil_m;
      double m_bkg_Y1S_4_pip_pim_recoil_m;
      double m_bkg_Y1S_5_pip_pim_recoil_m;
      double m_bkg_Y1S_6_pip_pim_recoil_m;
      double m_bkg_Y1S_7_pip_pim_recoil_m;

      double m_bkg_Y2S_1_pip_pim_recoil_m;
      double m_bkg_Y2S_2_pip_pim_recoil_m;

      double m_bkg_Y1S_1_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_2_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_3_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_4_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_5_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_6_pip_pim_gam_recoil_m;
      double m_bkg_Y1S_7_pip_pim_gam_recoil_m;

      double m_bkg_Y2S_1_pip_pim_gam_recoil_m;
      double m_bkg_Y2S_2_pip_pim_gam_recoil_m;

      int m_Y3S_to_X2P;

      int m_Y1S;

      int m_Y2S;   
      int m_Y2S_to_Y1S;

      int m_Y3S; 
      int m_Y3S_to_Y1S;
      int m_Y3S_to_Y2S;
      int m_Y3S_to_Y2S_to_Y1S;

      int m_Y2S_to_pi0pi0;
      int m_Y3S_to_pi0pi0;   
      int m_Y5S_to_pi0pi0; 

      int m_Y2S_to_pippim;
      int m_Y3S_to_pippim;   
      int m_Y5S_to_pippim; 
      
      int m_NChargedTracks;

      int m_trigbits;
      int m_istrig;

      double m_eBeam;
      double m_eBeamCM; // beam energy in CM frame
      Hep3Vector m_cmBoost;// vector to boost to CM frame: -CMBoost()
      Hep3Vector m_ip;
      HepLorentzVector m_pBeamCM;
      HepLorentzVector m_pBeam;

      HepLorentzVector m_pExpCM;
      HepLorentzVector m_pExp;

      //---------------------------------------------------------------------------------------------------------------
      // Declaring member functions  

      //      void prepareEventShapeVariables( Particle signal_part );
      void analyzeHepEvtRecordSignal();
      void analyzeBackgroundSignal();

      void writeTreesOut();

      void printCompleteMCHistory();

      vector<int> getPhotonChain( int id );
      int getIDHEPMCParentOfPhoton( int id );
      int getIDHEPMCParentOfMuons(int id);
      int getIDHEPMCParentOfPions(int id);
      int getIDHEPMCRecoilOfPions(int id);

      void getRelevantChargedIDHEPs(int id);
      void getRelevantPhotonIDHEPs(int id);
      //int getRelevantPhotonIDHEPs_2(int id);

      void selectChargedPions();
      void selectChargedMuons();

      void countChargedTracks();

      void selectPhotons();
      void selectRhos();
      void selectYnS();
      void selectWbj();
      void selectY5S();

      double getNeutralEnergy();
      double getChargedEnergy();
      
      int vetoPi0(int gam_id);

      void prepareCandidates();
      void prepareTreeRec();

      void diagnoseTreeRec();

      Helix getChargedTrackHelixAtIP ( const Mdst_charged &chrg, int mass_hypothesis);

      bool  goodChargedPion( const Mdst_charged& );
      bool  goodChargedMuon( const Mdst_charged& );

      bool goodGamma( const Mdst_gamma& );

      void photonSelectionEfficiencies();
      void chargedSelectionEfficiencies();
      void overallEfficiencies();

      double massFit2(Particle &p, int numDaughters);

      bool useMissingPhoton(double &gam_costh);

      bool goodTrackQuality( const Mdst_charged& chrg , string particle_type);

      void binary2char(unsigned);

      void reweight(double mc_gam_isr_boost_e);

      void clearValuesMC();

      //      double getCosThetaT(std::vector<Hep3Vector>& signalMomenta, std::vector<Hep3Vector>& otherMomenta);
     

      //---------------------------------------------------------------------------------------------------------------
      // Algorithms for best candidate selection
      
      struct bestPID {
        inline bool operator()(candidateRec const &cand1, candidateRec const &cand2) { 
	      double i=abs(cand1.Y5S_m-10.891);
	      double j=abs(cand2.Y5S_m-10.891);
	      return i < j; 
	      }
      };
      
      
      struct bestY5S {
	      inline bool operator()(candidateRec const &cand1, candidateRec const &cand2) { 
	      double i=abs(cand1.Y5S_m-10.891);
	      double j=abs(cand2.Y5S_m-10.891);
	      return i < j; 
	      }
     };

      
      struct bestEnergyY5S {
	      inline bool operator()(candidateRec const &cand1, candidateRec const &cand2) { 
        double i=abs(cand1.deltaE);
        double j=abs(cand2.deltaE);
        return i < j; 
	      }
      };

      struct bestGamWidth {
        inline bool operator()(candidateRec const &cand1, candidateRec const &cand2) { 
        double i=abs(cand1.gam_width);
        double j=abs(cand2.gam_width);
        return i < j; 
        }
      };
      
      struct bestPhotonEnergy {
        inline bool operator()(candidateRec const &cand1, candidateRec const &cand2) { 
	      double i = abs(cand1.gam_e - cand1.pip_pim_mup_mum_recoil_e);
	      double j = abs(cand2.gam_e - cand2.pip_pim_mup_mum_recoil_e);
	      return i < j; 
	      }
      };
    };

  //---------------------------------------------------------------------------------------------------------------
  // this must be done, because of the rules of our framework (basf)
    extern "C" Module_descr *mdcl_Analysis () { /* main */
      Analysis *module = new Analysis;
      Module_descr *dscr = new Module_descr ( "Analysis", module );
      
      //    dscr->define_param("r2max", "Keep events with R2 less than r2max", &module->r2max);
      
      dscr->define_param("file_name_out", "Output file name", "S", 4096, &module->outpath);
      dscr->define_param("inpath", "Input file name", "S", 4096, &module->inpath);
      
      IpProfile::define_global(dscr);
      return dscr;
    }
  
    void Analysis::hist_def( void ) { }
  
#if defined(BELLE_NAMESPACE)
}//namespace Belle
#endif
     
