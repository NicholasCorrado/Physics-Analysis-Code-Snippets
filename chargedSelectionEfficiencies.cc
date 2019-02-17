
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

    // hardcoded cuts--see goodTrackQuality function
      double dr_cut = 0.3;
      double dz_cut = 2.0;
      double pt_cut = 0.1;

      if ( mdst_chrg.charge() > 0 ) {

	if ( hep_chrg.get_ID() == m_mc_pip_id ) {
	  m_count_pip_rec++;
	  if ( fabs(dr_pi) <= dr_cut )               m_count_pip_good_dr++;
	  if ( fabs(dz_pi) <= dz_cut )               m_count_pip_good_dz++;
	  if ( pt > pt_cut )                         m_count_pip_good_pt++;
	  if ( goodTrackQuality( mdst_chrg, "pi" ) ) m_count_pip_good_trk++;
	}
 
	else if ( hep_chrg.get_ID() == m_mc_mup_id ) {
	  m_count_mup_rec++;
	  if ( fabs(dr_mu) <= dr_cut )               m_count_mup_good_dr++;
	  if ( fabs(dz_mu) <= dz_cut )               m_count_mup_good_dz++;
	  if ( pt > pt_cut )                         m_count_mup_good_pt++
	  if ( goodTrackQuality( mdst_chrg, "mu" ) ) m_count_mup_good_trk++;
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
