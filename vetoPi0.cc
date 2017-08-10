  void Analysis::vetoPi0(int gam_id) {

    Mdst_pi0_Manager& Pi0Mgr = Mdst_pi0_Manager::get_manager();
    for ( std::vector<Mdst_pi0>::iterator it_pi0 = Pi0Mgr.begin(); it_pi0 != Pi0Mgr.end(); it_pi0++ ) { 

      Mdst_pi0 &mdst_pi0 = (*it_pi0); 
      Particle pi0(mdst_pi0);
      
      Mdst_gamma const &mdst_gam_1 = pi0.relation().child(0).mdstGamma();
      Mdst_gamma const &mdst_gam_2 = pi0.relation().child(1).mdstGamma();

      if ( gam_id == mdst_gam_1.get_ID() || gam_id == mdst_gam_2.get_ID() ) 
        m_count_gam_from_pi0++;
	      return 1;   
      }
      else {
	      return 0;
      }
    }
  }
