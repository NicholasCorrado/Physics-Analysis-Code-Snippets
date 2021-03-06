/*
Prints the entire decay tree of a single simulated event, including all particle IDs and MDST information.
We also print additional information that is useful for my analysis in particular:
	1) Whether or not the photon is a noise photon
	2) Distance of a particle's origin to the collision vertex

static variables not defined in this snippet:
	isMC - true if the sample is Monte Carlo Simulation, false otherwise. Initialized 
 */

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
				
			if ( id != 0 ) {
				idhep = mc_particle.idhep();
			}
	
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
