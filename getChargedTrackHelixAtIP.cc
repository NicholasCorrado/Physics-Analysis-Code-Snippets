Helix getChargedTrackHelixAtIP( const Mdst_charged &chrg, int mass_hypothesis) {
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
