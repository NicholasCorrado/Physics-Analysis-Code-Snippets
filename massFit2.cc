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
		cout << "Non-zero RC from the fit! kmf_error = " << kmf_error << endl; 
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
	cout << "chi2 = " << chi2 << endl;
	*/  
	kmakemother kmm;
	makeMother(kmm,kmf,p,1);
	kmm.vertex(origin);
	unsigned kmm_error = kmm.make();

	if ( kmm_error != 0 ) {
		cout << "Non-zero RC from kmm,make()! kmm_error = " << kmm_error << endl; 
	}
	/*
	cout << "----------------Before Fit---------------------" << endl;
	cout << "pi: px = " << p.p().x() << endl;
	cout << "pi: py = " << p.p().y() << endl;
	cout << "pi: pz = " << p.p().z() << endl;
	*/

	p.momentum().momentumPosition(kmm.momentum(), // set "PI0" information.
				kmm.position(), // 4-momentum, position and these error.
				kmm.error()); 
	/*
	cout << "----------------After Fit---------------------" << endl;
	cout << "pi: px = " << p.p().x() << endl;
	cout << "pi: py = " << p.p().y() << endl;
	cout << "pi: pz = " << p.p().z() << endl;
	*/
	HepSymMatrix errVtx(3,0);
	p.momentum().decayVertex(kmf.vertex(),errVtx); // set decay point. (error matrix is meaningless.)

return chi2;

}
