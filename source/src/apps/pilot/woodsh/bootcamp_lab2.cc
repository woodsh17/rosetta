// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <core/pose/Pose.hh>

int main( int argc, char * argv [] )
{
	std::cout << "Hello World!" << std::endl;
	devel::init( argc, argv );
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	} else {
		std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}
	
	//construct a pose object from pdb file
	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
	
	//initialize a ScorFunction object
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function (); ;	
	//score pose
	core::Real score = sfxn->score( *mypose );
	
	std::cout << "Initial Score: " << score << std::endl;

	core::Real pert1 = numeric::random::gaussian();
	core::Real pert2 = numeric::random::gaussian();
	core::Real uniform_random_number = numeric::random::uniform();
	//grab number of residues in pose
	core::Size N = mypose->size();
	core::Size randres = static_cast< core::Size >( uniform_random_number * N +1);
	core::Real orig_phi = mypose->phi( randres );
	core::Real orig_psi = mypose->psi( randres );
	mypose->set_phi( randres, orig_phi + pert1 );
	mypose->set_psi( randres, orig_psi + pert2 );
	
	return 0;
}

