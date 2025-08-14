// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <iostream>

// Rosetta headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

#include <numeric/random/random.hh>

int main(int argc, char* argv[])
{
	std::cout << "Hello World!" << std::endl;
	
	// Initialize Rosetta with command-line arguments
	devel::init( argc, argv );

	// Read input PDB filenames
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	} else {
		std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}
	
	// construct a pose object from pdb file
	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

	// update FoldTree based on secondary structure
	mypose->fold_tree( protocols::bootcamp::fold_tree_from_ss( *mypose ) );
	if ( mypose->fold_tree().check_fold_tree() ) {
		std::cout << "FOLDABLE" << std::endl;
	} else {
		std::cout << "NOT FOLDABLE" << std::endl;
	}

	//Add cutpoints to pose
	core::pose::correctly_add_cutpoint_variants( *mypose );
	
	// initialize a ScorFunction object
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function (); ;	
	core::Real score = sfxn->score( *mypose );
	std::cout << "Initial Score: " << score << std::endl;

	// add linear_chainbreak term to scorefunction
	sfxn->set_weight(core::scoring::linear_chainbreak, 1.0);

	// Setup MoveMap, allow the backbone and sidechain to move 	
	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );

	// Setup Minimization 
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;	

	// Create MonteCarlo 
	protocols::moves::MonteCarlo montecarlo = protocols::moves::MonteCarlo( *mypose, *sfxn, 0.5); 

	// If you want to visulize in PyMol
	protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver ( *mypose, true, 0 );
	the_observer->pymol().apply( *mypose );

	// For speed make a working copy of the pose
	core::pose::Pose copy_pose;

	// Random perturbation and Monte Carlo loop
	for ( int i=0; i<500; i++ ) {
		// Random perturbation values	
		core::Real pert1 = numeric::random::gaussian();
		core::Real pert2 = numeric::random::gaussian();
		core::Real uniform_random_number = numeric::random::uniform();
		
		// Choose random residue to perturb
		core::Size N = mypose->size();
		core::Size randres = static_cast< core::Size >( uniform_random_number * N +1);
		
		// Apply random backbone torsion changes
		core::Real orig_phi = mypose->phi( randres );
		core::Real orig_psi = mypose->psi( randres );
		mypose->set_phi( randres, orig_phi + pert1 );
		mypose->set_psi( randres, orig_psi + pert2 );
		
		// Repack side chains
		core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
		repack_task->restrict_to_repacking();
		core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

		// Minimize pose 
		copy_pose = *mypose;
		atm.run( copy_pose, mm, *sfxn, min_opts );
		*mypose = copy_pose;

		// Accept pose based on MonteCarlo
		montecarlo.boltzmann( *mypose );

		// Output score
		std::cout << "New Score: " << sfxn->score( *mypose ) << std::endl;

		// print out acceptance rate every 100 interations 
		if (i % 100 == 0) {
			montecarlo.show_counters();
			std::cout << "Pose's Energy: " << mypose->energies().total_energy() << std::endl; 
		}
	}
	
	std::cout << "Final Score: " << sfxn->score(*mypose) << std::endl;
	return 0;
}
