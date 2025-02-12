// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file backrub.cc
/// @brief run backrub Monte Carlo
/// @author Colin A. Smith (colin.smith@mpibpc.mpg.de)
/// @details
/// Currently a work in progress. The goal is to match the features of rosetta++ -backrub_mc


#include <devel/init.hh>


// Protocols Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/backrub/BackrubProtocol.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/moves/MoverContainer.hh>

// Core Headers


// Utility Headers
#include <utility/excn/Exceptions.hh>

// Numeric Headers

// Platform Headers

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

void *
my_main( void* );

int
main( int argc, char * argv [] )
{
	try {

		OPT(in::path::database);
		OPT(in::file::s);
		OPT(in::file::l);
		OPT(in::file::movemap);
		OPT(in::ignore_unrecognized_res);
		OPT(out::nstruct);
		OPT(packing::resfile);
		OPT(constraints::cst_fa_weight);
		OPT(constraints::cst_fa_file);
		OPT(backrub::pivot_residues);
		OPT(backrub::pivot_atoms);
		OPT(backrub::min_atoms);
		OPT(backrub::max_atoms);
		OPT(backrub::ntrials);
		OPT(backrub::sc_prob);
		OPT(backrub::sm_prob);
		OPT(backrub::sc_prob_uniform);
		OPT(backrub::sc_prob_withinrot);
		OPT(backrub::mc_kt);
		OPT(backrub::mm_bend_weight);
		OPT(backrub::initial_pack);
		OPT(backrub::minimize_movemap);
		OPT(backrub::trajectory);
		OPT(backrub::trajectory_gz);
		OPT(backrub::trajectory_stride);

		// initialize Rosetta
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}


void *
my_main( void* )
{

	using namespace protocols::moves;
	SequenceMoverOP seqmov( new SequenceMover() );

	using namespace basic::options;
	if ( option[ OptionKeys::in::membrane ].user() ) {
		using namespace protocols::membrane;
		AddMembraneMoverOP add_memb( new AddMembraneMover() );
		seqmov->add_mover( add_memb );
	}

	protocols::backrub::BackrubProtocolOP backrub_protocol( new protocols::backrub::BackrubProtocol() );
	seqmov->add_mover( backrub_protocol );
	protocols::jd2::JobDistributor::get_instance()->go( seqmov );

	// write parameters for any sets of branching atoms for which there were not optimization coefficients
	backrub_protocol->write_database();
	return nullptr;
}
