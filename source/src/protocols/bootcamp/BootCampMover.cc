// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief Bootcamp mover
/// @author woodsh17 (hmwoods14@gmail.com)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/variant_util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& mypose){
	std::cout << "Hello World!" << std::endl;

	// update FoldTree based on secondary structure
	mypose.fold_tree( protocols::bootcamp::fold_tree_from_ss( mypose ) );
	if ( mypose.fold_tree().check_fold_tree() ) {
		std::cout << "FOLDABLE" << std::endl;
	} else {
		std::cout << "NOT FOLDABLE" << std::endl;
	}

	//Add cutpoints to pose
	core::pose::correctly_add_cutpoint_variants( mypose );
	
	// initialize a ScorFunction object
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function (); ;	
	core::Real score = sfxn->score( mypose );
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
	protocols::moves::MonteCarlo montecarlo = protocols::moves::MonteCarlo( mypose, *sfxn, 0.5); 

	// For speed make a working copy of the pose
	core::pose::Pose copy_pose;

	// Random perturbation and Monte Carlo loop
	for ( int i=0; i<500; i++ ) {
		// Random perturbation values	
		core::Real pert1 = numeric::random::gaussian();
		core::Real pert2 = numeric::random::gaussian();
		core::Real uniform_random_number = numeric::random::uniform();
		
		// Choose random residue to perturb
		core::Size N = mypose.size();
		core::Size randres = static_cast< core::Size >( uniform_random_number * N +1);
		
		// Apply random backbone torsion changes
		core::Real orig_phi = mypose.phi( randres );
		core::Real orig_psi = mypose.psi( randres );
		mypose.set_phi( randres, orig_phi + pert1 );
		mypose.set_psi( randres, orig_psi + pert2 );
		
		// Repack side chains
		core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( mypose );
		repack_task->restrict_to_repacking();
		core::pack::pack_rotamers( mypose, *sfxn, repack_task );

		// Minimize pose 
		copy_pose = mypose;
		atm.run( copy_pose, mm, *sfxn, min_opts );
		mypose = copy_pose;

		// Accept pose based on MonteCarlo
		montecarlo.boltzmann( mypose );

		// Output score
		std::cout << "New Score: " << sfxn->score( mypose ) << std::endl;

		// print out acceptance rate every 100 interations 
		if (i % 100 == 0) {
			montecarlo.show_counters();
			std::cout << "Pose's Energy: " << mypose.energies().total_energy() << std::endl; 
		}
	}
	
	std::cout << "Final Score: " << sfxn->score(mypose) << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Bootcamp mover", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns woodsh17 as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"woodsh17",
		"TODO: institution",
		"hmwoods14@gmail.com",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
