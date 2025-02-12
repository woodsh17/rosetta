// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_filters/ConstraintScoreCutoffFilter.cc
/// @brief
/// @details
///   Contains currently:
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/constraint_filters/ConstraintScoreCutoffFilter.hh>
#include <protocols/constraint_filters/ConstraintScoreCutoffFilterCreator.hh>

// Package Headers


// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

#include <protocols/jd2/util.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/tag/Tag.hh>


//Auto Headers


//// C++ headers
static basic::Tracer tr( "protocols.filters.ConstraintScoreCutoffFilter" );

namespace protocols {
namespace constraint_filters {


ConstraintScoreCutoffFilter::ConstraintScoreCutoffFilter() :
	parent("ConstraintScoreCutoffFilter"),
	cutoff_( 0.0 )
{}

ConstraintScoreCutoffFilter::ConstraintScoreCutoffFilter( core::Real cutoff_in ) :
	parent("ConstraintScoreCutoffFilter"),
	cutoff_(cutoff_in)
{}

void
ConstraintScoreCutoffFilter::set_score_type( core::scoring::ScoreType setting ) {
	score_type_ = setting;
}

void
ConstraintScoreCutoffFilter::set_constraints( core::scoring::constraints::ConstraintCOPs cst_in ) {
	constraints_ = cst_in;
}


bool
ConstraintScoreCutoffFilter::apply( core::pose::Pose const & pose ) const {
	core::Real cur_score = get_score( pose );
	if ( protocols::jd2::jd2_used() ) protocols::jd2::add_string_real_pair_to_current_job( get_user_defined_name(), cur_score );
	if ( cur_score <= cutoff() ) return true;
	return false;
} // apply_filter


core::Real
ConstraintScoreCutoffFilter::get_score( core::pose::Pose const & pose_in ) const {
	return get_score( pose_in, constraints_ );
}

core::Real
ConstraintScoreCutoffFilter::get_score( core::pose::Pose const & pose_in, core::scoring::constraints::ConstraintCOPs csts ) const {
	using namespace core::scoring;
	core::pose::Pose pose( pose_in );
	pose.constraint_set( nullptr );
	pose.add_constraints( csts );
	ScoreFunction scorefxn;
	scorefxn.set_weight( score_type_, 1.0 );
	return scorefxn( pose );
}

void
ConstraintScoreCutoffFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {
	if ( tag->hasOption("cutoff") ) {
		cutoff_ = tag->getOption<core::Real>("cutoff", 10000.0 );
	}
	if ( tag->hasOption("report_name") ) {
		set_user_defined_name( tag->getOption<std::string>("report_name") );
	} else {
		set_user_defined_name( "cst_cutoff_filter" );
	}
}

void
ConstraintScoreCutoffFilter::report( std::ostream & /*ostr*/, core::pose::Pose const & /*pose*/ ) const
{ }


filters::FilterOP
ConstraintScoreCutoffFilterCreator::create_filter() const { return utility::pointer::make_shared< ConstraintScoreCutoffFilter >(); }

std::string
ConstraintScoreCutoffFilterCreator::keyname() const { return "ConstraintScoreCutoffFilter"; }


} // filters
} // protocols
