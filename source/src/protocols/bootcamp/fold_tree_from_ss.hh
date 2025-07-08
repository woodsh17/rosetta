// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

namespace protocols {
namespace bootcamp {

// Given a string of secondary structure, returns the residues that span each structured region
utility::vector1< std::pair< core::Size, core::Size > > identify_secondary_structure_spans( std::string const & ss_string );

// Returns a fold tree based on bootcamp from a pose
core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose const & pose);

// when adding an edge to a foldtree, checks that it is not trying to add an edge where the first and last residue are equal
void add_edge_with_check( core::kinematics::FoldTree & foldtree, core::Size first, core::Size last, core::Size label);

// Returns a fold tree based on bootcamp from a secondary structure string
core::kinematics::FoldTree fold_tree_from_dssp_string(std::string const & ss_string);
	
}
}
