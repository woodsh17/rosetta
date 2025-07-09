// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit Headers
#include <protocols/bootcamp/fold_tree_from_ss.hh>

/// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace bootcamp {

// Given a string of secondary structure, returns the residues that span each structured region
utility::vector1< std::pair< core::Size, core::Size > > identify_secondary_structure_spans( std::string const & ss_string ) 
{
  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
  core::Size strand_start = -1;
  for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
    if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
      if ( int( strand_start ) == -1 ) {
        strand_start = ii;
      } else if ( ss_string[ii] != ss_string[strand_start] ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = ii;
      }
    } else {
      if ( int( strand_start ) != -1 ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = -1;
      }
    }
  } 
  if ( int( strand_start ) != -1 ) {
    // last residue was part of a ss-eleemnt                                                                                           
    ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
  }
  for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
    std::cout << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}

// Returns a fold tree based on bootcamp from a pose
core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose const & pose)
{
  core::scoring::dssp::Dssp dssp(pose);
  std::string ss_string = dssp.get_dssp_secstruct();
  std::cout << "Secondary structure string of pose from DSSP: " << ss_string << std::endl;
  return fold_tree_from_dssp_string(ss_string);
}

// Returns a fold tree based on bootcamp from a secondary structure string
core::kinematics::FoldTree fold_tree_from_dssp_string(std::string const & ss_string)
{
  //grab secondary structure boundaries from string
  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans(ss_string);
  core::Size num_ss_elements = ss_boundaries.size();

  core::Size jump_number(1);

  //create fold tree
  core::kinematics::FoldTree foldtree;

  //find middle residue of first ss element
  core::Size first_ss_mid_res = ss_boundaries[1].first + ( (ss_boundaries[1].second - ss_boundaries[1].first) / 2);

  foldtree.add_edge( first_ss_mid_res, 1, core::kinematics::Edge::PEPTIDE  );

  //loop over ss elements
  for ( core::Size i = 1; i <= num_ss_elements; i++ ) {

    core::Size current_ss_first_res = ss_boundaries[i].first;
    core::Size current_ss_last_res = ss_boundaries[i].second;
    core::Size current_ss_mid_res = current_ss_first_res + ( (current_ss_last_res - current_ss_first_res) / 2);


    if (i != 1 ) {

      //add jump between first ss middle residue and current ss middle residus
      foldtree.add_edge( first_ss_mid_res, current_ss_mid_res, jump_number++  );

      //add jump between first ss middle residue and middle residue of every loop
      core::Size mid_loop_res = ss_boundaries[i-1].second + ( (current_ss_first_res - ss_boundaries[i-1].second) / 2);
      //only add_edge if mid_loop_res is not in last ss_boundaries ie there really is a loop between
      if (mid_loop_res > ss_boundaries[i-1].second ){
        foldtree.add_edge( first_ss_mid_res, mid_loop_res, jump_number++ );
        //once jumps to middle residue of loop has been added, can add peptide edges from middle of loop
        foldtree.add_edge( mid_loop_res, ss_boundaries[i-1].second+1, core::kinematics::Edge::PEPTIDE );
        foldtree.add_edge( mid_loop_res, current_ss_first_res-1, core::kinematics::Edge::PEPTIDE );
      }

      //add peptide edge from middle of current ss to first of current ss
      foldtree.add_edge( current_ss_mid_res, current_ss_first_res, core::kinematics::Edge::PEPTIDE );
    }

    if (i != num_ss_elements) {

      //add peptide edge from middle of current ss to end of current ss
      foldtree.add_edge( current_ss_mid_res, current_ss_last_res, core::kinematics::Edge::PEPTIDE );

    } else {

      //add peptide edge from middle of last ss to end of chain
      foldtree.add_edge( current_ss_mid_res, ss_string.size(), core::kinematics::Edge::PEPTIDE );

    }
  }
  foldtree.delete_self_edges();
  return foldtree;
}

	
}
}
