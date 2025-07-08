// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/bootcamp/FoldTreeFromSS.cxxtest.hh
/// @brief

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers

/// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_hello_world() {
		TS_ASSERT( true );
	}
	
	utility::vector1< std::pair< core::Size, core::Size > >
	identify_secondary_structure_spans( std::string const & ss_string )
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

	void test_ss_1() {
	  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries_1 = identify_secondary_structure_spans("   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ");
	  utility::vector1< std::pair < core::Size, core::Size > > spanning_residues = {{4, 8}, {12, 19}, {22, 26}, {36, 41}, {45, 55}, {58, 62}, {65, 68}};
	  for ( core:: Size i = 1; i <= ss_boundaries_1.size(); i++ ) {
	    TS_ASSERT_EQUALS(ss_boundaries_1[i].first, spanning_residues[i].first);
	    TS_ASSERT_EQUALS(ss_boundaries_1[i].second, spanning_residues[i].second);
          }
	} 
	
	void test_ss_2() {
          utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans("HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH  ");
          utility::vector1< std::pair < core::Size, core::Size > > spanning_residues = {{1, 7}, {11, 22}, {29, 40}, {41, 50}, {51, 57}, {59, 62}, {63, 65}};
          for ( core:: Size i = 1; i <= ss_boundaries.size(); i++ ) { 
            TS_ASSERT_EQUALS(ss_boundaries[i].first, spanning_residues[i].first);
            TS_ASSERT_EQUALS(ss_boundaries[i].second, spanning_residues[i].second);
          }   
        }

	void test_ss_3() {
          utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans("EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE");
          utility::vector1< std::pair < core::Size, core::Size > > spanning_residues = {{1, 9}, {11, 18}, {20, 28}, {30, 30}, {32, 36}, {38, 38}, {40, 40}, {42, 42}, {44, 51}};
	  core::Size num_ss_elements(9);
          for ( core::Size i = 1; i <= num_ss_elements; i++ ) {  
            TS_ASSERT_EQUALS(ss_boundaries[i].first, spanning_residues[i].first);
            TS_ASSERT_EQUALS(ss_boundaries[i].second, spanning_residues[i].second);
          }
        }

	core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose const & pose) {
	  core::scoring::dssp::Dssp dssp(pose);
	  std::string ss_string = dssp.get_dssp_secstruct();
	  return fold_tree_from_dssp_string(ss_string); 
	}
	
	core::kinematics::FoldTree fold_tree_from_dssp_string(std::string const & ss_string) {
	  //grab secondary structure boundaries from string
	  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries = identify_secondary_structure_spans(ss_string);
	  core::Size num_ss_elements = ss_boundaries.size();
	
	  core::Size jump_number(1); 

	  //create fold tree
	  core::kinematics::FoldTree foldtree;

	  //find middle residue of first ss element
	  core::Size first_ss_mid_res = ss_boundaries[1].first + ( (ss_boundaries[1].second - ss_boundaries[1].first) / 2);

	 foldtree.add_edge( first_ss_mid_res, 1, core::kinematics::Edge::PEPTIDE );

	  //loop over ss elements 
	  for ( core::Size i = 1; i <= num_ss_elements; i++ ) {

	    core::Size current_ss_first_res = ss_boundaries[i].first;
            core::Size current_ss_last_res = ss_boundaries[i].second;
	    core::Size current_ss_mid_res = current_ss_first_res + ( (current_ss_last_res - current_ss_first_res) / 2);
	
	    if (i != num_ss_elements) {

	      //add peptide edge from middle of current ss to end of current ss
	      foldtree.add_edge( current_ss_mid_res, current_ss_last_res, core::kinematics::Edge::PEPTIDE ); 

	    } else {

	      //add peptide edge from middle of last ss to end of chain 
	      foldtree.add_edge ( current_ss_mid_res, ss_string.size(), core::kinematics::Edge::PEPTIDE );

	    }

	    
	    if (i != 1 ) {

	      //add peptide edge from middle of current ss to first of current ss
	      foldtree.add_edge( current_ss_mid_res, current_ss_first_res, core::kinematics::Edge::PEPTIDE );

	      //add jump between first ss middle residue and current ss middle residus
	      foldtree.add_edge( first_ss_mid_res, current_ss_mid_res, jump_number++ ); 
	    
	      //add jump between first ss middle residue and middle residue of every loop
	      core::Size mid_loop_res = ss_boundaries[i-1].second + ( (current_ss_first_res - ss_boundaries[i-1].second) / 2);
	      foldtree.add_edge( first_ss_mid_res, mid_loop_res, jump_number++ );

	      //once jumps to middle residue of loop has been added, can add peptide edges from middle of loop 
	      foldtree.add_edge( mid_loop_res, ss_boundaries[i-1].second+1, core::kinematics::Edge::PEPTIDE );
	      foldtree.add_edge ( mid_loop_res, current_ss_first_res-1, core::kinematics::Edge::PEPTIDE );
	 
 	    }
	  }
	  return foldtree;
	}

	void test_fold_tree_from_dssp_string() {
	  core::kinematics::FoldTree foldtree = fold_tree_from_dssp_string( "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ");
	  core::Size num_edges(38); 
	  TS_ASSERT_EQUALS( foldtree.size(), num_edges );
	  std::cout << foldtree << std::endl;

	  //Define expected edges: {start, stop, is_peptide}
	  std::vector< std::tuple<core::Size, core::Size, bool> > expected_edges = {
	    {7, 1, true},
	    {7, 10, true},
	    {7, 12, false},
	    {12, 11, true},
	    {12, 14, true},
	    {7, 18, false},
	    {18, 15, true},
	    {18, 21, true},
	    {7, 26, false},
	    {26, 22, true},
	    {26, 30, true},
	    {7, 35, false},
	    {35, 31, true},
	    {35, 39, true},
	    {7, 41, false},
	    {41, 40, true},
	    {41, 43, true},
	    {7, 48, false},
	    {48, 44, true},
	    {48, 53, true},
	    {7, 55, false},
	    {55, 54, true},
	    {55, 56, true},
	    {7, 59, false},
	    {59, 57, true},
	    {59, 62, true},
	    {7, 67, false},
	    {67, 63, true},
	    {67, 71, true},
	    {7, 76, false},
	    {76, 72, true},
	    {76, 80, true},
	    {7, 85, false},
	    {85, 81, true},
	    {85, 89, true},
	    {7, 92, false},
	    {92, 90, true},
	    {92, 99, true}
	  };

	  for (auto const& edge_tuple : expected_edges) {
	    core::Size start, stop;
	    bool is_pep; 
	    std::tie(start, stop, is_pep) = edge_tuple;
	
	    core::kinematics::Edge edge = foldtree.get_residue_edge(stop);
	    TS_ASSERT_EQUALS( edge.start(), start );
	    TS_ASSERT_EQUALS( edge.stop(), stop );
	    TS_ASSERT_EQUALS( edge.is_peptide(), is_pep ); 
	  }
	}
		
	// Shared finalization goes here.
	void tearDown() {
	}


};
