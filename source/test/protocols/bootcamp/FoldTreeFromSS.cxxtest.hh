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
		
	// Shared finalization goes here.
	void tearDown() {
	}


};
