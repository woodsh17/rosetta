// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
// libRosetta headers
#include <basic/options/option.hh>
// C++ headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/relax_main.hh>
#include <protocols/jd2/internal_util.hh>



int
main( int argc, char * argv [] )
{
	try {
		using namespace core;
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using namespace basic::options;
		// // initialize core
		// devel::init(argc, argv);
		//
		// utility::Inline_File_Provider *provider = utility::Inline_File_Provider::get_instance();
		//
		// provider->show_contents();
		//
		// std::istream *mystream;
		//
		// if(!provider->get_istream( "scoring/myfile.test", &mystream )){
		//   std::cout << "Cannot find file" << std::endl;
		// }else{
		//  std::cout << "Here" << std::endl;
		//  if( mystream->good() ){
		//   std::cout << "Here" << std::endl;
		//   std::string read_me;
		//   (*mystream) >> read_me;
		//   std::cout << "DATA:" << read_me << std::endl;
		//  }
		// }

		relax::ClassicRelax::register_options();
		jd2::register_options();
		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::relax::fast );
		devel::init(argc, argv);

		return relax::Relax_main( false );


		// utility::io::izstream data( "./scoring/weights/standard.wts" ); // may or may not include .wts extension
		// std::cout << __FILE__ << __LINE__ << std::endl;
		// if( data.good() ){
		// std::cout << __FILE__ << __LINE__ << std::endl;
		//  std::string test;
		//  std::cout << __FILE__ << __LINE__ << std::endl;
		//  data >> test;
		//  std::cout << __FILE__ << __LINE__ << std::endl;
		//  std::cout << test << std::endl;
		//  std::cout << __FILE__ << __LINE__ << std::endl;
		//  data.close();
		//  std::cout << __FILE__ << __LINE__ << std::endl;
		// }
		//std::cout << __FILE__ << __LINE__ << std::endl;
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}


