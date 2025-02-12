// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPAsymEzCB.hh
///
/// @brief  Fullatom implementation of asymetric EZ potential from
/// @details See Schramm et al 2012 at 10.1016/j.str.2012.03.016 for specific details.
///    Implemented in bins of 1A along the Z-axis; depth of 0 is middle of membrane. Positions
///    more than 30A from center are assigned a score of 0. Assigned scores are based on residue identity
///    and on Z-coordinate.
///    Last Modified: 8/21/2020
///
/// @author  Meghan Franklin (meghanwfranklin@gmail.com)

// Unit headers
#include <core/energy_methods/FaMPAsymEzCBEnergy.hh>
#include <core/energy_methods/FaMPAsymEzCBEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>

#include <core/scoring/Energies.hh>


#include <core/scoring/membrane/MembraneData.hh>

#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

static basic::Tracer TR( "core.energy_methods.FaMPAsymEzCBEnergy" );


using namespace core::scoring::methods;

namespace core {
namespace energy_methods {

// Membrane distance limits
static const core::Real MAX_Z_POSITION = 30.5;
static const core::Real Z_BIN_SHIFT = 31;
// e table dimensions
static const core::Real MAX_AA = 20;
static const core::Real ASYMEZ_TABLE_BINS = 61;

// Creator Methods //////////////////////////////////////////

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
FaMPAsymEzCBEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< FaMPAsymEzCBEnergy >();
}

/// @brief Return relevant score types
scoring::ScoreTypes
FaMPAsymEzCBEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( scoring::FaMPAsymEzCB );
	return sts;
}

// Constructors /////////////////////////////////////////////

FaMPAsymEzCBEnergy::FaMPAsymEzCBEnergy() :
	parent( utility::pointer::make_shared< FaMPAsymEzCBEnergyCreator >() )
{
	Size const max_aa( MAX_AA );
	Size const asymEZ_table_bins( ASYMEZ_TABLE_BINS );

	std::string tag,line;
	chemical::AA aa;
	asymEZ_CB_.dimension( max_aa, asymEZ_table_bins);

	// Read in updated table
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/MembranePotential/AsymEZ_CB.txt" );

	Size i=1;
	while ( stream  ) {
		getline( stream, line );
		if ( line[0] == '#' ) continue;
		std::istringstream l(line);
		l >> aa;
		if ( line.size() == 0 ) break;
		for ( Size j=1; j<=asymEZ_table_bins; ++j ) {
			if ( l.fail() ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/AsymEZ_CB.txt (FaMPAsymEzCBEnergy)");
			l >> asymEZ_CB_(aa,j);
		}
		if ( i++ > max_aa ) break;
	}
}

/// @brief Create a clone of this energy method
scoring::methods::EnergyMethodOP
FaMPAsymEzCBEnergy::clone() const
{
	return utility::pointer::make_shared< FaMPAsymEzCBEnergy >( *this );
}

// Scoring Methods ////////////////////////////////////////////////

/// @details looks up score for depth of each CB in membrane.
void
FaMPAsymEzCBEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	scoring::EnergyMap & emap
) const
{
	// asymEZ was only developed from/for proteins
	if ( ! ( rsd.is_protein() &&  rsd.type().is_canonical_aa() ) ) return;

	Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

	Real score = 0;

	Real const z_position( pose.conformation().membrane_info()->atom_z_position( pose.conformation(), rsd.seqpos(), atomindex_i ) );

	if ( (z_position <= (-1*MAX_Z_POSITION)) || (z_position >= MAX_Z_POSITION) ) {
		score = 0; //outside the membrane
	} else {
		//z-bin is the rounded z-coord + 30 (because Rosetta is indexed from 1, not -30)
		Size z_bin = static_cast< Size > ( round(z_position) + Z_BIN_SHIFT);
		score = asymEZ_CB_(rsd.aa(), z_bin);
		//std::cout << rsd.aa() << rsd.seqpos() << "\tz = " << z_position << "\tz_bin = " << z_bin << "\tFaMPAsymEzCB = " << score << std::endl;
	}

	emap[ scoring::FaMPAsymEzCB ] += score;
}


/// @details deals with the special cases of Gly - calculates position based on the CA position instead
std::string const &
FaMPAsymEzCBEnergy::representative_atom_name( chemical::AA const aa ) const
{
	// intend to error in this case
	debug_assert( aa >= 1 && aa <= chemical::num_canonical_aas );

	static std::string const cbeta_string(  "CB"  );
	static std::string const calpha_string( "CA"  );

	switch ( aa ) {
	case ( chemical::aa_ala ) : return cbeta_string;
	case ( chemical::aa_cys ) : return cbeta_string;
	case ( chemical::aa_asp ) : return cbeta_string;
	case ( chemical::aa_glu ) : return cbeta_string;
	case ( chemical::aa_phe ) : return cbeta_string;
	case ( chemical::aa_gly ) : return calpha_string;
	case ( chemical::aa_his ) : return cbeta_string;
	case ( chemical::aa_ile ) : return cbeta_string;
	case ( chemical::aa_lys ) : return cbeta_string;
	case ( chemical::aa_leu ) : return cbeta_string;
	case ( chemical::aa_met ) : return cbeta_string;
	case ( chemical::aa_asn ) : return cbeta_string;
	case ( chemical::aa_pro ) : return cbeta_string;
	case ( chemical::aa_gln ) : return cbeta_string;
	case ( chemical::aa_arg ) : return cbeta_string;
	case ( chemical::aa_ser ) : return cbeta_string;
	case ( chemical::aa_thr ) : return cbeta_string;
	case ( chemical::aa_val ) : return cbeta_string;
	case ( chemical::aa_trp ) : return cbeta_string;
	case ( chemical::aa_tyr ) : return cbeta_string;
	default :
		utility_exit_with_message( "ERROR: Failed to find amino acid " + chemical::name_from_aa( aa ) + " in FAMPAsymEZCB::representative_atom_name" );
		break;
	}
	// unreachable
	return calpha_string;
}

void
FaMPAsymEzCBEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}

core::Size
FaMPAsymEzCBEnergy::version() const
{
	return 1; // Initial versioning
}

} // energy_methods
} // core
