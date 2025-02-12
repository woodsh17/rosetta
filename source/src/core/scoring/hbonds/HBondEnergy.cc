// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HBondEnergy.fwd.hh
/// @brief  Hydrogen bond energy method forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit Headers
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/hbonds/HBondEnergyCreator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/hbonds/hbtrie/HBCPData.hh>
#include <core/scoring/hbonds/hbtrie/HBondTrie.fwd.hh>
#include <core/scoring/hbonds/hbtrie/HBCountPairFunction.hh>
#include <core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.hh>

#include <core/scoring/MinimizationData.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/trie/TrieCollection.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>

#include <core/scoring/dssp/Dssp.hh>  // for helix-len-dep hbond strength

#include <core/scoring/methods/EnergyMethodOptions.hh>

// membrane specific initialization
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneInfo.hh>

#include <basic/datacache/BasicDataCache.fwd.hh>

#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableData.hh>

#include <core/scoring/func/SmoothStepFunc.hh>

// Numeric Headers

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <utility/vector1.hh>
#include <boost/unordered_map.hpp>


#ifdef SERIALIZATION
// Project serialization headers
#include <core/scoring/trie/RotamerTrie.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.scoring.hbonds.HBondEnergy" );

namespace core {
namespace scoring {
namespace hbonds {

class HBondResidueMinData;
class HBondResPairMinData;

using HBondResidueMinDataOP = utility::pointer::shared_ptr<HBondResidueMinData>;
using HBondResidueMinDataCOP = utility::pointer::shared_ptr<const HBondResidueMinData>;

using HBondResPairMinDataOP = utility::pointer::shared_ptr<HBondResPairMinData>;
using HBondResPairMinDataCOP = utility::pointer::shared_ptr<const HBondResPairMinData>;

/// @brief A class to hold data for the HBondEnergy class used in
/// score and derivative evaluation.
class HBondResidueMinData : public basic::datacache::CacheableData
{
public:
	HBondResidueMinData()= default;
	~HBondResidueMinData() override = default;

	basic::datacache::CacheableDataOP clone() const override
	{ return utility::pointer::make_shared< HBondResidueMinData >( *this ); }

	void set_bb_don_avail( bool setting ) { bb_don_avail_ = setting; }
	void set_bb_acc_avail( bool setting ) { bb_acc_avail_ = setting; }

	bool bb_don_avail() const { return bb_don_avail_; }
	bool bb_acc_avail() const { return bb_acc_avail_; }

	void set_natoms( Size setting ) { natoms_ = setting; }
	Size natoms() const { return natoms_; }

	void set_nneighbors( Size setting ) { nneighbors_ = setting; }
	Size nneighbors() const { return nneighbors_; }

private:
	Size natoms_;
	Size nneighbors_;

	bool bb_don_avail_{ true };
	bool bb_acc_avail_{ true };
};

class HBondResPairMinData : public basic::datacache::CacheableData
{
public:
	HBondResPairMinData();
	~HBondResPairMinData() override = default;
	basic::datacache::CacheableDataOP clone() const override
	{ return utility::pointer::make_shared< HBondResPairMinData >( *this ); }

	void set_res1_data( HBondResidueMinDataCOP );
	void set_res2_data( HBondResidueMinDataCOP );

	HBondResidueMinData const & res1_data() const { return *res1_dat_; }
	HBondResidueMinData const & res2_data() const { return *res2_dat_; }

	//void update_natoms();

	void clear_hbonds();
	void add_hbond( HBond const & );

private:

	HBondResidueMinDataCOP res1_dat_;
	HBondResidueMinDataCOP res2_dat_;

};

HBondResPairMinData::HBondResPairMinData() = default;
void HBondResPairMinData::set_res1_data( HBondResidueMinDataCOP dat ) { res1_dat_ = dat; }
void HBondResPairMinData::set_res2_data( HBondResidueMinDataCOP dat ) { res2_dat_ = dat; }

static basic::Tracer tr( "core.scoring.hbonds.HbondEnergy" );

/// @details This must return a fresh instance of the HBondEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
HBondEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< HBondEnergy >( options.hbond_options() );
}

ScoreTypes
HBondEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hbond_lr_bb );
	sts.push_back( hbond_sr_bb );
	sts.push_back( hbond_bb_sc );
	sts.push_back( hbond_sr_bb_sc );
	sts.push_back( hbond_lr_bb_sc );
	sts.push_back( hbond_sc );
	sts.push_back( hbond_wat ); // hydrate/SPaDES protocol
	sts.push_back( wat_entropy ); // hydrate/SPaDES protocol
	sts.push_back( hbond_intra ); //Currently affects only RNA.
	sts.push_back( hbond );
	return sts;
}


/// ctor
HBondEnergy::HBondEnergy( HBondOptions const & opts ):
	parent( utility::pointer::make_shared< HBondEnergyCreator >() ),
	options_( utility::pointer::make_shared< HBondOptions >( opts )),
	database_( HBondDatabase::get_database(opts.params_database_tag()) )
{}

/// copy ctor
HBondEnergy::HBondEnergy( HBondEnergy const & src ):
	parent( src ),
	options_( utility::pointer::make_shared< HBondOptions >( *src.options_)) ,
	database_( src.database_)
{}

HBondEnergy::~HBondEnergy() = default;

/// clone
methods::EnergyMethodOP
HBondEnergy::clone() const
{
	return utility::pointer::make_shared< HBondEnergy >( *this );
}


void
HBondEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	using EnergiesCacheableDataType::HBOND_SET;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	pose.update_residue_neighbors();
	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( *options_ ) );

	// membrane object initialization
	if ( options_->Mbhbond() ) {
		Membrane_FAPotential const & memb_potential_ = ScoringManager::get_instance()->get_Membrane_FAPotential();
		memb_potential_.compute_fa_projection( pose );
		normal_ = MembraneEmbed_from_pose( pose ).normal();
		center_ = MembraneEmbed_from_pose( pose ).center();
		thickness_ = Membrane_FAEmbed_from_pose( pose ).thickness();
		steepness_ = Membrane_FAEmbed_from_pose( pose ).steepness();

		// membrane framework object initialization
	} else if ( options_->mphbond() ) {

		// Initialize membrane specific parameters
		normal_ = pose.conformation().membrane_info()->membrane_normal(pose.conformation());
		center_ = pose.conformation().membrane_info()->membrane_center(pose.conformation());
		thickness_ = pose.conformation().membrane_info()->membrane_thickness();
		steepness_ = pose.conformation().membrane_info()->membrane_steepness();
		membrane_core_ = pose.conformation().membrane_info()->membrane_core();

	} else {
		thickness_ = 0; // No membrane hydrogen bonding correction
	}

	hbond_set->setup_for_residue_pair_energies( pose );
	pose.energies().data().set( HBOND_SET, hbond_set );

	using namespace trie;
	using namespace hbtrie;

	TrieCollectionOP tries( new TrieCollection );
	tries->total_residue( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Do not compute energy for virtual residues.
		if ( pose.residue(ii).aa() == core::chemical::aa_vrt ) continue;

		HBondRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( ii ), pose );
		tries->trie( ii, one_rotamer_trie );
	}
	pose.energies().data().set( HBOND_TRIE_COLLECTION, tries );
}

void
HBondEnergy::prepare_rotamers_for_packing(
	pose::Pose const & pose,
	conformation::RotamerSetBase & set
) const
{
	using namespace hbtrie;

	HBondRotamerTrieOP rottrie = create_rotamer_trie( set, pose );
	//std::cout << "--------------------------------------------------" << std::endl << " HBROTAMER TRIE: " << set.resid() << std::endl;
	//rottrie->print();
	set.store_trie( methods::hbond_method, rottrie );
}

// Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
HBondEnergy::update_residue_for_packing( pose::Pose & pose, Size resid ) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	HBondRotamerTrieOP one_rotamer_trie = create_rotamer_trie( pose.residue( resid ), pose );

	// grab non-const & of the cached tries and replace resid's trie with a new one.
	auto & trie_collection( static_cast< TrieCollection & > (pose.energies().data().get( HBOND_TRIE_COLLECTION )));
	trie_collection.trie( resid, one_rotamer_trie );
}


void
HBondEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	pose.update_residue_neighbors();
	HBondSetOP hbond_set( new hbonds::HBondSet( *options_, pose.size() ) );

	// membrane object initialization
	if ( options_->Mbhbond() ) {
		Membrane_FAPotential const & memb_potential_ = ScoringManager::get_instance()->get_Membrane_FAPotential();
		memb_potential_.compute_fa_projection( pose );
		normal_ = MembraneEmbed_from_pose( pose ).normal();
		center_ = MembraneEmbed_from_pose( pose ).center();
		thickness_ = Membrane_FAEmbed_from_pose( pose ).thickness();
		steepness_ = Membrane_FAEmbed_from_pose( pose ).steepness();

		// membrane framework supported membrane object initialization
	} else if ( options_->mphbond() ) {

		// Initialize membrane parameters
		normal_ = pose.conformation().membrane_info()->membrane_normal(pose.conformation());
		center_ = pose.conformation().membrane_info()->membrane_center(pose.conformation());
		thickness_ = pose.conformation().membrane_info()->membrane_thickness();
		steepness_ = pose.conformation().membrane_info()->membrane_steepness();

	} else {
		thickness_ = 0; // no membrane hydrogen bonding correction
		membrane_core_ = 0;
	}

	//fpd  we need secstruct info in some cases ... don't change while minimizing
	//fpd      MUST be called before setup_for_residue_pair_energies
	if ( options_->length_dependent_srbb() && ! pose.energies().use_nblist() ) {
		dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	}

	hbond_set->setup_for_residue_pair_energies( pose );

	/// During minimization, keep the set of bb/bb hbonds "fixed" by using the old boolean values.
	if ( pose.energies().use_nblist() && pose.energies().data().has( HBOND_SET ) ) {
		auto const & existing_set = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
		hbond_set->copy_bb_donor_acceptor_arrays( existing_set );
	}
	pose.energies().data().set( HBOND_SET, hbond_set );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @details Note that this only evaluates sc-sc and sc-bb energies unless options_->decompose_bb_hb_into_pair_energies
/// is set to true, in which case, this function also evaluates bb-bb energies.
/// Note also that this function enforces the bb/sc hbond exclusion rule.
void
HBondEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET )));

	// hydrate/SPaDES protocol
	bool bond_near_wat = false;
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	// this only works because we have already called
	// hbond_set->setup_for_residue_pair_energies( pose )

	// Non-pairwise additive exclusion rules:
	// exclude backbone-backbone hbond if set in options_ (if, say, they were pre-computed)
	// exclude backbone-sidechain hbond if backbone-backbone hbond already in hbond_set*
	// exclude sidechain-backbone hbond if backbone-backbone hbond already in hbond_set*
	// * these two rules are only enforced as long as bb_donor_acceptor_check is "true"

	// mjo Historically, if this exclusion rule is not
	// enforced--accoring to Brian Kuhlman--"Serines are put up and down
	// helices".  According to John Karanicolas, amide acceptors have a
	// local energy minima where one lone pair moves in line with the
	// base-acceptor, and the other lone pair is delocalized in between
	// the base and acceptor atoms.  In this configuration it is
	// energetically disfavorable to make multiple hbonds with the
	// acceptor.

	// NOTE: "bsc" -> acc=bb don=sc
	//       "scb" -> don=sc don=bb

	bool exclude_bsc = false, exclude_scb = false;
	if ( rsd1.is_protein() ) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos());
	if ( rsd2.is_protein() ) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd2.seqpos());

	// Adjust hydrogen bonding potential to accomodate stronger hbonding in the
	// membrane hydrophobic core. Incorporating automatic deteciton for membrane poses (mpframework)
	if ( options_->Mbhbond() || options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap,
			pose,
			bond_near_wat);

		exclude_bsc = exclude_scb = false;
		if ( rsd2.is_protein() ) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd2.seqpos());
		if ( rsd1.is_protein() ) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos());

		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap,
			pose,
			bond_near_wat);

	} else {

		//ss-dep weights
		SSWeightParameters ssdep;
		ssdep.ssdep_ = options_->length_dependent_srbb();
		ssdep.l_ = options_->length_dependent_srbb_lowscale();
		ssdep.h_ = options_->length_dependent_srbb_highscale();
		ssdep.len_l_ = options_->length_dependent_srbb_minlength();
		ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
		Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

		boost::unordered_map< core::Size, core::Size > num_hbonds; //Not actually used, but required for identify_hbonds_1way function signature.

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap, num_hbonds, ssdep_weight_factor, bond_near_wat );

		exclude_bsc = exclude_scb = false;
		if ( rsd2.is_protein() ) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(rsd2.seqpos());
		if ( rsd1.is_protein() ) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos());

		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			!options_->decompose_bb_hb_into_pair_energies(), exclude_bsc, exclude_scb, false,
			*options_,
			emap, num_hbonds, ssdep_weight_factor, bond_near_wat );
	}
}

bool
HBondEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool res_moving_wrt_eachother
) const
{
	return res_moving_wrt_eachother;
}

bool
HBondEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

bool
HBondEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

/// @details computes the residue pair energy during minimization; this includes bb/bb energies
/// as opposed to the standard residue_pair_energy interface, which does not include bb/bb energies.
/// On the other hand, this interface presumes that no new bb/bb hydrogen bonds are formed during
/// the course of minimization, and no existing bb/bb hydrogen bonds are lost.
/// Note: this function does not directly enforce the bb/sc exclusion rule logic, but rather,
/// takes the boolean "bb_don_avail" and "bb_acc_avail" data stored in the pairdata object
void
HBondEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & pairdata,
	pose::Pose const & pose, //pba
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) )
			> std::pow( rsd1.nbr_radius() + rsd2.nbr_radius() + atomic_interaction_cutoff(), 2 ) ) return;

	debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData const > ( pairdata.get_data( hbond_respair_data ) ));
	auto const & hb_pair_dat( static_cast< HBondResPairMinData const & > ( pairdata.get_data_ref( hbond_respair_data ) ));

	using EnergiesCacheableDataType::HBOND_SET; // jklai
	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ))); // jklai

	// hydrate/SPaDES protocol
	bool bond_near_wat = false;
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	// Adjust hydrogen bonding potential to accomodate stronger hbonding in the
	// membrane hydrophobic core. Incorporating automatic deteciton for membrane poses (mpframework)
	if ( options_->Mbhbond() || options_->mphbond() ) {

		{ // scope        /// 1st == evaluate hbonds with donor atoms on rsd1
			/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
			bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
			/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
			bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

			identify_hbonds_1way_membrane(
				*database_,
				rsd1, rsd2,
				hb_pair_dat.res1_data().nneighbors(), hb_pair_dat.res2_data().nneighbors(),
				false /*calculate_derivative*/,
				false, exclude_bsc, exclude_scb, false,
				*options_,
				emap, pose, bond_near_wat);
		}

		{ // scope
			/// 2nd == evaluate hbonds with donor atoms on rsd2
			/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
			bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
			/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
			bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );

			identify_hbonds_1way_membrane(
				*database_,
				rsd2, rsd1,
				hb_pair_dat.res2_data().nneighbors(), hb_pair_dat.res1_data().nneighbors(),
				false /*calculate_derivative*/,
				false, exclude_bsc, exclude_scb, false,
				*options_,
				emap, pose, bond_near_wat);
		}

	} else {

		//ss-dep weights
		SSWeightParameters ssdep;
		ssdep.ssdep_ = options_->length_dependent_srbb();
		ssdep.l_ = options_->length_dependent_srbb_lowscale();
		ssdep.h_ = options_->length_dependent_srbb_highscale();
		ssdep.len_l_ = options_->length_dependent_srbb_minlength();
		ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
		Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

		{ // scope
			/// 1st == evaluate hbonds with donor atoms on rsd1
			/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
			bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
			/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
			bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

			identify_hbonds_1way(
				*database_,
				rsd1, rsd2,
				hb_pair_dat.res1_data().nneighbors(), hb_pair_dat.res2_data().nneighbors(),
				false /*calculate_derivative*/,
				false, exclude_bsc, exclude_scb, false,
				*options_,
				emap, ssdep_weight_factor, bond_near_wat);
		}

		{ // scope
			/// 2nd == evaluate hbonds with donor atoms on rsd2
			/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
			bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
			/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
			bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );

			identify_hbonds_1way(
				*database_,
				rsd2, rsd1,
				hb_pair_dat.res2_data().nneighbors(), hb_pair_dat.res1_data().nneighbors(),
				false /*calculate_derivative*/,
				false, exclude_bsc, exclude_scb, false,
				*options_,
				emap, ssdep_weight_factor, bond_near_wat);
		}
	}
}

/// @details Note that this function helps enforce the bb/sc exclusion rule by setting the donor and acceptor availability
/// for backbone donors and acceptors.  If the backbone-sidechain-exclusion rule is not being enforced, then this function
/// marks all donors and acceptors as being available.  If it is being enforced, then is uses the hbondset functions
/// don_bbg_in_bb_bb_hbond and acc_bbg_in_bb_bb_hbond.  The decisions made in this function impact the evaluation
/// of energies in the above residue_pair_energy_ext method.
void
HBondEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	basic::datacache::BasicDataCache &,
	ResSingleMinimizationData & res_data_cache
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	auto const & hbondset = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
	HBondResidueMinDataOP hbresdata( nullptr );
	if ( res_data_cache.get_data( hbond_res_data ) ) {
		hbresdata = utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData > ( res_data_cache.get_data( hbond_res_data ) );
		// assume that bb-don-avail and bb-acc-avail are already initialized
	} else {
		hbresdata = utility::pointer::make_shared< HBondResidueMinData >();
		hbresdata->set_nneighbors( hbondset.nbrs( rsd.seqpos() ) );
		if ( rsd.is_protein() ) {
			hbresdata->set_bb_don_avail( options_->bb_donor_acceptor_check() ? ! hbondset.don_bbg_in_bb_bb_hbond( rsd.seqpos() ) : true );
			hbresdata->set_bb_acc_avail( options_->bb_donor_acceptor_check() ? ! hbondset.acc_bbg_in_bb_bb_hbond( rsd.seqpos() ) : true );
		}
		res_data_cache.set_data( hbond_res_data, hbresdata );
	}
	hbresdata->set_natoms( rsd.natoms() );
}


void
HBondEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const & res1_data_cache,
	ResSingleMinimizationData const & res2_data_cache,
	ResPairMinimizationData & data_cache
) const
{
	HBondResPairMinDataOP hbpairdat;
	if ( data_cache.get_data( hbond_respair_data ) ) {
		// assume that hbpairdat has already been pointed at its two residues, and that it may need to update
		// its size based on a change to the number of atoms from a previous initialization.
		debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData > ( data_cache.get_data( hbond_respair_data ) ));
		hbpairdat = utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResPairMinData > ( data_cache.get_data( hbond_respair_data ) );
	} else {
		debug_assert( utility::pointer::dynamic_pointer_cast< HBondResidueMinData const > ( res1_data_cache.get_data( hbond_res_data ) ));
		debug_assert( utility::pointer::dynamic_pointer_cast< HBondResidueMinData const > ( res2_data_cache.get_data( hbond_res_data ) ));

		hbpairdat = utility::pointer::make_shared< HBondResPairMinData >();
		hbpairdat->set_res1_data( utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData const > ( res1_data_cache.get_data( hbond_res_data ) ));
		hbpairdat->set_res2_data( utility::pointer::static_pointer_cast< core::scoring::hbonds::HBondResidueMinData const > ( res2_data_cache.get_data( hbond_res_data ) ));
		data_cache.set_data( hbond_respair_data, hbpairdat );
	}
	//hbpairdat->update_natoms();
}

bool
HBondEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}


/// @details Triplication of the loops that iterate across hbond donors and acceptors
/// Find all hbonds for a pair of residues and add those found hbonds to the
/// hb_pair_dat object; these hbonds will be used for derivative evaluation, so evaluate
/// the F1/F2 derivative vectors now.  This function respects the exclude bsc and scb variables
/// to avoid hbonds.  Called by setup_for_derivatives_for_residue_pair
/// pba modified
void
HBondEnergy::hbond_derivs_1way(
	EnergyMap const & weights,
	HBondSet const & hbond_set,
	HBondDatabaseCOP database,
	pose::Pose const & pose,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	Real const ssdep_weight_factor,
	// output
	utility::vector1< DerivVectorPair > & don_atom_derivs,
	utility::vector1< DerivVectorPair > & acc_atom_derivs,
	bool bond_near_wat
) const
{
	EnergyMap emap;

	bool is_intra_res = (don_rsd.seqpos() == acc_rsd.seqpos());
	if ( is_intra_res && !calculate_intra_res_hbonds( don_rsd, hbond_set.hbond_options() ) ) return;

	// <f1,f2> -- derivative vectors
	HBondDerivs deriv;

	for ( Size const hatm : don_rsd.Hpos_polar() ) {
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		Vector const & hatm_xyz( don_rsd.atom(hatm).xyz() );
		Vector const & datm_xyz( don_rsd.atom(datm).xyz() );

		for ( Size const aatm : acc_rsd.accpt_pos() ) {

			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( ! datm_is_bb && exclude_bsc ) continue; // if the donor is sc, the acceptor bb, and exclude_b(a)sc(d)
			} else {
				if ( datm_is_bb && exclude_scb ) continue; // if the donor is bb, the acceptor sc, and exclude_sc(a)b(d)
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			debug_assert( base2 > 0 && base != base2 );

			hb_energy_deriv( *database, *options_, hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, true /*eval deriv*/, deriv);

			if ( unweighted_energy >= options_->max_hb_energy() ) continue;

			//pba buggy?
			Real weighted_energy = // evn weight * weight-set[ hbtype ] weight
				(! hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options() )) *
				hb_eval_type_weight( hbe_type.eval_type(), weights, is_intra_res, hbond_set.hbond_options().put_intra_into_total() );
			weighted_energy *= ssdep_weight_factor;

			// hydrate/SPaDES protocol
			// don't consider hb env dependency if hybrid hb env dependency and hb is near wat
			if ( hbond_set.hbond_options().water_hybrid_sf() && bond_near_wat ) {
				weighted_energy = hb_eval_type_weight( hbe_type.eval_type(), weights, is_intra_res);
			}

			// Readjust hydrogen bonding depth dependent weight based on z positions
			// Relying on nonzero thickness which should really be true here!!!
			if ( thickness_ != 0 || options_->Mbhbond() || options_->mphbond()  ) {
				weighted_energy = get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, don_rsd.seqpos(), acc_rsd.seqpos(), hatm, aatm, hatm_xyz, acc_rsd.atom(aatm ).xyz()) *
					hb_eval_type_weight( hbe_type.eval_type(), weights, is_intra_res, hbond_set.hbond_options().put_intra_into_total() );
			}

			HBDerivAssigner assigner( *options_, hbe_type, don_rsd, hatm, acc_rsd, aatm );
			for ( Size ii = 1; ii <= n_hb_atoms; ++ii ) {
				auto ii_which = which_atom_in_hbond(ii);
				if ( assigner.ind( ii_which ) == 0 ) continue;
				AssignmentScaleAndDerivVectID ii_asadvi = assigner.assignment( ii_which );
				if ( ii_asadvi.dvect_id_ == which_hb_unassigned ) continue;
				DerivVectorPair ii_deriv = ii_asadvi.scale_ * weighted_energy * deriv.deriv( ii_asadvi.dvect_id_ );
				if ( ii <= which_last_donor_atm ) {
					don_atom_derivs[ assigner.ind( ii_which ) ] += ii_deriv;
				} else {
					acc_atom_derivs[ assigner.ind( ii_which ) ] += ii_deriv;
				}
			}

		} // loop over acceptors
	} // loop over donors
}

void
HBondEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{

	if ( ! calculate_intra_res_hbonds( rsd, *options_ ) ) return;

	using EnergiesCacheableDataType::HBOND_SET;
	auto const & hbond_set = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));
	bool const exclude_scb( false );

	// hydrate/SPaDES protocol
	bool bond_near_wat = false; // used for water hybrid sf
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd.seqpos() )  ) bond_near_wat = true;
	}

	hbond_derivs_1way( weights, hbond_set, database_, pose, rsd, rsd, 1, 1, exclude_scb, exclude_scb, 1, atom_derivs, atom_derivs, bond_near_wat );
}

void
HBondEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // stores the hbond database in the cached hbond set
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) )
			> std::pow( rsd1.nbr_radius() + rsd2.nbr_radius() + atomic_interaction_cutoff(), 2 ) ) return;

	/// Iterate across all acceptor and donor atom pairs for these two residues, and write down the hydrogen bonds
	/// that are formed.

	using EnergiesCacheableDataType::HBOND_SET;

	auto const & hbondset = static_cast< HBondSet const & > (pose.energies().data().get( HBOND_SET ));

	debug_assert( utility::pointer::dynamic_pointer_cast< HBondResPairMinData const > ( min_data.get_data( hbond_respair_data ) ));
	auto const & hb_pair_dat = static_cast< HBondResPairMinData const & > ( min_data.get_data_ref( hbond_respair_data ) );

	Size const rsd1nneighbs( hb_pair_dat.res1_data().nneighbors() );
	Size const rsd2nneighbs( hb_pair_dat.res2_data().nneighbors() );

	//ss-dep weights
	SSWeightParameters ssdep;
	ssdep.ssdep_ = options_->length_dependent_srbb();
	ssdep.l_ = options_->length_dependent_srbb_lowscale();
	ssdep.h_ = options_->length_dependent_srbb_highscale();
	ssdep.len_l_ = options_->length_dependent_srbb_minlength();
	ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
	Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

	// hydrate/SPaDES protocol
	bool bond_near_wat = false; // water hybrid sf
	if ( hbondset.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	{ // scope
		/// 1st == find hbonds with donor atoms on rsd1
		/// case A: sc is acceptor, bb is donor && res2 is the acceptor residue -> look at the donor availability of residue 1
		bool exclude_scb( ! hb_pair_dat.res1_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res2 is the acceptor residue -> look at the acceptor availability of residue 2
		bool exclude_bsc( ! hb_pair_dat.res2_data().bb_acc_avail() );

		hbond_derivs_1way( weights, hbondset, database_, pose, rsd1, rsd2, rsd1nneighbs, rsd2nneighbs, exclude_bsc, exclude_scb, ssdep_weight_factor, r1_atom_derivs, r2_atom_derivs, bond_near_wat );
	}

	{ // scope
		/// 2nd == evaluate hbonds with donor atoms on rsd2
		/// case A: sc is acceptor, bb is donor && res1 is the acceptor residue -> look at the donor availability of residue 2
		bool exclude_scb( ! hb_pair_dat.res2_data().bb_don_avail() );
		/// case B: bb is acceptor, sc is donor && res1 is the acceptor residue -> look at the acceptor availability of residue 1
		bool exclude_bsc( ! hb_pair_dat.res1_data().bb_acc_avail() );

		hbond_derivs_1way( weights, hbondset, database_, pose, rsd2, rsd1, rsd2nneighbs, rsd1nneighbs, exclude_bsc, exclude_scb, ssdep_weight_factor, r2_atom_derivs, r1_atom_derivs, bond_near_wat );
	}
}

void
HBondEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	/// if we're including bb/bb energies in the energy maps, then they need to be calculated
	/// in the backbone_backbone_energy so that:
	/// residue_pair_energy = backbone_backbone_energy(r1,r2) + backbone_sidechain_energy(r1,r2) +
	///                       backbone_sidechain_energy(r2,r1) + sidechain_sidechain_energy(r1,r2)

	if ( ! options_->decompose_bb_hb_into_pair_energies() ) return;

	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >( pose.energies().data().get( HBOND_SET )));

	// hydrate/SPaDES protocol
	bool bond_near_wat = false; // water hybrid sf
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	// Adjust hydrogen bonding potential to accomodate stronger hbonding in the
	// membrane hydrophobic core. Incorporating automatic deteciton for membrane poses (mpframework)
	if ( options_->Mbhbond() || options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, pose, bond_near_wat);

		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, pose, bond_near_wat);

	} else {

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, 1.0, bond_near_wat);

		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
			false /*calculate_derivative*/,
			false, true, true, true, /* calc bb_bb, don't calc bb_sc, sc_bb, or sc_sc */
			*options_,
			emap, 1.0, bond_near_wat);

	}
}

/// @details Note: this function enforces the bb/sc hbond exclusion rule
void
HBondEnergy::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;
	if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & > ( pose.energies().data().get( HBOND_SET )));

	// this only works because we have already called
	// hbond_set->setup_for_residue_pair_energies( pose )

	// hydrate/SPaDES protocol
	bool bond_near_wat = false; // water hybrid sf
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	// Adjust hydrogen bonding potential to accomodate stronger hbonding in the
	// membrane hydrophobic core. Incorporating automatic deteciton for membrane poses (mpframework)
	if ( options_->Mbhbond() || options_->mphbond() ) {

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-donor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way_membrane function
		if ( !options_->bb_donor_acceptor_check() || ! rsd1.is_protein() || ! hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos()) ) {
			identify_hbonds_1way_membrane(
				*database_,
				rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
				false /*calculate_derivative*/,
				true, true, false, true,
				*options_,
				emap, pose, bond_near_wat);
		}

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-acceptor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way_membrane function
		if ( !options_->bb_donor_acceptor_check() || ! rsd1.is_protein() || ! hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos()) ) {
			identify_hbonds_1way_membrane(
				*database_,
				rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
				false /*calculate_derivative*/,
				true, false, true, true,
				*options_,
				emap, pose, bond_near_wat);
		}

	} else {
		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-donor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way function
		if ( !options_->bb_donor_acceptor_check() || !rsd1.is_protein() || !hbond_set.don_bbg_in_bb_bb_hbond(rsd1.seqpos()) ) {
			identify_hbonds_1way(
				*database_,
				rsd1, rsd2, hbond_set.nbrs(rsd1.seqpos()), hbond_set.nbrs(rsd2.seqpos()),
				false /*calculate_derivative*/,
				true, true, false, true,
				*options_,
				emap, 1.0, bond_near_wat);
		}

		/// If we're enforcing the bb/sc exclusion rule, and residue1 is a protein residue, and if residue 1's backbone-acceptor group
		/// is already participating in a bb/bb hbond, then do not evaluate the identify_hbonds_1way function
		if ( !options_->bb_donor_acceptor_check() || !rsd1.is_protein() || !hbond_set.acc_bbg_in_bb_bb_hbond(rsd1.seqpos()) ) {
			identify_hbonds_1way(
				*database_,
				rsd2, rsd1, hbond_set.nbrs(rsd2.seqpos()), hbond_set.nbrs(rsd1.seqpos()),
				false /*calculate_derivative*/,
				true, false, true, true,
				*options_,
				emap, 1.0, bond_near_wat);
		}
	}
}

void
HBondEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap & emap
) const
{
	Size nbrs1 = pose.energies().tenA_neighbor_graph().get_node( rsd1.seqpos() )->num_neighbors_counting_self_static();
	Size nbrs2 = pose.energies().tenA_neighbor_graph().get_node( rsd2.seqpos() )->num_neighbors_counting_self_static();

	using EnergiesCacheableDataType::HBOND_SET; // jklai
	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ))); // jklai

	// hydrate/SPaDES protocol
	bool bond_near_wat = false; // water hybrid sf
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) bond_near_wat = true;
	}

	// Adjust hydrogen bonding potential to accomodate stronger hbonding in the
	// membrane hydrophobic core. Incorporating automatic deteciton for membrane poses (mpframework)
	if ( options_->Mbhbond() || options_->mphbond() ) {

		identify_hbonds_1way_membrane(
			*database_,
			rsd1, rsd2, nbrs1, nbrs2,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap,
			pose, bond_near_wat);
		identify_hbonds_1way_membrane(
			*database_,
			rsd2, rsd1, nbrs2, nbrs1,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap,
			pose, bond_near_wat);

	} else {

		identify_hbonds_1way(
			*database_,
			rsd1, rsd2, nbrs1, nbrs2,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap, 1.0, bond_near_wat);
		identify_hbonds_1way(
			*database_,
			rsd2, rsd1, nbrs2, nbrs1,
			false /*calculate_derivative*/,
			true, true, true, false, *options_, emap, 1.0, bond_near_wat);
	}
}


void
HBondEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	debug_assert( set1.resid() != set2.resid() );

	if ( options_->exclude_DNA_DNA() && pose.residue( set1.resid() ).is_DNA() && pose.residue( set2.resid() ).is_DNA() ) return;

	using namespace methods;
	using namespace hbtrie;
	using namespace trie;
	using EnergiesCacheableDataType::HBOND_SET;

	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table1( energy_table, 0 ); // Filled copy
	ObjexxFCL::FArray2D< core::PackerEnergy > temp_table2( energy_table, 0 );

	HBondsTrieVsTrieCachedDataContainer container( weights ); //Note that weights MUST persist until container is destroyed.  It does not own weights, but stores a raw pointer to them!

	// save weight information so that its available during tvt execution
	// and also the neighbor counts for the two residues.
	container.set_res1( set1.resid() );
	container.set_res2( set2.resid() );

	container.set_rotamer_seq_sep( pose.residue( set2.resid() ).polymeric_oriented_sequence_distance( pose.residue( set1.resid() ) ) );

	if ( true ) { // super_hacky
		auto const & hbond_set
			( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( HBOND_SET )));
		container.set_res1_nb( hbond_set.nbrs( set1.resid() ) );
		container.set_res2_nb( hbond_set.nbrs( set2.resid() ) );
	} else {
		container.set_res1_nb( pose.energies().tenA_neighbor_graph().get_node( set1.resid() )->num_neighbors_counting_self() );
		container.set_res2_nb( pose.energies().tenA_neighbor_graph().get_node( set2.resid() )->num_neighbors_counting_self() );
	}
	HBondRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set1.get_trie( hbond_method ) ));
	HBondRotamerTrieCOP trie2( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set2.get_trie( hbond_method ) ));

	// figure out which trie countPairFunction needs to be used for this set
	TrieCountPairBaseOP cp( new HBCountPairFunction );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	//std::cout << "BEGIN TVT " << set1.resid() << " and " << set2.resid() << std::endl;
	trie1->trie_vs_trie( *trie2, *cp, *this, temp_table1, temp_table2, &container );
	//std::cout << "END TVT " << set1.resid() << " and " << set2.resid() << std::endl;

	/// add in the energies calculated by the tvt alg.
	energy_table += temp_table1;
}


//@brief overrides default rotamer/background energy calculation and uses
// the trie-vs-trie algorithm instead
void
HBondEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	using namespace methods;
	using namespace hbtrie;
	using namespace trie;
	using EnergiesCacheableDataType::HBOND_SET;
	using EnergiesCacheableDataType::HBOND_TRIE_COLLECTION;

	if ( options_->exclude_DNA_DNA() && residue.is_DNA() && pose.residue( set.resid() ).is_DNA() ) return;

	// allocate space for the trie-vs-trie algorithm
	utility::vector1< core::PackerEnergy > temp_vector1( set.num_rotamers(), 0.0 );
	utility::vector1< core::PackerEnergy > temp_vector2( set.num_rotamers(), 0.0 );

	// save weight information so that its available during tvt execution
	HBondsTrieVsTrieCachedDataContainer container(weights); //Note that weights MUST persist until container is destroyed.  It does not own weights, but stores a raw pointer to them!
	container.set_res1( set.resid() );      //Added by Parin Sripakdeevong (sripakpa@stanford.edu)
	container.set_res2( residue.seqpos() ); //Added by Parin Sripakdeevong (sripakpa@stanford.edu)
	container.set_rotamer_seq_sep( pose.residue( residue.seqpos() ).polymeric_oriented_sequence_distance( pose.residue( set.resid() ) ) );

	if ( true ) { // super_hacky
		auto const & hbond_set
			( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( HBOND_SET )));
		container.set_res1_nb( hbond_set.nbrs( set.resid() ) );
		container.set_res2_nb( hbond_set.nbrs( residue.seqpos() ) );
	} else {
		container.set_res1_nb( pose.energies().tenA_neighbor_graph().get_node( set.resid() )->num_neighbors_counting_self() );
		container.set_res2_nb( pose.energies().tenA_neighbor_graph().get_node( residue.seqpos() )->num_neighbors_counting_self() );
	}

	HBondRotamerTrieCOP trie1( utility::pointer::static_pointer_cast< trie::RotamerTrieBase const > ( set.get_trie( hbond_method ) ) );
	HBondRotamerTrieCOP trie2 = ( static_cast< TrieCollection const & >
		( pose.energies().data().get( HBOND_TRIE_COLLECTION )) ).trie( residue.seqpos() );

	//fpd
	if ( trie2 == nullptr ) return;

	TrieCountPairBaseOP cp( new HBCountPairFunction );

	/// now execute the trie vs trie algorithm.
	/// this steps through three rounds of type resolution before finally arriving at the
	/// actual trie_vs_trie method.  The type resolution calls allow the trie-vs-trie algorithm
	/// to be templated with full type knowledge (and therefore be optimized by the compiler for
	/// each variation on the count pair data used and the count pair funtions invoked.
	trie1->trie_vs_path( *trie2, *cp, *this, temp_vector1, temp_vector2, &container );

	/// add in the energies calculated by the tvt alg.
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		energy_vector[ ii ] += temp_vector1[ ii ];
	}
}


void
HBondEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	/// Don't add in bb/bb hbond energies during minimization
	if ( pose.energies().use_nblist() ) return;

	if ( options_->decompose_bb_hb_into_pair_energies() ) return;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >( pose.energies().data().get( HBOND_SET )));

	//EnergyMap hbond_emap(totals);

	// the current logic is that we fill the hbond set with backbone
	// hbonds only at the beginning of scoring. this is done to setup
	// the bb-bb hbond exclusion logic. so the hbondset should only
	// include bb-bb hbonds.
	// but see get_hb_don_chem_type in hbonds_geom.cc -- that only
	// classifies protein backbone donors as backbone, and the energy
	// accumulation by type is influenced by that via HBeval_lookup
	//
	// the important thing is that there's no double counting, which
	// is I think true since both fill_hbond_set and get_rsd-rsd-energy
	// use atom_is_backbone to check...
	//debug_assert( std::abs( bb_scE ) < 1e-3 && std::abs( scE ) < 1e-3 );


	// this is to replicate buggy behavior regarding protein-backbone -- dna-backbone hbonds
	Real original_bb_sc    = totals[ hbond_bb_sc ];
	Real original_sr_bb_sc = totals[ hbond_sr_bb_sc ];
	Real original_lr_bb_sc = totals[ hbond_lr_bb_sc ];
	Real original_sc       = totals[ hbond_sc ];
	Real original_wat      = totals[ hbond_wat ]; // hydrate/SPaDES scoring function
	Real original_ent      = totals[ wat_entropy ]; // hydrate/SPaDES scoring function
	Real original_intra    = totals[ hbond_intra ];
	// end replicate

	get_hbond_energies( hbond_set, totals );

	// begin replicate
	totals[ hbond_bb_sc ]    = original_bb_sc;
	totals[ hbond_sr_bb_sc ] = original_sr_bb_sc;
	totals[ hbond_lr_bb_sc ] = original_lr_bb_sc;
	totals[ hbond_sc ]       = original_sc;
	totals[ hbond_wat ]      = original_wat; // hydrate/SPaDES scoring function
	totals[ wat_entropy ]    = original_ent; // hydrate/SPaDES scoring function
	totals[ hbond_intra ]    = original_intra;
	// end replicate

}

/// @brief HACK!  MAX_R defines the maximum donorH to acceptor distance.
// The atomic_interaction_cutoff method is meant to return the maximum distance
// between two *heavy atoms* for them to have a zero interaction energy.
// I am currently assuming a 1.35 A maximum distance between a hydrogen and the
// heavy atom it is bound to, stealing this number from the CYS.params file since
// the HG in CYS is much further from it's SG than aliphatic hydrogens are from their carbons.
// This is a bad idea.  Someone come up with a way to fix this!
//
// At 4.35 A interaction cutoff, the hbond energy function is incredibly short ranged!
Distance
HBondEnergy::atomic_interaction_cutoff() const
{
	return MAX_R + 1.35; // MAGIC NUMBER
}

/// @brief the atomic interaction cutoff and the hydrogen interaction cutoff are the same.
Real
HBondEnergy::hydrogen_interaction_cutoff2() const
{
	return (MAX_R + 1.35) * ( MAX_R + 1.35 );
}


/// @brief HBondEnergy is context sensitive
void
HBondEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

bool
HBondEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	return ( weights[hbond_intra] > 0.0 || weights[hbond] > 0.0 );
}


void
HBondEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	if ( calculate_intra_res_hbonds( rsd, *options_ ) ) {
		TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
		Size const rsd_nb = tenA_neighbor_graph.get_node( rsd.seqpos() )->num_neighbors_counting_self_static();
		identify_intra_res_hbonds( *database_, rsd, rsd_nb, *options_, emap);
	}
}

void
create_rotamer_descriptor(
	conformation::Residue const & res,
	hbonds::HBondOptions const & options,
	hbonds::HBondSet const & hbond_set,
	trie::RotamerDescriptor< hbtrie::HBAtom, hbtrie::HBCPData > & rotamer_descriptor,
	bool near_wat = false,
	bool is_wat = false
)
{
	using namespace trie;
	using namespace hbtrie;

	Size const resid( res.seqpos() );

	Size n_to_add(0);
	utility::vector1< int > add_to_trie( res.natoms(), 0 );
	for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
		if ( res.atom_type_set()[ res.atom( jj ).type() ].is_acceptor() ) {
			//std::cout << "acc=" << jj << " ";
			add_to_trie[ jj ] = 1;
		} else if ( res.atom_is_hydrogen( jj ) && //else if ( res.atom_type_set()[ res.atom( jj ).type() ].is_hydrogen() &&
				res.atom_type_set()[ res.atom( res.type().atom_base( jj ) ).type() ].is_donor() ) {
			add_to_trie[ jj ] = 1;
			add_to_trie[ res.type().atom_base( jj ) ] = 1;
			//std::cout << "don=" << jj << " donb= " << res.type().atom_base( jj ) << " ";
		}
	}
	//std::cout << "rotamer trie for residue: " << res.type().name() << std::endl;
	for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
		if ( add_to_trie[ jj ] == 1 ) {
			++n_to_add;
			//std::cout << jj << " ";
		}
	}
	//std::cout << std::endl;

	if ( n_to_add == 0 ) {
		// What happens if there are NO hydrogen bonding atoms?  It would be nice if we didn't have to create
		// a trie at all, but I'm pretty sure the indexing logic requires that we have a place-holder rotamer.
		add_to_trie[ 1 ] = 1; ++n_to_add;
	}
	rotamer_descriptor.natoms( n_to_add );

	Size count_added_atoms = 0;
	for ( Size jj = 1; jj <= res.nheavyatoms(); ++jj ) {
		if ( add_to_trie[ jj ] == 0 ) continue;

		HBAtom newatom;
		HBCPData cpdata;

		newatom.xyz( res.atom(jj).xyz() );
		newatom.base_xyz( res.xyz( res.atom_base( jj )) );
		newatom.is_hydrogen( false );
		newatom.is_backbone( res.atom_is_backbone( jj ) ) ;

		//Following preserves hbond_sc, hbond_bb_sc, as preferred by other developers for protein/DNA.
		newatom.is_protein( res.is_protein() );
		newatom.is_dna( res.is_DNA() );

		// hydrate/SPaDES protocol
		newatom.near_wat( near_wat );
		newatom.is_wat( is_wat );

		if ( res.atom_type_set()[ res.atom( jj ).type() ].is_acceptor() ) {

			newatom.hb_chem_type( get_hb_acc_chem_type( jj, res ));
			newatom.base2_xyz( res.xyz( res.abase2( jj )) );

			cpdata.is_sc( ! res.type().atom_is_backbone( jj ) );

			/// Count-pair data is responsible for enforcing the sc/bb hbond exclusion rule.
			/// If we're not using the rule, set "avoid_sc_hbonds" to false.
			cpdata.avoid_sc_hbonds( options.bb_donor_acceptor_check() &&
				! cpdata.is_sc() &&
				res.type().is_protein() &&
				hbond_set.acc_bbg_in_bb_bb_hbond( resid ) );
		}


		RotamerDescriptorAtom< HBAtom, HBCPData > rdatom( newatom, cpdata );
		rotamer_descriptor.atom( ++count_added_atoms, rdatom );

		for ( Size kk = res.attached_H_begin( jj ),
				kk_end = res.attached_H_end( jj );
				kk <= kk_end; ++kk ) {
			if ( add_to_trie[ kk ] == 0 ) continue;

			HBAtom newhatom;
			newhatom.xyz( res.atom(kk).xyz() );
			newhatom.base_xyz( res.xyz( res.atom_base( kk )) );
			newhatom.base2_xyz( Vector( 0.0, 0.0, 0.0 ) );
			newhatom.hb_chem_type( get_hb_don_chem_type( res.atom_base(kk), res ));
			newhatom.is_hydrogen( true );
			newhatom.is_backbone( res.atom_is_backbone( kk ) ) ;

			//Following preserves hbond_sc, hbond_bb_sc, as preferred by other developers for protein/DNA.
			newhatom.is_protein( res.is_protein() );
			newhatom.is_dna( res.is_DNA() );

			// hydrate/SPaDES protocol
			newhatom.near_wat( near_wat );
			newhatom.is_wat( is_wat );

			HBCPData hcpdata;
			hcpdata.is_sc( ! res.type().atom_is_backbone( kk ) );

			/// Count-pair data is responsible for enforcing the sc/bb hbond exclusion rule.
			/// If we're not using the rule, set "avoid_sc_hbonds" to false.
			hcpdata.avoid_sc_hbonds( options.bb_donor_acceptor_check() &&
				! hcpdata.is_sc() &&
				res.type().is_protein() &&
				hbond_set.don_bbg_in_bb_bb_hbond( resid ) );

			RotamerDescriptorAtom< HBAtom, HBCPData > hrdatom( newhatom, hcpdata );
			rotamer_descriptor.atom( ++count_added_atoms, hrdatom );

		}
	}
}


void
HBondEnergy::atomistic_pair_energy(
	core::Size atm1, // Which atom in residue 1?
	conformation::Residue const & rsd1, // Residue 1
	core::Size atm2, // Which atom in residue 2?
	conformation::Residue const & rsd2, // Residue 2
	pose::Pose const & pose, // pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		// Residue internal energy: a crib of eval_intrares_energy, expanded to account for the atomistic evaluation.
		if ( ! calculate_intra_res_hbonds( rsd1, *options_ ) ) { return; }

		TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
		Size const rsd_nb = tenA_neighbor_graph.get_node( rsd1.seqpos() )->num_neighbors_counting_self_static();

		HBondSet hbond_set( *options_ );
		identify_intra_res_hbonds( *database_, rsd1, rsd_nb, false /*evaluate_derivative*/, hbond_set );

		for ( HBondOP const & hb: hbond_set.hbonds() ) {
			debug_assert( hb->don_res() == rsd1.seqpos() ); // Should be, as we're intrares.
			debug_assert( hb->acc_res() == rsd1.seqpos() );
			// Skip out if we have an hbond not between one of the atoms involved.
			if ( hb->acc_atm() != atm1 && hb->don_hatm() != atm1 ) { continue; }
			if ( hb->acc_atm() != atm2 && hb->don_hatm() != atm2 ) { continue; }

			Real const weighted_energy = hb->energy() * hb->weight();
			if ( options_->put_intra_into_total() ) {
				emap[ hbond ]       += weighted_energy;
			} else {
				emap[ hbond_intra ] += weighted_energy;
			}
		}
	} else {
		// Residue pair energy: a crib of residue_pair_energy(), expanded to account for the atomistic evaluation.
		if ( options_->exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) return;

		bool rsd1_is_donor = false;

		if ( rsd1.Hpos_polar().has_value( atm1 ) ) {
			if ( rsd2.accpt_pos().has_value( atm2 ) ) {
				rsd1_is_donor = true;
			} else { return; } // Not compatible
		} else if ( rsd1.accpt_pos().has_value( atm1 ) ) {
			if ( rsd2.Hpos_polar().has_value( atm2 ) ) {
				rsd1_is_donor = false;
			} else { return; } // Not compatible
		} else { return; } // Not an hbonder

		conformation::Residue const & don_rsd = rsd1_is_donor ? rsd1 : rsd2;
		conformation::Residue const & acc_rsd = rsd1_is_donor ? rsd2 : rsd1;
		Size const hatm = rsd1_is_donor ? atm1 : atm2;
		Size const aatm = rsd1_is_donor ? atm2 : atm1;

		auto const & hbond_set
			( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));

		bool const exclude_bb = !options_->decompose_bb_hb_into_pair_energies();
		bool const exclude_sc = false;
		bool exclude_bsc = false, exclude_scb = false;
		if ( don_rsd.is_protein() ) exclude_scb = options_->bb_donor_acceptor_check() && hbond_set.don_bbg_in_bb_bb_hbond(don_rsd.seqpos());
		if ( acc_rsd.is_protein() ) exclude_bsc = options_->bb_donor_acceptor_check() && hbond_set.acc_bbg_in_bb_bb_hbond(acc_rsd.seqpos());

		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);
		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) return;
		} else {
			if ( exclude_sc && exclude_bsc ) return;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		if ( acc_rsd.atom_is_backbone(aatm) ) {
			if ( datm_is_bb ) {
				if ( exclude_bb ) return;
			} else {
				if ( exclude_bsc ) return;
			}
		} else {
			if ( datm_is_bb ) {
				if ( exclude_scb ) return;
			} else {
				if ( exclude_sc ) return;
			}
		}

		// rough filter for existence of hydrogen bond
		if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) return;

		HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

		int const base ( acc_rsd.atom_base( aatm ) );
		int const base2( acc_rsd.abase2( aatm ) );
		debug_assert( base2 > 0 && base != base2 );

		Real unweighted_energy( 0.0 );

		hb_energy_deriv( *database_, *options_, hbe_type, datm_xyz, hatm_xyz,
			acc_rsd.atom(aatm ).xyz(),
			acc_rsd.atom(base ).xyz(),
			acc_rsd.atom(base2).xyz(),
			unweighted_energy, false, DUMMY_DERIVS);

		if ( unweighted_energy >= options_->max_hb_energy() ) return;

		core::Size const don_nb = hbond_set.nbrs(don_rsd.seqpos());
		core::Size const acc_nb = hbond_set.nbrs(rsd2.seqpos());

		core::Real environmental_weight;
		if ( options_->Mbhbond() || options_->mphbond() ) {
			environmental_weight =
				get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, don_rsd.seqpos(), acc_rsd.seqpos(), hatm, aatm, hatm_xyz,
				acc_rsd.atom(aatm ).xyz());

			// hydrate/SPaDES protocol for when bond is near water
			if ( options_->water_hybrid_sf() ) {
				if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) {
					environmental_weight = 1;
				}
			}
		} else {
			environmental_weight =
				(!options_->use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, *options_));

			// hydrate/SPaDES protocol for when bond is near water
			if ( options_->water_hybrid_sf() ) {
				if ( residue_near_water( pose, rsd1.seqpos() ) || residue_near_water( pose, rsd2.seqpos() ) ) {
					environmental_weight = 1;
				}
			}

			if ( get_hbond_weight_type(hbe_type.eval_type()) == hbw_SR_BB ) {
				SSWeightParameters ssdep;
				ssdep.ssdep_ = options_->length_dependent_srbb();
				ssdep.l_ = options_->length_dependent_srbb_lowscale();
				ssdep.h_ = options_->length_dependent_srbb_highscale();
				ssdep.len_l_ = options_->length_dependent_srbb_minlength();
				ssdep.len_h_ = options_->length_dependent_srbb_maxlength();
				Real ssdep_weight_factor = get_ssdep_weight(rsd1, rsd2, pose, ssdep);

				environmental_weight *= ssdep_weight_factor;
			}
		}

		Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;

		// hydrate/SPaDES protocol scoring function
		if ( options_->water_hybrid_sf() ) {
			if ( ( don_rsd.name() == "TP3" && acc_rsd.name() != "TP3") || ( acc_rsd.name() == "TP3" && don_rsd.name() != "TP3" ) ) {
				static core::scoring::func::FuncOP smoothed_step ( new core::scoring::func::SmoothStepFunc( -0.55,-0.45 ) );
				emap[ wat_entropy ] += 1.0 - smoothed_step->func( unweighted_energy );
			}
			if ( (don_rsd.name() == "TP3" || acc_rsd.name() == "TP3") ) {
				emap[ hbond_wat ] += hbE;
				return;
			}
		}

		increment_hbond_energy( hbe_type.eval_type(), emap, hbE );
	}
}


hbtrie::HBondRotamerTrieOP
HBondEnergy::create_rotamer_trie(
	conformation::RotamerSetBase const & rotset,
	pose::Pose const & pose
) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_SET;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ) ) );

	utility::vector1< RotamerDescriptor< HBAtom, HBCPData > > rotamer_descriptors( rotset.num_rotamers() );

	// hydrate/SPaDES protocol
	bool near_wat = false; // hybrid dependency on hb env weight
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, rotset.resid() ) ) {
			near_wat = true;
		}
	}

	// hydrate/SPaDES protocol
	bool is_wat = false; // hybrid water specific scoring
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( rotset.rotamer( 1 )->name() == "TP3" ) {
			is_wat = true;
		}
	}

	for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
		conformation::ResidueCOP ii_rotamer( rotset.rotamer( ii ) );
		create_rotamer_descriptor( *ii_rotamer, *options_, hbond_set, rotamer_descriptors[ ii ] , near_wat, is_wat );
		rotamer_descriptors[ ii ].rotamer_id( ii );
	}

	sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );

	return utility::pointer::make_shared< RotamerTrie< HBAtom, HBCPData > >( rotamer_descriptors, atomic_interaction_cutoff());
}

hbtrie::HBondRotamerTrieOP
HBondEnergy::create_rotamer_trie(
	conformation::Residue const & res,
	pose::Pose const & pose
) const
{
	using namespace trie;
	using namespace hbtrie;
	using EnergiesCacheableDataType::HBOND_SET;

	auto const & hbond_set
		( static_cast< hbonds::HBondSet const & >
		( pose.energies().data().get( HBOND_SET ) ) );

	utility::vector1< RotamerDescriptor< HBAtom, HBCPData > > rotamer_descriptors( 1 );

	// hydrate/SPaDES protocol
	bool near_wat = false; // hybrid dependency on hb env weight
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( residue_near_water( pose, res.seqpos() ) ) {
			near_wat = true;
		}
	}

	// hydrate/SPaDES protocol
	bool is_wat = false; // hybrid water specific scoring
	if ( hbond_set.hbond_options().water_hybrid_sf() ) {
		if ( res.name() == "TP3" ) {
			is_wat = true;
		}
	}

	create_rotamer_descriptor( res, *options_, hbond_set, rotamer_descriptors[ 1 ], near_wat, is_wat );
	rotamer_descriptors[ 1 ].rotamer_id( 1 );

	return utility::pointer::make_shared< RotamerTrie< HBAtom, HBCPData > >( rotamer_descriptors, atomic_interaction_cutoff());
}


/// @brief code to evaluate a hydrogen bond energy for the trie that
/// didn't belong in the header itself -- it certainly does enough work
/// such that inlining it would not likely produce a speedup.
Energy
HBondEnergy::drawn_out_heavyatom_hydrogenatom_energy(
	hbtrie::HBAtom const & at1, // atom 1 is the heavy atom, the acceptor
	hbtrie::HBAtom const & at2, // atom 2 is the hydrogen atom, the donor
	bool flipped, // is at1 from residue 1?
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
) const
{
#ifdef NDEBUG
	core::scoring::hbonds::hbtrie::HBondsTrieVsTrieCachedDataContainer const * const container( static_cast< core::scoring::hbonds::hbtrie::HBondsTrieVsTrieCachedDataContainer const * >(cached_data) );
#else
	core::scoring::hbonds::hbtrie::HBondsTrieVsTrieCachedDataContainer const * const container( dynamic_cast< core::scoring::hbonds::hbtrie::HBondsTrieVsTrieCachedDataContainer const * const >(cached_data) );
	debug_assert( container != nullptr );
#endif

	// When acc and don are both polymers and on the same chain:
	// ss = acc.seqpos - don.seqpos
	int ss = (flipped ? -1.0*container->rotamer_seq_sep() : container->rotamer_seq_sep());
	HBEvalTuple hbe_type = hbond_evaluation_type( at2, 0,   // donor atom
		at1, ss); // acceptor atom
	Energy hbenergy;

	//std::cout << "about to evaluate " << at1 << " " << at2 << std::endl;
	hb_energy_deriv( *database_, *options_, hbe_type,
		at2.base_xyz(), at2.xyz(), // donor heavy atom, donor hydrogen,
		at1.xyz(), at1.base_xyz(), at1.base2_xyz(), // acceptor, acceptor base, acceptor base2
		hbenergy, false, DUMMY_DERIVS );
	//std::cout << "drawn_out_heavyatom_hydrogenatom_energy " << hbenergy << " " << at1 << " " << at2 << std::endl;

	if ( hbenergy >= options_->max_hb_energy() ) return 0.0; // no hbond

	//std::cout << "drawn_out_heavyatom_hydrogenatom_energy " << hbenergy << " " << at1 << " " << at2 << std::endl;

	Real envweight( 1.0 );
	if ( options_->use_hb_env_dep() ) {
		envweight = ( flipped ? get_environment_dependent_weight( hbe_type, container->res2_nb(), container->res1_nb(), *options_ ) :
			get_environment_dependent_weight( hbe_type, container->res1_nb(), container->res2_nb(), *options_ ));
	}

	// hydrate/SPaDES protocol
	if ( options_->water_hybrid_sf() && ( at1.near_wat() || at2.near_wat() ) ) {
		envweight = 1.0;
	}

	//pba membrane specific correction
	if ( options_->Mbhbond() && options_->mphbond() ) {

		Real membrane_depth_dependent_weight( 1.0 );

		membrane_depth_dependent_weight = ( flipped ? get_membrane_depth_dependent_weight(normal_, center_, thickness_,
			steepness_, container->res2_nb(), container->res1_nb(), at2.xyz(), at1.xyz()) : get_membrane_depth_dependent_weight(normal_, center_, thickness_,
			steepness_, container->res2_nb(), container->res1_nb(), at2.xyz(), at1.xyz()) );

		envweight = membrane_depth_dependent_weight;
	}

	Real weighted_energy(0.0);

	// hydrate/SPaDES protocol
	if ( ( at1.is_wat() || at2.is_wat() ) && options_->water_hybrid_sf() ) {
		if ( container->res1() == container->res2() ) weighted_energy = container->weights()[ hbond_intra ] * hbenergy * envweight;
		else  weighted_energy = container->weights()[ hbond_wat ] * hbenergy * envweight;
	} else {
		weighted_energy = hb_eval_type_weight(hbe_type.eval_type(), container->weights(), container->res1()==container->res2()) * hbenergy * envweight;
	}

	// hydrate/SPaDES protocol
	if ( options_->water_hybrid_sf() ) {
		if ( ( at1.is_wat() && !at2.is_wat() ) || ( !at1.is_wat() && at2.is_wat() ) ) {
			static core::scoring::func::FuncOP smoothed_step ( new core::scoring::func::SmoothStepFunc( -0.55,-0.45 ) );
			weighted_energy += (1.0 - smoothed_step->func( hbenergy )) * container->weights()[ wat_entropy ];
		}
	}

	//std::cout << "weighted energy: " << weighted_energy << " " << hbenergy << " " << envweight << std::endl;

	return weighted_energy;
}


/// @details
/// Version 2: 2011-06-27 Fixing chi2 SER/THR and chi3 TYR derivatives when they act as acceptors.
core::Size
HBondEnergy::version() const
{
	//return 1; // Initial versioning
	// return 2;
	return 3; // fixing sp3-acceptor discontinuitiy w/ commit e42010b198a
}

} // hbonds
} // scoring
} // core

#ifdef    SERIALIZATION

typedef core::scoring::trie::RotamerTrie< core::scoring::hbonds::hbtrie::HBAtom, core::scoring::hbonds::hbtrie::HBCPData > HBRotTrie; SAVE_AND_LOAD_SERIALIZABLE( HBRotTrie );
CEREAL_REGISTER_TYPE( HBRotTrie )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_hbonds_HBondEnergy )
#endif
