// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/Abego.hh
/// @brief  ABEGO energy method class declaration
/// @author imv@uw.edu

#ifndef INCLUDED_core_energy_methods_Abego_hh
#define INCLUDED_core_energy_methods_Abego_hh

// Unit headers
#include <core/energy_methods/AbegoEnergy.fwd.hh>
#include <core/scoring/P_AA_ABEGO3.fwd.hh>
//#include <core/scoring/P_AA_ABEGO3.hh>

// Package headers
//#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/sequence/ABEGOManager.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class Abego : public core::scoring::methods::WholeStructureEnergy {
public:

	/// ctor
	Abego();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	// @brief Abego Energy is context independent and thus indicates that no context graphs need to
	// be maintained by class Energies
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const override;

	void setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void finalize_total_energy( pose::Pose &, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;

	// data
private:
	core::Size version() const override;
	core::scoring::P_AA_ABEGO3 const & paa_abego3_;
	mutable core::sequence::ABEGOManager abegoManager_;

	// @brief Sum of only the positive (unfavorable) ABEGO energies in the pose.
	mutable core::Real energy_positive_sum_;
	mutable core::Size energy_positive_sum_count_;

	// @brief Sum of all ABEGO energies across the pose.
	mutable core::Real energy_sum_;
	mutable core::Size energy_sum_count_;
};

} // scoring
} // core


#endif // INCLUDED_core_scoring_AbegoEnergy_HH
