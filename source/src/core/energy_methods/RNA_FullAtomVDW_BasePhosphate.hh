// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RNA_FullAtomVDW_BasePhosphate.hh
/// @brief  FullAtom VDW energy between the base and phosphate group in the same (intra) nucleotide
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_RNA_FullAtomVDW_BasePhosphate_hh
#define INCLUDED_core_scoring_rna_RNA_FullAtomVDW_BasePhosphate_hh

// Unit headers
#include <core/energy_methods/RNA_FullAtomVDW_BasePhosphate.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//typedef utility::pointer::shared_ptr< core::scoring::etable::TableLookupEtableEnergy > TableLookupEtableEnergyOP;

namespace core {
namespace energy_methods {


class RNA_FullAtomVDW_BasePhosphate : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief ctor
	RNA_FullAtomVDW_BasePhosphate( core::scoring::methods::EnergyMethodOptions const & options );

	/// @brief dtor
	~RNA_FullAtomVDW_BasePhosphate() override;

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_fast_pair_energy_attached_H( //copy from atom_pair_energy_inline.hh
		core::conformation::Residue const & res1,
		int const atomno1,
		core::conformation::Residue const & res2,
		Size const atomno2,
		Size const at1hbegin, //at1hbegin and at1hend define a range of hydrogen atom indices -- those h's bound to at1
		Size const at1hend,
		Size const at2hbegin,
		Size const at2hend,
		core::scoring::EnergyMap & emap
	) const;


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		core::scoring::EnergyMap & emap
	) const;


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		core::scoring::ScoreFunction const & /*sfxn*/, // needed for non-nblist minimization
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

private:

	/// @brief RNA_FullAtomVDW_BasePhosphate( is context independent; indicates that no context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	core::Size version() const override;

private:

	//methods::EnergyMethodOptions const & options_;
	core::scoring::etable::EtableEvaluatorOP etable_evaluator_;

	/*


	///////////////////////////////////////////////////////////////////////////////


	*/

};

} //scoring
} //core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
