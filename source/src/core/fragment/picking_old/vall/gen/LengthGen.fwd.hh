// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/LengthGen.fwd.hh
/// @brief  forward declaration for core::fragment::picking_old::vall::gen::LengthGen
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_fwd_hh
#define INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_fwd_hh

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief forward declaration for core::fragment::picking_old::vall::gen::LengthGen
class LengthGen;


/// @brief LengthGen owning pointer
typedef utility::pointer::shared_ptr< LengthGen > LengthGenOP;


/// @brief LengthGen const owning pointer
typedef utility::pointer::shared_ptr< LengthGen const > LengthGenCOP;


/// @brief LengthGen access pointer
typedef utility::pointer::weak_ptr< LengthGen > LengthGenAP;


/// @brief LengthGen const access pointer
typedef utility::pointer::weak_ptr< LengthGen const > LengthGenCAP;


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

#endif /* INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_FWD_HH */
