// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/ncbb/oop/OopMover.hh
/// @brief
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_oop_OopMover_hh
#define INCLUDED_protocols_simple_moves_oop_OopMover_hh
// Unit Headers
#include <protocols/ncbb/oop/OopMover.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace oop {

/// @details
class OopMover : public protocols::moves::Mover {

public:

	/// @brief
	OopMover( core::Size oop_seq_position );
	OopMover( core::Size oop_seq_position, core::Real phi_angle, core::Real psi_angle );

	~OopMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	virtual void set_phi( core::Real angle ) { phi_angle_ = angle; }
	virtual void set_psi( core::Real angle ) { psi_angle_ = angle; }

	virtual void update_hydrogens( core::pose::Pose & pose ) { update_hydrogens_( pose ); }

private:

	core::Size const oop_pre_pos_;
	core::Size const oop_post_pos_;
	core::Real phi_angle_;
	core::Real psi_angle_;

private:

	void update_hydrogens_( core::pose::Pose & pose );

};//end OopMover


}//namespace oop
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_oop_OopMover_hh
