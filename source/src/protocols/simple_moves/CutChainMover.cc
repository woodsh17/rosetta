// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ForceDisulfidesMover.cc
///
/// @brief
/// @author Sarel Fleishman

// unit headers
#include <protocols/simple_moves/CutChainMover.hh>
#include <protocols/simple_moves/CutChainMoverCreator.hh>

// type headers
#include <core/types.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

// package headers

// Project headers
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

// utility header
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.simple_moves.CutChainMover" );

std::string
CutChainMoverCreator::keyname() const
{
	return CutChainMoverCreator::mover_name();
}

protocols::moves::MoverOP
CutChainMoverCreator::create_mover() const {
	return utility::pointer::make_shared< CutChainMover >();
}

std::string
CutChainMoverCreator::mover_name()
{
	return "CutChain";
}

//Default Constructor
CutChainMover::CutChainMover() :
	protocols::moves::Mover("CutChain"),
	bond_length_(4.0), //define covalent bond length cut-off
	chain_id_(1)//define default main chain
{
}


CutChainMover::~CutChainMover() = default;

protocols::moves::MoverOP
CutChainMover::clone() const
{
	return utility::pointer::make_shared< CutChainMover >( *this );
}

protocols::moves::MoverOP
CutChainMover::fresh_instance() const
{
	return utility::pointer::make_shared< CutChainMover >();
}

//getters
std::string
CutChainMover::get_name() const {
	return "CutChain";
}

core::Real
CutChainMover::bond_length() const {
	return bond_length_;
}

core::Size
CutChainMover::chain_id() const {
	return chain_id_;
}

//setters
void
CutChainMover::bond_length(core::Real const length){
	bond_length_ = length;
}
void
CutChainMover::chain_id(core::Size const ID){
	chain_id_ = ID;
}
void
CutChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {
	bond_length( tag->getOption< core::Real >( "bond_length", 4.0 ) ); // EH!? Covalent bond length is more like 1.33 angstrom
	chain_id( tag->getOption< core::Size >( "chain_id", 1 ) );
	TR<<" bond_length: "<<bond_length();
	TR<<"Chain id: "<<chain_id();
	TR<<std::endl;
}

void CutChainMover::apply( core::pose::Pose & pose )
{
	create_subpose(pose);
	foldTree(pose);

}

core::Size
CutChainMover::chain_cut( core::pose::Pose & pose)
{
	core::Size cut_pos = 0;
	for ( core::Size resj = pose.conformation().chain_begin( chain_id_ ); resj <= pose.conformation().chain_end( chain_id_ )-1; ++resj ) {
		core::Real const distance = pose.residue( resj+1 ).xyz( "N" ).distance(pose.residue( resj ).xyz( "C" ));
		//TR<<"distance is: "<<distance<<std::endl;
		//TR<<"residue name is : "<<pose.residue(resj).name1()<<std::endl;
		if ( distance > bond_length() ) {
			cut_pos = resj;
			TR<<"Found cut at: "<<resj<<std::endl;
			break;
		}
	}
	return( cut_pos );//cut_pos is the amino acid position BEFORE the cut
}
core::Size
CutChainMover::chain_cut(core::pose::Pose & pose, core::Size start_res,core::Size end_res)
{
	core::Size cut_pos = 0;
	for ( core::Size resj = start_res; resj <end_res; ++resj ) {
		core::Real const distance = pose.residue( resj+1 ).xyz( "N" ).distance(pose.residue( resj ).xyz( "C" ));
		//   TR<<"distance is: "<<distance<<std::endl;
		//   TR<<"residue name is : "<<pose.residue(resj).name1()<<std::endl;
		if ( distance > bond_length() ) {
			cut_pos = resj;
			TR<<"Found cut at: "<<resj<<std::endl;
			break;
		}
	}
	return( cut_pos );//cut_pos is the amino acid position BEFORE the cut
}

void
CutChainMover::create_subpose(core::pose::Pose & pose)
{
	core::pose::Pose copy_pose( pose );
	pose.clear();

	//add to pose only residues from main chain
	for ( core::Size resj = copy_pose.conformation().chain_begin( chain_id_); resj <= copy_pose.conformation().chain_end( chain_id_); ++resj ) {
		core::conformation::Residue const & rsd( copy_pose.residue( resj) );
		pose.append_residue_by_bond( rsd );
	}

}

void
CutChainMover::foldTree (core::pose::Pose & pose){
	core::Size const s1 = chain_cut(pose);
	core::kinematics::FoldTree ft;
	ft.clear();
	ft.add_edge( 1, s1, -1 );
	ft.add_edge( s1, s1+1, 1 );
	ft.add_edge( s1+1, pose.conformation().chain_end( chain_id()), -1 );
	TR<<"old foldtree: "<<pose.fold_tree()<<std::endl;
	pose.fold_tree(ft);
	TR<<"new_foldtree: "<<pose.fold_tree()<<std::endl;
	pose.conformation().detect_disulfides();
}


} // simple_moves
} // protocols
