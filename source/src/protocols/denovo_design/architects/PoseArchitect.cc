// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/PoseArchitect.cc
/// @brief Design segments based on a pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitectCreator.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.denovo_design.architects.PoseArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

PoseArchitect::PoseArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	secstruct_( "" ),
	add_padding_( true )
{}

PoseArchitect::~PoseArchitect() = default;

DeNovoArchitectOP
PoseArchitect::clone() const
{
	return utility::pointer::make_shared< PoseArchitect >( *this );
}

std::string
PoseArchitect::type() const
{
	return architect_name();
}

void
PoseArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	PoseArchitect::provide_xml_schema( xsd );
}

void
PoseArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "add_padding", xsct_rosetta_bool, "Add padding to segments?" )
		+ XMLSchemaAttribute( "secstruct", xsct_dssp_string, "Desired secondary structure for the pose" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist );
	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );
	DeNovoArchitectFactory::xsd_architect_type_definition_w_attributes( xsd, architect_name(), "Design segments based on a pose", attlist );
}



void
PoseArchitect::parse_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	set_add_padding( tag->getOption< bool >( "add_padding", add_padding_ ) );
	set_secstruct( tag->getOption< std::string >( "secstruct", secstruct_ ) );
}

void
PoseArchitect::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

components::StructureDataOP
PoseArchitect::design( core::pose::Pose const & pose_in, core::Real & ) const
{
	using core::select::residue_selector::NotResidueSelector;
	using core::select::residue_selector::ResidueRanges;

	core::pose::PoseOP pose = pose_in.clone();
	if ( selector_ ) {
		// select regions to delete
		NotResidueSelector delete_selector;
		delete_selector.set_residue_selector( selector_ );
		ResidueRanges ranges( delete_selector.apply( *pose ) );
		for ( ResidueRanges::const_reverse_iterator it=ranges.rbegin(); it!=ranges.rend(); ++it ) {
			TR.Debug << "Deleting residue range " << it->start() << "-" << it->stop() << std::endl;
			pose->delete_residue_range_slow(it->start(), it->stop());
		}
	}

	components::StructureDataOP const sd( utility::pointer::make_shared< components::StructureData >( components::StructureDataFactory::get_instance()->create_from_pose( *pose, id() ) ) );
	if ( !secstruct_.empty() ) {
		if ( sd->pose_length() != secstruct_.size() ) {
			std::stringstream msg;
			msg << "PoseArchitect::design(): Length of user-provided secstruct (" << secstruct_.size()
				<< ") does not match input pose length (" << pose->size() << ")" << std::endl;
			utility_exit_with_message( msg.str() );
		}
		core::Size resid = 1;
		for ( std::string::const_iterator ss=secstruct_.begin(); ss!=secstruct_.end(); ++ss, ++resid ) {
			sd->set_ss( resid, *ss );
		}
	}

	// add template segments
	for ( auto s=sd->segments_begin(); s!=sd->segments_end(); ++s ) {
		core::Size const start = sd->segment( *s ).start();
		core::Size const stop = sd->segment( *s ).stop();
		sd->set_template_pose( *s, *pose, start, stop );
		if ( !add_padding_ ) {
			components::Segment seg = sd->segment( *s );
			seg.delete_lower_padding();
			seg.delete_upper_padding();
			sd->replace_segment( *s, seg );
		}
	}

	return sd;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
PoseArchitectCreator::keyname() const
{
	return PoseArchitect::architect_name();
}

DeNovoArchitectOP
PoseArchitectCreator::create_architect( std::string const & id ) const
{
	return utility::pointer::make_shared< PoseArchitect >( id );
}

} //architects
} //denovo_design
} //protocols

