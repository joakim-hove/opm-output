/*
  Copyright (c) 2016 Statoil ASA
  Copyright (c) 2013-2015 Andreas Lauser
  Copyright (c) 2013 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 Uni Research AS
  Copyright (c) 2015 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <unordered_map>

#include "config.h"

#include "EclipseIO.hpp"

#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Units/Dimension.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/EclipseState/Eclipse3DProperties.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/Utility/Functional.hpp>
#include <opm/output/eclipse/Summary.hpp>
#include <opm/output/eclipse/Tables.hpp>
#include <opm/output/eclipse/RegionCache.hpp>
#include <opm/output/eclipse/RestartIO.hpp>

#include <cstdlib>
#include <memory>     // unique_ptr
#include <utility>    // move

#include <ert/ecl/EclKW.hpp>
#include <ert/ecl/FortIO.hpp>
#include <ert/ecl/EclFilename.hpp>

#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_rft_file.h>
#include <ert/ecl/ecl_rst_file.h>
#include <ert/ecl_well/well_const.h>
#include <ert/ecl/ecl_rsthead.h>
#include <ert/util/util.h>
#define OPM_XWEL      "OPM_XWEL"
#define OPM_IWEL      "OPM_IWEL"

// namespace start here since we don't want the ERT headers in it
namespace Opm {
namespace {




void writeKeyword( ERT::FortIO& fortio ,
                   const std::string& keywordName,
                   const std::vector<int> &data ) {
    ERT::EclKW< int > kw( keywordName, data );
    kw.fwrite( fortio );
}

/*
  This overload hardcodes the common assumption that properties which
  are stored internally as double values in OPM should be stored as
  float values in the ECLIPSE formatted binary files.
*/

void writeKeyword( ERT::FortIO& fortio ,
                   const std::string& keywordName,
                   const std::vector<double> &data) {

    ERT::EclKW< float > kw( keywordName, data );
    kw.fwrite( fortio );

}


using restart_file = ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >;

restart_file open_rst( const std::string& filename,
                       bool truncate,
                       bool unifout,
                       int report_step ) {

    if( !unifout )
        return restart_file{ ecl_rst_file_open_write( filename.c_str() ) };

    if( truncate )
        return restart_file{ ecl_rst_file_open_write_seek( filename.c_str(), report_step ) };

    return restart_file{ ecl_rst_file_open_append( filename.c_str() ) };
}



/// Convert OPM phase usage to ERT bitmask
inline int ertPhaseMask( const Phases& phase ) {
    return ( phase.active( Phase::WATER ) ? ECL_WATER_PHASE : 0 )
         | ( phase.active( Phase::OIL ) ? ECL_OIL_PHASE : 0 )
         | ( phase.active( Phase::GAS ) ? ECL_GAS_PHASE : 0 );
}

class RFT {
    public:
    RFT( const std::string&  output_dir,
         const std::string&  basename,
         bool format );

        void writeTimeStep( std::vector< const Well* >,
                            const EclipseGrid& grid,
                            int report_step,
                            time_t current_time,
                            double days,
                            ert_ecl_unit_enum,
                            const data::Solution& cells);
    private:
        std::string filename;
};


    RFT::RFT( const std::string& output_dir,
              const std::string& basename,
              bool format ) :
        filename( ERT::EclFilename( output_dir, basename, ECL_RFT_FILE, format ) )
{}

inline ert_ecl_unit_enum to_ert_unit( UnitSystem::UnitType t ) {
    using ut = UnitSystem::UnitType;
    switch ( t ) {
        case ut::UNIT_TYPE_METRIC: return ECL_METRIC_UNITS;
        case ut::UNIT_TYPE_FIELD:  return ECL_FIELD_UNITS;
        case ut::UNIT_TYPE_LAB:    return ECL_LAB_UNITS;
    }

    throw std::invalid_argument("unhandled enum value");
}

void RFT::writeTimeStep( std::vector< const Well* > wells,
                         const EclipseGrid& grid,
                         int report_step,
                         time_t current_time,
                         double days,
                         ert_ecl_unit_enum unitsystem,
                         const data::Solution& cells) {

    using rft = ERT::ert_unique_ptr< ecl_rft_node_type, ecl_rft_node_free >;
    const std::vector<double>& pressure = cells.data("PRESSURE");
    const std::vector<double>& swat = cells.data("SWAT");
    const std::vector<double>& sgas = cells.has("SGAS") ? cells.data("SGAS") : std::vector<double>( pressure.size() , 0 );
    ERT::FortIO fortio(filename, std::ios_base::out);

    for ( const auto& well : wells ) {
        if( !( well->getRFTActive( report_step )
            || well->getPLTActive( report_step ) ) )
            continue;

        auto* rft_node = ecl_rft_node_alloc_new( well->name().c_str(), "RFT",
                current_time, days );

        for( const auto& completion : well->getCompletions( report_step ) ) {
            const size_t i = size_t( completion.getI() );
            const size_t j = size_t( completion.getJ() );
            const size_t k = size_t( completion.getK() );

            if( !grid.cellActive( i, j, k ) ) continue;

            const auto index = grid.activeIndex( i, j, k );
            const double depth = grid.getCellDepth( i, j, k );
            const double press = pressure[ index ];
            const double satwat = swat[ index ];
            const double satgas = sgas[ index ];

            auto* cell = ecl_rft_cell_alloc_RFT(
                            i, j, k, depth, press, satwat, satgas );

            ecl_rft_node_append_cell( rft_node, cell );
        }

        rft ecl_node( rft_node );
        ecl_rft_node_fwrite( ecl_node.get(), fortio.get(), unitsystem );
    }

    fortio.close();
}

inline std::string uppercase( std::string x ) {
    std::transform( x.begin(), x.end(), x.begin(),
        []( char c ) { return std::toupper( c ); } );

    return x;
}

}

class EclipseIO::Impl {
    public:
        Impl( const EclipseState& es, EclipseGrid grid );
        void writeINITFile( const data::Solution& simProps, const NNC& nnc) const;
        void writeEGRIDFile( const NNC& nnc ) const;

        const EclipseState& es;
        EclipseGrid grid;
        out::RegionCache regionCache;
        std::string outputDir;
        std::string baseName;
        out::Summary summary;
        RFT rft;
        time_t sim_start_time;
        std::array< int, 3 > cartesianSize;
        bool output_enabled;
        int ert_phase_mask;
        bool first_restart = true;
};

EclipseIO::Impl::Impl( const EclipseState& eclipseState,
                           EclipseGrid grid_)
    : es( eclipseState )
    , grid( std::move( grid_ ) )
    , regionCache( es , grid )
    , outputDir( eclipseState.getIOConfig().getOutputDir() )
    , baseName( uppercase( eclipseState.getIOConfig().getBaseName() ) )
    , summary( eclipseState, eclipseState.getSummaryConfig() , grid )
    , rft( outputDir.c_str(), baseName.c_str(), es.getIOConfig().getFMTOUT() )
    , sim_start_time( es.getSchedule().posixStartTime() )
    , output_enabled( eclipseState.getIOConfig().getOutputEnabled() )
    , ert_phase_mask( ertPhaseMask( eclipseState.runspec().phases() ) )
{}



void EclipseIO::Impl::writeINITFile( const data::Solution& simProps, const NNC& nnc) const {
    const auto& units = this->es.getUnits();
    const IOConfig& ioConfig = this->es.cfg().io();

    std::string  initFile( ERT::EclFilename( this->outputDir,
                                             this->baseName,
                                             ECL_INIT_FILE,
                                             ioConfig.getFMTOUT() ));

    ERT::FortIO fortio( initFile,
                        std::ios_base::out,
                        ioConfig.getFMTOUT(),
                        ECL_ENDIAN_FLIP );


    // Write INIT header. Observe that the PORV vector is treated
    // specially; that is because for this particulat vector we write
    // a total of nx*ny*nz values, where the PORV vector has been
    // explicitly set to zero for inactive cells. The convention is
    // that the active/inactive cell mapping can be inferred by
    // reading the PORV vector.
    {

        const auto& opm_data = this->es.get3DProperties().getDoubleGridProperty("PORV").getData();
        auto ecl_data = opm_data;

        for (size_t global_index = 0; global_index < opm_data.size(); global_index++)
            if (!this->grid.cellActive( global_index ))
                ecl_data[global_index] = 0;


        ecl_init_file_fwrite_header( fortio.get(),
                                     this->grid.c_ptr(),
                                     NULL,
                                     to_ert_unit( units.getType( ) ),
                                     this->ert_phase_mask,
                                     this->sim_start_time);

        units.from_si( UnitSystem::measure::volume, ecl_data );
        writeKeyword( fortio, "PORV" , ecl_data );
    }

    // Writing quantities which are calculated by the grid to the INIT file.
    ecl_grid_fwrite_depth( this->grid.c_ptr() , fortio.get() , to_ert_unit( units.getType( )) );
    ecl_grid_fwrite_dims( this->grid.c_ptr() , fortio.get() , to_ert_unit( units.getType( )) );

    // Write properties from the input deck.
    {
        const auto& properties = this->es.get3DProperties();
        using double_kw = std::pair<std::string, UnitSystem::measure>;
        std::vector<double_kw> doubleKeywords = {{"PORO"  , UnitSystem::measure::identity },
                                                 {"PERMX" , UnitSystem::measure::permeability },
                                                 {"PERMY" , UnitSystem::measure::permeability },
                                                 {"PERMZ" , UnitSystem::measure::permeability }};

        for (const auto& kw_pair : doubleKeywords) {
            const auto& opm_property = properties.getDoubleGridProperty(kw_pair.first);
            auto ecl_data = opm_property.compressedCopy( this->grid );

            units.from_si( kw_pair.second, ecl_data );
            writeKeyword( fortio, kw_pair.first, ecl_data );
        }
    }


    // Write properties which have been initialized by the simulator.
    {
        for (const auto& prop : simProps) {
            auto ecl_data = this->grid.compressedVector( prop.second.data );
            writeKeyword( fortio, prop.first, ecl_data );
        }
    }

    // Write tables
    {
        Tables tables( this->es.getUnits() );
        tables.addPVTO( this->es.getTableManager().getPvtoTables() );
        tables.addPVTG( this->es.getTableManager().getPvtgTables() );
        tables.fwrite( fortio );
    }

    // Write all integer field properties from the input deck.
    {
        const auto& properties = this->es.get3DProperties().getIntProperties();

        // It seems that the INIT file should always contain these
        // keywords, we therefor call getKeyword() here to invoke the
        // autocreation property, and ensure that the keywords exist
        // in the properties container.
        properties.getKeyword("PVTNUM");
        properties.getKeyword("SATNUM");
        properties.getKeyword("EQLNUM");
        properties.getKeyword("FIPNUM");

        for (const auto& property : properties) {
            auto ecl_data = property.compressedCopy( this->grid );
            writeKeyword( fortio , property.getKeywordName() , ecl_data );
        }
    }


    // Write NNC transmissibilities
    {
        std::vector<double> tran;
        for( const NNCdata& nd : nnc.nncdata() )
            tran.push_back( nd.trans );

        units.from_si( UnitSystem::measure::transmissibility , tran );
        writeKeyword( fortio, "TRANNNC" , tran );
    }
}


void EclipseIO::Impl::writeEGRIDFile( const NNC& nnc ) const {
    const auto& ioConfig = this->es.getIOConfig();

    std::string  egridFile( ERT::EclFilename( this->outputDir,
                                              this->baseName,
                                              ECL_EGRID_FILE,
                                              ioConfig.getFMTOUT() ));

    {
        int idx = 0;
        auto* ecl_grid = const_cast< ecl_grid_type* >( this->grid.c_ptr() );
        for (const NNCdata& n : nnc.nncdata())
            ecl_grid_add_self_nnc( ecl_grid, n.cell1, n.cell2, idx++);

        ecl_grid_fwrite_EGRID2(ecl_grid, egridFile.c_str(), to_ert_unit( this->es.getDeckUnitSystem().getType()));
    }
}


void EclipseIO::writeInitial( data::Solution simProps, const NNC& nnc) {
    if( !this->impl->output_enabled )
        return;

    {
        const auto& es = this->impl->es;
        const IOConfig& ioConfig = es.cfg().io();

        simProps.convertFromSI( es.getUnits() );
        if( ioConfig.getWriteINITFile() )
            this->impl->writeINITFile( simProps , nnc );

        if( ioConfig.getWriteEGRIDFile( ) )
            this->impl->writeEGRIDFile( nnc );
    }

    this->impl->summary.set_initial( simProps );
}



// implementation of the writeTimeStep method
void EclipseIO::writeTimeStep(int report_step,
                              bool  isSubstep,
                              double secs_elapsed,
                              data::Solution cells,
                              data::Wells wells)
{

    if( !this->impl->output_enabled )
        return;


    time_t current_posix_time = this->impl->sim_start_time + secs_elapsed;
    const auto& es = this->impl->es;
    const auto& grid = this->impl->grid;
    const auto& units = es.getUnits();
    const auto& ioConfig = es.getIOConfig();
    const auto& restart = es.cfg().restart();


    const auto days = units.from_si( UnitSystem::measure::time, secs_elapsed );
    const auto& schedule = es.getSchedule();

    /*
       This routine can optionally write RFT and/or restart file; to
       be certain that the data correctly converted to output units
       the conversion is done once here - and not closer to the actual
       use-site.
    */
    cells.convertFromSI( units );
    // Write restart file
    if(!isSubstep && restart.getWriteRestartFile(report_step))
    {
        const size_t ncwmax     = schedule.getMaxNumCompletionsForWells(report_step);
        const size_t numWells   = schedule.numWells(report_step);
        const auto sched_wells  = schedule.getWells(report_step);

        std::string filename = ERT::EclFilename( this->impl->outputDir,
                                                 this->impl->baseName,
                                                 ioConfig.getUNIFOUT() ? ECL_UNIFIED_RESTART_FILE : ECL_RESTART_FILE,
                                                 report_step,
                                                 ioConfig.getFMTOUT() );

        restart_file rst_file = open_rst( filename ,
                                          this->impl->first_restart,
                                          ioConfig.getUNIFOUT(),
                                          report_step );

        RestartIO::save( );
        this->impl->first_restart = false;


        {
            ecl_rsthead_type rsthead_data = {};
            rsthead_data.sim_time   = current_posix_time;
            rsthead_data.nactive    = grid.getNumActive(),
            rsthead_data.nx         = grid.getNX();
            rsthead_data.ny         = grid.getNY();
            rsthead_data.nz         = grid.getNZ();
            rsthead_data.nwells     = numWells;
            rsthead_data.niwelz     = RestartIO::NIWELZ;
            rsthead_data.nzwelz     = RestartIO::NZWELZ;
            rsthead_data.niconz     = RestartIO::NICONZ;
            rsthead_data.ncwmax     = ncwmax;
            rsthead_data.phase_sum  = this->impl->ert_phase_mask;
            rsthead_data.sim_days   = days;
            rsthead_data.unit_system= to_ert_unit( units.getType() );

            RestartIO::writeHeader( rst_file , report_step, &rsthead_data);
        }

        {
            const auto& phases = es.runspec().phases();
            const auto opm_xwel = RestartIO::serialize_OPM_XWEL( wells, report_step, sched_wells, phases, grid );
            const auto opm_iwel = RestartIO::serialize_OPM_IWEL( wells, sched_wells );
            const auto iwel_data = RestartIO::serialize_IWEL(report_step, sched_wells);
            const auto icon_data = RestartIO::serialize_ICON(report_step , ncwmax, sched_wells);
            const auto zwel_data = RestartIO::serialize_ZWEL( sched_wells );

            RestartIO::write_kw( rst_file, ERT::EclKW< int >( IWEL_KW, iwel_data) );
            RestartIO::write_kw( rst_file, ERT::EclKW< const char* >(ZWEL_KW, zwel_data ) );
            RestartIO::write_kw( rst_file, ERT::EclKW< double >( OPM_XWEL, opm_xwel ) );
            RestartIO::write_kw( rst_file, ERT::EclKW< int >( OPM_IWEL, opm_iwel ) );
            RestartIO::write_kw( rst_file, ERT::EclKW< int >( ICON_KW, icon_data ) );
        }

        RestartIO::writeSolution( rst_file , cells );
    }


    const auto unit_type = es.getDeckUnitSystem().getType();
    this->impl->rft.writeTimeStep( schedule.getWells( report_step ),
                                   grid,
                                   report_step,
                                   current_posix_time,
                                   days,
                                   to_ert_unit( unit_type ),
                                   cells );

    if( isSubstep ) return;

    this->impl->summary.add_timestep( report_step,
                                      secs_elapsed,
                                      this->impl->es,
                                      this->impl->regionCache,
                                      wells ,
                                      cells_si );
    this->impl->summary.write();
 }



EclipseIO::EclipseIO( const EclipseState& es, EclipseGrid grid)
    : impl( new Impl( es, std::move( grid ) ) )
{
    if( !this->impl->output_enabled )
        return;
    {
        const auto& outputDir = this->impl->outputDir;

        // make sure that the output directory exists, if not try to create it
        if ( !util_entry_exists( outputDir.c_str() ) ) {
            util_make_path( outputDir.c_str() );
        }

        if( !util_is_directory( outputDir.c_str() ) ) {
            throw std::runtime_error( "The path specified as output directory '"
                                      + outputDir + "' is not a directory");
        }
    }
}

std::pair< data::Solution, data::Wells >
EclipseIO::load_from_restart_file( const EclipseState& es, const std::map<std::string, UnitSystem::measure>& keys, int numcells ) const {
    return RestartIO::load( es , keys , numcells );
}

EclipseIO::~EclipseIO() {}

} // namespace Opm
