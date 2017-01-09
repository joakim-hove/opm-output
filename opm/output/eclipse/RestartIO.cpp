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
#include <string>
#include <vector>

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/output/eclipse/RestartIO.hpp>

#include <ert/ecl/EclKW.hpp>
#include <ert/ecl/FortIO.hpp>
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

namespace Opm {
namespace RestartIO  {

void save() {

}


namespace {

    inline int to_ert_welltype( const Well& well, size_t timestep ) {
        if( well.isProducer( timestep ) ) return IWEL_PRODUCER;

        switch( well.getInjectionProperties( timestep ).injectorType ) {
        case WellInjector::WATER:
            return IWEL_WATER_INJECTOR;
        case WellInjector::GAS:
            return IWEL_GAS_INJECTOR;
        case WellInjector::OIL:
            return IWEL_OIL_INJECTOR;
        default:
            return IWEL_UNDOCUMENTED_ZERO;
        }
    }

    std::vector<double> double_vector( const ecl_kw_type * ecl_kw ) {
        size_t size = ecl_kw_get_size( ecl_kw );

        if (ecl_kw_get_type( ecl_kw ) == ECL_DOUBLE_TYPE ) {
            const double * ecl_data = ecl_kw_get_double_ptr( ecl_kw );
            return { ecl_data , ecl_data + size };
        } else {
            const float * ecl_data = ecl_kw_get_float_ptr( ecl_kw );
            return { ecl_data , ecl_data + size };
        }

    }


    inline data::Solution restoreSOLUTION( ecl_file_view_type* file_view,
                                           const std::map<std::string, UnitSystem::measure>& keys,
                                           int numcells,
                                           const UnitSystem& units ) {

        data::Solution sol;
        for (const auto& pair : keys) {
            const std::string& key = pair.first;
            UnitSystem::measure dim = pair.second;
            if( !ecl_file_view_has_kw( file_view, key.c_str() ) )
                throw std::runtime_error("Read of restart file: "
                                         "File does not contain "
                                         + key
                                         + " data" );

            const ecl_kw_type * ecl_kw = ecl_file_view_iget_named_kw( file_view , key.c_str() , 0 );
            if( ecl_kw_get_size(ecl_kw) != numcells)
                throw std::runtime_error("Restart file: Could not restore "
                                         + std::string( ecl_kw_get_header( ecl_kw ) )
                                         + ", mismatched number of cells" );

            std::vector<double> data = double_vector( ecl_kw );
            units.to_si( dim , data );

            sol.insert( key, dim, data , data::TargetType::RESTART_SOLUTION );
        }

        return sol;
    }


using rt = data::Rates::opt;
    data::Wells restore_wells( const ecl_kw_type * opm_xwel,
                               const ecl_kw_type * opm_iwel,
                               int restart_step,
                               const EclipseState& es ) {

    const auto& sched_wells = es.getSchedule().getWells( restart_step );
    const EclipseGrid& grid = es.getInputGrid( );
    std::vector< rt > phases;
    const auto& phase = es.runspec().phases();
    if( phase.active( Phase::WATER ) ) phases.push_back( rt::wat );
    if( phase.active( Phase::OIL ) )   phases.push_back( rt::oil );
    if( phase.active( Phase::GAS ) )   phases.push_back( rt::gas );


    const auto well_size = [&]( size_t acc, const Well* w ) {
        return acc
            + 2 + phases.size()
            + (  w->getCompletions( restart_step ).size()
              * (phases.size() + data::Completion::restart_size) );
    };

    const auto expected_xwel_size = std::accumulate( sched_wells.begin(),
                                                     sched_wells.end(),
                                                     size_t( 0 ),
                                                     well_size );

    if( ecl_kw_get_size( opm_xwel ) != expected_xwel_size ) {
        throw std::runtime_error(
                "Mismatch between OPM_XWEL and deck; "
                "OPM_XWEL size was " + std::to_string( ecl_kw_get_size( opm_xwel ) ) +
                ", expected " + std::to_string( expected_xwel_size ) );
    }

    if( ecl_kw_get_size( opm_iwel ) != sched_wells.size() )
        throw std::runtime_error(
                "Mismatch between OPM_IWEL and deck; "
                "OPM_IWEL size was " + std::to_string( ecl_kw_get_size( opm_iwel ) ) +
                ", expected " + std::to_string( sched_wells.size() ) );

    data::Wells wells;
    const double * opm_xwel_data = ecl_kw_get_double_ptr( opm_xwel );
    const int * opm_iwel_data = ecl_kw_get_int_ptr( opm_iwel );
    for( const auto* sched_well : sched_wells ) {
        data::Well& well = wells[ sched_well->name() ];

        well.bhp = *opm_xwel_data++;
        well.temperature = *opm_xwel_data++;
        well.control = *opm_iwel_data++;

        for( auto phase : phases )
            well.rates.set( phase, *opm_xwel_data++ );

        for( const auto& sc : sched_well->getCompletions( restart_step ) ) {
            const auto i = sc.getI(), j = sc.getJ(), k = sc.getK();
            if( !grid.cellActive( i, j, k ) || sc.getState() == WellCompletion::SHUT ) {
                opm_xwel_data += data::Completion::restart_size + phases.size();
                continue;
            }

            const auto active_index = grid.activeIndex( i, j, k );

            well.completions.emplace_back();
            auto& completion = well.completions.back();
            completion.index = active_index;
            completion.pressure = *opm_xwel_data++;
            completion.reservoir_rate = *opm_xwel_data++;
            for( auto phase : phases )
                completion.rates.set( phase, *opm_xwel_data++ );
        }
    }

    return wells;
}
}

/* should take grid as argument because it may be modified from the simulator */
std::pair< data::Solution, data::Wells > load( const EclipseState& es, const std::map<std::string, UnitSystem::measure>& keys, int numcells ) {

    const InitConfig& initConfig         = es.getInitConfig();
    const auto& ioConfig                 = es.getIOConfig();
    int restart_step                     = initConfig.getRestartStep();
    const std::string& restart_file_root = initConfig.getRestartRootName();
    bool output                          = false;
    const std::string filename           = ioConfig.getRestartFileName(restart_file_root,
                                                                       restart_step,
                                                                       output);
    const bool unified                   = ioConfig.getUNIFIN();
    ERT::ert_unique_ptr< ecl_file_type, ecl_file_close > file(ecl_file_open( filename.c_str(), 0 ));
    ecl_file_view_type * file_view;

    if( !file )
        throw std::runtime_error( "Restart file " + filename + " not found!" );

    if( unified ) {
        file_view = ecl_file_get_restart_view( file.get() , -1 , restart_step , -1 , -1 );
        if (!file_view)
            throw std::runtime_error( "Restart file " + filename
                                      + " does not contain data for report step "
                                      + std::to_string( restart_step ) + "!" );
    } else
        file_view = ecl_file_get_global_view( file.get() );

    const ecl_kw_type* opm_xwel = ecl_file_view_iget_named_kw( file_view , "OPM_XWEL", 0 );
    const ecl_kw_type* opm_iwel = ecl_file_view_iget_named_kw( file_view, "OPM_IWEL", 0 );
    return {
        restoreSOLUTION( file_view, keys, numcells, es.getUnits() ),
        restore_wells( opm_xwel, opm_iwel,
                       restart_step,
                       es )
    };
}



std::vector<int> serialize_ICON( int report_step,
                                 int ncwmax,
                                 const std::vector<const Well*>& sched_wells) {

    size_t well_offset = 0;
    std::vector<int> data( sched_wells.size() * ncwmax * RestartIO::NICONZ , 0 );
    for (const Well* well : sched_wells) {
        const auto& completions = well->getCompletions( report_step );
        size_t completion_offset = 0;
        for( const auto& completion : completions) {
            size_t offset = well_offset + completion_offset;
            data[ offset + ICON_IC_INDEX ] = 1;

            data[ offset + ICON_I_INDEX ] = completion.getI() + 1;
            data[ offset + ICON_J_INDEX ] = completion.getJ() + 1;
            data[ offset + ICON_K_INDEX ] = completion.getK() + 1;
            data[ offset + ICON_DIRECTION_INDEX ] = completion.getDirection();
            {
                const auto open = WellCompletion::StateEnum::OPEN;
                data[ offset + ICON_STATUS_INDEX ] = completion.getState() == open
                    ? 1
                    : 0;
            }

            completion_offset += RestartIO::NICONZ;
        }
        well_offset += ncwmax * RestartIO::NICONZ;
    }
    return data;
}

std::vector<int> serialize_IWEL( size_t step,
                                 const std::vector<const Well *>& wells) {

    std::vector<int> data( wells.size() * RestartIO::NIWELZ , 0 );
    size_t offset = 0;
    for (const auto well : wells) {
        const auto& completions = well->getCompletions( step );

        data[ offset + IWEL_HEADI_INDEX ] = well->getHeadI( step ) + 1;
        data[ offset + IWEL_HEADJ_INDEX ] = well->getHeadJ( step ) + 1;
        data[ offset + IWEL_CONNECTIONS_INDEX ] = completions.size();
        data[ offset + IWEL_GROUP_INDEX ] = 1;

        data[ offset + IWEL_TYPE_INDEX ] = to_ert_welltype( *well, step );
        data[ offset + IWEL_STATUS_INDEX ] =
            well->getStatus( step ) == WellCommon::OPEN ? 1 : 0;

        offset += RestartIO::NIWELZ;
    }
    return data;
}

std::vector< int > serialize_OPM_IWEL( const data::Wells& wells,
                                       const std::vector< const Well* > sched_wells ) {

    const auto getctrl = [&]( const Well* w ) {
        const auto itr = wells.find( w->name() );
        return itr == wells.end() ? 0 : itr->second.control;
    };

    std::vector< int > iwel( sched_wells.size(), 0.0 );
    std::transform( sched_wells.begin(), sched_wells.end(), iwel.begin(), getctrl );

    return iwel;
}

std::vector< double > serialize_OPM_XWEL( const data::Wells& wells,
                                      int report_step,
                                      const std::vector< const Well* > sched_wells,
                                      const Phases& phase_spec,
                                      const EclipseGrid& grid ) {

    using rt = data::Rates::opt;

    std::vector< rt > phases;
    if( phase_spec.active( Phase::WATER ) ) phases.push_back( rt::wat );
    if( phase_spec.active( Phase::OIL ) )   phases.push_back( rt::oil );
    if( phase_spec.active( Phase::GAS ) )   phases.push_back( rt::gas );

    std::vector< double > xwel;
    for( const auto* sched_well : sched_wells ) {

        if( wells.count( sched_well->name() ) == 0 ) {
            const auto elems = (sched_well->getCompletions( report_step ).size()
                               * (phases.size() + data::Completion::restart_size))
                + 2 /* bhp, temperature */
                + phases.size();

            // write zeros if no well data is provided
            xwel.insert( xwel.end(), elems, 0.0 );
            continue;
        }

        const auto& well = wells.at( sched_well->name() );

        xwel.push_back( well.bhp );
        xwel.push_back( well.temperature );
        for( auto phase : phases )
            xwel.push_back( well.rates.get( phase ) );

        for( const auto& sc : sched_well->getCompletions( report_step ) ) {
            const auto i = sc.getI(), j = sc.getJ(), k = sc.getK();

            const auto rs_size = phases.size() + data::Completion::restart_size;
            if( !grid.cellActive( i, j, k ) || sc.getState() == WellCompletion::SHUT ) {
                xwel.insert( xwel.end(), rs_size, 0.0 );
                continue;
            }

            const auto active_index = grid.activeIndex( i, j, k );
            const auto at_index = [=]( const data::Completion& c ) {
                return c.index == active_index;
            };
            const auto& completion = std::find_if( well.completions.begin(),
                                                   well.completions.end(),
                                                   at_index );

            if( completion == well.completions.end() ) {
                xwel.insert( xwel.end(), rs_size, 0.0 );
                continue;
            }

            xwel.push_back( completion->pressure );
            xwel.push_back( completion->reservoir_rate );
            for( auto phase : phases )
                xwel.push_back( completion->rates.get( phase ) );
        }
    }

    return xwel;
};


std::vector<const char*> serialize_ZWEL( const std::vector<const Well *>& wells) {
    std::vector<const char*> data( wells.size( ) * RestartIO::NZWELZ , "");
    size_t offset = 0;

    for (const auto& well : wells) {
        data[ offset ] = well->name().c_str();
        offset += RestartIO::NZWELZ;
    }
    return data;
}



void writeHeader(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , int stepIdx, ecl_rsthead_type* rsthead_data ) {
  ecl_util_set_date_values( rsthead_data->sim_time,
                            &rsthead_data->day,
                            &rsthead_data->month,
                            &rsthead_data->year );
  ecl_rst_file_fwrite_header( rst_file.get() , stepIdx, rsthead_data );

}

void writeSolution(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , const data::Solution& solution) {
    ecl_rst_file_start_solution( rst_file.get() );
    for (const auto& elm: solution) {
        if (elm.second.target == data::TargetType::RESTART_SOLUTION)
            ecl_rst_file_add_kw( rst_file.get() , ERT::EclKW<float>(elm.first, elm.second.data).get());
     }
     ecl_rst_file_end_solution( rst_file.get() );

     for (const auto& elm: solution) {
        if (elm.second.target == data::TargetType::RESTART_AUXILLARY)
            ecl_rst_file_add_kw( rst_file.get() , ERT::EclKW<float>(elm.first, elm.second.data).get());
     }
}


}
}

