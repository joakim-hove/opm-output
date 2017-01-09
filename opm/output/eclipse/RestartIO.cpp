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

    inline data::Solution restoreSOLUTION( ecl_file_type* file,
                                           const std::map<std::string, UnitSystem::measure>& keys,
                                           int numcells,
                                           const UnitSystem& units ) {

        data::Solution sol;
        for (const auto& pair : keys) {
            const std::string& key = pair.first;
            UnitSystem::measure dim = pair.second;
            if( !ecl_file_has_kw( file, key.c_str() ) )
                throw std::runtime_error("Read of restart file: "
                                         "File does not contain "
                                         + key
                                         + " data" );

            const ecl_kw_type * ecl_kw = ecl_file_iget_named_kw( file , key.c_str() , 0 );
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
data::Wells restore_wells( const double* xwel_data,
                           size_t xwel_data_size,
                           const int* iwel_data,
                           size_t iwel_data_size,
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

    if( xwel_data_size != expected_xwel_size ) {
        throw std::runtime_error(
                "Mismatch between OPM_XWEL and deck; "
                "OPM_XWEL size was " + std::to_string( xwel_data_size ) +
                ", expected " + std::to_string( expected_xwel_size ) );
    }

    if( iwel_data_size != sched_wells.size() )
        throw std::runtime_error(
                "Mismatch between OPM_IWEL and deck; "
                "OPM_IWEL size was " + std::to_string( iwel_data_size ) +
                ", expected " + std::to_string( sched_wells.size() ) );

    data::Wells wells;
    for( const auto* sched_well : sched_wells ) {
        data::Well& well = wells[ sched_well->name() ];

        well.bhp = *xwel_data++;
        well.temperature = *xwel_data++;
        well.control = *iwel_data++;

        for( auto phase : phases )
            well.rates.set( phase, *xwel_data++ );

        for( const auto& sc : sched_well->getCompletions( restart_step ) ) {
            const auto i = sc.getI(), j = sc.getJ(), k = sc.getK();
            if( !grid.cellActive( i, j, k ) || sc.getState() == WellCompletion::SHUT ) {
                xwel_data += data::Completion::restart_size + phases.size();
                continue;
            }

            const auto active_index = grid.activeIndex( i, j, k );

            well.completions.emplace_back();
            auto& completion = well.completions.back();
            completion.index = active_index;
            completion.pressure = *xwel_data++;
            completion.reservoir_rate = *xwel_data++;
            for( auto phase : phases )
                completion.rates.set( phase, *xwel_data++ );
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
    using ft = ERT::ert_unique_ptr< ecl_file_type, ecl_file_close >;
    ft file( ecl_file_open( filename.c_str(), 0 ) );

    if( !file )
        throw std::runtime_error( "Restart file " + filename + " not found!" );

    if( unified &&
        !ecl_file_select_rstblock_report_step( file.get(), restart_step ) ) {
        throw std::runtime_error( "Restart file " + filename
                + " does not contain data for report step "
                + std::to_string( restart_step ) + "!" );
    }

    ecl_kw_type* xwel = ecl_file_iget_named_kw( file.get(), "OPM_XWEL", 0 );
    const double* xwel_data = ecl_kw_get_double_ptr( xwel );
    const auto xwel_size = ecl_kw_get_size( xwel );

    ecl_kw_type* iwel = ecl_file_iget_named_kw( file.get(), "OPM_IWEL", 0 );
    const int* iwel_data = ecl_kw_get_int_ptr( iwel );
    const auto iwel_size = ecl_kw_get_size( iwel );

    return {
        restoreSOLUTION( file.get(), keys, numcells, es.getUnits() ),
        restore_wells( xwel_data, xwel_size,
                       iwel_data, iwel_size,
                       restart_step,
                       es )
    };
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

