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
#ifndef RESTART_IO_HPP
#define RESTART_IO_HPP

#include <vector>

#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>

#include <ert/ecl/ecl_rsthead.h>
#include <ert/ecl/ecl_rst_file.h>
#include <ert/util/util.h>

namespace Opm {

class EclipseGrid;
class EclipseState;
class Phases;

namespace RestartIO {
    static const int NIWELZ = 11; //Number of data elements per well in IWEL array in restart file
    static const int NZWELZ = 3;  //Number of 8-character words per well in ZWEL array restart file
    static const int NICONZ = 15; //Number of data elements per completion in ICON array restart file

    /**
     * The constants NIWELZ and NZWELZ referes to the number of
     * elements per well that we write to the IWEL and ZWEL eclipse
     * restart file data arrays. The constant NICONZ refers to the
     * number of elements per completion in the eclipse restart file
     * ICON data array.These numbers are written to the INTEHEAD
     * header.
     *
     * Observe that all of these values are our "current-best-guess"
     * for how many numbers are needed; there might very well be third
     * party applications out there which have a hard expectation for
     * these values.
     */




// Should be private:
void writeHeader(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , int stepIdx, ecl_rsthead_type* rsthead_data );
void writeSolution(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , const data::Solution& solution);

std::vector< double > serialize_OPM_XWEL( const data::Wells& wells,
                                      int report_step,
                                      const std::vector< const Well* > sched_wells,
                                      const Phases& phase_spec,
                                      const EclipseGrid& grid );

std::vector< int > serialize_OPM_IWEL( const data::Wells& wells,
                                       const std::vector< const Well* > sched_wells );

void save();



std::pair< data::Solution, data::Wells > load( const EclipseState& es, const std::map<std::string, UnitSystem::measure>& keys, int numcells );

}
}
#endif
