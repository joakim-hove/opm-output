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

#include <ert/ecl/EclKW.hpp>
#include <ert/ecl/ecl_rsthead.h>
#include <ert/ecl/ecl_rst_file.h>
#include <ert/util/util.h>

namespace Opm {

class EclipseGrid;
class EclipseState;
class Phases;
class Schedule;

namespace RestartIO {




// Should be private:
void writeHeader(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file ,
                 int report_step,
                 time_t posix_time,
                 double sim_days,
                 int ert_phase_mask,
                 const Schedule& schedule,
                 const EclipseGrid& grid);
void writeSolution(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , const data::Solution& solution);

void writeWell(ERT::ert_unique_ptr< ecl_rst_file_type, ecl_rst_file_close >& rst_file , int report_step, const EclipseState& es , const EclipseGrid& grid, const data::Wells& wells);

void save(const std::string& filename,
          int report_step,
          double seconds_elapsed,
          data::Solution cells,
          data::Wells wells,
          const EclipseState& es,
          const EclipseGrid& grid);



std::pair< data::Solution, data::Wells > load( const EclipseState& es, const std::map<std::string, UnitSystem::measure>& keys, int numcells );

}
}
#endif
