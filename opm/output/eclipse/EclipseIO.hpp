/*
  Copyright (c) 2013 Andreas Lauser
  Copyright (c) 2013 Uni Research AS
  Copyright (c) 2014 IRIS AS

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

#ifndef OPM_ECLIPSE_WRITER_HPP
#define OPM_ECLIPSE_WRITER_HPP

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>

#include <string>
#include <vector>
#include <array>
#include <memory>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>

namespace Opm {

class EclipseState;

/*!
 * \brief A class to write the reservoir state and the well state of a
 *        blackoil simulation to disk using the Eclipse binary format.
 */
class EclipseIO {
public:
    /*!
     * \brief Sets the common attributes required to write eclipse
     *        binary files using ERT.
     */
    EclipseIO( const EclipseState&, EclipseGrid );




    /// Write the static eclipse data (grid, PVT curves, etc) to disk, and set
    /// up additional initial properties
    ///
    /// - simProps contains a list of properties which must be
    ///   calculated by the simulator, e.g. the transmissibility
    ///   vectors TRANX, TRANY and TRANZ.
    ///
    /// - The NNC argument is distributed between the EGRID and INIT
    ///   files.
    ///
    /// If you want FOE in the SUMMARY section, you must pass the initial
    /// oil-in-place "OIP" key in Solution
    void writeInitial( data::Solution = data::Solution(), const NNC& = NNC());

    /*!
     * \brief Write a reservoir state and summary information to disk.
     *
     *
     * The reservoir state can be inspected with visualization tools like
     * ResInsight.
     *
     * The summary information can then be visualized using tools from
     * ERT or ECLIPSE. Note that calling this method is only
     * meaningful after the first time step has been completed.
     *
     * The optional simProps vector contains fields which have been
     * calculated by the simulator and are written to the restart
     * file. Examples of such fields would be the relative
     * permeabilities KRO, KRW and KRG and fluxes. The keywords which
     * can be added here are represented with mnenonics in the RPTRST
     * keyword.
     */

    void writeTimeStep( int report_step,
                        bool isSubstep,
                        double seconds_elapsed,
                        data::Solution,
                        data::Wells);


    /*
      Will load solution data and wellstate from the restart
      file. This method will consult the IOConfig object to get
      filename and report step to restart from.

      The map keys should be a map of keyword names and their
      corresponding dimension object, i.e. to load the state from a
      simple two phase simulation you would pass:

         keys = {{"PRESSURE" , UnitSystem::measure::pressure},
                 {"SWAT" , UnitSystem::measure::identity }}

      For a three phase black oil simulation you would add pairs for
      SGAS, RS and RV. If you ask for keys which are not found in the
      restart file an exception will be raised, the happens if the
      size of a vector is wrong.

      The function will consult the InitConfig object in the
      EclipseState object to determine which file and report step to
      load.

      The return value is of type 'data::Solution', which is the same
      container type which is used by the EclipseIO, but observe
      that the dim and target elements carry no information:

         - The returned double data has been converted to SI.
         . The target is unconditionally set to 'RESTART_SOLUTION'
    */
    std::pair< data::Solution, data::Wells >
    loadRestart(const std::map<std::string, UnitSystem::measure>& keys) const;


    EclipseIO( const EclipseIO& ) = delete;
    ~EclipseIO();

private:
    class Impl;
    std::unique_ptr< Impl > impl;

};

} // namespace Opm


#endif // OPM_ECLIPSE_WRITER_HPP
