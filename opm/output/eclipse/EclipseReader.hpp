#ifndef ECLIPSEREADER_HPP
#define ECLIPSEREADER_HPP

#include <map>
#include <utility>

#include <opm/parser/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/data/Solution.hpp>

namespace Opm {

    class EclipseState;

    /*
      Will load solution data and wellstate from the restart file. The
      map keys should be a map of keyword names and their
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
      container type which is used by the EclipseWriter, but observe
      that the dim and target elements carry no information:

         - The returned double data has been converted to SI.
         . The target is unconditionally set to 'RESTART_SOLUTION'
    */

    std::pair< data::Solution, data::Wells >
    load_from_restart_file( const EclipseState&, const std::map<std::string, UnitSystem::measure>& keys,  int numcells );


}

#endif // ECLIPSEREADER_HPP
