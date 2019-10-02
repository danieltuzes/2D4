/*
 * SDDDST Simple Discrete Dislocation Dynamics Toolkit
 * Copyright (C) 2015-2019  Gábor Péterffy <peterffy95@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifndef SDDDST_CORE_STRESS_PROTOCOL_H
#define SDDDST_CORE_STRESS_PROTOCOL_H

#include "dislocation.h"

#include <memory>
#include <string>
#include <vector>

namespace sdddstCore {

    enum class StressProtocolStepType
    {
        Original,
        EndOfBigStep,
        EndOfFirstSmallStep,
        EndOfSecondSmallStep
    };

    /**
     * @brief The StressProtocol class is the base class for all external stress protocol classes,
     * all parameters for calculation is used from a SimulationData object or should be preset
     */
    class StressProtocol
    {
    public:
        StressProtocol(double initExtStress);
        StressProtocol();

        /**
         * @brief calculateStress calculates the stress value for the given situation
         * @param simulationTime
         * @param dislocations
         * @param type
         */
        virtual void calcExtStress(double simulationTime, StressProtocolStepType type);

        /**
         * @brief getStress returns with the stress value at the given situation
         * @param type
         * @return
         */
        virtual double getExtStress(StressProtocolStepType type) const;

        /**
         * @brief getType
         * @return returns with the type of the applied stress field
         */
        virtual std::string getType() const;

    protected:
        double initExtStress;
    };
}

#endif
