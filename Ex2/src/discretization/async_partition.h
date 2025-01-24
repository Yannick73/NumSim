#pragma once

#include <cassert>
#include <vector>
#include <set>
#include "settings.h"
#include "timekeeper.h"
#include "discretization/discretization.h"
#include "discretization/partition_shell.h"
#include "discretization/partition_information.h"
#include "boundary/boundary.h"
#include "boundary/dirichlet.h"
#include "boundary/async_neighbour_boundary.h"

//! Class to encapsulate the discretization and its boundaries
class AsyncPartition : public PartitionShell
{
public:
    AsyncPartition(const std::shared_ptr<Discretization> discretization,
              const Settings &settings,
              const PartitionInformation &pi);

    //! set applies all boundary conditions on fix boundary an on neighbours
    void setBoundaryUV() override;
    void setBoundaryFG() override;
    void setBoundaryP() override;
    //! only exchange the pressure values with the respective neighbours w/o setting dirichlet
    void exchangeP() override;
    //! also only used in pressure solver, but specifically only for CG
    void exchangeRA() override;
    //! used before paraview output
    void exchangeUV() override;

    //! sets up the appropriate send and start-receive for all neighbours and pushes the 
    //! neighbours into the neighbourRecvQueue for later use
    inline void setupExchange(std::vector<std::shared_ptr<AsyncNeighbourBoundary>> &neighbourRecvQueue, 
                              void (AsyncNeighbourBoundary::*setupFun)())
    {
        #ifdef TIMER
        timer_.setT0();
        #endif
        for(std::shared_ptr<AsyncNeighbourBoundary> neighbour : asyncNeighbours_)
        {
            (neighbour.get()->*setupFun)();
            neighbourRecvQueue.push_back(neighbour);
        }
        #ifdef TIMER
        timer_.addTimeSinceT0();
        #endif
    }

    //! this little function takes sets the data of the first incoming data 
    inline void setFirstIncomingData(std::vector<std::shared_ptr<AsyncNeighbourBoundary>> &neighbourRecvQueue,
                                    void (AsyncNeighbourBoundary::*setFun)())
    {
        // Although adding/removing elements is the domain of lists, the vector 
        // datatype is so lightweight&in cache, that it is very likely to outpeform lists (for N very small).
        // However it migt be interesting to test this hypothesis some time.
        #ifdef TIMER
        timer_.setT0();
        #endif
        while(neighbourRecvQueue.size() > 0)
        {
            for(int n = 0; n < neighbourRecvQueue.size(); n++)
            {
                if(neighbourRecvQueue[n]->mpiHandler_.queryRecvComplete())
                {
                    (neighbourRecvQueue[n].get()->*setFun)();
                    neighbourRecvQueue.erase(neighbourRecvQueue.begin() + n);
                    break;
                }
            }
        }
        #ifdef TIMER
        timer_.addTimeSinceT0();
        #endif
    }

protected:
    //! Extra async neighbours vector, making them neighbour-derived migt solve this,
    //! but would remove neighbour specific virtual functions. One may solve this again by making
    //! a virtual SyncNeighbour class, but well...
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> asyncNeighbours_;
};