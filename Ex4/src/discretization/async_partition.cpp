#include "discretization/async_partition.h"

AsyncPartition::AsyncPartition(const std::shared_ptr<Discretization> discretization,
                               const Settings &settings,
                               const PartitionInformation &pi) :
                               PartitionShell(discretization, pi)
{
    assert(discretization != nullptr);

    // generate the boundaries, horizontal boundaries take priority, hence why they come first
    // odd y position and even y positions are flipped with regards of priority
    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletBottom>(discretization, settings.dirichletBcBottom));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourBottom>(discretization, pi.bottomRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));

        if(pi.ownTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletTop>(discretization, settings.dirichletBcTop));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourTop>(discretization, pi.topRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));
    }
    else
    {
        if(pi.ownTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletTop>(discretization, settings.dirichletBcTop));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourTop>(discretization, pi.topRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));

        if(pi.ownBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletBottom>(discretization, settings.dirichletBcBottom));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourBottom>(discretization, pi.bottomRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));
    }

    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletLeft> (discretization, settings.dirichletBcLeft));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourLeft> (discretization, pi.leftRank()));

        if(pi.ownRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletRight> (discretization, settings.dirichletBcRight));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourRight> (discretization, pi.rightRank()));
    }
    else
    {
        if(pi.ownRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletRight> (discretization, settings.dirichletBcRight));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourRight> (discretization, pi.rightRank()));

        if(pi.ownLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletLeft> (discretization, settings.dirichletBcLeft));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourLeft> (discretization, pi.leftRank()));
    }

    if(pi.getPartPosZ() & 0b1)
    {
        if(pi.ownFrontBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletFront> (discretization, settings.dirichletBcFront));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourFront> (discretization, pi.frontRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));

        if(pi.ownHindBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletHind> (discretization, settings.dirichletBcHind));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourHind> (discretization, pi.hindRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));
    }
    else
    {
        if(pi.ownHindBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletHind> (discretization, settings.dirichletBcHind));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourHind> (discretization, pi.hindRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));

        if(pi.ownFrontBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletFront> (discretization, settings.dirichletBcFront));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourFront> (discretization, pi.frontRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));
    }

    // checks, that each direction exist, while the dirichlet constraints take priority
    // also puts every neighbour in a common boundaries object
    std::set<BoundaryEdge> directions;
    for(std::shared_ptr<Dirichlet> boundary : fixBoundaries_)
    {
        assert(boundary != nullptr);
        directions.insert(boundary->edge_);
    }
    for(std::shared_ptr<AsyncNeighbourBoundary> neigbour : asyncNeighbours_)
    {
        assert(neigbour != nullptr);
        directions.insert(neigbour->edge_);
    }
    // each direction should be unique and as such should have each an entry in the set
    assert(directions.size() == 6);
}

// First this method sets up async receive, sends its data
// then it sets the dirichlet boundary, finally setting the first incoming ghost data
void AsyncPartition::setBoundaryUVW()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeUVW);
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setUVW();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvUVW);
}

void AsyncPartition::setBoundaryFGH()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeFGH);
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setFGH();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvFGH);
}

void AsyncPartition::setBoundaryP()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeP);
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setP();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvP);
}

void AsyncPartition::exchangeP()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeP);
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvP);
}

void AsyncPartition::exchangeUVW()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeUVW);
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvUVW);
}