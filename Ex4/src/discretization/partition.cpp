#include "discretization/partition.h"

Partition::Partition(const std::shared_ptr<Discretization> discretization,
                     const Settings &settings,
                     const PartitionInformation &pi) :
                     PartitionShell(discretization, pi)
{
    assert(discretization != nullptr);

    // generate the boundaries, horizontal boundaries take priority, hence why they come first
    // odd y position and even y positions are flipped with regards of priority
    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownPartitionContainsBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletBottom>(discretization, settings.dirichletBcBottom));
        else
            neighbours_.push_back(std::make_shared<NeighbourSouth>(discretization, pi.bottomNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));

        if(pi.ownPartitionContainsTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletTop>(discretization, settings.dirichletBcTop));
        else
            neighbours_.push_back(std::make_shared<NeighbourNorth>(discretization, pi.topNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));
    }
    else
    {
        if(pi.ownPartitionContainsTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletTop>(discretization, settings.dirichletBcTop));
        else
            neighbours_.push_back(std::make_shared<NeighbourNorth>(discretization, pi.topNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));

        if(pi.ownPartitionContainsBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletBottom>(discretization, settings.dirichletBcBottom));
        else
            neighbours_.push_back(std::make_shared<NeighbourSouth>(discretization, pi.bottomNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));
    }

    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownPartitionContainsLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletLeft> (discretization, settings.dirichletBcLeft));
        else
            neighbours_.push_back(std::make_shared<NeighbourWest> (discretization, pi.leftNeighbourRankNo()));

        if(pi.ownPartitionContainsRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletRight> (discretization, settings.dirichletBcRight));
        else
            neighbours_.push_back(std::make_shared<NeighbourEast> (discretization, pi.rightNeighbourRankNo()));
    }
    else
    {
        if(pi.ownPartitionContainsRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletRight> (discretization, settings.dirichletBcRight));
        else
            neighbours_.push_back(std::make_shared<NeighbourEast> (discretization, pi.rightNeighbourRankNo()));

        if(pi.ownPartitionContainsLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletLeft> (discretization, settings.dirichletBcLeft));
        else
            neighbours_.push_back(std::make_shared<NeighbourWest> (discretization, pi.leftNeighbourRankNo()));
    }

    // checks, that each direction exist, while the dirichlet constraints take priority
    // also puts every neighbour in a common boundaries object
    std::set<BoundaryEdge> directions;
    for(std::shared_ptr<Dirichlet> boundary : fixBoundaries_)
    {
        assert(boundary != nullptr);
        directions.insert(boundary->edge_);
    }
    for(std::shared_ptr<NeighbourBoundary> neigbour : neighbours_)
    {
        assert(neigbour != nullptr);
        directions.insert(neigbour->edge_);
    }
    // each direction should be unique and as such should have each an entry in the set
    assert(directions.size()  == 4);
}

void Partition::setBoundaryUV()
{
    fixDirichletBoundary();
    exchangeUV();
}

void Partition::setBoundaryFG()
{
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setFGH();
    }

    #ifdef TIMER
    timer_.setT0();
    #endif
    for(std::shared_ptr<NeighbourBoundary> neighbour : neighbours_)
    {
        neighbour->setFG();
    }
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif
}

void Partition::setBoundaryP()
{
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setP();
    }
    exchangeP();
}

void Partition::exchangeP()
{
    #ifdef TIMER
    timer_.setT0();
    #endif
    for(std::shared_ptr<NeighbourBoundary> neighbour : neighbours_)
    {
        neighbour->setP();
    }
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif
}

void Partition::exchangeUV()
{
    #ifdef TIMER
    timer_.setT0();
    #endif
    for(std::shared_ptr<NeighbourBoundary> neighbour : neighbours_)
    {
        neighbour->setUV();
    }
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif
}

void Partition::exchangeRA()
{
    #ifdef TIMER
    timer_.setT0();
    #endif
    for(std::shared_ptr<NeighbourBoundary> neighbour : neighbours_)
    {
        neighbour->setRA();
    }
    #ifdef TIMER
    timer_.addTimeSinceT0();
    #endif
}