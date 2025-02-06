#include "discretization/async_partition.h"

std::shared_ptr<DomainBoundary> createTopBoundary(const std::shared_ptr<Discretization> discretization,
                                                  const Settings &settings)
{
    if(settings.boundaryTop == "inflow")
        return std::make_shared<DirichletTop>(discretization, settings.dirichletBcTop);
    else if(settings.boundaryTop == "pressure")
        return std::make_shared<PressureTop>(discretization, settings.pressureTop);
    else if(settings.boundaryTop == "slip")
        return std::make_shared<SlipTop>(discretization);
    else if(settings.boundaryTop == "outflow")
        return std::make_shared<OutflowTop>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryTop + "\" for top boundary");
}

std::shared_ptr<DomainBoundary> createBottomBoundary(const std::shared_ptr<Discretization> discretization,
                                                     const Settings &settings)
{
    if(settings.boundaryBottom == "inflow")
        return std::make_shared<DirichletBottom>(discretization, settings.dirichletBcBottom);
    else if(settings.boundaryBottom == "pressure")
        return std::make_shared<PressureBottom>(discretization, settings.pressureBottom);
    else if(settings.boundaryBottom == "slip")
        return std::make_shared<SlipBottom>(discretization);
    else if(settings.boundaryBottom == "outflow")
        return std::make_shared<OutflowBottom>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryBottom + "\" for bottom boundary");
}

std::shared_ptr<DomainBoundary> createLeftBoundary(const std::shared_ptr<Discretization> discretization,
                                                   const Settings &settings)
{
    if(settings.boundaryLeft == "inflow")
        return std::make_shared<DirichletLeft>(discretization, settings.dirichletBcLeft);
    else if(settings.boundaryLeft == "pressure")
        return std::make_shared<PressureLeft>(discretization, settings.pressureLeft);
    else if(settings.boundaryLeft == "slip")
        return std::make_shared<SlipLeft>(discretization);
    else if(settings.boundaryLeft == "outflow")
        return std::make_shared<OutflowLeft>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryLeft + "\" for left boundary");
}

std::shared_ptr<DomainBoundary> createRightBoundary(const std::shared_ptr<Discretization> discretization,
                                                    const Settings &settings)
{
    if(settings.boundaryRight == "inflow")
        return std::make_shared<DirichletRight>(discretization, settings.dirichletBcRight);
    else if(settings.boundaryRight == "pressure")
        return std::make_shared<PressureRight>(discretization, settings.pressureRight);
    else if(settings.boundaryRight == "slip")
        return std::make_shared<SlipRight>(discretization);
    else if(settings.boundaryRight == "outflow")
        return std::make_shared<OutflowRight>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryRight + "\" for right boundary");
}

std::shared_ptr<DomainBoundary> createHindBoundary(const std::shared_ptr<Discretization> discretization,
                                                   const Settings &settings)
{
    if(settings.boundaryHind == "inflow")
        return std::make_shared<DirichletHind>(discretization, settings.dirichletBcHind);
    else if(settings.boundaryHind == "pressure")
        return std::make_shared<PressureHind>(discretization, settings.pressureHind);
    else if(settings.boundaryHind == "slip")
        return std::make_shared<SlipHind>(discretization);
    else if(settings.boundaryHind == "outflow")
        return std::make_shared<OutflowHind>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryHind + "\" for hind boundary");
}

std::shared_ptr<DomainBoundary> createFrontBoundary(const std::shared_ptr<Discretization> discretization,
                                                    const Settings &settings)
{
    if(settings.boundaryFront == "inflow")
        return std::make_shared<DirichletFront>(discretization, settings.dirichletBcFront);
    else if(settings.boundaryFront == "pressure")
        return std::make_shared<PressureFront>(discretization, settings.pressureFront);
    else if(settings.boundaryFront == "slip")
        return std::make_shared<SlipFront>(discretization);
    else if(settings.boundaryFront == "outflow")
        return std::make_shared<OutflowFront>(discretization);
    else
        throw std::runtime_error("Unknown boundary type \"" + settings.boundaryFront + "\" for front boundary");
}

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
            domainBoundaries_.push_back(createBottomBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourBottom>(discretization, pi.bottomRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));

        if(pi.ownTopBoundary())
            domainBoundaries_.push_back(createTopBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourTop>(discretization, pi.topRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));
    }
    else
    {
        if(pi.ownTopBoundary())
            domainBoundaries_.push_back(createTopBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourTop>(discretization, pi.topRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));

        if(pi.ownBottomBoundary())
            domainBoundaries_.push_back(createBottomBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourBottom>(discretization, pi.bottomRank(),
                                !pi.ownLeftBoundary(), !pi.ownRightBoundary(), !pi.ownHindBoundary(), !pi.ownFrontBoundary()));
    }

    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownLeftBoundary())
            domainBoundaries_.push_back(createLeftBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourLeft> (discretization, pi.leftRank()));

        if(pi.ownRightBoundary())
            domainBoundaries_.push_back(createRightBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourRight> (discretization, pi.rightRank()));
    }
    else
    {
        if(pi.ownRightBoundary())
            domainBoundaries_.push_back(createRightBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourRight> (discretization, pi.rightRank()));

        if(pi.ownLeftBoundary())
            domainBoundaries_.push_back(createLeftBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourLeft> (discretization, pi.leftRank()));
    }

    if(pi.getPartPosZ() & 0b1)
    {
        if(pi.ownFrontBoundary())
            domainBoundaries_.push_back(createFrontBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourFront> (discretization, pi.frontRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));

        if(pi.ownHindBoundary())
            domainBoundaries_.push_back(createHindBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourHind> (discretization, pi.hindRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));
    }
    else
    {
        if(pi.ownHindBoundary())
            domainBoundaries_.push_back(createHindBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourHind> (discretization, pi.hindRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));

        if(pi.ownFrontBoundary())
            domainBoundaries_.push_back(createFrontBoundary(discretization, settings));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourFront> (discretization, pi.frontRank(), 
                                        !pi.ownLeftBoundary(), !pi.ownRightBoundary()));
    }

    // checks, that each direction exist, while the dirichlet constraints take priority
    // also puts every neighbour in a common boundaries object
    std::set<BoundaryEdge> directions;
    for(std::shared_ptr<DomainBoundary> boundary : domainBoundaries_)
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
    for(std::shared_ptr<DomainBoundary> domainBoundary : domainBoundaries_)
    {
        domainBoundary->setUVW();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvUVW);
}

void AsyncPartition::setBoundaryFGH()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeFGH);
    for(std::shared_ptr<DomainBoundary> domainBoundary : domainBoundaries_)
    {
        domainBoundary->setFGH();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvFGH);
}

void AsyncPartition::setBoundaryP()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeP);
    for(std::shared_ptr<DomainBoundary> domainBoundary : domainBoundaries_)
    {
        domainBoundary->setP();
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