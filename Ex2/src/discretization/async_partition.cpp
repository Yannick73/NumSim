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
        if(pi.ownPartitionContainsBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletSouth>(discretization, settings.dirichletBcBottom));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourSouth>(discretization, pi.bottomNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));

        if(pi.ownPartitionContainsTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletNorth>(discretization, settings.dirichletBcTop));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourNorth>(discretization, pi.topNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));
    }
    else
    {
        if(pi.ownPartitionContainsTopBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletNorth>(discretization, settings.dirichletBcTop));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourNorth>(discretization, pi.topNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));

        if(pi.ownPartitionContainsBottomBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletSouth>(discretization, settings.dirichletBcBottom));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourSouth>(discretization, pi.bottomNeighbourRankNo(),
                                !pi.ownPartitionContainsLeftBoundary(), !pi.ownPartitionContainsRightBoundary()));
    }

    if(pi.getPartPosY() & 0b1)
    {
        if(pi.ownPartitionContainsLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletWest> (discretization, settings.dirichletBcLeft));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourWest> (discretization, pi.leftNeighbourRankNo()));

        if(pi.ownPartitionContainsRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletEast> (discretization, settings.dirichletBcRight));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourEast> (discretization, pi.rightNeighbourRankNo()));
    }
    else
    {
        if(pi.ownPartitionContainsRightBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletEast> (discretization, settings.dirichletBcRight));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourEast> (discretization, pi.rightNeighbourRankNo()));

        if(pi.ownPartitionContainsLeftBoundary())
            fixBoundaries_.push_back(std::make_shared<DirichletWest> (discretization, settings.dirichletBcLeft));
        else
            asyncNeighbours_.push_back(std::make_shared<AsyncNeighbourWest> (discretization, pi.leftNeighbourRankNo()));
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
    assert(directions.size()  == 4);
}

// First this method sets up async receive, sends its data
// then it sets the dirichlet boundary, finally setting the first incoming ghost data
void AsyncPartition::setBoundaryUV()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeUV);
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setUV();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvUV);
}

void AsyncPartition::setBoundaryFG()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeFG);
    for(std::shared_ptr<Dirichlet> fixBoundary : fixBoundaries_)
    {
        fixBoundary->setFG();
    }
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvFG);
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

void AsyncPartition::exchangeUV()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeUV);
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvUV);
}

void AsyncPartition::exchangeRA()
{
    std::vector<std::shared_ptr<AsyncNeighbourBoundary>> neighbourRecvQueue;
    setupExchange(neighbourRecvQueue, &AsyncNeighbourBoundary::exchangeRA);
    setFirstIncomingData(neighbourRecvQueue, &AsyncNeighbourBoundary::setRecvRA);
}