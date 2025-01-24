#pragma once

#include <vector>
#include <cassert>
#include <mpi.h>

// Maybe use something like the data-enum in neighbour as MPI-tag?

//! minimalistic class to wrap the MPI communication
class MPI_Wrapper
{
public:
    MPI_Wrapper(int neighbourRank) : neighbourRank_(neighbourRank) { }

    //! setups the receive request, it requires a valid receive buffer
    inline void startReceive(std::vector<double> &recvBuf, int len)
    {
        // std::cout << "startReceive handler " << this << std::endl;
        #ifndef NDEBUG
        assert(recvBuf.size() >= len);
        #endif
        MPI_Irecv(recvBuf.data(), len, MPI_DOUBLE, neighbourRank_, 0, MPI_COMM_WORLD, &recvRequest_);
        expectedRecvLen_ = len;
    }

    //! wrapper for sending
    inline void send(std::vector<double> &sendBuf, int len)
    {
        #ifndef NDEBUG
        assert(sendBuf.size() >= len);
        #endif
        MPI_Request sendReq;    // request handle is ignored
        MPI_Isend(sendBuf.data(), len, MPI_DOUBLE, neighbourRank_, 0, MPI_COMM_WORLD, &sendReq);
    }

    //! checks, if the asynchronoues receive is complete 
    //! w/o destroying the request handle
    inline bool queryRecvComplete()
    {
        // std::cout << "queryRecvComplete handler " << this << std::endl;
        MPI_Status recvStatus;
        int requestStatusComplete;
        MPI_Request_get_status(recvRequest_, &requestStatusComplete, &recvStatus);
        int count;
        MPI_Get_count(&recvStatus, MPI_DOUBLE, &count);
        return requestStatusComplete && count == expectedRecvLen_;
    }

    //! blocking wait for the receive to finish, destroys the handle
    inline void waitForRecvComplete()
    {
        MPI_Wait(&recvRequest_, MPI_STATUS_IGNORE);
    }

protected:
    //! Rank to communicate with
    const int neighbourRank_;

    //! used to check, if a transaction is complete
    MPI_Request recvRequest_;
    //! length of last transaction
    int expectedRecvLen_;
};