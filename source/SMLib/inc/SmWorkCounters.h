// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmWorkCounters.h
* PURPOSE: (SolidCraft fork) Process-wide counts of the work done by the
*          steps of a boolean whose cost once grew with the square of the
*          body's size, so a test can check how that work grows on any
*          machine, where a ratio of times depends on the machine.
**********************************************************************/

#ifndef __SMWORKCOUNTERS_H__
#define __SMWORKCOUNTERS_H__

#include <atomic>
#include <cstdint>

/*******************************************************************//**
PURPOSE: Three counters, each added to where its work is done:
         - eClosureFaceuses: faceuses a new face's closure test visits
           (SmBrep::TestClosureSplitRegion, through
           SmTopologyTraverser::ClosureTraversalFromBothSides);
         - eFacePairs: pairs of faces of the two bodies examined for
           intersection (SmTopologyIntersector::IntersectTopology);
         - eShellFaceusesMoved: faceuses moved from one shell to another when
           a new edge joins two shells (SmBrep::MakeEdgeTopology).

NOTES: Process-wide, not per thread or per Brep. Atomic and relaxed, so
       concurrent calls never race, but the counts of concurrent calls mix.
***********************************************************************/
class SmWorkCounters
{
public:
  enum Counter
  {
    eClosureFaceuses    = 0,
    eFacePairs          = 1,
    eShellFaceusesMoved = 2,
    eNumCounters        = 3
  };

  static void Add(Counter eCounter, int64_t lCount = 1)
  {
    s_aCounts[eCounter].fetch_add(lCount, std::memory_order_relaxed);
  }

  static int64_t Get(Counter eCounter)
  {
    return s_aCounts[eCounter].load(std::memory_order_relaxed);
  }

  static void Reset()
  {
    for (int ii = 0; ii < eNumCounters; ii++)
      { s_aCounts[ii].store(0, std::memory_order_relaxed); }
  }

private:
  static inline std::atomic<int64_t> s_aCounts[eNumCounters] = {};
};

#endif // __SMWORKCOUNTERS_H__
