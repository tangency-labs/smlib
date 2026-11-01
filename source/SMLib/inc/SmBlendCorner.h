// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCorner.h
* PURPOSE: SolidCraft change. The blend corner: one construction for the
*   corner events at a filleted vertex, replacing the SmFilletCorner
*   hierarchy event by event. First, the corner helpers both share.
**********************************************************************/

#ifndef __SMBLENDCORNER_H__
#define __SMBLENDCORNER_H__

#ifndef __SMFILLETCORNER_H__
#include <SmFilletCorner.h>
#endif

/*******************************************************************//**
    The corner helpers the blend corner and the corner hierarchy share
    (moved from SmFilletCorner.cpp, where they were file statics).
***********************************************************************/

// Of a side edgeuse at a corner, itself or its radial: the one whose face also holds
// pFilEdgeuse's edge (the one the solver's rail/side-edge solve takes). NULL when neither.
SmEdgeuse * sm_FindSideEdgeuse(SmEdgeuse * pFilEdgeuse, SmEdgeuse * pSideEU, int iDebugLevel);

// Whether two edges at a vertex are tangent there (end tangents within dAngTolDeg); the angle
// between them (degrees, folded to [0, 90]) in *pdAngDeg when not NULL.
SmBoolean sm_TestG1Continuity(const SmVertex * pVert, const SmEdge * pEdge1, const SmEdge * pEdge2,
                              double dAngTolDeg, double * pdAngDeg, int iDebugLevel);

// The face the sectors of two edgeuses share, not pEU3's face (when not NULL); NULL when none.
SmFace * sm_GetFaceByTwoEdgeuses(SmEdgeuse * pEU1, SmEdgeuse * pEU2, SmEdgeuse * pEU3);

#endif // __SMBLENDCORNER_H__
