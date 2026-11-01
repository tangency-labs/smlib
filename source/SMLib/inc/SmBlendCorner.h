// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCorner.h
* PURPOSE: SolidCraft change. The blend corner: one construction for the
*   corner events at a filleted vertex, replacing the SmFilletCorner
*   hierarchy event by event. Its first event is the Turn: two blends that
*   share a face, with one sharp unblended side edge between their other
*   faces (a boss's base at its sharp corner, a pocket's rim). A Turn is
*   Mitred (the default; still SmFilletNx2Corner's) or Swept (here).
**********************************************************************/

#ifndef __SMBLENDCORNER_H__
#define __SMBLENDCORNER_H__

#ifndef __SMFILLETCORNER_H__
#include <SmFilletCorner.h>
#endif

class SmTorus;

/*******************************************************************//**
PURPOSE: The corner events SmBlendCorner builds, and a Turn's two modes.
***********************************************************************/
enum SmBlendCornerEvent
{
  SM_BCE_NONE, // not an event of the blend corner (yet): the corner hierarchy builds it
  SM_BCE_TURN  // two blends sharing a face turn at a sharp side edge between their other faces
};

enum SmBlendTurn
{
  SM_BT_MITRE, // each blend run on to their intersection (the default; G0)
  SM_BT_SWEPT  // the blends' section revolved about the side edge (G1)
};

/*******************************************************************//**
PURPOSE: A Swept turn at a sharp convex side edge E.

NOTES: Two rounds, R1 and R2, share the common face F0; E lies between
   their other faces. Each round ends on its own section in the plane
   through E square to its edge: the side vertex S, where both rounds' rails
   on their other faces meet E, and one plate vertex per round (A, B), where
   its rail on F0 meets that section (S's mates). Corner edges: the two end
   sections S-A and B-S (each round's own cross-section) and the arc A-B,
   where the corner's face meets F0. That face is the end section revolved
   about E through the corner's turn: where E is a straight line in both
   end sections' planes and the radii are equal, an exact SmTorus (centre on
   E, major radius the ball centre's distance from E, minor radius r; a horn
   torus at S when the rounds' walls hold E), tangent to both rounds along
   their end sections and to F0 along A-B. A sharp concave side edge is
   refused (the ball cannot pass it), naming the radius it must be rounded
   at first; other side edges, drafted walls and unequal radii are not
   built yet (refused, SM_FILERR_CORNER_UNSUPPORTED).
***********************************************************************/
class SM_EXPORT SmBlendCorner : public SmFilletCorner
{
protected:
  SmBlendTurn               m_eTurn;        // the turn's mode (SM_BT_SWEPT: only that is built here so far)
  SmEdgeuse               * m_pSideEU;      // the side edge E's edgeuse the rails of the first round meet
  SmTorus                 * m_pTurnSurface; // the corner's face's surface, from CalcCornerEdgeGeom until CreateCornerPatch hands it to its face
  SmTArray<SmCurve*>        m_vTurn3d;      // its loop's 3D curves (natural iso-curves, the pole at S left out)
  SmTArray<SmBSplineCurve*> m_vTurnUV;      // and their uv lines
  SmTArray<SmOrientType>    m_vTurnOrients; // and their orientations in the loop

  SmStatus BuildSweptSurface();  // m_pTurnSurface and its loop, from the processed vertices
  void     DeleteSweptSurface();

public:
  SmBlendCorner(const SmContext      & crContext,
                const SmVertex       * cpVertex,
                SmTArray<SmEdgeuse*> * pSolverEUs,
                SmTArray<SmEdgeuse*> * pAllEUs,
                SmFilletExecutive    * pExec,
                double                 dApproxTol3d,
                double                 dTangencyTolRadians,
                SmBlendTurn            eTurn);
  virtual ~SmBlendCorner();

  // The event at a filleted vertex with two blended edgeuses among three (pAllEUs, the corner's
  // ring), the blends of one convexity: a Turn when their edges do not continue each other. With
  // SmFilletExecutive::GetSweptTurns() the Turn is Swept, and a sharp concave side edge is refused
  // (SM_ERR, SM_FILERR_NO_ROOM noted with the radius that edge must be rounded at first).
  static SmStatus Classify(const SmVertex             * cpVertex,
                           const SmTArray<SmEdgeuse*> & crSolverEUs,
                           const SmTArray<SmEdgeuse*> & crAllEUs,
                           SmFilletExecutive          * pExec,
                           double                       dTangencyTolDeg,
                           SmBlendCornerEvent         & reEvent,
                           SmBlendTurn                & reTurn);

  virtual SmStatus MakeCornerTopology();
  virtual SmStatus CalcCornerEdgeGeom();
  virtual SmStatus CreateCornerPatch();

  SmBlendTurn GetTurn() const { return m_eTurn; }

  virtual SmFilletCornerType GetCornerType() const { return SM_FCR_BLEND; }
  virtual void GetName(TCHAR * pcMyName, size_t lMyNameAllocLen) const { smos_WStrCpy(pcMyName, lMyNameAllocLen, _T("SmBlendCorner")); }

} ; // end class SmBlendCorner

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
