// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/**********************************************************************//**
* FILE NAME --- SmBlendCorner.h
* PURPOSE: SolidCraft change. The blend corner: one construction for the
*   corner events at a filleted vertex, replacing the SmFilletCorner
*   hierarchy event by event. Its first event is the Turn: two blends that
*   share a face, with one sharp unblended side edge between their other
*   faces (a boss's base at its sharp corner, a pocket's rim), Mitred (the
*   default) or Swept; and the Join, two blends whose edges continue each
*   other. It builds every corner of two blends of one convexity among three
*   or four edges.
**********************************************************************/

#ifndef __SMBLENDCORNER_H__
#define __SMBLENDCORNER_H__

#ifndef __SMFILLETCORNER_H__
#include <SmFilletCorner.h>
#endif

#ifndef __SMBLENDSETBACK_H__
#include <SmBlendSetback.h>
#endif

#include <vector>

class SmTorus;

/*******************************************************************//**
PURPOSE: The corner events SmBlendCorner builds, and a Turn's two modes.
***********************************************************************/
enum SmBlendCornerEvent
{
  SM_BCE_NONE, // not an event of the blend corner (yet): the corner hierarchy builds it
  SM_BCE_TURN, // two blends turn at the unblended side edge(s) between their faces
  SM_BCE_JOIN, // two blends whose edges continue each other (G1) at the vertex
  SM_BCE_BALL, // three rounds of one radius and convexity, the common ball's sphere between them (SmBlendCornerBall.cpp)
  SM_BCE_END,    // one blend ends at the vertex, its side faces there not blended
  SM_BCE_SETBACK // three or more blends, every edge there blended, each set back and the hole filled (SmBlendSetback)
};

// How an End ends (SmBlendCornerEnd.cpp).
enum SmBlendEnd
{
  SM_BEND_ON_END_FACE, // on the end face(s) between its side edges (an end rib), capped at a mixed convexity
  SM_BEND_OPEN,        // at a lamina edge (a sheet): its end section
  SM_BEND_RUNS_OUT,    // its side faces continue each other through the vertex: its rails meet there
  SM_BEND_AT_SMOOTH    // its edge runs on tangentially into an edge not blended: cut square to the edge there
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
  SmBlendCornerEvent        m_eEvent;       // the corner's event (Classify)
  SmBlendTurn               m_eTurn;        // a Turn's mode
  SmBlendEnd                m_eEnd;         // an End's kind
  SmSurface               * m_pEndCap;      // an End at a smooth vertex: the plane square to its edge there
  SmBoolean                 m_bG1Case;      // the two blended edges meet G1 (MakeTwoBlendTopology)
  SmBoolean                 m_bMitre;       // the two blends are mitred at the common face's concave corner
  SmEdgeuse               * m_pSideEU;      // Swept: the side edge E's edgeuse the rails of the first round meet
  SmTorus                 * m_pTurnSurface; // Swept: the corner's face's surface, from CalcCornerEdgeGeom until CreateCornerPatch hands it to its face
  SmTArray<SmCurve*>        m_vTurn3d;      // Swept: its loop's 3D curves (natural iso-curves, the pole at S left out)
  SmTArray<SmBSplineCurve*> m_vTurnUV;      // and their uv lines
  SmTArray<SmOrientType>    m_vTurnOrients; // and their orientations in the loop
  // Setback: each blended edge's setback (in m_pAllEUs's order), the fill's order (1: G1, 2: G2),
  // whether the corner's edges are convex, and the fill (from CalcCornerEdgeGeom until
  // CreateCornerPatch hands its patches to faces)
  std::vector<double>       m_vSetBacks;
  ULONG                     m_lFillOrder;
  SmBoolean                 m_bConvexSetBack;
  SmBoolean                 m_bSetBackOverflows; // the rule's setbacks do not fit: the blends overflow (refused if the corner is built)
  SmBlendSetback::Fill    * m_pFill;

  // the Ball (SmBlendCornerBall.cpp)
  SmSurface               * m_pBallSurface; // the corner's face's sphere, from CalcCornerEdgeGeom until CreateCornerPatch hands it to its face
  SmTArray<SmCurve*>        m_vBall3d;      // its loop's 3D curves (natural iso-curves, the pole left out)
  SmTArray<SmBSplineCurve*> m_vBallUV;      // and their uv lines
  SmTArray<SmOrientType>    m_vBallOrients; // and their orientations in the loop

  SmBoolean IsBall() const { return ( m_eEvent == SM_BCE_BALL ) ? TRUE : FALSE; }
  SmStatus  MakeBallTopology();
  SmStatus  BuildBallSurface();  // m_pBallSurface and its loop, from the faces and the processed contacts
  void      DeleteBallSurface();
  SmStatus  BallCentre(SmPoint3d & rCentre, double & rdRadius, SmPoint3d asContact[3]) const; // the faces' triple offset point, the contacts
  SmStatus  PlaceBallContacts(); // the contacts and the rounds' sections there in closed form
  SmStatus  CreateBallPatch();

  SmBoolean IsSetBack() const { return ( m_eEvent == SM_BCE_SETBACK ) ? TRUE : FALSE; }
  SmBoolean IsSwept() const { return ( m_eEvent == SM_BCE_TURN && m_eTurn == SM_BT_SWEPT ) ? TRUE : FALSE; }
  SmBoolean IsMitreTurn() const { return ( m_eEvent == SM_BCE_TURN && m_eTurn == SM_BT_MITRE && m_pAllEUs != NULL && m_pAllEUs->GetSize() == 3 ) ? TRUE : FALSE; }

  // the Mitre turn at a corner of three edges: each round runs on past the vertex to where the two
  // meet, both cut along their crossing (the mitre), from where their rails cross on the common
  // face (X) to where they meet the side edge (S). Where the rails meet the side edge at two
  // places (unequal radii, or a slope), the larger round meets the smaller one's side face
  // between them: on the face itself at a convex face corner, on its continuation, where a patch
  // of it fills the gap, at a concave one. The mitre is the rounds' exact plane section where
  // both are cylinders of one radius cut by the plane they are symmetric about (an ellipse arc),
  // else their certified intersection.
  SmStatus MakeMitreTopology();
  SmStatus AdjustMitreTopology();
  SmStatus SplitMitreSideVertex(SmFilletVertex * pSide, SmBoolean bGapOnContinuation);
  SmStatus ExactMitreCurve(SmFilletEdge * pMitre, SmBoolean & rbExact);
  SmStatus CreateMitrePatch();

  // the Swept turn
  SmStatus MakeSweptTopology();
  SmStatus BuildSweptSurface();  // m_pTurnSurface and its loop, from the processed vertices
  void     DeleteSweptSurface();
  SmStatus CreateSweptPatch();

  // the Join: the blends end on one shared section (SmBlendCornerJoin.cpp)
  SmStatus  MakeJoinTopology();

  // the End (SmBlendCornerEnd.cpp)
  SmStatus  MakeEndOnEndFace();
  SmStatus  MakeEndOpen();
  SmStatus  MakeEndRunningOut();
  SmStatus  MakeEndAtSmoothVertex();
  SmStatus  CalcEndCapEdges();
  SmStatus  RefuseEndPastSideEdge(SmEdgeuse * pFilEdgeuse, SmEdgeuse * pSideEU, SmFace * pOtherFace);
  SmStatus  RadiusReachingSideEdgeEnd(SmEdgeuse * pFilEdgeuse, SmEdgeuse * pSideEU, SmFace * pOtherFace,
                                      SmBoolean & rbPast, double & rdFit);
  SmStatus  CalcEndEdgeGeom();
  SmStatus  CreateEndPatchOfEnd();

  // the Setback (SmBlendCornerSetback.cpp): per face at the vertex, the rail ends of the two
  // rounds beside it (the one before's SM_FV_MATE, the one after's SM_FV_SETBACK) and the
  // fill's base point between them (SM_FV_PRECOMPUTED), the two halves of the face's spring
  // curve (SM_FE_PRECOMPUTED) and each round's end section (SM_FE_CROSS_SECTION)
  SmStatus MakeSetBackTopology();
  SmStatus BuildSetBackFill();     // the fill, its base points and the spring curves, from the processed rail ends and the rounds' surfaces
  SmStatus CreateSetBackPatches(); // the fill's patches as faces
  SmStatus SetBackNotBuilt(const TCHAR * cpWhy);

  // the Join, and a corner of four edges: the two blends end on one common section (blending), or
  // run on to their intersection (bevel); moved from SmFilletNx2Corner.
  SmStatus MakeTwoBlendTopology();
  SmStatus Make4x2ConcaveCornerTopology();
  SmStatus AdjustConcaveCornerTopology();
  SmStatus AdjustTwoBlendTopology();
  SmStatus CreateTwoBlendPatch();

public:
  SmBlendCorner(const SmContext      & crContext,
                const SmVertex       * cpVertex,
                SmTArray<SmEdgeuse*> * pSolverEUs,
                SmTArray<SmEdgeuse*> * pAllEUs,
                SmFilletExecutive    * pExec,
                double                 dApproxTol3d,
                double                 dTangencyTolRadians,
                SmBlendCornerEvent     eEvent,
                SmBlendTurn            eTurn);
  virtual ~SmBlendCorner();

  // The event at a filleted vertex with two blended edgeuses of one convexity among three or
  // four (pAllEUs, the corner's ring): a Join when their edges continue each other (G1), else a
  // Turn. With SmFilletExecutive::GetSweptTurns() a Turn is Swept, refused (SM_ERR, the reason
  // noted) at a sharp concave side edge (SM_FILERR_NO_ROOM, with the radius that edge must be
  // rounded at first) and where the Swept build does not reach yet (SM_FILERR_CORNER_UNSUPPORTED).
  static SmStatus Classify(const SmVertex             * cpVertex,
                           const SmTArray<SmEdgeuse*> & crSolverEUs,
                           const SmTArray<SmEdgeuse*> & crAllEUs,
                           SmFilletExecutive          * pExec,
                           double                       dTangencyTolDeg,
                           SmBlendCornerEvent         & reEvent,
                           SmBlendTurn                & reTurn);

  // The Ball at a filleted vertex where all three edges are blended (pAllEUs, the corner's ring):
  // constant circular rounds of one radius, all convex or all concave, between three faces each a
  // plane or a natural quadric (a sphere, a cylinder or a cone; a cylinder or sphere wider than the
  // ball) whose normals at the vertex are independent, the vertex not set back
  // (SmFilletExecutive::IsSetBackVertex).
  // reEvent SM_BCE_BALL, else SM_BCE_NONE (the corner hierarchy builds it).
  static SmStatus ClassifyBall(const SmVertex             * cpVertex,
                               const SmTArray<SmEdgeuse*> & crSolverEUs,
                               const SmTArray<SmEdgeuse*> & crAllEUs,
                               SmFilletExecutive          * pExec,
                               SmBlendCornerEvent         & reEvent);
  // Whether the corner of a filleted vertex where every edge is blended (crSolverEUs, crAllEUs,
  // three or more) is a Setback, which SmBlendSetback builds: every face there a plane, every edge
  // a straight line, no two edges tangent, all of one convexity, each blend a constant round of
  // one radius (SM_FSG_CIRCULAR) or the curvature-continuous section (SM_FSG_BLEND_CURVE), and
  // either set back on request (SmFilletExecutive::GetSetBackFill), or G2 (curvature-continuous
  // sections), or rounds of unequal radii, which no one ball touches (owner, 2026-11-01). Equal
  // circular rounds stay where they were (no event: SM_BCE_NONE). The fill's order in rlOrder:
  // 2 for curvature-continuous sections or a setback with an ease, else 1.
  static SmStatus ClassifySetBack(const SmVertex             * cpVertex,
                                  const SmTArray<SmEdgeuse*> & crSolverEUs,
                                  const SmTArray<SmEdgeuse*> & crAllEUs,
                                  SmFilletExecutive          * pExec,
                                  double                       dTangencyTolDeg,
                                  SmBlendCornerEvent         & reEvent,
                                  ULONG                      & rlOrder);

  virtual SmStatus MakeCornerTopology();
  virtual SmStatus AdjustTopology();
  virtual SmStatus CalcCornerVertGeom();
  virtual SmStatus CalcCornerEdgeGeom();
  virtual SmStatus CreateCornerPatch();

  virtual SmFilletSolver * GetTangentFilletSolver( const SmFilletSolver *pFS );
  virtual SmBoolean        IsTangentCorner() const { return m_bG1Case; }
  virtual SmBoolean        IsClosedJoin() const;
  virtual SmBoolean        IsEndOnEndFace() const { return ( m_eEvent == SM_BCE_END && m_eEnd == SM_BEND_ON_END_FACE ) ? TRUE : FALSE; }

  // the End at a vertex where one blended edgeuse ends: its kind
  void SetEnd( SmBlendEnd eEnd ) { m_eEnd = eEnd; }
  // whether a round ending at cpCornerVert runs out there: its side faces continue each other
  static SmBoolean EndRunsOut( const SmVertex * pCornerVert, SmEdgeuse * pFilEdgeuse );
  // whether the blended edge runs on at cpCornerVert, tangentially, into an edge of its convexity
  // that is not blended (crAllEUs: the edgeuses round the vertex)
  static SmBoolean EndAtSmoothVertex( const SmVertex * cpCornerVert, SmEdgeuse * pFilEdgeuse,
                                      const SmTArray<SmEdgeuse*> & crAllEUs, double dTanTolDeg );
  virtual SmSurface * GetEndCapSurface() const { return m_pEndCap; }
  // an End on its end face, before it is built: the least radius at which a rail reaches its
  // side edge's far end, where the rail runs past it (the refusal's reason: DescribeEndPastSideEdge)
  SmStatus FitOfEndPastSideEdges( SmBoolean & rbFound, double & rdFit, SmEdge *& rpEdge );
  static SmStatus DescribeEndPastSideEdge( double dFit, TCHAR * cWhy );

  virtual double   GetSetBackDistAlong(const SmEdge * pEdge);
  void             SetFillOrder(ULONG lOrder) { m_lFillOrder = lOrder; }

  SmBlendCornerEvent GetEvent() const { return m_eEvent; }
  SmBlendTurn        GetTurn()  const { return m_eTurn; }

  // whether the two blends are mitred at the common face's concave corner (MakeTwoBlendTopology),
  // and the two filleted edgeuses, outward from the corner's vertex
  SmBoolean IsMitre() const { return m_bMitre; }
  SmStatus  GetFilletedEdgeuses(SmEdgeuse *& rpFirst, SmEdgeuse *& rpSecond) const;

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

// Whether the face corner at pCornerVert between two edgeuses (through the tangent edge pOptRefEU
// when not NULL) has the convexity opposite to theirs (SmFilletCorner.cpp).
// The corner helpers of the End (SmFilletCorner.cpp).
SmEdgeuse * sm_FindClosedSideEdgeuse(ULONG lRailIndex, SmFilletSolver * pFilSolver, SmEdgeuse * pEU, int iDebugLevel);
SmBoolean   sm_IsTangentSeamAtVertex(SmEdgeuse * pFilEdgeuse, SmEdgeuse * pSideEU, const SmVertex * pVertex, double dTanTolDeg);
SmBoolean   same_convexity_edges(SmEdgeuse * pEU1, SmEdgeuse * pEU2, const SmVertex * pVertex, int iDebugLevel);
SmBoolean   same_convexity_edges_round_vertex(SmEdgeuse * pFilEdgeuse, SmEdge * pEdge);
double      sm_FindSurfaceExtensionDistAtCorner(SmEdgeuse * pFilEdgeuse, SmFilletSolver * pFilSolver, const SmVertex * pVertex);

SmBoolean sm_CheckFaceComplimentaryCornerConcavity(const SmVertex * pCornerVert, SmEdgeuse * pTestEU1, SmEdgeuse * pTestEU2,
                                                   SmEdgeuse * pOptRefEU, int iDebugLevel);

#endif // __SMBLENDCORNER_H__
