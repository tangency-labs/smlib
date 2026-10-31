// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmOffsetGeometryCreateion.h
* PURPOSE: Header file for the class. 
**********************************************************************/

#ifndef __SMOFFSETGEOMETRYCREATION_H__
#define __SMOFFSETGEOMETRYCREATION_H__


#ifndef __SMOS_TYPES_H__
#include <SmTypes.h>
#endif

#ifndef __SMTOPO_TYPES_H__
#include <SmTopoTypes.h>
#endif

#ifndef __SMCURV_TYPES_H__
#include <SmCurveTypes.h>
#endif

#ifndef __SMSURF_TYPES_H__
#include <SmSurfTypes.h>
#endif

#ifndef __SMTARRAY_H__
#include <SmTArray.h>
#endif

#ifndef __SMFACE_H__
#include <SmFace.h>
#endif

#ifndef __SMSURFACE_H__
#include <SmSurface.h>
#endif

#include <SmSSIData.h>

class SmOffsetExecutive;


/*******************************************************************//**
PURPOSE: The following class is used by the offset classs to
            generate geometry for the offset.  This class can be
            subclassed to create different kinds of offsets.

NOTES: 
***********************************************************************/
class SM_EXPORT SmOffsetGeometryCreation
{
protected:
    double    m_dOffsetDistance;       //
    double    m_dThisApproxTol3d;     //
    // SolidCraft change: half a degree, not two. Faces meeting at a crease up to two degrees were
    // taken as tangent, their offsets left to meet by extending their trims instead of being
    // intersected: the four facets of a floor falling 1 degree to a drain, which meet at 1.2
    // degrees, then had no edge where their offsets cross. Half a degree is SolidCraft's one
    // threshold for a smooth edge (its shaded view, tangent chains, blends).
    double    m_dTangencyTolRadians;   //
    SmBoolean m_bMakeAnalytics;        // TRUE =
                                       // FALSE=

public:

    SmOffsetGeometryCreation
      (double dOffsetDistance, 
       double dThisApproxTol3d) 
     : m_dOffsetDistance(dOffsetDistance), 
       m_dThisApproxTol3d(dThisApproxTol3d), 
       m_dTangencyTolRadians(0.5*SM_PI/180.0),
       m_bMakeAnalytics(TRUE) 
       { }

    virtual ~SmOffsetGeometryCreation();

    double GetThisApproxTol3d()    const { SM_ASSERT_TOL(m_dThisApproxTol3d) ; return m_dThisApproxTol3d; }
    double GetTangencyTolRadians() const { SM_ASSERT_TOL(m_dTangencyTolRadians) ; return m_dTangencyTolRadians; }
    double GetTangencyTolDegrees() const { SM_ASSERT_TOL(m_dTangencyTolRadians) ; return m_dTangencyTolRadians * 180.0 / SM_PI; }
    double GetOffsetDistance()     const { return m_dOffsetDistance; }

    SmStatus ReplaceEdgeCurves(SmEdgeuse        * pEdgeuse,          // in : when rp3DCurve is a new shape, 
                                                                     //           - assigned rp3DCurve and rpUVCurve deleting any 
                                                                     //               preExisting Curve and UVTrimCurves.
                                                                     //      else - not changed.
                               SmCurve         *& rp3DCurve,         // i/o: candidate 3dCurve for Edgeuse,     deleted and set to NULL when not stored in pEdgeuse
                               SmBSplineCurve  *& rpUVCurve,         // i/o: candidate UVTrimCurve for Edgeuse, deleted and set to NULL when not stored in pEdgeuse
                               const SmEdgeuse  * cpOptOrigEdgeuse,  // in : optional Edgeuse used to compute orientation of rp3DCurve to Edgeuse, NULL to ignore
                               SmBoolean & rbDelete3DCurve,          // o: if true, then caller should delete 3DCurve.
                               SmBoolean & rbDeleteUVCurve);         // o: if true, then caller should delete UVCurve.

    SmStatus MateEdges(const SmContext   & crContext,
                       SmOffsetExecutive * cpOffsetExecutive,
                       const SmEdgeuse   * pEdgeuse,
                       const SmEdgeuse   * pOtherEdgeuse,
                       SmEdgeuse         * pOffsetEdgeuse,
                       SmEdgeuse         * pOtherOffsetEdgeuse,
                       SmMarkType          eMarkType) ;     // in : uses without incrementing eMarkType value

    void SetMakeAnalytics(SmBoolean bMakeAnalytics) 
    { m_bMakeAnalytics = bMakeAnalytics; }

    // SolidCraft: each edge of pBrep whose curve is a B-spline between two faces on analytic
    // surfaces made the exact line or conic those surfaces meet on, between the same ends.
    static SmStatus MakeEdgesOnAnalyticsExact(const SmContext & crContext,  // in : context for new obj construction
                                              SmBrep          * pBrep);     // in : target Brep to modify

    SmStatus UpdateFaceSurfaceDomain(SmFace * pFace,
                                     SmSurface * pExtendedSurface);

    virtual SmStatus VertexOffset(const SmContext         & crContext,
                                  const SmVertex          * pOriginalVertex,
                                  const SmTArray<SmEdge*> & rBoundingEdges,
                                  SmTArray<SmCurve*>      & r3DTrimmingCurves,
                                  SmTArray<SmOrientType>  & rOrients,
                                  SmSurface              *& rpNewSurface);

    virtual SmStatus EdgeuseOffset(const SmContext        & crContext,
                                  SmOffsetExecutive       * cpOffsetExecutive,
                                  const SmEdgeuse         * cpOriginalEdgeuse,
                                  const SmEdgeuse         * cpOtherOriginalEdgeuse,
                                  SmEdgeuse               * cpOffsetEdgeuse,
                                  SmEdgeuse               * cpRadialOffsetEdgeuse,
                                  SmSurface              *& rpNewSurface,
                                  SmTArray<SmCurve*>      & r3DTrimmingCurves,
                                  SmTArray<SmOrientType>  & rTrimOrientations,
                                  SmSSIData               & rSelfIntersections,
                                  SmMarkType                eMarkType,               // in : uses without increment eMarkType Value
                                  SmBoolean               * pbOptClosingWall = NULL) ; // out: TRUE = rpNewSurface is a step wall: it spans the gap between
                                                                                       //      two edges and is never extended (as a lamina side wall)

    // create offset surface - Approximate When necessary - Out Surf->Domain(s) may be trimmed but not scaled
    virtual SmStatus FaceuseOffset
    (
      const SmContext       & crContext,                  // in : context for new obj construction
      const SmFaceuse       * pOriginalFaceuse,           // in : Faceuse->Surface to copy and offset
      SmBoolean               bSkipSelfIntersection,      // in : TRUE =
                                                          //      FALSE=
      SmSurface * &           rNewSurface,                // out: New Surface; NULL when, and only when, rbCollapses
      SmSSIData             & rSelfIntersections,         // out: List of SelfXSect curves (when bSkipSelfIntersection == FALSE)
      SmBoolean             & rbCollapses                 // out: SolidCraft: TRUE = the offset shrinks the face's surface onto its spine
                                                          //      or past it (SmSurface::OffsetCollapses), and there is no offset surface
    ) const;
                                                         
    SmStatus TrimExtendedEdge(const SmContext   & crContext,            // in : context for new object construction 
                              SmOffsetExecutive * cpOffsetExecutive,    // in : offset operation context data
                              const SmEdgeuse   * pBaseEdgeuse,         // in : 1st target edgeuse
                              const SmEdgeuse   * pOtherBaseEdgeuse,    // in : 2nd target edgeuse
                              SmEdgeuse         * pOffsetEdgeuse,       // in : Offset of 1st target edgeuse
                              SmEdgeuse         * pOtherOffsetEdgeuse,  // NotUsed: in : Offset of 2nd target edgeuse
                              SmCurve           * p3DCurve,             // in : XSectCurve between edgeuse->Face->ExtendedSurfaces
                              SmBSplineCurve    * pUVCurveBase,         // in : UVTrimCurve on 1st target edgeuse->Face->ExtendedSurface
                              SmBSplineCurve    * pUVCurveOtherBase,    // in : UVTrimCurve on 2nd target edgeuse->Face->ExtendedSurface 
                              SmBoolean         & bGoodEdge,            // out: TRUE = TrimPoints were found for the input 3DCurve
                              SmExtent1d        & rTrimInterval,        // out: TrimInterval for the 3DCurve
                              SmMarkType          eMarkType) ;          // in : uses without increment eMarkType value

    SmStatus ExtendTrimEdges(const SmContext   & crContext,
                             SmOffsetExecutive * cpOffsetExecutive,
                             const SmEdgeuse   * pEdgeuse,
                             const SmEdgeuse   * pOtherEdgeuse,
                             SmEdgeuse         * pOffsetEdgeuse,
                             SmEdgeuse         * pOtherOffsetEdgeuse,
                             SmMarkType          eMarkType) ;       // in : uses without incrementing eMarkType value

    // TRUE when the faces of two offset edgeuses of a tangent edge (only one of them moved)
    // stayed parallel across it but apart: everywhere along the edge they lie apart along both
    // faces' normals, so no extension makes them meet and a wall between the edges closes them.
    SmStatus IsStepAcrossEdge(SmEdgeuse * pOffsetEdgeuse,
                              SmEdgeuse * pOtherOffsetEdgeuse,
                              SmBoolean & rbStep) const;

    // The ruled wall between two offset edgeuses that lie apart (IsStepAcrossEdge).
    SmStatus CreateStepWall(const SmContext & crContext,
                            SmEdgeuse       * pOffsetEdgeuse,
                            SmEdgeuse       * pOtherOffsetEdgeuse,
                            SmSurface      *& rpWall) const;

    SmBoolean HaveTangentGaps(ULONG lNumSamples,
                              SmEdgeuse * cpOffsetEdgeuse,
                              SmEdgeuse * cpOtherOffsetEdgeuse);


} ; // end class SmOffsetGeometryCreation

/*******************************************************************//**
PURPOSE: 

NOTES: 
***********************************************************************/
class SM_EXPORT SmLocalOperationGeometryCreation : public SmOffsetGeometryCreation
{
public:
  SmTArray<SmFace*>    m_vFaces;         //      
  SmTArray<SmSurface*> m_vSurfaces;      //      

public:
  SmLocalOperationGeometryCreation
    (double dOffsetDistance, 
     double dThisApproxTol3d,
     SmTArray<SmFace*> & rFaces, 
     SmTArray<SmSurface*> & rSurfaces) 
   : SmOffsetGeometryCreation(dOffsetDistance,dThisApproxTol3d)
    {
      m_bMakeAnalytics = FALSE ;
      m_vFaces.Append(rFaces) ; 
      m_vSurfaces.Append(rSurfaces) ;
    }
  
  // create offset surface - Approximate When necessary - Out Surf->Domain(s) may be trimmed but not scaled  
  virtual SmStatus FaceuseOffset(const SmContext       & crContext,                   // in : context for new obj construction
                                 const SmFaceuse       * pOriginalFaceuse,            // in : Faceuse->Surface to copy and offset
                                 SmBoolean               bSkipSelfIntersection,       // NotUsed: in : TRUE =
                                                                                      //      FALSE=
                                 SmSurface*            & rNewSurface,                 // out: New Surface
                                 SmSSIData             & rSelfIntersections,          // out: List of SelfXSect curves (when bSkipSelfIntersection == FALSE)
                                 SmBoolean             & rbCollapses) const ;         // out: always FALSE: the surfaces are given

} ; // end class SmLocalOperationGeometryCreation


#endif // !__SMOFFSETGEOMETRYCREATION_H__
