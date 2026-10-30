// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Sheets (include/smlib.h, "Sheets"): sewing a body's faces along the edges that coincide, and a
// sheet's offset along its faces' normals.
//
// smlib_stitch glues coincident vertices and edges with SmBrep::StitchFaces at the caller's
// tolerance, once, and then orients the faces with SmBrep::OrientTrimmedSurfaces: outward where
// the faces close into a solid, else only made to agree with each other across the edges they
// share (bJustUnifyNormals, which needs no ray firing). SmBrep::StitchAndOrient is not used: it
// caps the tolerance at 100 zone tolerances and stitches a second time at a larger one, so the
// caller's tolerance would not be the one that decides.
//
// StitchFaces glues an edge only to an edge one face uses (SmStitch::m_bMakingManifoldSolid): a
// sheet laid on edges two faces of a solid already share (a cover over a box's cut top, its rim
// on the cut's edges) was left beside the solid, unglued. Edges still lamina after that pass are
// glued at the same tolerance to the edges they lie on whatever uses those (bAllowSpineEdges), and
// a gluing that closes a set of faces makes the space they close a region of its own
// (SmBrep::GlueEdgesGeneral's region nesting: a new region is material). A face then lying
// between two material regions (the cut face under the cover) is inside the part and is removed,
// the two regions joined: the sewn faces enclose the space as one solid. A sheet glued so that a
// solid's inside and outside become one region (one that runs on along the solid's own faces
// from those edges, or overlaps them: no sector between the faces to put it in) is refused.
//
// smlib_offset_sheet is SM_API's offset (SmApiOffsetBrepFull, through
// SmOffsetExecutive::OffsetBrep) of a sheet, or of faces that close on themselves (a sphere's one
// face): every face onto its surface's offset along its normal, as smlib_thicken moves a sheet
// back by half before a symmetric thicken.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiGeneral.h>
#include <SmApiQueries.h>
#include <SmRegion.h>
#include <SmShell.h>

#include <memory>

using namespace smlib_shim;

namespace
{

using Owned = std::unique_ptr<SmBrep>;

Owned copy_of(const smlib_brep *b, int32_t &s)
{
    SmBrep *r = nullptr;
    s = status(SmApiBrepCopy(b->p, r));
    if (s != SMLIB_OK)
    {
        delete r;
        return nullptr;
    }
    if (r == nullptr)
    {
        s = SMLIB_ERR_EMPTY;
    }
    return Owned(r);
}

// How many regions of `b` are material: neither the infinite region nor a void.
ULONG material_regions(SmBrep &b)
{
    SmTArray<SmRegion *> regions;
    b.GetRegions(regions);
    ULONG n = 0;
    for (ULONG i = 0; i < regions.GetSize(); ++i)
    {
        n += regions[i] != b.GetInfiniteRegion() && !regions[i]->IsVoid() ? 1 : 0;
    }
    return n;
}

// Remove every face of `b` lying between two material regions (neither the infinite region nor a
// void): it is inside the part. Its regions are joined (SmBrep::DeleteFace's region nesting), and
// its edges and vertices nothing else uses go with it.
int32_t remove_faces_inside(SmBrep &b)
{
    SmTArray<SmFace *> faces;
    b.GetFaces(faces);
    const SmRegion *infinite = b.GetInfiniteRegion();
    std::vector<SmFace *> inside;
    for (ULONG i = 0; i < faces.GetSize(); ++i)
    {
        SmFaceuse *u1 = nullptr, *u2 = nullptr;
        faces[i]->GetFaceuses(u1, u2);
        if (u1 == nullptr || u2 == nullptr || u1->GetShell() == nullptr || u2->GetShell() == nullptr)
        {
            continue;
        }
        const SmRegion *r1 = u1->GetShell()->GetRegion(), *r2 = u2->GetShell()->GetRegion();
        if (r1 != nullptr && r2 != nullptr && r1 != r2 && r1 != infinite && r2 != infinite && !r1->IsVoid() && !r2->IsVoid())
        {
            inside.push_back(faces[i]);
        }
    }
    // A face is removed only with the Brep open for editing (SmBrep::DeleteFace refuses
    // otherwise).
    SmTemporaryChangeValue<SmBoolean> editing(b.m_bEditingEnabled, TRUE);
    for (SmFace *f : inside)
    {
        int32_t s = status(b.DeleteFace(f, TRUE, TRUE));
        if (s != SMLIB_OK)
        {
            return s;
        }
    }
    return SMLIB_OK;
}

// Whether an edge of `b` is used by more than two faces.
bool has_spine_edge(SmBrep &b)
{
    SmTArray<SmEdge *> edges;
    b.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i)
    {
        if (edges[i]->IsSpine())
        {
            return true;
        }
    }
    return false;
}

// `w`, stitched as a manifold, with `lamina` edges left free: those glued at `tolerance` to the
// edges they lie on whatever uses those, on a copy, which replaces `w` when the space they close
// is the solid's and no edge is left with more than two faces (the faces inside removed). Sheets
// that close no space with the solid's faces (a sheet hung on a solid's edge) stay as the manifold
// stitch left them. A gluing that leaves fewer material regions than it found is refused (see
// the notes above).
int32_t close_onto_shared_edges(Owned &w, double tolerance, ULONG &lamina)
{
    SmBrep *raw = nullptr;
    int32_t s = status(SmApiBrepCopy(w.get(), raw));
    Owned c(raw);
    if (s != SMLIB_OK || c == nullptr)
    {
        return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
    }
    const ULONG material = material_regions(*c);
    ULONG more = 0, left = 0;
    double vertex_gap = 0.0, edge_gap = 0.0;
    s = status(c->StitchFaces(tolerance, more, left, vertex_gap, edge_gap, nullptr, TRUE));
    if (s != SMLIB_OK)
    {
        return s;
    }
    if (more == 0)
    {
        return SMLIB_OK;
    }
    // Gluing an edge onto an edge two faces share puts its faces in the sector between those
    // faces that their own region is in (SmEdge::FindRadialSector). A sheet that leaves a
    // solid's edges along the solid's own faces, or overlaps them, has no such sector to find;
    // glued elsewhere, it joins the solid's inside to the outside, and the solid's material is
    // lost: that sewing is refused.
    if (material_regions(*c) < material)
    {
        return static_cast<int32_t>(SMLIB_ERR_UNSUPPORTED);
    }
    s = remove_faces_inside(*c);
    if (s != SMLIB_OK)
    {
        return s;
    }
    if (has_spine_edge(*c))
    {
        return SMLIB_OK;
    }
    w = std::move(c);
    lamina = left;
    return SMLIB_OK;
}

} // namespace

extern "C" int32_t smlib_stitch(const smlib_brep *brep, double tolerance, int32_t *out_closed, smlib_brep **out)
{
    if (!valid(brep) || out == nullptr || out_closed == nullptr || !std::isfinite(tolerance) || !(tolerance > 0.0) || tolerance > 1e3)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_stitch", handle(brep) + ", tolerance=" + num(tolerance));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(brep, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            ULONG stitched = 0;
            ULONG lamina = 0;
            double vertex_gap = 0.0;
            double edge_gap = 0.0;
            s = status(w->StitchFaces(tolerance, stitched, lamina, vertex_gap, edge_gap));
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (lamina != 0)
            {
                s = close_onto_shared_edges(w, tolerance, lamina);
                if (s != SMLIB_OK)
                {
                    return s;
                }
            }
            SmBoolean maybe_open = FALSE;
            s = status(w->OrientTrimmedSurfaces(TRUE, maybe_open, lamina != 0 ? TRUE : FALSE));
            if (s != SMLIB_OK)
            {
                return s;
            }
            *out_closed = lamina == 0 && w->IsManifoldSolid() ? 1 : 0;
            return take(w.release(), out);
        }
    );
}

extern "C" int32_t smlib_offset_sheet(const smlib_brep *sheet, double distance, smlib_brep **out)
{
    if (!valid(sheet) || out == nullptr || !std::isfinite(distance) || distance == 0.0 || std::fabs(distance) > 1e6)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_offset_sheet", handle(sheet) + ", distance=" + num(distance));
    c.input("sheet", sheet);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned w = copy_of(sheet, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            SmBrep *moved = nullptr;
            s = status(SmApiOffsetBrepFull(w.get(), distance, TRUE, TRUE, moved));
            if (s != SMLIB_OK || moved == nullptr)
            {
                if (moved != w.get())
                {
                    delete moved;
                }
                return s != SMLIB_OK ? s : SMLIB_ERR_EMPTY;
            }
            // OffsetBrep leaves its input allocated (possibly changed) beside the result.
            if (moved != w.get())
            {
                w.reset(moved);
            }
            return take(w.release(), out);
        }
    );
}
