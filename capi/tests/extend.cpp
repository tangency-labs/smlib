// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_extend_sheet: a square plane sheet continued past one side in each mode, a bilinear
// B-spline sheet and a parabolic one continued naturally (by the distance along the surface), a
// cylinder's strip continued round its axis as the cylinder and by run-offs the distance long
// across it, and the shapes and arguments it refuses.

#include "../src/smlib_shim.h"

#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>
#include <SmApiSurfaces.h>
#include <SmBrep.h>
#include <SmCone.h>
#include <SmFace.h>
#include <SmSurface.h>

#include <cmath>
#include <cstdio>

static int failures = 0;

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);          \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

static smlib_brep *wrap(SmBrep *p)
{
    return p == nullptr ? nullptr : new smlib_brep{p};
}

static double area_of(const smlib_brep *b)
{
    double a = -1.0;
    return smlib_brep_area(b, 1e-7, &a) == SMLIB_OK ? a : -1.0;
}

int main()
{
    CHECK(smlib_init() == SMLIB_OK);
    const SmContext &ctx = *SmApiGetOrCreateContext();

    // The square [0, 10] x [0, 10] at z = 0, continued 5 past its side x = 10: 15 x 10 in each
    // mode (a plane's G1 and G2 run-offs are the plane).
    SmBrep *p = nullptr;
    CHECK(SmApiCreatePlane(SmVector3d(0.0, 0.0, 0.0), 10.0, 10.0, p) == SM_SUCCESS);
    smlib_brep *square = wrap(p);
    const double side[3] = {10.0, 5.0, 0.0};
    for (int32_t mode : {SMLIB_EXTEND_NATURAL, SMLIB_EXTEND_TANGENT, SMLIB_EXTEND_CURVATURE})
    {
        smlib_brep *r = nullptr;
        CHECK(smlib_extend_sheet(square, side, 1e-3, 5.0, mode, &r) == SMLIB_OK);
        CHECK(std::fabs(area_of(r) - 150.0) < 1e-6);
        int64_t n[3] = {0, 0, 0};
        CHECK(smlib_brep_counts(r, n) == SMLIB_OK && n[0] == 1 && n[1] == 4);
        smlib_brep_free(r);
    }
    // The input is unchanged.
    CHECK(std::fabs(area_of(square) - 100.0) < 1e-6);

    // A bilinear B-spline sheet (a flat square 10 x 10) continued 4 past its side y = 10.
    SmSurface *s = nullptr;
    CHECK(SmApiCreateSurfaceFromCornerPoints(SmPoint3d(0, 0, 0), SmPoint3d(10, 0, 0), SmPoint3d(0, 10, 0), SmPoint3d(10, 10, 0), s) ==
          SM_SUCCESS);
    SmBrep *bs = new (ctx) SmBrep();
    SmFace *f = nullptr;
    CHECK(s != nullptr && bs->CreateFaceFromSurface(s, s->GetNaturalUVDomain(), f) == SM_SUCCESS);
    smlib_brep *spline = wrap(bs);
    const double top[3] = {5.0, 10.0, 0.0};
    smlib_brep *grown = nullptr;
    CHECK(smlib_extend_sheet(spline, top, 1e-3, 4.0, SMLIB_EXTEND_NATURAL, &grown) == SMLIB_OK);
    CHECK(std::fabs(area_of(grown) - 140.0) < 1e-4);

    // The parabolic sheet z = x²/20 (0 <= x, y <= 10: a quadratic Bézier extruded 10 along y)
    // continued naturally 5 past its side x = 10: the parabola itself, 5 further along its arc
    // (CreateExtendedSurface alone stops about 2.5% short), area 10 · (s(10) + 5), s the arc
    // length from x = 0.
    {
        const int64_t kinds[3] = {SMLIB_CURVE_NURBS, 2, 3};
        const double reals[18] = {0, 0, 0, 1, 5, 0, 0, 1, 10, 0, 5, 1, 0, 0, 0, 1, 1, 1};
        const double along[3] = {0.0, 1.0, 0.0};
        smlib_brep *parabola = nullptr;
        CHECK(smlib_extrude(kinds, 3, reals, 18, along, 10.0, 0.0, 0, 0, &parabola, nullptr) == SMLIB_OK);
        auto arc = [](double x) { return 0.5 * x * std::sqrt(1.0 + x * x / 100.0) + 5.0 * std::asinh(x / 10.0); };
        CHECK(std::fabs(area_of(parabola) - 10.0 * arc(10.0)) < 1e-5);
        const double edge[3] = {10.0, 5.0, 5.0};
        smlib_brep *longer = nullptr;
        CHECK(smlib_extend_sheet(parabola, edge, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &longer) == SMLIB_OK);
        const double want = 10.0 * (arc(10.0) + 5.0);
        std::printf("parabola continued 5 along its arc: area %.9f, want %.9f\n", area_of(longer), want);
        CHECK(std::fabs(area_of(longer) - want) < 1e-5);
        smlib_brep_free(parabola);
        smlib_brep_free(longer);
    }

    // A quarter of a cylinder's side (radius 10 about z, heights 0..20: an arc extruded)
    // continued naturally 5 round its axis past its side at angle 90 degrees: the cylinder itself,
    // a cone face of SMLib's again (the extended B-spline's parameters are past what a cone holds),
    // 10 * (pi/2 + 1/2) * 20.
    {
        const int64_t kinds[1] = {SMLIB_CURVE_ARC};
        const double arc[12] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 10, 0, 1.5707963267948966};
        const double up[3] = {0.0, 0.0, 1.0};
        smlib_brep *strip = nullptr;
        CHECK(smlib_extrude(kinds, 1, arc, 12, up, 20.0, 0.0, 0, 0, &strip, nullptr) == SMLIB_OK);
        const double straight[3] = {0.0, 10.0, 10.0};
        smlib_brep *round = nullptr;
        CHECK(smlib_extend_sheet(strip, straight, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &round) == SMLIB_OK);
        const double want = 10.0 * (1.5707963267948966 + 0.5) * 20.0;
        std::printf("cylinder strip continued 5 round its axis: area %.9f, want %.9f\n", area_of(round), want);
        CHECK(std::fabs(area_of(round) - want) < 1e-5);
        SmTArray<SmFace *> faces;
        if (round != nullptr)
        {
            round->p->GetFaces(faces);
        }
        CHECK(faces.GetSize() == 1 && faces[0]->GetSurface()->IsKindOf(SmCone_TYPE));
        // Its G1 and G2 run-offs past the same side, 5 long across the surface whatever their
        // curvature: the strip stays its section swept 20 along z, so each adds 5 x 20, 100 pi +
        // 100. (CreateExtendedSurface's G2 run-off of the circle runs about 2.3 times the distance
        // asked and more than doubles with it, which an ask scaled by distance over length alone
        // left 0.28% short after 12 rounds.)
        for (int32_t mode : {SMLIB_EXTEND_TANGENT, SMLIB_EXTEND_CURVATURE})
        {
            smlib_brep *runoff = nullptr;
            CHECK(smlib_extend_sheet(strip, straight, 1e-3, 5.0, mode, &runoff) == SMLIB_OK);
            const double want_runoff = 100.0 * 3.14159265358979323846 + 100.0;
            std::printf("cylinder strip run-off (mode %d) 5 long: area %.9f, want %.9f\n", mode, area_of(runoff), want_runoff);
            CHECK(std::fabs(area_of(runoff) - want_runoff) < 1e-5);
            int64_t n[3] = {0, 0, 0};
            CHECK(smlib_brep_counts(runoff, n) == SMLIB_OK && n[0] == 1);
            smlib_brep_free(runoff);
        }
        smlib_brep_free(strip);
        smlib_brep_free(round);
    }

    // Refused: a box (not one face), a pick on no edge, bad arguments.
    const double o[3] = {0, 0, 0}, sz[3] = {10, 10, 10};
    smlib_brep *box = nullptr;
    CHECK(smlib_make_box(o, sz, &box) == SMLIB_OK);
    smlib_brep *none = nullptr;
    CHECK(smlib_extend_sheet(box, side, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_UNSUPPORTED && none == nullptr);
    const double mid[3] = {5.0, 5.0, 0.0};
    CHECK(smlib_extend_sheet(square, mid, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_PICK_NONE);
    CHECK(smlib_extend_sheet(square, side, 1e-3, -1.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_extend_sheet(square, side, 1e-3, 5.0, 7, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_extend_sheet(nullptr, side, 1e-3, 5.0, SMLIB_EXTEND_NATURAL, &none) == SMLIB_ERR_ARG);

    for (smlib_brep *b : {square, spline, grown, box})
    {
        smlib_brep_free(b);
    }
    if (failures == 0)
    {
        std::printf("extend: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
