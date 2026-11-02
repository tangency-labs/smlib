// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

/*******************************************************************//**
* FILE NAME --- SmBlendNetScope.h
* PURPOSE: SolidCraft change. The scope in which a blend by its definition
*   (SmProceduralBlend) evaluates as its net, the certified B-spline made for
*   the uses that take an approximation: meshing and display.
**********************************************************************/

#ifndef __SMBLENDNETSCOPE_H__
#define __SMBLENDNETSCOPE_H__

#ifndef __SMOS_TYPES_H__
#include <SmTypes.h>
#endif

/*******************************************************************//**
PURPOSE: While one lives in a thread, every SmProceduralBlend's Evaluate,
    EvaluateSimple and EvaluatePoint in that thread are its net's
    (SmBSplineSurface's), within the net's certified bound of the
    definition. SmTess holds one: a mesh is made on the nets.

NOTES: Scopes nest; the thread leaves the net's evaluation when the last
    one ends.
***********************************************************************/
class SM_EXPORT SmBlendNetScope
{
public:
  SmBlendNetScope();
  ~SmBlendNetScope();
  SmBlendNetScope(const SmBlendNetScope &) = delete;
  SmBlendNetScope & operator=(const SmBlendNetScope &) = delete;

  // whether the thread is inside one
  static SmBoolean IsActive();
};

#endif // __SMBLENDNETSCOPE_H__
