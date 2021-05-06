
/** \file
 * \ingroup modifiers
 */

#include "BKE_modifier.h"
#include "DNA_mesh_types.h"
#include "DNA_modifier_types.h"

// not sure, which of those is really used
#include "BLI_utildefines.h"

#include "DNA_defaults.h"
#include "DNA_object_types.h"
#include "DNA_screen_types.h"

#include "BKE_context.h"
#include "BKE_screen.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"

#include <stdio.h>

#include "BKE_mesh.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include <math.h>
#include <stdlib.h>

const float PI = 3.1415926535897932384626433832795;

static void setCoord(float co[], float x, float y, float z)
{
  co[0] = x;
  co[1] = y;
  co[2] = z;
}
static void setCoord2(float coTarget[], float coSource[])
{
  coTarget[0] = coSource[0];
  coTarget[1] = coSource[1];
  coTarget[2] = coSource[2];
}

static float *subtract(float from[], float amount[])
{
  float *value = malloc(3 * sizeof(float));
  value[0] = from[0] - amount[0];
  value[1] = from[1] - amount[1];
  value[2] = from[2] - amount[2];
  return value;
}

static float vLength(float v[])
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
static float vLength2(float v[])
{
  return sqrt(v[0] * v[0] + v[1] * v[1]);
}
static float *vNormalize2(float v[])
{
  float length = vLength2(v);
  float *out = malloc(2 * sizeof(float));
  out[0] = v[0] / length;
  out[1] = v[1] / length;
  return out;
}

static Mesh *modifyMesh(struct ModifierData *md,
                        const struct ModifierEvalContext *ctx,
                        struct Mesh *mesh)
{  // Convert the generic ModifierData to our modifier's DNA data.
   // This is ensured to be valid by the architecture.

  if (mesh->totpoly == 0 || mesh->mpoly[0].totloop < 4)
    return mesh;

  QuarterPipeModifierData *pmd = (QuarterPipeModifierData *)md;
  int steps = pmd->num_olives;

  Mesh *result = BKE_mesh_new_nomain(
      2 + 2 * steps /* vertices */, 0, 0, 4 * steps /* loops */, steps /* face */);
  MVert *mvert = result->mvert;
  MLoop *mloop = result->mloop;

  MVert *origin = mesh->mvert;

  const int faceVs = 4;

  float anglePlus = PI / 2.f / steps;
  float angle = anglePlus;
  MVert *top = (MVert[]){origin[0], origin[2]};

  float size = 1.f;
  float height[2];
  float *pipeStart[2] = {origin[0].co, origin[2].co};
  float *lengthDir = subtract(pipeStart[1], pipeStart[0]);
  float heightDir[2][3];  // points upwards
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      heightDir[i][j] = origin[i * 2].co[j] - origin[i * 2 + 1].co[j];
    }
    height[i] = vLength(heightDir[i]);
    printf("h% ", height[i]);
  }

  float widthFactor = 1.f;
  float width[] = {height[0] * widthFactor, height[1] * widthFactor};  // _

  float lengthDir2D[2] = {lengthDir[0], lengthDir[1]};
  float *lengthDir2D_1N = vNormalize2((float[]){-lengthDir[1], lengthDir[0]});
  float widthDir[2][3];  // points in width direction
  for (int i = 0; i < 2; i++) {
    widthDir[i][0] = lengthDir2D_1N[0] * width[i];
    widthDir[i][1] = lengthDir2D_1N[1] * width[i];
    widthDir[i][2] = 0.f;
  }

  float pipeOrigin[2][3];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      pipeOrigin[i][j] = pipeStart[i][j] + widthDir[i][j];
    }
  }


  setCoord2(mvert[0].co, pipeStart[0]);  // top left
  setCoord2(mvert[1].co, pipeStart[1]);  // top right

  mvert += 2;
  for (int step = 0; step < steps; step++, angle += anglePlus, mvert += 2, mloop += faceVs) {

    float w = -cosf(angle);
    float h = -sinf(angle);

    for (int i = 0; i < 2; i++) {
      int vI = i;
      if (step % 2 == 0)
        vI = 1 - vI;

      for (int j = 0; j < 3; j++) {
        mvert[vI].co[j] = pipeOrigin[i][j] + w * widthDir[i][j] + h * heightDir[i][j];
      }
    }

    for (int i = 0; i < faceVs; i++) {
      mloop[i].v = step * 2 + i;
    }
    result->mpoly[step].loopstart = step * 4;
    result->mpoly[step].totloop = 4;
  }

  // Fill edge data automatically
  BKE_mesh_calc_edges(result, true, false);
  return result;
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA ob_ptr;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, &ob_ptr);

  uiLayoutSetPropSep(layout, true);

  uiItemR(layout, ptr, "num_olives", 0, NULL, ICON_NONE);

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  modifier_panel_register(region_type, eModifierType_QuarterPipe, panel_draw);
}

ModifierTypeInfo modifierType_QuarterPipe = {
    /* name */ "QuarterPipe",
    /* structName */ "QuarterPipeModifierData",
    /* structSize */ sizeof(QuarterPipeModifierData),
    /* srna */ &RNA_QuarterPipeModifier,
    /* type */ eModifierTypeType_Constructive,
    /* flags */ eModifierTypeFlag_AcceptsMesh | eModifierTypeFlag_SupportsEditmode |
        eModifierTypeFlag_EnableInEditmode,  // eModifierTypeFlag_AcceptsCVs |
                                             // eModifierTypeFlag_SupportsMapping |
                                             // eModifierTypeFlag_SupportsEditmode |
                                             // eModifierTypeFlag_EnableInEditmode //
                                             // fromMOD_solidify.c

    /* icon */ ICON_MOD_SOLIDIFY,

    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ NULL,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ NULL,
    /* deformMatricesEM */ NULL,
    /* modifyMesh */ modifyMesh,
    /* modifyHair */ NULL,
    /* modifyGeometrySet */ NULL,

    /* initData */ NULL,
    /* requiredDataMask */ NULL,
    /* freeData */ NULL,
    /* isDisabled */ NULL,
    /* updateDepsgraph */ NULL,
    /* dependsOnTime */ NULL,
    /* dependsOnNormals */ NULL,
    /* foreachIDLink */ NULL,
    /* foreachTexLink */ NULL,
    /* freeRuntimeData */ NULL,

    /* panelRegister */ panelRegister,
    /* blendWrite */ NULL,
    /* blendRead */ NULL,
};
