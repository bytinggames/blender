
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

static Mesh *modifyMesh(struct ModifierData *md,
                        const struct ModifierEvalContext *ctx,
                        struct Mesh *mesh)
{// Convert the generic ModifierData to our modifier's DNA data.
// This is ensured to be valid by the architecture.
  QuarterPipeModifierData *pmd = (QuarterPipeModifierData *)md;
  int steps = pmd->num_olives;

  float size = 1.f;

  Mesh *result = BKE_mesh_new_nomain(
      2 + 2 * steps /* vertices */, 0, 0, 4 * steps /* loops */, steps /* face */);
  MVert *mvert = result->mvert;
  MLoop *mloop = result->mloop;

  MVert *origin = mesh->mvert;

  const int faceVs = 4;

  float anglePlus = steps / PI / 2.f;
  float angle = 0.f;
  MVert *top = (MVert[]){origin[0], origin[2]};

  setCoord2(mvert[0].co, top[0].co);  // top left
  setCoord2(mvert[1].co, top[1].co);  // top right

  mvert += 2;
  for (int step = 0; step < steps; step++, angle += anglePlus, mvert += 2, mloop += faceVs) {
    setCoord(mvert[0].co, top[1].co[0], top[1].co[1], top[1].co[2] - size);  // bottom right
    setCoord(mvert[1].co, top[0].co[0], top[0].co[1], top[0].co[2] - size);  // bottom left

    for (int i = 0; i < faceVs; i++) {
      mloop[i].v = step * 2 + i;
    }
    result->mpoly[step].loopstart = step * 4;
    result->mpoly[step].totloop = 4;

    top = mvert;
  }


  // get first two vertices
  // extrude X times, each time changing the angle a bit until 90° are reached
    //

  ////mvert[0].co = float[3]{1, 1, 1};
  //setCoord(mvert[0].co, -1.f, -pmd->num_olives, 0.f);
  //// Fill coordinates of the 4 vertices
  ///*mvert[0].co[0] = -1.f;
  //mvert[0].co[1] = -pmd->num_olives;
  //mvert[0].co[2] = 0.f;*/

  //mvert[1].co[0] = -1.f;
  //mvert[1].co[1] = pmd->num_olives;
  //mvert[1].co[2] = 0.f;

  //mvert[2].co[0] = 1.f;
  //mvert[2].co[1] = pmd->num_olives;
  //mvert[2].co[2] = 0.f;

  //mvert[3].co[0] = 1.f;
  //mvert[3].co[1] = -pmd->num_olives;
  //mvert[3].co[2] = 0.f;

  //// Fill the loops
  //result->mloop[0].v = 0;
  //result->mloop[1].v = 1;
  //result->mloop[2].v = 2;
  //result->mloop[3].v = 3;

  //// Fill the face info, i.e. its first loop and total number of loops
  //result->mpoly[0].loopstart = 0;
  //result->mpoly[0].totloop = 4;

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
                                                // eModifierTypeFlag_EnableInEditmode // fromMOD_solidify.c

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
