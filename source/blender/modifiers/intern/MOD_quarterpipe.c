
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

#include "BLI_math.h"
#include "bmesh.h"

#include "DNA_customdata_types.h"


typedef struct nodeEdge {
  BMEdge *val;
  struct nodeEdge *next;
  bool v1AfterV2;
} nodeEdge_t;

nodeEdge_t *nodeEdge_t_new()
{
  nodeEdge_t *n = malloc(sizeof(nodeEdge_t));
  if (n == NULL)
    printf("well... rails == NULL, malloc didn't work");  // TODO: ask bernie
  n->v1AfterV2 = false;
  return n;
}

typedef struct nodeRail {
  nodeEdge_t *val;
  struct nodeRail *next;
} nodeRail_t;

/* Does crappy fan triangulation of poly, may not be so accurate for
 * concave faces */
static int isect_ray_poly(const float ray_start[3],
                          const float ray_dir[3],
                          BMFace *f,
                          float *r_lambda,
                          float r_normal[3])
{
  BMVert *v, *v_first = NULL, *v_prev = NULL;
  BMIter iter;
  float best_dist = FLT_MAX;
  bool hit = false;

  BM_ITER_ELEM (v, &iter, f, BM_VERTS_OF_FACE) {
    if (!v_first) {
      v_first = v;
    }
    else if (v_prev != v_first) {
      float dist;
      bool curhit;
      float uv[2];
      curhit = isect_ray_tri_v3(ray_start, ray_dir, v_first->co, v_prev->co, v->co, &dist, &uv);
      if (curhit && dist < best_dist) {
        hit = true;
        best_dist = dist;

        float a[3], b[3];
        sub_v3_v3v3(a, v_prev->co, v_first->co);
        sub_v3_v3v3(b, v->co, v_first->co);

        cross_v3_v3v3(r_normal, a, b);
      }
    }

    v_prev = v;
  }

  *r_lambda = best_dist;

  if (best_dist != FLT_MAX)
    normalize_v3(r_normal);

  return hit;
}

static int isect_ray_bmesh(const float ray_start[3], const float ray_dir[3], BMesh *bm, float *r_lambda, float r_normal[3])
{
  BMIter iter;
  float best_dist = FLT_MAX;
  bool hit = false;

  BMFace *f;

  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    float dist;
    bool curhit;

    float normal[3];
    curhit = isect_ray_poly(ray_start, ray_dir, f, &dist, normal);
    if (curhit && dist < best_dist) {
      hit = true;
      best_dist = dist;
      copy_v3_v3(r_normal, normal);
    }
  }

  *r_lambda = best_dist;
  return hit;
}

static BMVert *extrude(BMesh *bm, BMVert *vert, const float co[3], BMEdge **r_edge)
{
  BMVert *v = BM_vert_create(bm, co, NULL, BM_CREATE_NOP);
  BMEdge* e = BM_edge_create(bm, vert, v, NULL, BM_CREATE_NOP);
  if (r_edge != NULL)
    *r_edge = e;
  return v;
}
static BMVert *extrudeRelative(BMesh *bm, BMVert *vert, const float offset[3], BMEdge *r_edge)
{
  float co[3];
  copy_v3_v3(co, vert->co);
  add_v3_v3(co, offset);
  return extrude(bm, vert, co, r_edge);
}

static bool isThereANeighbourLineOnVertex(BMEdge* edgeWithoutFace, bool v1)
{
  BMEdge *e = edgeWithoutFace;

  BMEdge *eBefore = e;
  BMEdge *eIter;
  if (v1)
    eIter = e->v1_disk_link.next;
  else
    eIter = e->v2_disk_link.next;

  while (eIter != e) {
    if (eIter->l == NULL) {
      return true;
    }

    // check which disk link contains the last edge, then continue on that disk link (see
    // https://wiki.blender.org/wiki/Source/Modeling/BMesh/Design)
    bool _v1 = eIter->v1_disk_link.prev == eBefore;

    eBefore = eIter;

    if (_v1)
      eIter = eIter->v1_disk_link.next;
    else
      eIter = eIter->v2_disk_link.next;
  }
  return false;
}

static bool isEndOfRail(nodeRail_t *rails, BMEdge *edge)
{
  nodeRail_t *it = rails;
  while (it != NULL) {
    nodeEdge_t *itEdge = it->val;
    if (itEdge == NULL)
      continue;
    while (itEdge->next != NULL) {
        itEdge = itEdge->next;
    }

    if (itEdge->val == edge)
      return true;

    it = it->next;
  }
  return false;
}

// returns first vreated edge
static BMEdge *extrudeMain(BMesh *bm, BMVert *v, float rayDir[3], float slopeDir[3], const int steps)
{
  float bestNormal[3];
  float distance;
  if (!isect_ray_bmesh(v->co, rayDir, bm, &distance, bestNormal)) {
    distance = 1.f;
    bestNormal[0] = 0.f;
    bestNormal[1] = 0.f;
    bestNormal[2] = 1.f;
  }

  float tangentU[3];
  float tangent[3];

  cross_v3_v3v3(tangentU, slopeDir, bestNormal);
  cross_v3_v3v3(tangent, tangentU, bestNormal);
  normalize_v3(tangent);

  if (dot_v3v3(tangent, slopeDir) < 0) {
    negate_v3(tangent);
  }


  float extrudeDir[3];
  mul_v3_v3fl(extrudeDir, rayDir, distance);
  mul_v3_fl(tangent, distance);
  //add_v3_v3(rayDir, tangent);

  float pipeCenter[3];
  copy_v3_v3(pipeCenter, v->co);
  add_v3_v3(pipeCenter, tangent);
  float anglePlus = M_PI_2 / steps;
  float angle = anglePlus;

  float vGenerate[3];
  float t[3];
  float r[3];
  BMVert *vLast = v;
  BMEdge *firstGenEdge = NULL;
  for (int step = 0; step < steps; step++, angle += anglePlus) {
    // calc v from pipe center +
    float s = sinf(angle);
    float c = cosf(angle);
    // pipeCenter - tangent * c + rayDir * s;
    mul_v3_v3fl(t, tangent, c);
    mul_v3_v3fl(r, extrudeDir, s);
    copy_v3_v3(vGenerate, pipeCenter);
    sub_v3_v3(vGenerate, t);
    add_v3_v3(vGenerate, r);
    if (firstGenEdge == NULL)
      vLast = extrude(bm, vLast, vGenerate, &firstGenEdge);
    else
      vLast = extrude(bm, vLast, vGenerate, NULL);

  }
  return firstGenEdge;
}

static void getSlopeDir(float slopeDir[3], BMEdge *edge, float rayDir[3], bool v1AfterV2)
{
  float edgeDir[3];
  if (v1AfterV2)
    sub_v3_v3v3(edgeDir, edge->v1->co, edge->v2->co);
  else
    sub_v3_v3v3(edgeDir, edge->v2->co, edge->v1->co);
  cross_v3_v3v3(slopeDir, rayDir, edgeDir);
  normalize_v3(slopeDir);
}

static void fillMain(BMesh *bm, BMEdge *baseEdge, BMEdge *edge1, BMEdge *edge2, int steps)
{
  for (int step = 0; step < steps; step++) {

    BMEdge *newBaseEdge = BM_edge_create(bm, edge1->v2, edge2->v2, NULL, BM_CREATE_NOP);
    BMEdge *edges[4] = {edge1, newBaseEdge, edge2, baseEdge};
    BMVert *verts[4] = {edge1->v1, edge1->v2, edge2->v2, edge2->v1};
    BM_face_create(bm, verts, edges, 4, NULL, BM_CREATE_NOP);

    edge1 = edge1->v2_disk_link.next;
    edge2 = edge2->v2_disk_link.next;
    baseEdge = newBaseEdge;
  }
}

static Mesh *modifyMesh(struct ModifierData *md,
                         const struct ModifierEvalContext *ctx,
                         struct Mesh *mesh)
{
  QuarterPipeModifierData *pmd = (QuarterPipeModifierData *)md;
  int steps = pmd->num_olives;

  Mesh *result;
  BMesh *bm;
  CustomData_MeshMasks cd_mask_extra = {
      .vmask = CD_MASK_ORIGINDEX, .emask = CD_MASK_ORIGINDEX, .pmask = CD_MASK_ORIGINDEX};

  bm = BKE_mesh_to_bmesh_ex(mesh,
                            &((struct BMeshCreateParams){0}),
                            &((struct BMeshFromMeshParams){
                                .calc_face_normal = true,
                                .cd_mask_extra = cd_mask_extra,
                            }));


  nodeRail_t *rails = NULL;
  nodeRail_t *lastRail = NULL;

  // for each edge
  //  that has no faces (referred as line)
  //   that has max 1 neighbouring edge
  //    that is not already the end of a rail (continuous lines)
  //     create a rail recursively
  //      iterate until you reach a vertex that has not two neighbours
  

  BMEdge *e;
  // for each edge
  BMIter iter;
  BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
    if (e->l != NULL)  // ...that has no faces (check if no loops attached -> no faces on the edge)
      continue;

    // ...that has max 1 neighbouring edge
    // check for an end point
    bool v1EndPoint = e->v1_disk_link.next == e;
    bool v2EndPoint = e->v2_disk_link.next == e;
    if (!v1EndPoint && !v2EndPoint)
      continue;

    //  int neighbourLinesCount = 0;
    // bool neighbourOnV1 = false;
    //  if (isThereANeighbourLineOnVertex(e, true)) {
    //  neighbourLinesCount++;
    //    neighbourOnV1 = true;
    //}
    //  if (isThereANeighbourLineOnVertex(e, false))
    //    neighbourLinesCount++;

    //// ...that has max 1 neighbouring line
    // if (neighbourLinesCount > 1)
    //  continue;

    // ...that is not already the end of a rail (continuous lines)
    if (isEndOfRail(rails, e))
      continue;

    // ...create rail recursively
    // ...iterate until you reach a vertex that has not two neighbours

    // add a new rail
    if (rails == NULL) {
      rails = malloc(sizeof(nodeRail_t));
      if (rails == NULL)
        printf("well... rails == NULL, malloc didn't work");  // TODO: ask bernie
      lastRail = rails;
    }
    else {
      lastRail->next = malloc(sizeof(nodeRail_t));
      lastRail = lastRail->next;
    }

    lastRail->next = NULL;

    nodeEdge_t *currentEdgeList = nodeEdge_t_new();
    lastRail->val = currentEdgeList;

    currentEdgeList->next = NULL;
    currentEdgeList->val = e;

    // iterate through next edges
    if (v1EndPoint && v2EndPoint)  // but first check if there are any edges
      continue;

    BMDiskLink *diskLink;
    if (v1EndPoint)
      diskLink = &e->v2_disk_link;
    else {
      diskLink = &e->v1_disk_link;
      currentEdgeList->v1AfterV2 = true;
    }

      BMEdge *prev = e;
      BMEdge *current = NULL;
      while (true) {
        if (diskLink->next == current)
          break;
          current = diskLink->next;
        if (current->l != NULL)  // check if next edge has no faces
          break;

        // add edge to the rail
        currentEdgeList->next = nodeEdge_t_new();
        currentEdgeList = currentEdgeList->next;
        currentEdgeList->val = current;
        currentEdgeList->next = NULL;

        if (current->v1_disk_link.next == prev)
          diskLink = &current->v2_disk_link;
        else {
          diskLink = &current->v1_disk_link;
          currentEdgeList->v1AfterV2 = true;
        }
        prev = current;
      }
  }

  float rayDir[3] = {0.f, 0.f, -1.f};
  float slopeDirNext[3];  // slope dir for prev v
  float slopeDir1[3]; // slope dir for prev v
  float slopeDir2[3]; // slope dir for next v
  float *slopeDirv1;  // slope dir for v1
  float *slopeDirv2;  // slope dir for v2

  float store[3];

  nodeRail_t *railIter = rails;
  while (railIter != NULL) {

    nodeEdge_t *edgeIter = railIter->val;
    nodeEdge_t *edgePrev = NULL;

    // slopeDir =  normalize(rayDir x edgeDir)
    getSlopeDir(slopeDir1, edgeIter->val, rayDir, edgeIter->v1AfterV2);
    copy_v3_v3(slopeDir2, slopeDir1);
    if (edgeIter->next != NULL) {
      // slopeDir =  normalize(rayDir x edgeDir + rayDir x next.edgeDir)
      getSlopeDir(slopeDirNext, edgeIter->next->val, rayDir, edgeIter->next->v1AfterV2);
      add_v3_v3(slopeDir2, slopeDirNext);
      normalize_v3(slopeDir2);
    }

    bool invert = edgeIter->v1AfterV2;
    if (invert) {
      slopeDirv1 = slopeDir2;
      slopeDirv2 = slopeDir1;
    }
    else {
      slopeDirv1 = slopeDir1;
      slopeDirv2 = slopeDir2;
    }

    BMEdge *edge1 = extrudeMain(bm, edgeIter->val->v1, rayDir, slopeDirv1, steps);
    BMEdge *edge2 = extrudeMain(bm, edgeIter->val->v2, rayDir, slopeDirv2, steps);

    if (invert) {
      BMEdge *store = edge1;
      edge1 = edge2;
      edge2 = store;
    }

    fillMain(bm, edgeIter->val, edge1, edge2, steps);

    edgePrev = edgeIter;
    edgeIter = edgeIter->next;
    while (edgeIter != NULL) {

      BMVert *v;
      if (edgeIter->val->v1 == edgePrev->val->v1 || edgeIter->val->v1 == edgePrev->val->v2)
        v = edgeIter->val->v2;
      else
        v = edgeIter->val->v1;

      copy_v3_v3(slopeDir2, slopeDirNext);
      if (edgeIter->next != NULL) {
        // slopeDir =  normalize(rayDir x edgeDir + rayDir x next.edgeDir)
        getSlopeDir(slopeDirNext, edgeIter->next->val, rayDir, edgeIter->next->v1AfterV2);
        add_v3_v3(slopeDir2, slopeDirNext);
        normalize_v3(slopeDir2);
      }

      edge1 = edge2;
      edge2 = extrudeMain(bm, v, rayDir, slopeDir2, steps);

      fillMain(bm, edgeIter->val, edge1, edge2, steps);

      edgePrev = edgeIter;
      edgeIter = edgeIter->next;
    }

    railIter = railIter->next;
  }

  result = BKE_mesh_from_bmesh_for_eval_nomain(bm, &cd_mask_extra, mesh);
  BM_mesh_free(bm);

  if (result->mloopuv != NULL) {
    for (int i = 0; i < result->totloop; i++) {
      if (result->mloopuv[i].uv[0] == 0 && result->mloopuv[i].uv[1] == 0) {
        result->mloopuv[i].uv[0] = 6.5f / 8.f;
        result->mloopuv[i].uv[1] = 3.5f / 8.f;
      }
    }
  }

  return result;
}

static Mesh *modifyMesh2(struct ModifierData *md,
                        const struct ModifierEvalContext *ctx,
                        struct Mesh *mesh)
{  // Convert the generic ModifierData to our modifier's DNA data.
   // This is ensured to be valid by the architecture.

  if (mesh->totpoly == 0 || mesh->mpoly[0].totloop < 4)
    return mesh;

  QuarterPipeModifierData *pmd = (QuarterPipeModifierData *)md;
  int steps = pmd->num_olives;

  Mesh *result = BKE_mesh_new_nomain_from_template(mesh,
      2 + 2 * steps /* vertices */, 0, 0, 4 * steps /* loops */, steps /* face */);
  BKE_mesh_copy_settings(result, mesh);
  MVert *mvert = result->mvert;
  MLoop *mloop = result->mloop;
  MLoopUV *mloopuv = NULL;
  if (mesh->mloopuv != NULL) {
    mloopuv = malloc(result->totloop * sizeof(MLoopUV));
  }

  const int faceVs = 4;

  float anglePlus = M_PI_2 / steps;
  float angle = anglePlus;

  float size = 1.f;
  float height[2];

  // get pipeStart vertices

  // get two heighest vertices
  int heighestV = 0;
  int secondHeighestV = -1;
  int heighestVChild = -1;
  int secondHeighestVChild = -1;
  for (int i = 1; i < 4; i++) {
    float h = mesh->mvert[i].co[2];
    if (h > mesh->mvert[heighestV].co[2]) {
      heighestV = i;
    }
  }
  for (int i = 0; i < 4; i++) {
    if (heighestV == i)
      continue;
    float h = mesh->mvert[i].co[2];
    if (secondHeighestV == -1 || h > mesh->mvert[secondHeighestV].co[2]) {
      secondHeighestV = i;
    }
  }

  // swap two heighest indices, if the vertex order says so
  if (heighestV > secondHeighestV) {
    SWAP(int, heighestV, secondHeighestV);
  }
  // get children
  bool foundFristChild = false;
  for (int i = 0; i < 4; i++) {
    if (i != heighestV && i != secondHeighestV) {

      if (!foundFristChild) {
        if (len_squared_v2v2(mesh->mvert[i].co, mesh->mvert[heighestV].co) <
            len_squared_v2v2(mesh->mvert[i].co, mesh->mvert[secondHeighestV].co))
          heighestVChild = i;
        else
          secondHeighestVChild = i;
        foundFristChild = true;
      }
      else {
        if (heighestVChild == -1)
          heighestVChild = i;
        else
          secondHeighestVChild = i;
      }
    }
  }

  float *originVertices[4] = {
      mesh->mvert[heighestV].co,
      mesh->mvert[secondHeighestV].co,
      mesh->mvert[heighestVChild].co,
      mesh->mvert[secondHeighestVChild].co,
  };

  
  float *pipeStart[2] = {originVertices[0], originVertices[1]};

  float lengthDir[3];
  sub_v3_v3v3(lengthDir, pipeStart[1], pipeStart[0]);
  float heightDir[2][3];  // points upwards
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      heightDir[i][j] = originVertices[i][j] - originVertices[i + 2][j];
    }
    height[i] = len_v3(heightDir[i]);
  }

  float widthFactor = 1.f;
  float width[] = {height[0] * widthFactor, height[1] * widthFactor};  // _

  float lengthDir2D[2] = {lengthDir[0], lengthDir[1]};
  float lengthDir2D_1N[2] = {-lengthDir[1], lengthDir[0]};
  normalize_v2(lengthDir2D_1N);
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

  // TODO: fix rotation on object

  copy_v3_v3(mvert[0].co, pipeStart[0]);
  copy_v3_v3(mvert[1].co, pipeStart[1]);

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

    if (step % 2 == 0) {
      for (int i = 0; i < faceVs; i++) {
        mloop[i].v = step * 2 + i;
      }
    }
    else {
      for (int i = 0; i < faceVs; i++) {
        mloop[i].v = step * 2 + 3 - i;
      }
    }

    if (mloopuv != NULL) {

      for (int i = 0; i < faceVs; i++) {
        mloopuv[i] = mesh->mloopuv[0];  // just copy the first uv
      }
      mloopuv += faceVs;
    }

    result->mpoly[step].loopstart = step * 4;
    result->mpoly[step].totloop = 4;
  }

  if (mloopuv != NULL) {
    result->mloopuv = mloopuv;
  }
  //result->mat = mesh->mat;
  //result->texflag = mesh->texflag;
  //result->flag = mesh->flag;
  //result->totcol = mesh->totcol;

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

static void initData(ModifierData *md)
{
  QuarterPipeModifierData *qmd = (QuarterPipeModifierData *)md;
  qmd->num_olives = 4;
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

    /* initData */ initData,
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
