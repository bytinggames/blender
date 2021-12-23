
const float growTopVertPlank = 0.1f; // grows the top plank by 10% of the distance to the ground (from the top)

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
  bool ring;
  nodeEdge_t *lastEdgeInRing;
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
  float rayStart[3];
  mul_v3_v3fl(rayStart, slopeDir, 0.1f);
  add_v3_v3(rayStart, v->co);
  float overGroundOffset = 0.1f;
  rayStart[2] += overGroundOffset;  // start the ray a bit higher, in case the vertex is directly
                                    // positioned on the ground
  if (!isect_ray_bmesh(rayStart, rayDir, bm, &distance, bestNormal)) {
    distance = 1.f;
    bestNormal[0] = 0.f;
    bestNormal[1] = 0.f;
    bestNormal[2] = 1.f;
  }
  else
    distance -= overGroundOffset;

  
  float loweredCo[3]; // for making the top plank higher (starting the quarter a little bit deeper)
  copy_v3_v3(loweredCo, v->co);
  if (growTopVertPlank != 0) {
    float growTopPlank = distance * growTopVertPlank;
    distance -= growTopPlank;
    loweredCo[2] -= growTopPlank;
  }

  float tangentU[3];
  float tangent[3];

  cross_v3_v3v3(tangentU, slopeDir, bestNormal);
  cross_v3_v3v3(tangent, tangentU, bestNormal);
  normalize_v3(tangent);

  if (dot_v3v3(tangent, slopeDir) < 0) {
    negate_v3(tangent);
  }

  float solidifyThickness[3];
  mul_v3_v3fl(solidifyThickness, tangent, 0.1f);

  float wallFraction = 0.331f;
  float wallDistance = distance * wallFraction;
  float slopeHeight = distance * (1.f - wallFraction);

  float wallDir[3];
  float extrudeDir[3];
  mul_v3_v3fl(wallDir, rayDir, wallDistance);
  mul_v3_v3fl(extrudeDir, rayDir, slopeHeight);
  mul_v3_fl(tangent, slopeHeight);
  //add_v3_v3(rayDir, tangent);

  float pipeCenter[3];
  copy_v3_v3(pipeCenter, loweredCo);
  add_v3_v3(pipeCenter, tangent);
  add_v3_v3(pipeCenter, solidifyThickness);
  int slopeSteps = steps - 1;
  float anglePlus = M_PI_2 / slopeSteps;
  float angle = anglePlus;

  float vGenerate[3];
  float t[3];
  float r[3];
  BMVert *vLast = v;
  BMEdge *firstGenEdge = NULL;

  // move base vertex by solidify thickness
  add_v3_v3(v->co, solidifyThickness);

  // generate wall vertex
  copy_v3_v3(vGenerate, pipeCenter);
  sub_v3_v3(vGenerate, tangent);
  add_v3_v3(vGenerate, wallDir);
  vLast = extrude(bm, vLast, vGenerate, &firstGenEdge);

  for (int step = 0; step < slopeSteps; step++, angle += anglePlus) {
    // calc v from pipe center +
    float s = sinf(angle);
    float c = cosf(angle);
    // pipeCenter - tangent * c + rayDir * s;
    mul_v3_v3fl(t, tangent, c);
    mul_v3_v3fl(r, extrudeDir, s);
    copy_v3_v3(vGenerate, pipeCenter);
    sub_v3_v3(vGenerate, t);
    add_v3_v3(vGenerate, r);
    add_v3_v3(vGenerate, wallDir);
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
    lastRail->ring = false;

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
        if (diskLink->next == current) // dead end
        {
          // check if last vertex has the same position as the first vertex (ring detection)
          // get last vertex
          float *lastVertex;
          if (currentEdgeList->v1AfterV2)
            lastVertex = current->v1->co;
          else
            lastVertex = current->v2->co;
          // get first vertex
          float *firstVertex;
          if (lastRail->val->v1AfterV2)
            firstVertex = lastRail->val->val->v2->co;
          else
            firstVertex = lastRail->val->val->v1->co;

          if (equals_v3v3(lastVertex, firstVertex)) {
            // then a ring is detected
            lastRail->ring = true;
            lastRail->lastEdgeInRing = currentEdgeList;
          }
          break;
        }
        if (!(diskLink->next->v1_disk_link.next == diskLink->next->v1_disk_link.prev &&
            diskLink->next->v2_disk_link.next == diskLink->next->v2_disk_link.prev)) // detected branch: not allowed
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
    BMEdge *firstExtrudeEdge = NULL;

    // slopeDir =  normalize(rayDir x edgeDir)
    getSlopeDir(slopeDir1, edgeIter->val, rayDir, edgeIter->v1AfterV2);
    copy_v3_v3(slopeDir2, slopeDir1);
    if (edgeIter->next != NULL) {
      // slopeDir =  normalize(rayDir x edgeDir + rayDir x next.edgeDir)
      getSlopeDir(slopeDirNext, edgeIter->next->val, rayDir, edgeIter->next->v1AfterV2);
      add_v3_v3(slopeDir2, slopeDirNext);
      normalize_v3(slopeDir2);
    }
    // if it is a ring, modify the slope dir by regarding the last edge too
    if (railIter->ring) {
      nodeEdge_t *last = railIter->lastEdgeInRing;
      float slopeDirLast[3];
      getSlopeDir(slopeDirLast, last->val, rayDir, last->v1AfterV2);
      add_v3_v3(slopeDir1, slopeDirLast);
      normalize_v3(slopeDir1);
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
    firstExtrudeEdge = edge1;

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
      bool fillRing = false;
      if (edgeIter->next != NULL) {
        // slopeDir =  normalize(rayDir x edgeDir + rayDir x next.edgeDir)
        getSlopeDir(slopeDirNext, edgeIter->next->val, rayDir, edgeIter->next->v1AfterV2);
        add_v3_v3(slopeDir2, slopeDirNext);
        normalize_v3(slopeDir2);
      }
      else if (railIter->ring) {
        fillRing = true;
      }

      edge1 = edge2;
      if (!fillRing)
        edge2 = extrudeMain(bm, v, rayDir, slopeDir2, steps);
      else {
        edge2 = firstExtrudeEdge;  // use existing first edge instead to close the ring
        // connect last vertex with first one (merge last to first)
        // instead of merging, delete the last edge, and draw a new one to the first vertex
        copy_v3_v3(v->co, edge2->v1->co);
        BM_edge_kill(bm, edgeIter->val);
        BM_vert_kill(bm, v);
        BMVert *vPrev;
        if (edgeIter->val->v1 == edgePrev->val->v1 || edgeIter->val->v1 == edgePrev->val->v2)
          vPrev = edgeIter->val->v1;
        else
          vPrev = edgeIter->val->v2;

        edgeIter->val = BM_edge_create(bm, vPrev, edge2->v1, NULL, BM_CREATE_NOP);
        edgeIter->v1AfterV2 = false;
      }

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
