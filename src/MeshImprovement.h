// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_MESHIMPROVEMENT_H
#define FLOATTETWILD_MESHIMPROVEMENT_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Types.hpp>

namespace floatTetWild {
    void init(Mesh &mesh, AABBWrapper& tree);
    void optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops = {{1, 1, 1, 1}});
    void cleanup_empty_slots(Mesh &mesh, double percentage = 0.7);
    void operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 5> &ops = {{1, 1, 1, 1, 1}});
    bool update_scaling_field(Mesh &mesh, Scalar max_energy);
    void compute_surface_based_sizing(Mesh &mesh);

    int get_max_p(const Mesh &mesh);

    void correct_tracked_surface_orientation(Mesh &mesh, AABBWrapper& tree);
    void get_tracked_surface(Mesh& mesh, Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &V, Eigen::Matrix<int, Eigen::Dynamic, 3> &F, int c_id = 0);
    void boolean_operation(Mesh& mesh, int op);
    void boolean_operation(Mesh& mesh, const json &csg_tree_with_ids);
    void filter_outside(Mesh& mesh, bool invert_faces = false);
    void filter_outside_floodfill(Mesh& mesh, bool invert_faces = false);
    void mark_outside(Mesh& mesh, bool invert_faces = false);

    // Per-input "is the tet inside this input mesh?" data, computed once via
    // winding number against each ORIGINAL input geometry (params.surface_sizing_Vs/Fs).
    // Consumed by output_tracked_surface_per_input.
    //   kept_t_ids[idx]  : tet IDs of non-removed tets, in iteration order
    //   Ws[i](idx)       : winding number of tet kept_t_ids[idx] against input i
    // empty() == true when per-input geometry is unavailable (single-input or
    // legacy callers).
    struct PerInputData {
        std::vector<int> kept_t_ids;
        std::vector<Eigen::VectorXd> Ws;
        bool empty() const { return Ws.empty() || kept_t_ids.empty(); }
        int n_inputs() const { return (int)Ws.size(); }
    };

    PerInputData compute_per_input_data(Mesh& mesh);

    // Drop-in replacement for filter_outside in the multi-input (--inputs) path.
    // Uses the per-input winding numbers in `data` to decide inside/outside,
    // voting per surface-bounded region (flood fill blocked at is_surface_fs
    // faces): a region that is >=90% inside some input is kept whole, <=10%
    // is dropped whole, anything in between falls back to the per-tet WN
    // test. Region voting fixes the per-tet 0.5-threshold flicker along
    // internal shared walls (simplify keeps only one of the two coincident
    // wall sheets, so WN near the open sheet is W +- ~0.5): no more interior
    // tets dropped at the wall, no more exterior tets kept past the boundary.
    // Bypasses:
    //   - get_tracked_surface (heavy bfs_orient + duplicate-vertex merge)
    //   - the merged-surface fast_winding_number
    //   - _sf.stl / _tracked_surface.stl writes inside filter_outside
    // Falls back to the original filter_outside when `data.empty()`.
    void filter_outside_per_input(Mesh& mesh, const PerInputData& data,
                                  bool invert_faces = false);

    // Overwrite _tracked_surface.stl with a per-input boundary surface that is
    // robust to is_surface_fs holes (FAIL subdivide_tets, untangle clears).
    // A tet face is emitted (from the kept side, oriented into that side's
    // volume) when either:
    //  (a) WN rule: some input mesh i contains this tet but not the neighbor
    //      (winding-number based; outer boundaries, robust to surface holes),
    //  (b) region rule: both tets are kept but belong to different regions of
    //      the surface-blocked flood fill (connected components of kept tets
    //      whose adjacency does not cross an is_surface_fs face). This
    //      restores internal walls between solids that live in the SAME input
    //      file (identical winding numbers on both sides, invisible to (a))
    //      and the intersection boundaries of overlapping solids.
    // Internal walls between two kept regions are emitted twice with opposite
    // normals. If the tracked face set has a hole in a wall, the two regions
    // merge and (b) degrades to (a)'s behavior — never worse than WN-only.
    void output_tracked_surface_per_input(Mesh& mesh, const std::string& path,
                                          const PerInputData& data);
    void smooth_open_boundary(Mesh& mesh, const AABBWrapper& tree);
    void smooth_open_boundary_aux(Mesh& mesh, const AABBWrapper& tree);
    void manifold_surface(Mesh& mesh);
    void manifold_edges(Mesh& mesh);
    void manifold_vertices(Mesh& mesh);

    void output_info(Mesh& mesh, const AABBWrapper& tree);
    void check_envelope(Mesh& mesh, const AABBWrapper& tree);
    void output_surface(Mesh& mesh, const std::string& filename);

    void untangle(Mesh &mesh);
}

#endif //FLOATTETWILD_MESHIMPROVEMENT_H
