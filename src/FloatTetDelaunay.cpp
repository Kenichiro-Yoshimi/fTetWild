// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <floattetwild/FloatTetDelaunay.h>

#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

#include <iterator>
#include <algorithm>
#include <bitset>
#include <cmath>

#include <floattetwild/Predicates.hpp>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/MeshIO.hpp>

namespace floatTetWild {
	namespace {
        void
        get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
            min = vertices.front();
            max = vertices.front();

            for (size_t j = 0; j < vertices.size(); j++) {
                for (int i = 0; i < 3; i++) {
                    min(i) = std::min(min(i), vertices[j](i));
                    max(i) = std::max(max(i), vertices[j](i));
                }
            }

//            const Scalar dis = std::max((max - min).minCoeff() * params.box_scale, params.eps_input * 2);
            const Scalar dis = std::max(params.ideal_edge_length, params.eps_input * 2);
            for (int j = 0; j < 3; j++) {
                min[j] -= dis;
                max[j] += dis;
            }

            logger().debug("min = {} {} {}", min[0], min[1], min[2]);
            logger().debug("max = {} {} {}", max[0], max[1], max[2]);
        }

        bool comp(const std::array<int, 4> &a, const std::array<int, 4> &b) {
            return std::tuple<int, int, int>(a[0], a[1], a[2]) < std::tuple<int, int, int>(b[0], b[1], b[2]);
        }

        void match_surface_fs(Mesh &mesh, const std::vector<Vector3> &input_vertices,
                              const std::vector<Vector3i> &input_faces, std::vector<bool> &is_face_inserted) {
            std::vector<std::array<int, 4>> input_fs(input_faces.size());
            for (int i = 0; i < input_faces.size(); i++) {
                input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
                std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
            }
            std::sort(input_fs.begin(), input_fs.end(), comp);

//            for(auto& f: input_fs){
//                cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<endl;
//            }
//            cout<<"/////"<<endl;

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    std::array<int, 3> f = {{t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4]}};
                    std::sort(f.begin(), f.end());
//                    cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
                    auto bounds = std::equal_range(input_fs.begin(), input_fs.end(),
                                                   std::array<int, 4>({{f[0], f[1], f[2], -1}}),
                                                   comp);
//                    bool is_matched = false;
//                    int total_ori = 0;
                    for (auto it = bounds.first; it != bounds.second; ++it) {
//                        is_matched = true;
                        int f_id = (*it)[3];
                        is_face_inserted[f_id] = true;
//                        int ori = Predicates::orient_3d(mesh.tet_vertices[t[j]].pos,
//                                                        input_vertices[input_faces[f_id][0]],
//                                                        input_vertices[input_faces[f_id][1]],
//                                                        input_vertices[input_faces[f_id][2]]);
//                        if (ori == Predicates::ORI_POSITIVE)
//                            total_ori++;
//                        else if (ori == Predicates::ORI_NEGATIVE)
//                            total_ori--;
                    }
//                    if (is_matched)
//                        t.is_surface_fs[j] = total_ori;
//                    else
//                        t.is_surface_fs[j] = NOT_SURFACE;

//                    if(is_matched)
//                        cout<<"matched: "<<total_ori<<endl;
                }
            }
        }

        void match_bbox_fs(Mesh &mesh, const Vector3 &min, const Vector3 &max) {
            auto get_bbox_fs = [&](const MeshTet &t, int j) {
                std::array<int, 6> cnts = {{0, 0, 0, 0, 0, 0}};
                for (int k = 0; k < 3; k++) {
                    Vector3 &pos = mesh.tet_vertices[t[(j + k + 1) % 4]].pos;
                    for (int n = 0; n < 3; n++) {
                        if (pos[n] == min[n])
                            cnts[n * 2]++;
                        else if (pos[n] == max[n])
                            cnts[n * 2 + 1]++;
                    }
                }
                for (int i = 0; i < cnts.size(); i++) {
                    if (cnts[i] == 3)
                        return i;
                }
                return NOT_BBOX;
            };

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    t.is_bbox_fs[j] = get_bbox_fs(t, j);
                }
            }
        }

        void
        compute_voxel_points(const Vector3 &min, const Vector3 &max, const Parameters &params, const AABBWrapper &tree,
                             std::vector<Vector3> &voxels) {
            const Vector3 diag = max - min;
            Vector3i n_voxels = (diag / (params.bbox_diag_length * params.box_scale)).cast<int>();

            for (int d = 0; d < 3; ++d)
                n_voxels(d) = std::max(n_voxels(d), 1);

            const Vector3 delta = diag.array() / n_voxels.array().cast<Scalar>();

            voxels.clear();
            voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

//            const double sq_distg = std::min(params.ideal_edge_length / 2, 10 * params.eps);
            const double sq_distg = 100 * params.eps_2;
            GEO::vec3 nearest_point;

            for (int i = 0; i <= n_voxels(0); ++i) {
                const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
                for (int j = 0; j <= n_voxels(1); ++j) {
                    const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
                    for (int k = 0; k <= n_voxels(2); ++k) {
                        const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);
                        if (tree.get_sq_dist_to_sf(Vector3(px, py, pz)) > sq_distg)
                            voxels.emplace_back(px, py, pz);
                    }
                }
            }
        }

        // =====================================================================
        // Octree-based adaptive background grid
        // =====================================================================

        struct OctreeNode {
            Vector3 center;
            Scalar half_size;   // half-length of the cube side
            int depth;
            int children[8];    // -1 = no child (leaf)
            bool is_leaf;

            OctreeNode() : depth(0), is_leaf(true) {
                for (int i = 0; i < 8; i++) children[i] = -1;
            }
        };

        // Child index from octant: bit 0 = x, bit 1 = y, bit 2 = z
        inline Vector3 child_center(const Vector3 &parent_center, Scalar parent_half, int child_idx) {
            Scalar quarter = parent_half * 0.5;
            return Vector3(
                parent_center[0] + ((child_idx & 1) ? quarter : -quarter),
                parent_center[1] + ((child_idx & 2) ? quarter : -quarter),
                parent_center[2] + ((child_idx & 4) ? quarter : -quarter)
            );
        }

        // Refine a leaf node into 8 children
        void refine_node(std::vector<OctreeNode> &nodes, int node_id) {
            Vector3 center = nodes[node_id].center;
            Scalar half = nodes[node_id].half_size;
            int depth = nodes[node_id].depth;

            nodes[node_id].is_leaf = false;
            Scalar child_half = half * 0.5;

            for (int c = 0; c < 8; c++) {
                int child_id = (int)nodes.size();
                nodes[node_id].children[c] = child_id;

                OctreeNode child;
                child.center = child_center(center, half, c);
                child.half_size = child_half;
                child.depth = depth + 1;
                child.is_leaf = true;
                nodes.push_back(child);
            }
        }

        // Check if cell should be refined: the cell is near the surface
        // and larger than the target grid spacing.
        //
        // The octree's purpose is to redistribute background grid points:
        //   - Near surfaces: same density as the uniform grid (grid_spacing)
        //   - Far from surfaces: coarser (saves initial element count)
        // The mesh improvement step handles final sizing via target_edge_length.
        bool should_refine_cell(const OctreeNode &node, const AABBWrapper &tree,
                                Scalar grid_spacing, int max_depth) {
            if (node.depth >= max_depth)
                return false;

            Scalar cell_size = node.half_size * 2.0;

            // Don't refine if cell is already at or below grid spacing
            if (cell_size <= grid_spacing)
                return false;

            // Refinement distance: refine cells whose interior or vicinity
            // intersects the surface. Use cell diagonal + grid_spacing as buffer.
            Scalar refine_dist = cell_size * std::sqrt(3.0) + grid_spacing;
            Scalar refine_dist_sq = refine_dist * refine_dist;

            // Check center
            if (tree.get_sq_dist_to_sf(node.center) < refine_dist_sq)
                return true;

            // Check 8 corners
            for (int i = 0; i < 8; i++) {
                Vector3 corner(
                    node.center[0] + ((i & 1) ? node.half_size : -node.half_size),
                    node.center[1] + ((i & 2) ? node.half_size : -node.half_size),
                    node.center[2] + ((i & 4) ? node.half_size : -node.half_size)
                );
                if (tree.get_sq_dist_to_sf(corner) < refine_dist_sq)
                    return true;
            }

            return false;
        }

        // Check if a cell overlaps with the bounding box
        bool cell_overlaps_bbox(const OctreeNode &node,
                                const Vector3 &bbox_min, const Vector3 &bbox_max) {
            Vector3 cell_min = node.center - Vector3(node.half_size, node.half_size, node.half_size);
            Vector3 cell_max = node.center + Vector3(node.half_size, node.half_size, node.half_size);
            for (int d = 0; d < 3; d++) {
                if (cell_min[d] > bbox_max[d] || cell_max[d] < bbox_min[d])
                    return false;
            }
            return true;
        }

        void build_octree(std::vector<OctreeNode> &nodes, const AABBWrapper &tree,
                          Scalar grid_spacing, int max_depth,
                          const Vector3 &bbox_min, const Vector3 &bbox_max) {
            // Process nodes level by level (BFS-like, but we just iterate since
            // new nodes are appended at the end)
            size_t i = 0;
            while (i < nodes.size()) {
                if (nodes[i].is_leaf &&
                    cell_overlaps_bbox(nodes[i], bbox_min, bbox_max) &&
                    should_refine_cell(nodes[i], tree, grid_spacing, max_depth)) {
                    refine_node(nodes, (int)i);
                }
                i++;
            }
        }

        // Find the leaf node containing a given point by traversing the tree from root.
        // Returns -1 if the point is outside the root cell.
        int find_leaf(const std::vector<OctreeNode> &nodes, const Vector3 &pt) {
            // Check if point is inside root
            const auto &root = nodes[0];
            for (int d = 0; d < 3; d++) {
                if (pt[d] < root.center[d] - root.half_size - 1e-10 ||
                    pt[d] > root.center[d] + root.half_size + 1e-10)
                    return -1;
            }

            int node_id = 0;
            while (!nodes[node_id].is_leaf) {
                int child_idx = 0;
                if (pt[0] > nodes[node_id].center[0]) child_idx |= 1;
                if (pt[1] > nodes[node_id].center[1]) child_idx |= 2;
                if (pt[2] > nodes[node_id].center[2]) child_idx |= 4;

                int next = nodes[node_id].children[child_idx];
                if (next < 0) return node_id; // shouldn't happen in a proper octree
                node_id = next;
            }
            return node_id;
        }

        // 2:1 balance: ensure no leaf cell is more than 2x bigger than its neighbor.
        // Uses tree-traversal for neighbor finding: O(n * depth) per iteration.
        void balance_octree(std::vector<OctreeNode> &nodes, int max_depth) {
            bool changed = true;
            while (changed) {
                changed = false;

                // Collect current leaves (indices may grow during iteration,
                // so snapshot the current set)
                std::vector<int> leaves;
                for (int i = 0; i < (int)nodes.size(); i++) {
                    if (nodes[i].is_leaf)
                        leaves.push_back(i);
                }

                for (int idx : leaves) {
                    if (!nodes[idx].is_leaf)
                        continue;
                    int my_depth = nodes[idx].depth;
                    if (my_depth <= 0)
                        continue;

                    Scalar step = nodes[idx].half_size * 2.0;
                    const Vector3 &my_center = nodes[idx].center;

                    // Check 6 face-neighbors
                    const Vector3 offsets[6] = {
                        { step, 0, 0}, {-step, 0, 0},
                        {0,  step, 0}, {0, -step, 0},
                        {0, 0,  step}, {0, 0, -step}
                    };

                    for (int n = 0; n < 6; n++) {
                        Vector3 neighbor_pt = my_center + offsets[n];
                        int neighbor_id = find_leaf(nodes, neighbor_pt);
                        if (neighbor_id < 0)
                            continue; // outside root

                        // If neighbor is more than 1 level coarser, refine it
                        if (nodes[neighbor_id].is_leaf &&
                            nodes[neighbor_id].depth < my_depth - 1 &&
                            nodes[neighbor_id].depth < max_depth) {
                            refine_node(nodes, neighbor_id);
                            changed = true;
                        }
                    }
                }
            }
        }

        // Octree-guided adaptive grid point generation.
        //
        // Parameters (from inputs.json or CLI):
        //   max_cell_size: far-field grid spacing (coarsest level)
        //
        // Strategy:
        //   1. Generate a regular grid at bbox_diag * box_scale spacing
        //      (same as compute_voxel_points)
        //   2. Build an octree that identifies local cell sizes
        //   3. For each grid point, query the octree to get local cell size
        //   4. Skip points where stride > 1 and indices not aligned
        //   5. Always keep bbox boundary points
        //
        // All kept points lie on the regular lattice → no Delaunay quality loss
        void compute_octree_points(const Vector3 &min, const Vector3 &max,
                                   const Parameters &params, const AABBWrapper &tree,
                                   std::vector<Vector3> &voxels) {
            const Vector3 diag = max - min;

            // Finest grid spacing = same as uniform grid (bbox_diag * box_scale).
            // This matches the original compute_voxel_points behavior exactly.
            Scalar finest_spacing = params.bbox_diag_length * params.box_scale;

            // Determine coarsest grid spacing (far-field resolution)
            Scalar coarsest_spacing = params.octree_max_cell_size;
            if (coarsest_spacing <= 0.0) {
                coarsest_spacing = params.bbox_diag_length * params.box_scale;
            }
            coarsest_spacing = std::max(coarsest_spacing, finest_spacing);

            Scalar grid_spacing = finest_spacing;

            logger().info("Octree: finest_spacing={}, coarsest_spacing={} (ratio={})",
                          finest_spacing, coarsest_spacing, coarsest_spacing / finest_spacing);

            // Compute fine grid dimensions (for index alignment)
            Vector3i n_voxels = (diag / grid_spacing).cast<int>();
            for (int d = 0; d < 3; ++d)
                n_voxels(d) = std::max(n_voxels(d), 1);
            const Vector3 delta = diag.array() / n_voxels.array().cast<Scalar>();

            // Build octree
            Scalar max_side = std::max({diag[0], diag[1], diag[2]});
            Vector3 center = (min + max) * 0.5;
            Scalar root_half = max_side * 0.5;

            int max_depth = params.octree_max_depth;
            if (max_depth <= 0) {
                max_depth = (int)std::ceil(std::log2(max_side / grid_spacing));
                max_depth = std::max(max_depth, 2);
                max_depth = std::min(max_depth, 10);
            }

            Scalar finest_cell = root_half * 2.0 / std::pow(2.0, max_depth);
            logger().info("Octree: root_half={}, max_depth={}, grid_spacing={}, finest_cell={}",
                          root_half, max_depth, grid_spacing, finest_cell);

            std::vector<OctreeNode> nodes;
            nodes.reserve(1024);

            OctreeNode root_node;
            root_node.center = center;
            root_node.half_size = root_half;
            root_node.depth = 0;
            root_node.is_leaf = true;
            nodes.push_back(root_node);

            build_octree(nodes, tree, grid_spacing, max_depth, min, max);
            logger().info("Octree after build: {} nodes", nodes.size());

            balance_octree(nodes, max_depth);

            int n_leaves = 0;
            for (const auto &n : nodes)
                if (n.is_leaf) n_leaves++;
            logger().info("Octree after balance: {} nodes, {} leaves", nodes.size(), n_leaves);

            // Generate grid points with octree-guided thinning (filter approach)
            // Iterate over the full fine grid, skip far-field points based on
            // octree cell size. Always keep bbox boundary points.
            voxels.clear();
            voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

            const double sq_distg = 100 * params.eps_2;
            int n_skipped = 0;

            for (int i = 0; i <= n_voxels(0); ++i) {
                const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
                for (int j = 0; j <= n_voxels(1); ++j) {
                    const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
                    for (int k = 0; k <= n_voxels(2); ++k) {
                        const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);

                        // Check octree stride FIRST (cheap) before surface distance (expensive)
                        Vector3 pt(px, py, pz);
                        int leaf_id = find_leaf(nodes, pt);
                        if (leaf_id >= 0 && nodes[leaf_id].is_leaf) {
                            Scalar local_cell_size = nodes[leaf_id].half_size * 2.0;
                            local_cell_size = std::min(local_cell_size, coarsest_spacing);
                            int stride = std::max(1, (int)std::round(local_cell_size / grid_spacing));
                            int stride_pow2 = 1;
                            while (stride_pow2 * 2 <= stride)
                                stride_pow2 *= 2;

                            if (stride_pow2 > 1) {
                                bool on_bbox_face = (i == 0 || i == n_voxels(0) ||
                                                     j == 0 || j == n_voxels(1) ||
                                                     k == 0 || k == n_voxels(2));
                                if (!on_bbox_face &&
                                    (i % stride_pow2 != 0 || j % stride_pow2 != 0 || k % stride_pow2 != 0)) {
                                    n_skipped++;
                                    continue;
                                }
                            }
                        }

                        // Surface distance check (expensive) only for kept points
                        if (tree.get_sq_dist_to_sf(pt) <= sq_distg)
                            continue;

                        voxels.emplace_back(px, py, pz);
                    }
                }
            }

            int uniform_total = (n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1);
            logger().info("Octree grid: {} background points, {} skipped (uniform grid: {} total)",
                          voxels.size(), n_skipped, uniform_total);
        }
    }

//#include <igl/unique_rows.h>
//#include <floattetwild/Predicates.hpp>
//    extern "C" floatTetWild::Scalar orient3d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc, const floatTetWild::Scalar *pd);

	void FloatTetDelaunay::tetrahedralize(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper &tree,
	        Mesh &mesh, std::vector<bool> &is_face_inserted) {
        const Parameters &params = mesh.params;
        auto &tet_vertices = mesh.tet_vertices;
        auto &tets = mesh.tets;

        is_face_inserted.resize(input_faces.size(), false);

        Vector3 min, max;
        get_bb_corners(params, input_vertices, min, max);
        mesh.params.bbox_min = min;
        mesh.params.bbox_max = max;

        std::vector<Vector3> boxpoints; //(8);
        // for (int i = 0; i < 8; i++) {
        //     auto &p = boxpoints[i];
        //     std::bitset<sizeof(int) * 8> flag(i);
        //     for (int j = 0; j < 3; j++) {
        //         if (flag.test(j))
        //             p[j] = max[j];
        //         else
        //             p[j] = min[j];
        //     }
        // }


        std::vector<Vector3> voxel_points;
        if (params.use_octree) {
            compute_octree_points(min, max, params, tree, voxel_points);
        } else {
            compute_voxel_points(min, max, params, tree, voxel_points);
        }

        const int n_pts = input_vertices.size() + boxpoints.size() + voxel_points.size();
        tet_vertices.resize(n_pts);
//        std::vector<double> V_d;
//        V_d.resize(n_pts * 3);

        size_t index = 0;
        int offset = 0;
        for (int i = 0; i < input_vertices.size(); i++) {
            tet_vertices[offset + i].pos = input_vertices[i];
            // tet_vertices[offset + i].is_on_surface = true;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = input_vertices[i](j);
        }
        offset += input_vertices.size();
        for (int i = 0; i < boxpoints.size(); i++) {
            tet_vertices[i + offset].pos = boxpoints[i];
            // tet_vertices[i + offset].is_on_bbox = true;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = boxpoints[i](j);
        }
        offset += boxpoints.size();
        for (int i = 0; i < voxel_points.size(); i++) {
            tet_vertices[i + offset].pos = voxel_points[i];
            // tet_vertices[i + offset].is_on_bbox = false;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = voxel_points[i](j);
        }

        std::vector<double> V_d;
        V_d.resize(n_pts * 3);
        for (int i = 0; i < tet_vertices.size(); i++) {
            for (int j = 0; j < 3; j++)
                V_d[i * 3 + j] = tet_vertices[i].pos[j];
        }

        GEO::Delaunay::initialize();
        GEO::Delaunay_var T = GEO::Delaunay::create(3, "BDEL");
        T->set_vertices(n_pts, V_d.data());
        //
        tets.resize(T->nb_cells());
        const auto &tet2v = T->cell_to_v();
        for (int i = 0; i < T->nb_cells(); i++) {
            for (int j = 0; j < 4; ++j) {
                const int v_id = tet2v[i * 4 + j];

                tets[i][j] = v_id;
                tet_vertices[v_id].conn_tets.push_back(i);
            }
            std::swap(tets[i][1], tets[i][3]);
        }

        for (int i = 0; i < mesh.tets.size(); i++) {
            auto &t = mesh.tets[i];
            if (is_inverted(mesh.tet_vertices[t[0]].pos, mesh.tet_vertices[t[1]].pos,
                            mesh.tet_vertices[t[2]].pos, mesh.tet_vertices[t[3]].pos)) {
                cout << "EXIT_INV" << endl;
                exit(0);
            }
        }

        //fortest
//        Eigen::MatrixXd VV(mesh.tet_vertices.size(), 3), VVo;
//        Eigen::VectorXi _1, _2;
//        for(int i=0;i<mesh.tet_vertices.size();i++){
//            VV.row(i) = mesh.tet_vertices[i].pos;
//        }
//        igl::unique_rows(VV, VVo, _1, _2);
//        cout<<VV.rows()<<" "<<VVo.rows()<<endl;
//        cout<<T->nb_vertices()<<endl;
//        cout<<mesh.tet_vertices.size()<<endl;
//
//        cout<<"T->nb_finite_cells() = "<<T->nb_finite_cells()<<endl;
//        cout<<"T->nb_cells() = "<<T->nb_cells()<<endl;
//        for (int i=0;i< mesh.tets.size();i++) {
//            auto &t = mesh.tets[i];
//            if (-GEO::PCK::orient_3d(mesh.tet_vertices[t[0]].pos.data(), mesh.tet_vertices[t[1]].pos.data(),
//                                     mesh.tet_vertices[t[2]].pos.data(), mesh.tet_vertices[t[3]].pos.data()) <= 0) {
//                cout << "inverted found!!!! 1" << endl;
//                cout<<i<<endl;
//            }
//        }
//        for (int i=0;i< mesh.tets.size();i++) {
//            auto &t = mesh.tets[i];
//            if (orient3d(mesh.tet_vertices[t[0]].pos.data(), mesh.tet_vertices[t[1]].pos.data(),
//                         mesh.tet_vertices[t[2]].pos.data(), mesh.tet_vertices[t[3]].pos.data()) <= 0) {
//                cout << "inverted found!!!! 2" << endl;
//                cout<<i<<endl;
//            }
//        }
//        for (int i=0;i< mesh.tets.size();i++) {
//            auto &t = mesh.tets[i];
//            if (is_inverted(mesh.tet_vertices[t[0]].pos, mesh.tet_vertices[t[1]].pos,
//                         mesh.tet_vertices[t[2]].pos, mesh.tet_vertices[t[3]].pos)) {
//                cout << "inverted found!!!! 3" << endl;
//                cout<<i<<endl;
//                t.print();
//
//                cout<<std::setprecision(17)<<tet_vertices[t[0]].pos.transpose()<<endl;
//                cout<<tet_vertices[t[1]].pos.transpose()<<endl;
//                cout<<tet_vertices[t[2]].pos.transpose()<<endl;
//                cout<<tet_vertices[t[3]].pos.transpose()<<endl;
//
//                cout<<(tet_vertices[t[0]].pos[0] == tet_vertices[t[1]].pos[0])<<endl;
//                cout<<(tet_vertices[t[1]].pos[0] == tet_vertices[t[2]].pos[0])<<endl;
//                cout<<(tet_vertices[t[2]].pos[0] == tet_vertices[t[3]].pos[0])<<endl;
//
//                cout<<(tet_vertices[t[0]].pos[1] == tet_vertices[t[3]].pos[1])<<endl;
//                cout<<(tet_vertices[t[1]].pos[1] == tet_vertices[t[2]].pos[1])<<endl;
//
//                cout<<(tet_vertices[t[0]].pos[2] == tet_vertices[t[2]].pos[2])<<endl;
//                cout<<(tet_vertices[t[1]].pos[2] == tet_vertices[t[3]].pos[2])<<endl;
//            }
//        }
//        pausee();
//        //fortest

//        //set opp_t_ids
//        for(int t_id = 0;t_id<mesh.tets.size();t_id++) {
//            auto &t = mesh.tets[t_id];
//            for (int j = 0; j < 4; j++) {
//                if (t.opp_t_ids[j] >= 0)
//                    continue;
//                std::vector<int> pair;
//                set_intersection(tet_vertices[t[(j + 1) % 4]].conn_tets,
//                                 tet_vertices[t[(j + 2) % 4]].conn_tets,
//                                 tet_vertices[t[(j + 3) % 4]].conn_tets, pair);
//                if (pair.size() == 2) {
//                    int opp_t_id = pair[0] == t_id ? pair[1] : pair[0];
//                    t.opp_t_ids[j] = opp_t_id;
//                    auto &opp_t = mesh.tets[opp_t_id];
//                    for (int k = 0; k < 4; k++) {
//                        if (opp_t[k] != t[(j + 1) % 4] && opp_t[k] != t[(j + 2) % 4] && opp_t[k] != t[(j + 3) % 4])
//                            opp_t.opp_t_ids[k] = t_id;
//                    }
//                }
//            }
//        }

        //match faces: should be integer with sign
        //match bbox 8 facets: should be -1 and 0~5
//        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted);
        match_bbox_fs(mesh, min, max);

//        MeshIO::write_mesh("delaunay.msh", mesh);
    }
}
