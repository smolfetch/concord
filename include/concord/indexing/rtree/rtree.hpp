#pragma once

#include "../../core/types/point.hpp"
#include "../../geometry/bounding.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace concord {
    namespace indexing {

        /**
         * @brief R-tree spatial index implementation
         *
         * Hierarchical spatial data structure for efficient range queries,
         * nearest neighbor searches, and spatial operations.
         *
         * @tparam T The type of data stored in the R-tree
         */
        template <typename T> class RTree {
          public:
            struct Entry {
                AABB bounds;
                T data;

                Entry(const AABB &bounds, const T &data) : bounds(bounds), data(data) {}
            };

            struct Node {
                AABB bounds;
                std::vector<Entry> entries;
                std::vector<std::shared_ptr<Node>> children;
                bool is_leaf = true;

                Node() = default;
            };

          private:
            std::shared_ptr<Node> root_;
            size_t max_entries_ = 4; // M parameter
            size_t min_entries_ = 2; // m parameter

          public:
            /**
             * @brief Construct R-tree with specified parameters
             * @param max_entries Maximum entries per node (M)
             * @param min_entries Minimum entries per node (m)
             */
            RTree(size_t max_entries = 4, size_t min_entries = 2)
                : max_entries_(max_entries), min_entries_(min_entries) {
                root_ = std::make_shared<Node>();
            }

            /**
             * @brief Insert an entry into the R-tree
             * @param bounds Bounding box of the entry
             * @param data Data associated with the entry
             */
            void insert(const AABB &bounds, const T &data) {
                Entry entry(bounds, data);
                insert_entry(root_, entry);
            }

            /**
             * @brief Search for entries that intersect with the given bounds
             * @param query_bounds The bounding box to search within
             * @return Vector of entries that intersect with query_bounds
             */
            std::vector<Entry> search(const AABB &query_bounds) const {
                std::vector<Entry> results;
                search_recursive(root_, query_bounds, results);
                return results;
            }

            /**
             * @brief Find all entries within a given distance from a point
             * @param point The query point
             * @param radius The search radius
             * @return Vector of entries within the radius
             */
            std::vector<Entry> search_radius(const Point &point, double radius) const {
                AABB query_bounds(Point(point.x - radius, point.y - radius, point.z - radius),
                                  Point(point.x + radius, point.y + radius, point.z + radius));

                auto candidates = search(query_bounds);
                std::vector<Entry> results;

                for (const auto &entry : candidates) {
                    // Check actual distance to entry bounds
                    if (entry.bounds.distance_to_point(point) <= radius) {
                        results.push_back(entry);
                    }
                }

                return results;
            }

            /**
             * @brief Get the total number of entries in the R-tree
             */
            size_t size() const { return count_entries(root_); }

            /**
             * @brief Clear all entries from the R-tree
             */
            void clear() { root_ = std::make_shared<Node>(); }

          private:
            void insert_entry(std::shared_ptr<Node> node, const Entry &entry) {
                if (node->is_leaf) {
                    node->entries.push_back(entry);
                    update_bounds(node);

                    if (node->entries.size() > max_entries_) {
                        split_node(node);
                    }
                } else {
                    // Find the child node with minimum area increase
                    auto best_child = choose_subtree(node, entry.bounds);
                    insert_entry(best_child, entry);
                    update_bounds(node);
                }
            }

            std::shared_ptr<Node> choose_subtree(std::shared_ptr<Node> node, const AABB &bounds) {
                double min_enlargement = std::numeric_limits<double>::max();
                std::shared_ptr<Node> best_child = nullptr;

                for (auto &child : node->children) {
                    double enlargement = child->bounds.union_with(bounds).area() - child->bounds.area();
                    if (enlargement < min_enlargement) {
                        min_enlargement = enlargement;
                        best_child = child;
                    }
                }

                return best_child;
            }

            void split_node(std::shared_ptr<Node> node) {
                // Simple split: divide entries roughly in half
                auto new_node = std::make_shared<Node>();
                new_node->is_leaf = node->is_leaf;

                size_t split_point = node->entries.size() / 2;

                // Move half the entries to the new node
                new_node->entries.insert(new_node->entries.end(), node->entries.begin() + split_point,
                                         node->entries.end());
                node->entries.erase(node->entries.begin() + split_point, node->entries.end());

                update_bounds(node);
                update_bounds(new_node);

                // If this was the root, create a new root
                if (node == root_) {
                    auto new_root = std::make_shared<Node>();
                    new_root->is_leaf = false;
                    new_root->children.push_back(node);
                    new_root->children.push_back(new_node);
                    update_bounds(new_root);
                    root_ = new_root;
                }
            }

            void update_bounds(std::shared_ptr<Node> node) {
                if (node->entries.empty() && node->children.empty()) {
                    return;
                }

                AABB bounds;
                bool first = true;

                for (const auto &entry : node->entries) {
                    if (first) {
                        bounds = entry.bounds;
                        first = false;
                    } else {
                        bounds = bounds.union_with(entry.bounds);
                    }
                }

                for (const auto &child : node->children) {
                    if (first) {
                        bounds = child->bounds;
                        first = false;
                    } else {
                        bounds = bounds.union_with(child->bounds);
                    }
                }

                node->bounds = bounds;
            }

            void search_recursive(std::shared_ptr<Node> node, const AABB &query_bounds,
                                  std::vector<Entry> &results) const {
                if (!node->bounds.intersects(query_bounds)) {
                    return;
                }

                if (node->is_leaf) {
                    for (const auto &entry : node->entries) {
                        if (entry.bounds.intersects(query_bounds)) {
                            results.push_back(entry);
                        }
                    }
                } else {
                    for (const auto &child : node->children) {
                        search_recursive(child, query_bounds, results);
                    }
                }
            }

            size_t count_entries(std::shared_ptr<Node> node) const {
                if (node->is_leaf) {
                    return node->entries.size();
                }

                size_t count = 0;
                for (const auto &child : node->children) {
                    count += count_entries(child);
                }
                return count;
            }
        };

    } // namespace indexing
} // namespace concord
