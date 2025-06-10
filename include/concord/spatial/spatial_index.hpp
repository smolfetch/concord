#pragma once

#include "../core/types/types.hpp"
#include "../geometry/types_bounding.hpp"
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>

namespace concord {

    // R-Tree for spatial indexing
    template<typename T>
    class RTree {
    public:
        struct Entry {
            AABB bbox;
            T data;
            
            Entry(const AABB& b, const T& d) : bbox(b), data(d) {}
        };
        
    private:
        static constexpr size_t MIN_ENTRIES = 4;
        static constexpr size_t MAX_ENTRIES = 8;
        
        struct Node {
            std::vector<Entry> entries;
            std::vector<std::unique_ptr<Node>> children;
            bool is_leaf = true;
            
            inline AABB getBoundingBox() const {
                if (entries.empty()) return AABB{};
                
                AABB result = entries[0].bbox;
                for (size_t i = 1; i < entries.size(); ++i) {
                    // Expand to include this entry's bbox
                    result.min_point.enu.x = std::min(result.min_point.enu.x, entries[i].bbox.min_point.enu.x);
                    result.min_point.enu.y = std::min(result.min_point.enu.y, entries[i].bbox.min_point.enu.y);
                    result.min_point.enu.z = std::min(result.min_point.enu.z, entries[i].bbox.min_point.enu.z);
                    
                    result.max_point.enu.x = std::max(result.max_point.enu.x, entries[i].bbox.max_point.enu.x);
                    result.max_point.enu.y = std::max(result.max_point.enu.y, entries[i].bbox.max_point.enu.y);
                    result.max_point.enu.z = std::max(result.max_point.enu.z, entries[i].bbox.max_point.enu.z);
                }
                return result;
            }
        };
        
        std::unique_ptr<Node> root_;
        
        inline void insert(Node* node, const Entry& entry) {
            if (node->is_leaf) {
                node->entries.push_back(entry);
                if (node->entries.size() > MAX_ENTRIES) {
                    splitNode(node);
                }
            } else {
                // Find best child to insert into
                Node* best_child = chooseBestChild(node, entry.bbox);
                insert(best_child, entry);
            }
        }
        
        inline Node* chooseBestChild(Node* node, const AABB& bbox) {
            Node* best = node->children[0].get();
            double best_enlargement = calculateEnlargement(best->getBoundingBox(), bbox);
            
            for (size_t i = 1; i < node->children.size(); ++i) {
                double enlargement = calculateEnlargement(node->children[i]->getBoundingBox(), bbox);
                if (enlargement < best_enlargement) {
                    best_enlargement = enlargement;
                    best = node->children[i].get();
                }
            }
            
            return best;
        }
        
        inline double calculateEnlargement(const AABB& bbox, const AABB& new_bbox) {
            AABB expanded = bbox;
            expanded.expand(new_bbox.min_point);
            expanded.expand(new_bbox.max_point);
            
            return expanded.volume() - bbox.volume();
        }
        
        inline void splitNode(Node* node) {
            // Simple split algorithm - in practice you'd want a more sophisticated approach
            if (node->entries.size() <= MAX_ENTRIES) return;
            
            // Create new node
            auto new_node = std::make_unique<Node>();
            new_node->is_leaf = node->is_leaf;
            
            // Move half the entries to the new node
            size_t split_point = node->entries.size() / 2;
            new_node->entries.assign(
                node->entries.begin() + split_point,
                node->entries.end()
            );
            node->entries.resize(split_point);
            
            // If this was the root, create a new root
            if (node == root_.get()) {
                auto new_root = std::make_unique<Node>();
                new_root->is_leaf = false;
                new_root->children.push_back(std::move(root_));
                new_root->children.push_back(std::move(new_node));
                root_ = std::move(new_root);
            }
        }
        
        inline void searchRecursive(Node* node, const AABB& query_bbox, std::vector<T>& results) const {
            if (!node) return;
            
            AABB node_bbox = node->getBoundingBox();
            if (!node_bbox.intersects(query_bbox)) {
                return;
            }
            
            if (node->is_leaf) {
                for (const auto& entry : node->entries) {
                    if (entry.bbox.intersects(query_bbox)) {
                        results.push_back(entry.data);
                    }
                }
            } else {
                for (const auto& child : node->children) {
                    searchRecursive(child.get(), query_bbox, results);
                }
            }
        }
        
    public:
        RTree() : root_(std::make_unique<Node>()) {}
        
        inline void insert(const AABB& bbox, const T& data) {
            Entry entry{bbox, data};
            insert(root_.get(), entry);
        }
        
        inline std::vector<T> search(const AABB& query_bbox) const {
            std::vector<T> results;
            searchRecursive(root_.get(), query_bbox, results);
            return results;
        }
        
        inline std::vector<T> searchPoint(const Point& point) const {
            AABB point_bbox{point, point};
            return search(point_bbox);
        }
        
        inline void clear() {
            root_ = std::make_unique<Node>();
        }
    };
    
    // Quadtree for 2D spatial indexing
    template<typename T>
    class QuadTree {
    private:
        struct Node {
            AABB boundary;
            std::vector<std::pair<Point, T>> points;
            std::unique_ptr<Node> children[4]; // NW, NE, SW, SE
            
            static constexpr size_t CAPACITY = 16;
            
            Node(const AABB& bounds) : boundary(bounds) {}
            
            inline bool isLeaf() const {
                return children[0] == nullptr;
            }
            
            inline void subdivide() {
                Point center = boundary.center();
                
                // Create quadrants
                children[0] = std::make_unique<Node>(AABB{ // NW
                    boundary.min_point,
                    Point{ENU{center.enu.x, boundary.max_point.enu.y, boundary.max_point.enu.z}, Datum{}}
                });
                children[1] = std::make_unique<Node>(AABB{ // NE
                    Point{ENU{center.enu.x, center.enu.y, boundary.min_point.enu.z}, Datum{}},
                    boundary.max_point
                });
                children[2] = std::make_unique<Node>(AABB{ // SW
                    Point{ENU{boundary.min_point.enu.x, boundary.min_point.enu.y, boundary.min_point.enu.z}, Datum{}},
                    Point{ENU{center.enu.x, center.enu.y, boundary.max_point.enu.z}, Datum{}}
                });
                children[3] = std::make_unique<Node>(AABB{ // SE
                    Point{ENU{center.enu.x, boundary.min_point.enu.y, boundary.min_point.enu.z}, Datum{}},
                    Point{ENU{boundary.max_point.enu.x, center.enu.y, boundary.max_point.enu.z}, Datum{}}
                });
            }
        };
        
        std::unique_ptr<Node> root_;
        
        inline bool insert(Node* node, const Point& point, const T& data) {
            if (!node->boundary.contains(point)) {
                return false;
            }
            
            if (node->points.size() < Node::CAPACITY && node->isLeaf()) {
                node->points.emplace_back(point, data);
                return true;
            }
            
            if (node->isLeaf()) {
                node->subdivide();
                
                // Redistribute existing points
                auto old_points = std::move(node->points);
                node->points.clear();
                
                for (const auto& [pt, dt] : old_points) {
                    bool inserted = false;
                    for (int i = 0; i < 4; ++i) {
                        if (insert(node->children[i].get(), pt, dt)) {
                            inserted = true;
                            break;
                        }
                    }
                    if (!inserted) {
                        node->points.emplace_back(pt, dt); // Keep in parent if doesn't fit children
                    }
                }
            }
            
            // Try to insert in children
            for (int i = 0; i < 4; ++i) {
                if (insert(node->children[i].get(), point, data)) {
                    return true;
                }
            }
            
            // If all children reject, keep in this node
            node->points.emplace_back(point, data);
            return true;
        }
        
        inline void query(Node* node, const AABB& range, std::vector<T>& results) const {
            if (!node || !node->boundary.intersects(range)) {
                return;
            }
            
            for (const auto& [point, data] : node->points) {
                if (range.contains(point)) {
                    results.push_back(data);
                }
            }
            
            if (!node->isLeaf()) {
                for (int i = 0; i < 4; ++i) {
                    query(node->children[i].get(), range, results);
                }
            }
        }
        
    public:
        QuadTree(const AABB& boundary) : root_(std::make_unique<Node>(boundary)) {}
        
        inline bool insert(const Point& point, const T& data) {
            return insert(root_.get(), point, data);
        }
        
        inline std::vector<T> query(const AABB& range) const {
            std::vector<T> results;
            query(root_.get(), range, results);
            return results;
        }
        
        inline std::vector<T> queryRadius(const Point& center, double radius) const {
            AABB range{
                Point{ENU{center.enu.x - radius, center.enu.y - radius, center.enu.z - radius}, Datum{}},
                Point{ENU{center.enu.x + radius, center.enu.y + radius, center.enu.z + radius}, Datum{}}
            };
            
            auto candidates = query(range);
            std::vector<T> results;
            
            for (const auto& candidate : candidates) {
                // Note: This assumes T has a way to get its position
                // In practice, you'd need to store points separately or have a position accessor
            }
            
            return results;
        }
        
        inline void clear() {
            root_ = std::make_unique<Node>(root_->boundary);
        }
    };
    
    // Spatial hash grid for fast approximate queries
    template<typename T>
    class SpatialHashGrid {
    private:
        struct Cell {
            std::vector<std::pair<Point, T>> items;
        };
        
        double cell_size_;
        std::unordered_map<int64_t, Cell> grid_;
        
        inline int64_t hash(int x, int y) const {
            // Simple hash combining x and y coordinates
            return (static_cast<int64_t>(x) << 32) | static_cast<int64_t>(y);
        }
        
        inline std::pair<int, int> getGridCoords(const Point& point) const {
            int x = static_cast<int>(std::floor(point.enu.x / cell_size_));
            int y = static_cast<int>(std::floor(point.enu.y / cell_size_));
            return {x, y};
        }
        
    public:
        SpatialHashGrid(double cell_size) : cell_size_(cell_size) {}
        
        inline void insert(const Point& point, const T& data) {
            auto [x, y] = getGridCoords(point);
            int64_t key = hash(x, y);
            grid_[key].items.emplace_back(point, data);
        }
        
        inline std::vector<T> query(const Point& center, double radius) const {
            std::vector<T> results;
            
            int min_x = static_cast<int>(std::floor((center.enu.x - radius) / cell_size_));
            int max_x = static_cast<int>(std::floor((center.enu.x + radius) / cell_size_));
            int min_y = static_cast<int>(std::floor((center.enu.y - radius) / cell_size_));
            int max_y = static_cast<int>(std::floor((center.enu.y + radius) / cell_size_));
            
            double radius_sq = radius * radius;
            
            for (int x = min_x; x <= max_x; ++x) {
                for (int y = min_y; y <= max_y; ++y) {
                    int64_t key = hash(x, y);
                    auto it = grid_.find(key);
                    if (it != grid_.end()) {
                        for (const auto& [point, data] : it->second.items) {
                            double dx = point.enu.x - center.enu.x;
                            double dy = point.enu.y - center.enu.y;
                            if (dx*dx + dy*dy <= radius_sq) {
                                results.push_back(data);
                            }
                        }
                    }
                }
            }
            
            return results;
        }
        
        void clear() {
            grid_.clear();
        }
        
        size_t size() const {
            size_t total = 0;
            for (const auto& [key, cell] : grid_) {
                total += cell.items.size();
            }
            return total;
        }
    };

} // namespace concord
