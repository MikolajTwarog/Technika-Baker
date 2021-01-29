#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <vector>
#include <climits>
#include <queue>

#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "problems.hpp"

using namespace boost;

namespace boost {
    enum edge_faces_t { edge_faces };

    BOOST_INSTALL_PROPERTY(edge, faces);
}

typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

typedef adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int, EdgeProperty>
        >
        Graph;

template <typename Edge>
struct my_visitor : public planar_face_traversal_visitor
{
    property_map<Graph, edge_faces_t>::type &faces;
    int current_face = 0;

    my_visitor(property_map<Graph, edge_faces_t>::type &f)
    : faces(f) { }

    void end_face() {
//        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
//        std::cout << e << " ";
        faces[e].push_back(current_face);
    }
};


template <typename Edge, typename Problem>
struct tree_builder : public planar_face_traversal_visitor
{
    const property_map<Graph, edge_faces_t>::type &faces;
    std::vector<Problem> &tree;
    int current_face = 0;
    int parent_it = -1;
    int edge_it = 0;

    tree_builder(property_map<Graph, edge_faces_t>::type &f, std::vector<Problem> &t)
            : faces(f), tree(t){ }

    void end_face() {
        current_face++;
    }

    void next_edge(Edge e)
    {
        if(current_face > 0) {
            int neighbor = faces[e][0];
            if (neighbor == current_face)
                neighbor = faces[e][1];

            if (neighbor == 0) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};

template <typename Problem>
void root_tree(std::vector<Problem> &tree, int root) {
    std::queue<int> q;
    q.push(root);

    while(!q.empty()) {
        int node = q.front();
        q.pop();

        int parent_it = 0;
        auto child = tree[node].children.begin();
        for(int i = 0; child != tree[node].children.end(); child++, i++) {
            if(*child == tree[node].parent) {
                tree[node].children.erase(child--);
                parent_it = i;
            } else {
                q.push(*child);
                tree[*child].parent = node;
            }
        }

        std::rotate(tree[node].children.begin(),
                    tree[node].children.begin()+parent_it,
                    tree[node].children.end());
    }
}

template <typename Problem>
void calculate(std::vector<Problem> &tree, int root){
    if(tree[root].children.empty()) {
        return;
    }

    calculate(tree, tree[root].children[0]);
    tree[root].val = tree[tree[root].children[0]].val;

    for(int i=1; i < tree[root].children.size(); i++) {
        calculate(tree, tree[root].children[i]);
        tree[root].merge(tree[tree[root].children[i]]);
    }

    tree[root].adjust();
}

template <typename Problem>
int baker(Graph &g, std::vector<std::vector< graph_traits<Graph>::edge_descriptor > > &embedding, int root) {
    property_map<Graph, edge_faces_t>::type faces = get(edge_faces_t(), g);
    my_visitor<graph_traits<Graph>::edge_descriptor> my_vis(faces);
    planar_face_traversal(g, &embedding[0], my_vis);
    std::vector<Problem> tree(my_vis.current_face);

    tree_builder<graph_traits<Graph>::edge_descriptor, Problem> tree_b(faces, tree);
    planar_face_traversal(g, &embedding[0], tree_b);
    root_tree(tree, root);

    calculate(tree, root);

//    int k=0;
//    for(auto i : tree) {
//        std::cout << k++ << " p " << i.parent << " c ";
//        for(auto j : i.children) {
//            std::cout << j << " ";
//        }
//        std::cout << "\n";
//    }
//
//    int y = 1;
//    calculate(tree, y);
//    std::cout << tree[y].val[0] << " " << tree[y].val[1] << " " << tree[y].val[2] << " " << tree[y].val[3] << "\n";

    return tree[root].result();
}