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
#include "../utils/level_face_traversal.hpp"

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
struct face_getter : public planar_face_traversal_visitor
{
    property_map<Graph, edge_faces_t>::type &faces;
    std::vector<std::vector<int> >& vertices_in_face;
    int current_face = 0;

    face_getter(property_map<Graph, edge_faces_t>::type& f, std::vector<std::vector<int> >& o)
        : faces(f), vertices_in_face(o){ }

    void begin_face() {
        vertices_in_face.emplace_back();
    }

    void end_face() {
        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
        std::cout << e << "\n";
        faces[e].push_back(current_face);
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        vertices_in_face[current_face].push_back(v);
    }
};


template <typename Edge, typename Problem, typename PlanarEmbedding>
struct tree_builder : public planar_face_traversal_visitor
{
    const property_map<Graph, edge_faces_t>::type &faces;
    std::vector<Problem> &tree;
    int current_face = 0;
    int parent_it = -1;
    int edge_it = 0;
    Graph graph;
    PlanarEmbedding embedding;

    tree_builder(property_map<Graph, edge_faces_t>::type& f, std::vector<Problem>& t, Graph& g, PlanarEmbedding& emb)
            : faces(f), tree(t), graph(g), embedding(emb){ }

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
                tree[last].label.first = e.m_source;
                tree[last].label.second = e.m_target;

                if (e != *embedding[target(e, graph)].begin())
                    std::swap(tree[last].label.first, tree[last].label.second);
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};



template <typename Problem>
void root_tree(std::vector<Problem> &tree, int node) {
    if (tree[node].children.empty())
        return;

    int parent_it = 0;
    auto child = tree[node].children.begin();
    for(int i = 0; child != tree[node].children.end(); child++, i++) {
        if(*child == tree[node].parent) {
            tree[node].children.erase(child--);
            parent_it = i;
        } else {
            tree[*child].parent = node;
            root_tree(tree, *child);
        }
    }

    int last = tree[node].children.size() - 1;

    std::rotate(tree[node].children.begin(),
                tree[node].children.begin()+parent_it,
                tree[node].children.end());

    tree[node].label.first = tree[tree[node].children[0]].label.first;
    tree[node].label.second = tree[tree[node].children[last]].label.second;
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

typedef std::vector<std::vector< graph_traits<Graph>::edge_descriptor > > PlanarEmbedding;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef std::vector<Edge> vec_t;

template <typename Problem>
int baker(Graph &g, PlanarEmbedding &embedding, int root) {
    property_map<Graph, edge_faces_t>::type faces = get(edge_faces_t(), g);
    std::vector<std::vector<int> > outer_face;
    face_getter<graph_traits<Graph>::edge_descriptor> my_vis(&faces, outer_face);
    planar_face_traversal(g, &embedding[0], my_vis);
    std::vector<Problem> tree(my_vis.current_face);

    tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, tree, g, embedding);
    planar_face_traversal(g, &embedding[0], tree_b);
    root_tree(tree, root);

    calculate(tree, root);

    return tree[root].result();
}