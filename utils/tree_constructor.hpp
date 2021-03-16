//
// Created by mikolajtwarog on 2021-03-16.
//

#ifndef TECHNIKA_BAKER_TREE_CONSTRUCTOR_HPP
#define TECHNIKA_BAKER_TREE_CONSTRUCTOR_HPP

using namespace boost;

template<typename Graph, typename Problem, typename PlanarEmbedding>
class tree_constructor {
    typedef typename graph_traits<Graph>::edge_descriptor Edge;

    template <typename Edge>
    struct my_visitor : public planar_face_traversal_visitor
    {
        std::map<Edge, std::vector<int> > &faces;
        std::vector<std::vector<int> >& vertices_in_face;
        int current_face = 0;

        my_visitor(std::map<Edge, std::vector<int> >& f, std::vector<std::vector<int> >& o)
                : faces(f), vertices_in_face(o){ }

        void begin_face() {
            vertices_in_face.emplace_back();
        }

        void end_face() {
            current_face++;
        }

        void next_edge(Edge e)
        {
            faces[e].push_back(current_face);
        }

        template <typename Vertex>
        void next_vertex(Vertex v) {
            vertices_in_face[current_face].push_back(v);
        }
    };


    template <typename Edge>
    struct tree_builder : public planar_face_traversal_visitor
    {
        std::map<Edge, std::vector<int> > &faces;
        ::tree<Problem> &tree;
        int current_face = 0;
        int last_vertex;
        Graph graph;
        PlanarEmbedding embedding;

        tree_builder(std::map<Edge, std::vector<int> >& f, ::tree<Problem>& t, Graph& g, PlanarEmbedding& emb)
                : faces(f), tree(t), graph(g), embedding(emb){ }

        void end_face() {
            current_face++;
        }

        template <typename Vertex>
        void next_vertex(Vertex v) {
            tree[current_face].face.push_back(v);
            last_vertex = v;
        }

        void next_edge(Edge e)
        {
            if(current_face != tree.outer_face) {
                int neighbor = faces[e][0] == current_face ? faces[e][1] : faces[e][0];

                if (neighbor == tree.outer_face) {
                    tree.emplace_back();
                    int last = tree.size() - 1;
                    tree[current_face].children.push_back(last);
                    tree[last].parent = current_face;
                    tree[last].label.first = last_vertex;
                    tree[last].label.second = last_vertex == e.m_source ? e.m_target : e.m_source;
                } else {
                    tree[current_face].children.push_back(neighbor);
                }
            }
        }
    };



    void find_outer_face(std::vector<int>& outer_face) {
        std::vector<int> ll;
        std::vector<int> component;
        for (int i = 0; i < embedding.size(); i++) {
            if (!embedding[i].empty()) {
                component.push_back(i);
            }
            ll.push_back(1);
        }
        std::map<typename graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        level_face_traversal(g, embedding, my_vis, 1, ll, component);

        for (const auto& face : vertices_in_face) {
            Edge current_e;
            Edge next_e(face[0], face[1], nullptr);

            bool res = true;
            for (int i = 1; i < face.size() - 1; i++) {
                current_e = next_e;
                next_e.m_source = face[i];
                next_e.m_target = face[i + 1];

                int dis = (get_edge_it(current_e, face[i], embedding) - get_edge_it(next_e, face[i], embedding)) % embedding[face[i]].size();
                if (dis != -1 && dis != embedding[face[i]].size() - 1) {
                    res = false;
                    break;
                }

            }

            if (res) {
                for (int v : face) {
                    outer_face.push_back(v);
                }
                return;
            }
        }
    }

    int name_levels() {
        std::vector<int> outer_face;
        find_outer_face(g, embedding, outer_face);

        for (int i = 0; i < embedding.size(); i++) {
            if (!embedding[i].empty()) {
                vertex_level[i] = -1;
            }
        }

        for (int v : outer_face) {
            vertex_level[v] = 1;
        }

        std::queue<Edge> next_level_edges;

        for (int i = 0; i < outer_face.size(); i++) {
            int v = outer_face[i];
            int w = outer_face[(i + 1) % outer_face.size()];
            for (Edge e : embedding[v]) {
                if (vertex_level[e.m_source] == -1 || vertex_level[e.m_target] == -1) {
                    next_level_edges.push(e);
                }

                if (e.m_source == w) {
                    std::swap(e.m_source, e.m_target);
                    int e_it = get_edge_it(e, w, embedding);
                    std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
                }
            }
        }

        Edge next_level_edge;
        Edge current_edge;
        int starting_v;
        int current_v;
        int current_egde_it;

        int level = 1;
        std::vector<int> current_level;

        while (!next_level_edges.empty()) {
            next_level_edge = next_level_edges.front();
            next_level_edges.pop();

            if (vertex_level[next_level_edge.m_source] > -1 && vertex_level[next_level_edge.m_target] > -1) {
                continue;
            }

            level = vertex_level[next_level_edge.m_source] == -1 ?
                    vertex_level[next_level_edge.m_target] : vertex_level[next_level_edge.m_source];
            level++;

            starting_v = vertex_level[next_level_edge.m_source] == -1 ?
                         next_level_edge.m_source : next_level_edge.m_target;

            current_level.push_back(starting_v);

            current_egde_it = get_edge_it(next_level_edge, starting_v, embedding);

            vertex_level[starting_v] = level;
            bool res = false;
            for (int j = (current_egde_it + 1) % embedding[starting_v].size(); j != current_egde_it;
                 j = (j + 1) % embedding[starting_v].size()) {
                Edge e_j = embedding[starting_v][j];
                if (vertex_level[e_j.m_source] == -1
                    || vertex_level[e_j.m_target] == -1) {
                    current_edge = e_j;
                    res = true;
                    break;
                }
            }
            if (!res) {
                continue;
            }

            vertex_level[starting_v] = -1;

            current_v = current_edge.m_source == starting_v ? current_edge.m_target : current_edge.m_source;
            current_egde_it = get_edge_it(current_edge, current_v, embedding);

            current_level.push_back(current_v);
            vertex_level[current_v] = level;

            while (current_v != starting_v) {
//                Edge e = embedding[current_v][current_egde_it];

                for (int j = (current_egde_it + 1) % embedding[current_v].size(); j != current_egde_it;
                     j = (j + 1) % embedding[current_v].size()) {
                    Edge e_j = embedding[current_v][j];
                    if (vertex_level[e_j.m_source] == -1
                        || vertex_level[e_j.m_target] == -1) {
                        current_edge = e_j;
                        break;
                    }
                }

                current_v = current_edge.m_source == current_v ? current_edge.m_target : current_edge.m_source;
                current_egde_it = get_edge_it(current_edge, current_v, embedding);

                current_level.push_back(current_v);
                vertex_level[current_v] = level;
            }

            vertex_level[starting_v] = level;

            current_level.pop_back();

            for (int i = 0; i < current_level.size(); i++) {
                int v = current_level[i];
                int w = current_level[(i+1) % current_level.size()];
                for (Edge& e : embedding[v]) {
                    if (vertex_level[e.m_source] == -1
                        || vertex_level[e.m_target] == -1) {
                        next_level_edges.push(e);
                    }

                    if (e.m_source == w) {
                        std::swap(e.m_source, e.m_target);
                        int e_it = get_edge_it(e, w, embedding);
                        std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
                    }
                }
            }
        }
        return level;
    }

    void triangulate(std::vector<int>& face, std::vector<int>& component, int turn) {
        int level = vertex_level[face[0]];
        Edge connecting_e;
        int starting_v;
        int starting_v_it = -1;

        for (int i = 0; i < face.size(); i++) {
            int v = face[i];
            for (Edge e : embedding[v]) {
                int neighbour = e.m_source == v ? e.m_target : e.m_source;
                if (std::find(component.begin(), component.end(), neighbour) != component.end()) {
                    connecting_e = e;
                    starting_v = v;
                    starting_v_it = i;
                    break;
                }
            }
            if (starting_v_it > -1)
                break;
        }

        std::rotate(face.begin(), face.begin() + starting_v_it, face.end());

        for (int i = 0; i < face.size(); i++) {
            int first = face[i];
            int second = face[(i + 1) % face.size()];

            Edge temp1(first, second, nullptr);

            int edge_it = get_edge_it(temp1, first, embedding);
            int edge_it2 = get_edge_it(temp1, second, embedding);


            edge_it = (edge_it + (1 * turn) + embedding[first].size()) % embedding[first].size();
            connecting_e = embedding[first][edge_it];

            int target = connecting_e.m_source == first ? connecting_e.m_target : connecting_e.m_source;
            Edge temp2(first, target, nullptr);

            if (!boost::edge(second, target, gr).second) {
                add_edge(second, target, gr);
                int target_e_it = get_edge_it(temp2, target, embedding);
                if (turn == -1) {
                    edge_it2++;
                } else {
                    target_e_it++;
                }
                embedding[second].emplace(embedding[second].begin() + edge_it2, second, target, nullptr);
                embedding[target].emplace(embedding[target].begin() + target_e_it, second, target, nullptr);
                added_edges.emplace(second, target);
            }
        }
    }

    int find_third(int one, int two) {
        auto& one_edges = embedding[one];
        auto& two_edges = embedding[two];

        int level = vertex_level[one];

        if (one == two) {
            for (auto e_i : one_edges) {
                int target_i = e_i.m_source == one ? e_i.m_target : e_i.m_source;

                if (vertex_level[target_i] == level + 1) {
                    return target_i;
                }
            }
        }

        for (auto e_i : one_edges) {
            int target_i = e_i.m_source == one ? e_i.m_target : e_i.m_source;

            if (vertex_level[target_i] == level -1) {
                auto& target_i_edges = embedding[target_i];

                for (auto e_j : target_i_edges) {
                    int target_j = e_j.m_source == target_i ? e_j.m_target : e_j.m_source;
                    if (target_j == two) {
                        return target_i;
                    }
                }
            }
        }

        return -1;
    }

    void get_leaves(::tree<Problem>& t, std::vector<int>& leaves, int v) {
        if (t[v].children.empty()) {
            leaves.push_back(v);
            return;
        }

        for (auto child : t[v].children) {
            get_leaves(t, leaves, child);
        }
    }

    void root_tree_2(::tree<Problem> &t, int node) {
        if (t[node].children.empty())
            return;

        int parent_it = 0;
        auto child = t[node].children.begin();
        for (int i = 0; child != t[node].children.end(); child++, i++) {
            if (*child == t[node].parent) {
                t[node].children.erase(child--);
                parent_it = i;
            } else {
                t[*child].parent = node;
                root_tree_2(t, *child);
            }
        }

        int last = t[node].children.size() - 1;

        std::rotate(t[node].children.begin(),
                    t[node].children.begin() + parent_it,
                    t[node].children.end());

        t[node].label.first = t[t[node].children[0]].label.first;
        t[node].label.second = t[t[node].children[last]].label.second;


    }

    void root_tree_with_root(::tree<Problem> &t, int root) {
        int node;

        for (int i = 0; i < t.size(); i++) {
            for (int v : t[i].face) {
                if (v == root) {
                    node = i;
                    break;
                }
            }
        }

        t.root = node;

        root_tree_2(t, node);

        auto &children = t[node].children;

        for (int i = 0; i < children.size(); i++) {
            int child = children[i];

            if (t[child].label.first == root) {
                std::rotate(t[node].children.begin(),
                            t[node].children.begin() + i,
                            t[node].children.end());
                break;
            }
        }

        t[node].label.first = root;
        t[node].label.second = root;

        if (!t[node].component_tree.empty()) {
            std::vector<int> face;
            for (int x : t[node].children) {
                face.push_back(t[x].label.first);
            }
            if (t[t[node].children.back()].label.second != face[0]) {
                face.push_back(t[t[node].children.back()].label.second);
            }
            std::vector<int> component;
            get_component(component, check_for_component(face), g, embedding);
//            std::reverse(component.begin(), component.end());
            triangulate(face, component, 1);
            triangulate(component, face, 1);
            int v = find_third(t[node].label.first, t[node].label.second);
            root_tree_with_root(t[node].component_tree, v);
        }
    }

    int check_for_component(const std::vector<int> &face) {
        int level = vertex_level[face[0]];

        for (int v : face) {
            for (auto e : embedding[v]) {
                int target = e.m_source == v ? e.m_target : e.m_source;
                if (vertex_level[target] == level + 1) {
                    return target;
                }
            }
        }

        return -1;
    }




    std::vector<bool> visited;

    void get_component(std::vector<int>& component, int v) {
        if (vertex_level[v] == 1) {
            find_outer_face(g, embedding, component);

            return;
        }

        std::queue<int> q;
        q.push(v);
        int level = vertex_level[v];
        int starting_v = -1;
        int connecting_e_it = -1;

        while (!q.empty()) {
            int current_v = q.front();
            q.pop();
//            component.push_back(current_v);

            for (int i = 0; i < embedding[current_v].size(); i++) {
                Edge e = embedding[current_v][i];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (!visited[neighbour] && vertex_level[neighbour] == level){
                    q.push(neighbour);
                    visited[neighbour] = true;
                }

                if (vertex_level[neighbour] == level - 1) {
                    starting_v = current_v;
                    connecting_e_it = i;
                    break;
                }
            }

            if (starting_v > -1) {
                break;
            }
        }

        component.push_back(starting_v);

        Edge current_e;
        int current_v;
        bool res = false;

        auto &edges = embedding[starting_v];
        for (int i = (connecting_e_it - 1 + + edges.size()) % edges.size(); i != connecting_e_it;
             i = (i - 1 + + edges.size()) % edges.size()) {
            Edge e = edges[i];
            int neighbour = e.m_source == starting_v ? e.m_target : e.m_source;

            if (vertex_level[neighbour] == level) {
                current_e = e;
                current_v = neighbour;
                res = true;
                break;
            }
        }
        if (!res) {
            return;
        }

        while (current_v != starting_v) {
            component.push_back(current_v);

            int current_e_it = get_edge_it(current_e, current_v, embedding);
            auto& edges2 = embedding[current_v];
            int i = (current_e_it - 1 + edges2.size()) % edges2.size();
            for (;i != current_e_it; i = (i - 1 + edges2.size()) % edges2.size()) {
                Edge e = edges2[i];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    break;
                }
            }
        }
    }

    ::tree<Problem> build_tree(int v) {
        int level = vertex_level[v];

        std::vector<int> component;
        get_component(component, v, g, embedding, vertex_level);

        std::reverse(component.begin(), component.end());

        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        level_face_traversal(g, embedding, my_vis, level, vertex_level, component);
        ::tree<Problem> t(my_vis.current_face, level);

        int outer = 0;
        for (auto& face : vertices_in_face) {
            int start = -1;
            for (int i = 0; i < face.size(); i++) {
                if (face[i] == component[0]) {
                    start = i;
                    break;
                }
            }
            if (start == -1) {
                continue;
            }

            std::rotate(face.begin(), face.begin() + start, face.end());

            if (face == component) {
                break;
            }
            outer++;
        }

        t.outer_face = outer;
        tree_builder<Edge> tree_b(faces, t, g, embedding);
        level_face_traversal(g, embedding, tree_b, level, vertex_level, component);

        for (int i = 1; i < vertices_in_face.size(); i++) {
            auto &face = vertices_in_face[i];

            int v_in_c = check_for_component(face);
            if (v_in_c > -1) {
                t[i].component_tree = build_tree(v_in_c, g, embedding, vertex_level);
                t[i].component_tree.enclosing_tree = &t;
                t[i].component_tree.enclosing_face = i;
            }
        }

        return t;
    }

    Graph g;
    PlanarEmbedding embedding;
    std::vector<int> vertex_level;
    ::tree<Problem> tree;

public:
    tree_constructor(Graph g, PlanarEmbedding emb): g(g), embedding(emb) {
        name_levels();

        int v = 0;
        graph_traits<Graph>::vertex_iterator vi, vi_end;
        for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
            if (vertex_level[*vi] == 1) {
                v = *vi;
                break;
            }
        }

        tree = build_tree(v);
    }

    ::tree<Problem> get_tree() {
        return tree;
    }
};


#endif //TECHNIKA_BAKER_TREE_CONSTRUCTOR_HPP
