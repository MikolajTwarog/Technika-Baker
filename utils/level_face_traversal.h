//
// Created by mikolajtwarog on 2021-02-04.
//

#ifndef TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_H
#define TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_H

#include <vector>
#include <set>
#include <map>
#include <boost/next_prior.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>


template < typename Graph, typename Visitor >
inline void level_face_traversal(
        const Graph& g, std::vector<std::vector<typename boost::graph_traits<Graph>::edge_descriptor> >& embedding,
        Visitor& visitor, int level, const std::vector<int>& vertex_level)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_t;

    std::map<edge_t, std::map<vertex_t, edge_t> > next_edge;
    std::map<edge_t, std::set<vertex_t> > visited;

    visitor.begin_traversal();

    for (int v = 0; v < embedding.size(); v++) {
        if (vertex_level[v] != level)
            continue;

        auto edges = embedding[v];

        for (int i = 0; i < edges.size(); i++) {
            edge_t e = edges[i];
            if (vertex_level[source(e, g)] != level
                && vertex_level[target(e, g)] != level)
                continue;

            for (int j = (i+1)%edges.size(); j != i; j = (j+1)%edges.size()) {
                edge_t e_j = edges[j];
                if (vertex_level[source(e_j, g)] == level
                    && vertex_level[target(e_j, g)] == level) {
                    next_edge[e][v] = e_j;
                }
            }

        }
    }

    std::vector<edge_t> edges_cache;
    std::vector<vertex_t> vertices_in_edge;
    typename boost::graph_traits< Graph >::edge_iterator fi, fi_end;

    for (boost::tie(fi, fi_end) = edges(g); fi != fi_end; ++fi) {
        edge_t e(*fi);
        if (vertex_level[source(e, g)] == level
            && vertex_level[target(e, g)] == level) {
            edges_cache.push_back(e);
        }
    }

    for (auto e : edges_cache) {
//            std::cout << e << "\n";

        vertices_in_edge.clear();
        vertices_in_edge.push_back(source(e, g));
        vertices_in_edge.push_back(target(e,g));

        for (auto v : vertices_in_edge) {
            std::set<vertex_t>& e_visited = visited[e];
            typename std::set<vertex_t>::iterator e_visited_found
                = e_visited.find(v);

            if (e_visited_found == visited[e].end())
                visitor.begin_face();

            while (visited[e].find(v) == visited[e].end())
            {
                visitor.next_vertex(v);
                visitor.next_edge(e);
                visited[e].insert(v);
                v = source(e, g) == v ? target(e, g) : source(e, g);
                e = next_edge[e][v];
            }

            if (e_visited_found == e_visited.end())
                visitor.end_face();
        }
    }
}

#endif //TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_H
