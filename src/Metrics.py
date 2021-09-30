#!/usr/bin/env python3
import heapq
import graph_tool.all as gt


def Leniency(graph, safe_rooms):
    total_rooms = len(graph.get_vertices())
    safe_rooms_count = 0
    for v in graph.vertices():
        if graph.vp.vertex_property[v].abbrev in safe_rooms:
            safe_rooms_count += 1
    score = safe_rooms_count / total_rooms
    return score


def MissionLinearity(graph, start, end):
    path_reverse = _dijkstra(graph, start, end)

    num_nodes_shortest_path = len(path_reverse)
    total_nodes = len(graph.get_vertices())

    score = num_nodes_shortest_path / total_nodes
    return score


def PathRedundancy(graph, non_critical_rooms):
    total_rooms = len(graph.get_vertices())
    non_critical_rooms_count = 0
    for v in graph.vertices():
        if graph.vp.vertex_property[v].abbrev in non_critical_rooms:
            out_n = graph.get_out_degrees([v])[0]
            if out_n == 0:
                non_critical_rooms_count = non_critical_rooms_count + 1
    score = non_critical_rooms_count / total_rooms
    return score


def MissionLinearityComplex(graph, start, end, lock_rooms):
    return


def MapLinearity(graph):
    total_rooms_with_exits = 0
    single_exits = 0
    double_exits = 0
    for v in graph.vertices():
        out_vertices = graph.get_out_neighbors(v)
        if len(out_vertices) > 0:
            total_rooms_with_exits = total_rooms_with_exits + 1
        if len(out_vertices) == 1:
            single_exits = single_exits + 1
        if len(out_vertices) == 2:
            double_exits = double_exits + 1

    score = (single_exits + (0.5 * double_exits)) / total_rooms_with_exits
    return score


def _dijkstra(graph, start, target):
    graph.vertex_properties["distance"] = graph.new_vertex_property("int")
    graph.vertex_properties["visted"] = graph.new_vertex_property("bool")
    graph.vertex_properties["previous"] = graph.new_vertex_property("int")
    for v in graph.vertices():
        graph.vp.distance[v] = 99999
        graph.vp.visted[v] = False
        graph.vp.previous[v] = -1

    graph.vp.distance[start] = 0
    unvistied_queue = [(0, start)]
    heapq.heapify(unvistied_queue)
    while len(unvistied_queue) > 0:
        u = heapq.heappop(unvistied_queue)
        current_distance, current_vertex = u
        graph.vp.visted[current_vertex] = True
        if current_distance > graph.vp.distance[current_vertex]:
            continue
        for nxt in graph.get_out_neighbors(current_vertex):
            new_dist = graph.vp.distance[current_vertex] + 1
            if new_dist < graph.vp.distance[nxt]:
                graph.vp.distance[nxt] = new_dist
                graph.vp.previous[nxt] = current_vertex
                heapq.heappush(unvistied_queue, (new_dist, nxt))

    path_reverse = [target]
    v = graph.vp.previous[target]
    path_reverse.append(v)
    while v != -1:
        v = graph.vp.previous[v]
        path_reverse.append(v)

    return path_reverse
