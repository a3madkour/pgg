#!/usr/bin/env python3

from Metrics import MapLinearity, MissionLinearity, Leniency, PathRedundancy
import graph_tool.all as gt
from Rule import Rule, RHS, LHS
from VertexProperty import VertexProperty
import pickle

def makeGraph(mongoGraph):
    g = gt.Graph()
    IDtoIndex = {}
    vertices = []
    g.vertex_properties["vertex_property"] = g.new_vertex_property("python::object")
    for i, node in enumerate(mongoGraph["nodes"]):
        v = g.add_vertex()
        IDtoIndex[node["id"]] = i
        abbrev = node["abbrev"]
        label = node["label"].strip()
        terminality = False
        mark = node["mark"]

        vertices.append(v)
        g.vp.vertex_property[v] = VertexProperty(label, abbrev, terminality, mark)

    for edge in mongoGraph["edges"]:
        u = vertices[IDtoIndex[edge["from"]]]
        v = vertices[IDtoIndex[edge["to"]]]
        g.add_edge(u, v)
    return g


def makeGrammar(grammar):
    axiom = makeGraph(grammar["axiom"])

    rules = []
    for i, rule in enumerate(grammar["rules"]):
        rule_name = rule["name"]
        rule_lhs = rule["lhs"]
        lhs_graph = makeGraph(rule_lhs)
        lhs = LHS(lhs_graph)
        rule_rhss = rule["rhs"]
        rhss = []
        for j, rule_rhs in enumerate(rule_rhss):
            rhs_graph = makeGraph(rule_rhs["graph"])
            rhs = RHS(rhs_graph, rule_rhs["probability"])
            rhss.append(rhs)
        rules.append(Rule(rule_name, lhs, rhss, 1))

    return (axiom, rules)


def hasMoreThanFour(graph):
    for v in graph.vertices():
        if len(graph.get_all_neighbors(v)) > 4:
            return True
    return False


def genSample(gen):
    config = {"max_applications": 30, "sampling_method": "grammar"}
    sample, sample_id, gen_tree = gen.applyRules(config)
    graph = makeMongoGraph(sample, sample_id, gen_tree)

    return graph


def genSampleRecipe(gen, recipe, cont):
    config = {"max_applications": 30, "sampling_method": "grammar"}
    sample, sample_id, gen_tree = gen.applyRecipe(config, recipe, cont)
    graph = makeMongoGraph(sample, sample_id, gen_tree)

    return graph


def containsKM(graph):
    for v in graph.vertices():
        if graph.vp.vertex_property[v].abbrev == "km":
            return True
    return False


def findExit(graph):
    exit = 0
    for v in graph.vertices():
        if graph.vp.vertex_property[v].abbrev == "g":
            exit = v
            break
    return exit


def makeMongoGraph(sample, sample_id, gen_tree):
    nodes = []
    safe_rooms = {"e", "g", "l", "lf", "lm", "n", "k", "kf", "km"}
    non_critical_rooms = ["n", "l", "t"]

    exitie = findExit(sample)
    leni = Leniency(sample, safe_rooms)
    mission_lin = MissionLinearity(sample, 0, exitie)
    map_lin = MapLinearity(sample)
    path_redun = PathRedundancy(sample, non_critical_rooms)

    for v in sample.vertices():
        node = {}
        node["index"] = int(v)
        node["abbrev"] = sample.vp.vertex_property[v].abbrev
        node["label"] = sample.vp.vertex_property[v].label
        nodes.append(node)

    edges = []
    for e in sample.edges():
        edge = {}
        u, v = e
        edge["from"] = int(u)
        edge["to"] = int(v)
        edges.append(edge)

    thebytes = pickle.dumps(gen_tree)
    graph = {
        "nodes": nodes,
        "edges": edges,
        "gen_chain": thebytes,
        "Sample ID": sample_id,
        "Leniency": leni,
        "Mission Linearity": mission_lin,
        "Map Linearity": map_lin,
        "Path Redundancy": path_redun,
    }

    return graph
