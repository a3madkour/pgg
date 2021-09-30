#!/usr/bin/python

import logging
import random
from GenTree import GenTree
import uuid
import sys
import Rule
import graph_tool.all as gt
import pickle
import copy
from os import listdir, mkdir
from os.path import isfile, join, isdir
from VertexProperty import VertexProperty


class Generator(object):
    def __init__(self, axiom, rules, name=""):
        """
        Constructor.
        Inputs:
            * graph - Graph object on the LHS.
            * abbrevs - Vertex labels
            * prob - Probablity
        Outputs: N/A
        """
        self._axiom = axiom
        self._rules = rules
        self._rules_map = {}
        self._rule_prob_map = {}
        self._default_config = {"max_applications": 20, "sampling_method": "grammar"}
        for rule in rules:
            self._rules_map[rule.name] = rule
        if name == "":
            self._name = str(uuid.uuid4())
        else:
            self._name = name

    @property
    def axiom(self):
        return self._axiom

    @axiom.setter
    def axiom(self, value):
        self._axiom = value

    @property
    def rules_map(self):
        return self._rules_map

    @rules_map.setter
    def rules_map(self, value):
        self._rules_map = value

    @property
    def rules(self):
        return self._rules

    @rules.setter
    def rules(self, value):
        self._rules = value
        for rule in rules:
            self.rules_maps[rule.name] = rule

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    def _loadTreesFromDir(self, directory):
        onlyfiles = [
            directory + f for f in listdir(directory) if isfile(join(directory, f))
        ]
        trees = []
        for filename in onlyfiles:
            input_file = open(filename, "rb")
            tree = pickle.load(input_file)
            trees.append(tree)

        return trees

    def learnParameters(self, directory):
        # TODO seperate the learning logic from the loading of the trees from files
        logging.debug("In learnParameters")
        onlyfiles = [
            directory + f for f in listdir(directory) if isfile(join(directory, f))
        ]

        total_rule_set_map = {}
        total_rule_map = {}
        trees = self._loadTreesFromDir(directory)
        # for filename in onlyfiles:
        #     input_file = open(filename, "rb")
        # tree = pickle.load(input_file)
        for tree in trees:
            for rule_set in tree.rule_set_map:
                if rule_set not in total_rule_set_map:
                    total_rule_set_map[rule_set] = {"rules": {}, "total_count": 0}

                for rule_id in rule_set:
                    if rule_id not in total_rule_set_map[rule_set]["rules"]:
                        total_rule_set_map[rule_set]["rules"][rule_id] = 0

                for rule_id in tree.rule_set_map[rule_set]["rules"]:
                    total_rule_set_map[rule_set]["rules"][rule_id] = (
                        total_rule_set_map[rule_set]["rules"][rule_id]
                        + tree.rule_set_map[rule_set]["rules"][rule_id]
                    )

                total_rule_set_map[rule_set]["total_count"] = (
                    total_rule_set_map[rule_set]["total_count"]
                    + tree.rule_set_map[rule_set]["total_count"]
                )

            for rule_name in tree.rule_map:
                if rule_name not in total_rule_map:
                    total_rule_map[rule_name] = {"rhs_index": {}, "total_count": 0}
                for rhs_index in tree.rule_map[rule_name]["rhs_index"]:
                    if rhs_index in total_rule_map[rule_name]["rhs_index"]:
                        total_rule_map[rule_name]["rhs_index"][rhs_index] = (
                            total_rule_map[rule_name]["rhs_index"][rhs_index]
                            + tree.rule_map[rule_name]["rhs_index"][rhs_index]
                        )
                    else:
                        total_rule_map[rule_name]["rhs_index"][
                            rhs_index
                        ] = tree.rule_map[rule_name]["rhs_index"][rhs_index]

                total_rule_map[rule_name]["total_count"] = (
                    total_rule_map[rule_name]["total_count"]
                    + tree.rule_map[rule_name]["total_count"]
                )

        for rule in self.rules:
            if rule.name in total_rule_map:
                logging.debug("Changing the rhs probabilities for: " + rule.name)
                logging.debug("------------------------------------------")
                for (i, rhs) in enumerate(rule.rhss):
                    logging.debug(
                        "index "
                        + str(i)
                        + " count: "
                        + str(total_rule_map[rule.name]["rhs_index"][i])
                    )
                    logging.debug(
                        "total counter:" + str(total_rule_map[rule.name]["total_count"])
                    )
                    rhs.prob = float(
                        total_rule_map[rule.name]["rhs_index"][i]
                        / total_rule_map[rule.name]["total_count"]
                    )
                    logging.debug("new " + str(i) + "th rhs prob: " + str(rhs.prob))
                logging.debug("------------------------------------------")

        for rule_set in total_rule_set_map:
            if rule_set not in self._rule_prob_map:
                self._rule_prob_map[rule_set] = {}

            logging.debug(
                "Checking the rule probabilites in rule set: " + str(rule_set)
            )
            logging.debug("------------------------------------------")
            for rule_name in total_rule_set_map[rule_set]["rules"]:
                logging.debug(
                    "rule "
                    + rule_name
                    + "'s count: "
                    + str(total_rule_set_map[rule_set]["rules"][rule_name])
                )
                logging.debug(
                    "total counter:" + str(total_rule_set_map[rule_set]["total_count"])
                )
                self._rule_prob_map[rule_set][rule_name] = float(
                    total_rule_set_map[rule_set]["rules"][rule_name]
                    / total_rule_set_map[rule_set]["total_count"]
                )
                logging.debug(
                    "new "
                    + rule_name
                    + " prob: "
                    + str(self._rule_prob_map[rule_set][rule_name])
                )
            logging.debug("------------------------------------------")

    # a way to apply partial rules/rule-set/recipes

    def learnParametersBin(self, trees_bin):
        total_rule_set_map = {}
        total_rule_map = {}
        for tree_bin in trees_bin:
            tree = pickle.loads(tree_bin)
            for rule_set in tree.rule_set_map:
                if rule_set not in total_rule_set_map:
                    total_rule_set_map[rule_set] = {"rules": {}, "total_count": 0}

                for rule_id in rule_set:
                    if rule_id not in total_rule_set_map[rule_set]["rules"]:
                        total_rule_set_map[rule_set]["rules"][rule_id] = 0

                for rule_id in tree.rule_set_map[rule_set]["rules"]:
                    total_rule_set_map[rule_set]["rules"][rule_id] = (
                        total_rule_set_map[rule_set]["rules"][rule_id]
                        + tree.rule_set_map[rule_set]["rules"][rule_id]
                    )

                total_rule_set_map[rule_set]["total_count"] = (
                    total_rule_set_map[rule_set]["total_count"]
                    + tree.rule_set_map[rule_set]["total_count"]
                )

            for rule_name in tree.rule_map:
                if rule_name not in total_rule_map:
                    total_rule_map[rule_name] = {"rhs_index": {}, "total_count": 0}
                for rhs_index in tree.rule_map[rule_name]["rhs_index"]:
                    if rhs_index in total_rule_map[rule_name]["rhs_index"]:
                        total_rule_map[rule_name]["rhs_index"][rhs_index] = (
                            total_rule_map[rule_name]["rhs_index"][rhs_index]
                            + tree.rule_map[rule_name]["rhs_index"][rhs_index]
                        )
                    else:
                        total_rule_map[rule_name]["rhs_index"][
                            rhs_index
                        ] = tree.rule_map[rule_name]["rhs_index"][rhs_index]

                total_rule_map[rule_name]["total_count"] = (
                    total_rule_map[rule_name]["total_count"]
                    + tree.rule_map[rule_name]["total_count"]
                )

        for rule in self.rules:
            if rule.name in total_rule_map:
                logging.debug("Changing the rhs probabilities for: " + rule.name)
                logging.debug("------------------------------------------")
                for (i, rhs) in enumerate(rule.rhss):
                    logging.debug(
                        "index "
                        + str(i)
                        + " count: "
                        + str(total_rule_map[rule.name]["rhs_index"][i])
                    )
                    logging.debug(
                        "total counter:" + str(total_rule_map[rule.name]["total_count"])
                    )
                    rhs.prob = float(
                        total_rule_map[rule.name]["rhs_index"][i]
                        / total_rule_map[rule.name]["total_count"]
                    )
                    logging.debug("new " + str(i) + "th rhs prob: " + str(rhs.prob))
                logging.debug("------------------------------------------")

        for rule_set in total_rule_set_map:
            if rule_set not in self._rule_prob_map:
                self._rule_prob_map[rule_set] = {}

            logging.debug(
                "Checking the rule probabilites in rule set: " + str(rule_set)
            )
            logging.debug("------------------------------------------")
            for rule_name in total_rule_set_map[rule_set]["rules"]:
                logging.debug(
                    "rule "
                    + rule_name
                    + "'s count: "
                    + str(total_rule_set_map[rule_set]["rules"][rule_name])
                )
                logging.debug(
                    "total counter:" + str(total_rule_set_map[rule_set]["total_count"])
                )
                self._rule_prob_map[rule_set][rule_name] = float(
                    total_rule_set_map[rule_set]["rules"][rule_name]
                    / total_rule_set_map[rule_set]["total_count"]
                )
                logging.debug(
                    "new "
                    + rule_name
                    + " prob: "
                    + str(self._rule_prob_map[rule_set][rule_name])
                )
            logging.debug("------------------------------------------")

    def applyRecipe(self, config, recipe, cont=True):

        # Check if the rule in recipe is in the rules of the grammar
        sample = self.axiom.copy()
        for v in self.axiom.vertices():
            sample.vp.vertex_property[v] = copy.deepcopy(
                self.axiom.vp.vertex_property[v]
            )
        sample_id = str(uuid.uuid4())
        gen_tree = GenTree(sample_id)
        for num_recipe, recipe_entry in enumerate(recipe):
            if recipe_entry.name not in self.rules_map:
                print(
                    "We don't have the recipe_entry ",
                    recipe_entry.name,
                    " in our rules",
                )

            # Apply it a random time the number of times that the recipe dictates.
            num_ap = 1
            if recipe_entry.min_ap != recipe_entry.max_ap:
                # This is another parameter we can tune
                num_ap = random.randint(recipe_entry.min_ap, recipe_entry.max_ap)
            # print("rule: ", recipe_entry.name)
            for i in range(num_ap):
                # print("rule: ", recipe_entry.name)
                # Check if rule is applicable at the current time
                matchingRules = self._findMatchingRules(
                    sample, [self.rules_map[recipe_entry.name]]
                )

                if len(matchingRules) < 1:
                    # print("did not apply: " + )
                    # print("whops")
                    continue

                # if recipe_entry.name == "moveLockBack":
                #     print("we are applying movelockback")
                (rule, mapping) = self._pickRule(matchingRules, config, gen_tree)
                logging.debug("Applying rule " + rule.name)
                self._applyRule(sample, rule, mapping, config, gen_tree)
                # sample.vertex_properties["abbrevs"] = sample.new_vertex_property(
                #     "string"
                # )

                # for v in sample.vertices():
                #     sample.vp.abbrevs[v] = sample.vp.vertex_property[v].abbrev

                # gt.graph_draw(
                #     sample,
                #     vertex_text=sample.vp.abbrevs,
                #     output="test-rule-"
                #     + str(num_recipe)
                #     + "-"
                #     + rule.name
                #     + "-"
                #     + str(i)
                #     + "-khalifa.png",
                # )

        # print(gen_tree)
        # print(gen_tree.prob)
        # self._saveSampleGenTree(gen_tree)
        # self._saveSampleGraph(sample, sample_id)
        # self._drawSampleGraph(sample, sample_id)
        if cont:
            for i in range(int(config["max_applications"])):

                # check if there are no non-terminals left. If so break the loop we are done.
                # if False not in self.axiom.vp.terminality.get_array():
                #     print("No more non-terminals")
                # break

                matchingRules = self._findMatchingRules(sample, self.rules)

                if len(matchingRules) == 0:
                    # print("no matching rules")
                    break
                    # raise RuntimeError('No Rule match the given graph.')

                # How can we specify which rule to pick? Probabilities associated with each rule not just rhs
                (rule, mapping) = self._pickRule(matchingRules, config, gen_tree)
                logging.debug("Applying rule " + rule.name)
                self._applyRule(sample, rule, mapping, config, gen_tree)

        logging.debug("Out applyRecipe")
        return sample, sample_id, gen_tree

    """

    Transformation engine for graphs. Given a set of Rules of the form
    lhs ==> rhs, and using a starting graph G, uses graph isomorphic searching
    to find instances of a lhs in G and replaces the lhs vertices with the rhs.
    The engine continues to apply these transformations until G contains a
    given number of vertices. This assumes that the rules generally
    increase the number of vertices.

    Usage: Either use the all-inclusive main() method which checks the
    command-line arguments, or call applyRules() yourself which
    takes a starting graph, list of rules, and a dictionary of
    configuration options.
    """

    # write a method that copys the graph tuple, stores the derication tree and returns it. make it the only public method
    # --------------------------------------------------------------------------
    def applyRules(self, config):
        """
        Randomly applies a Rule from the given list of Rules to the
        specified starting graph until the graph contains at least the number
        of vertices specified by the config option "min_vertices". This assumes
        that the rules generally increase the number of vertices.
        Inputs:
            * start_graph_tuple - Tuple with the graph, abbrevs, labels of the start graph to begin applying transformations
            * rules - list of Rule objects
            * config - dictionary of options
        Outputs: None
        """

        logging.debug("In applyRules")
        sample = self.axiom.copy()
        for v in self.axiom.vertices():
            sample.vp.vertex_property[v] = copy.deepcopy(
                self.axiom.vp.vertex_property[v]
            )
        sample_id = str(uuid.uuid4())
        gen_tree = GenTree(sample_id)
        # check max application or if there are no longer any non-terminals
        for i in range(int(config["max_applications"])):

            # check if there are no non-terminals left. If so break the loop we are done.
            # if False not in self.axiom.vp.terminality.get_array():
            #     print("No more non-terminals")
            # break

            matchingRules = self._findMatchingRules(sample, self.rules)

            if len(matchingRules) == 0:
                # print("no matching rules")
                break
                # raise RuntimeError('No Rule match the given graph.')

            # How can we specify which rule to pick? Probabilities associated with each rule not just rhs
            (rule, mapping) = self._pickRule(matchingRules, config, gen_tree)
            logging.debug("Applying rule " + rule.name)
            self._applyRule(sample, rule, mapping, config, gen_tree)

        # print(gen_tree)
        # print(gen_tree.prob)
        # self._saveSampleGenTree(gen_tree)
        # self._saveSampleGraph(sample, sample_id)
        # self._drawSampleGraph(sample, sample_id)
        logging.debug("Out applyRules")
        return sample, sample_id, gen_tree

    # --------------------------------------------------------------------------
    # PRIVATE METHODS - These aren't the methods you're looking for.
    # --------------------------------------------------------------------------
    def _saveSampleGenTree(self, gen_tree):
        dir_path = "../data/generation-trees/" + self.name
        if not isdir(dir_path):
            mkdir(dir_path)

        gen_tree.saveTree(dir_path + "/" + gen_tree.sample_id + ".pkl")

    def _saveSampleGraph(self, sample_graph, sample_id):
        dir_path = "../data/sample-graphs/" + self.name
        if not isdir(dir_path):
            mkdir(dir_path)

        sample_graph.save(dir_path + "/" + sample_id + ".gt")

    def _drawSampleGraph(self, sample_graph, sample_id):
        dir_path = "../data/imgs/" + self.name
        if not isdir(dir_path):
            mkdir(dir_path)
        sample_graph.vertex_properties["abbrevs"] = sample_graph.new_vertex_property(
            "string"
        )

        for v in sample_graph.vertices():
            sample_graph.vp.abbrevs[v] = sample_graph.vp.vertex_property[v].abbrev

        gt.graph_draw(
            sample_graph,
            vertex_text=sample_graph.vp.abbrevs,
            output=dir_path + "/" + sample_id + ".png",
        )

    def _findRuleProperties(self, matchingRules):
        rule_map_duplicate = {}
        rule_probs = []
        rule_names = []
        for (i, (rule, _)) in enumerate(matchingRules):
            # print("Rule:", rule.name)
            rule_names.append(rule.name)
            rule_probs.append(rule.prob)
            if rule.name in rule_map_duplicate:
                rule_map_duplicate[rule.name].append({i: rule.prob})
            else:
                rule_map_duplicate[rule.name] = [{i: rule.prob}]

        return rule_probs, rule_names, rule_map_duplicate

    def _pickRule(self, matchingRules, config, gen_tree):
        if config["sampling_method"] == "grammar":
            rule_probs, rule_names = self._computeProps(matchingRules)
            rule_mapping = random.choices(
                population=matchingRules, weights=rule_probs, k=1
            )
            (rule, mapping) = rule_mapping[0]
            gen_tree.addNode(rule.name, frozenset(rule_names))
            rule_index = matchingRules.index(rule_mapping[0])
        else:
            rule_probs, rule_names, _ = self._findRuleProperties(matchingRules)
            rule_mapping = random.choice(matchingRules)
            (rule, mapping) = rule_mapping
            gen_tree.addNode(rule.name, frozenset(rule_names))
            rule_index = matchingRules.index(rule_mapping)

        gen_tree.prob *= rule_probs[rule_index]
        return (rule, mapping)

    def _computeProps(self, matchingRules):

        # DONT USE RULE PROB MAP IF SAMPLING METHOD IS UNIFORM
        rule_probs, rule_names, rule_map_duplicate = self._findRuleProperties(
            matchingRules
        )
        if frozenset(rule_names) in self._rule_prob_map:
            rule_prob_set = self._rule_prob_map[frozenset(rule_names)]
            # print(rule_prob_set)
            for rule_name in rule_prob_set:
                for rule_index_prob_map in rule_map_duplicate[rule_name]:
                    for index in rule_index_prob_map:
                        rule_index_prob_map[index] = rule_prob_set[rule_name]

        for rule_name in rule_map_duplicate:
            num_rule_matches = len(rule_map_duplicate[rule_name])
            for rule_index_prob_map in rule_map_duplicate[rule_name]:
                # print("rule_index_prob_map: ",rule_index_prob_map)
                for index in rule_index_prob_map:
                    rule_probs[index] = float(
                        rule_index_prob_map[index] / num_rule_matches
                    )

        rule_probs = [float(i) / sum(rule_probs) for i in rule_probs]
        return rule_probs, rule_names

    def _findMatchingRules(self, big_graph, rules):
        """
        Finds all the rules whose LHS graph can be found in graph. A
        rule LHS matches if the text-only abbrevs (e.g., "A") and the
        edges match (i.e., searching doesn't use the vertex number).
        Inputs:
            * graph - Graph tuple for graph to which to apply the rule
            * rules - list of Rule objects to search
        Outputs: list of (Rule, mapping) tuples where Rule
            is a Rule whose LHS can be found in graph, and mapping is
            a {vid->vid} dictionary (LHS->graph) of where the LHS can be found.
        """
        logging.debug("In _findMatchingRules")
        solutions = []
        for rule in rules:
            # logging.debug("Checking rule LHS %s " % rule.lhs.graph)
            logging.debug("Checking rule: " + rule.name + "'s LHS %s " % rule.lhs.graph)

            # find the graph isomorphism
            listOfMatches = gt.subgraph_isomorphism(
                # define a new vertex map with vector<string> of both access level and vertex label
                rule.lhs.graph,
                big_graph,
                # vertex_label=(rule.lhs.graph.vp.labels, big_graph.vp.labels)
            )
            if len(listOfMatches) > 0:
                for match in listOfMatches:
                    # print("match---------------------------- ")
                    # TODO: Go through the list of matches and check if each of
                    # the vertices that map from the sub to the bigger graph
                    # are y'know equal with respect to vertex property.
                    # print(match.get_array())
                    # returns the values, the index is the index of the subgraph
                    # just loop through the size of vertices of subgraph and get it from the match you dumb dumb
                    itMatches = True
                    minAccessLevel = 0
                    maxAccessLevel = 0
                    for u in rule.lhs.graph.vertices():
                        v = match[u]
                        minAccessLevel = min(
                            minAccessLevel, big_graph.vp.vertex_property[v].access_level
                        )
                        maxAccessLevel = max(
                            maxAccessLevel, big_graph.vp.vertex_property[v].access_level
                        )
                        rule.lhs.graph.vp.vertex_property[u].relative_access_level = 0
                    if minAccessLevel > 0:
                        print("Min access level is not 0 anymore!")
                        print(minAccessLevel)
                    # print("min: %s, max: %s" % (minAccessLevel, maxAccessLevel))
                    for u in rule.lhs.graph.vertices():
                        v = match[u]
                        # if rule.lhs.graph.vp.vertex_property[u] == big_graph.vp.vertex_property[v]:
                        # print('okay?')
                        # if rule.lhs.graph.vp.labels[u] != big_graph.vp.labels[v]:
                        # TODO print the the entire graph of matches and let's see their access levels
                        # rule.lhs.graph.vp.vertex_property[u].clearMatch()
                        # print(
                        #     "Abbrev of lhs: ",
                        #     rule.lhs.graph.vp.vertex_property[u].abbrev,
                        # )
                        # print(
                        #     "Abbrev of graph: ", big_graph.vp.vertex_property[v].abbrev
                        # )
                        # print("whoohoo")
                        rule.lhs.graph.vp.vertex_property[
                            u
                        ].relative_access_level = minAccessLevel

                        # if (
                        #     rule.lhs.graph.vp.vertex_property[u]
                        #     == big_graph.vp.vertex_property[v]
                        # ):

                        #     if big_graph.vp.vertex_property[v].label == "Start":
                        #         print("maxAccessLevel: ", maxAccessLevel)
                        #         print(
                        #             "self.relative_access_level: ",
                        #             rule.lhs.graph.vp.vertex_property[
                        #                 u
                        #             ].relative_access_level,
                        #         )
                        #         print(
                        #             "other.relative_access_level: ",
                        #             big_graph.vp.vertex_property[
                        #                 v
                        #             ].relative_access_level,
                        #         )
                        #         print(
                        #             "self.access_level: ",
                        #             rule.lhs.graph.vp.vertex_property[u].access_level,
                        #         )
                        #         print(
                        #             "other.access_level: ",
                        #             big_graph.vp.vertex_property[v].access_level,
                        #         )

                        # TODO:if vertex has degree than is four or more then it is false too
                        if (
                            rule.lhs.graph.vp.vertex_property[u]
                            != big_graph.vp.vertex_property[v]
                        ):

                            itMatches = False
                            break

                    if itMatches:
                        new_rule = copy.deepcopy((rule))
                        solutions.append((new_rule, match))
                        # solutions.append((rule, match))

                    logging.debug("Rule " + rule.name + " matches")
            else:
                logging.debug("Rule " + rule.name + " does not match")
                # logging.debug("Rule %s does not match" % rule.lhs.graph)
        logging.debug("Out _findMatchingRules")

        return solutions

    # --------------------------------------------------------------------------
    def _applyRule(self, graph, rule, lhs_mapping, config, gen_tree):
        """
        Applies the given rule to the given graph. The general idea is to
        transform the portion of the graph identified by mapping (which
        corresponds to the rule's LHS) to look like the RHS of the
        rule, by adding and/or removing vertices and edges.
        Inputs:
            graph - Graph tuple for graph to which to apply the rule
            rule - Rule to apply
            lhs_mapping - {vid->vid} mapping from rule.lhs
                to graph
        Outputs: None
        """
        # Pick an rhs, for now pick the random, but should pick according to the config

        logging.debug("In _applyRule")
        if config["sampling_method"] == "grammar":
            weights = []
            for rhs_c in rule.rhss:
                weights.append(rhs_c.prob)
            rhs_ls = random.choices(rule.rhss, weights=weights)
            rhs = rhs_ls[0]
            gen_tree.addRuleApplication(rule, rule.rhss.index(rhs))
            gen_tree.prob *= rhs.prob
        else:
            rhs = random.choice(rule.rhss)
            gen_tree.addRuleApplication(rule, rule.rhss.index(rhs))
            gen_tree.prob *= rhs.prob
        rhs_mapping, lhs_rhs_mapping = self._mapRHSToGraph(
            graph, rule.lhs, rhs, lhs_mapping
        )
        self._deleteMissingVertices(graph, rule.lhs, rhs, lhs_mapping, lhs_rhs_mapping)
        self._deleteMissingEdges(graph, rule.lhs, rhs, lhs_mapping, rhs_mapping)
        self._addNewVertices(graph, rule.lhs, rhs, rhs_mapping)
        self._addNewEdges(graph, rule.lhs, rhs, rhs_mapping)
        logging.debug("Out _applyRule")

    # --------------------------------------------------------------------------
    def _mapRHSToGraph(self, big_graph, lhs, rhs, lhs_mapping):
        """
        Maps to rule's rhs vertices to graph. For rhs vertices that
        appear in the lhs, we use the lhs_mapping to determine which graph
        vertex the rhs vertex maps to. For rhs vertices that are new (i.e.,
        don't exist in the lhs), we ignore them for now and update the
        rhs mapping when we add the new vertices to graph. LHS and RHS
        vertices are considered the same if both their label and their
        number match.
        Inputs:
            * graph - Graph tuple for graph to which to apply the rule
            * rule - Rule to apply
            * lhs2GraphMapping - {vid->vid} mapping from rule.lhs
              vertices to graph
        Outputs: {vid->vid} mapping from rule.rhs vertices to graph
        """

        lhs_rhs_mapping = {}
        rhs_mapping = {}
        for rhs_v in rhs.graph.vertices():
            # clear match info for rhs
            rhs.graph.vp.vertex_property[rhs_v].clearMatch()
            rhs.graph.vp.vertex_property[rhs_v].relative_access_level = 0
            for lhs_v in lhs.graph.vertices():
                rhs.graph.vp.vertex_property[
                    rhs_v
                ].relative_access_level = lhs.graph.vp.vertex_property[
                    lhs_v
                ].relative_access_level
                if (
                    rhs.graph.vp.vertex_property[rhs_v].mark
                    == lhs.graph.vp.vertex_property[lhs_v].mark
                ):
                    match = big_graph.vp.vertex_property[lhs_mapping[lhs_v]]
                    # adjust the rhs' access level based on the lhs
                    # print("changing RHS relative accesss")
                    # print(
                    #     "rhs access rel: "
                    #     + str(rhs.graph.vp.vertex_property[rhs_v].relative_access_level)
                    # )
                    # print(
                    #     "lhs access rel: "
                    #     + str(lhs.graph.vp.vertex_property[lhs_v].relative_access_level)
                    # )
                    # print(
                    #     "Matching: "
                    #     + "lhs: "
                    #     + str(lhs.graph.vp.vertex_property[lhs_v])
                    #     + " "
                    #     + "rhs: "
                    #     + str(rhs.graph.vp.vertex_property[rhs_v])
                    #     + " "
                    #     + " with: "
                    #     + str(big_graph.vp.vertex_property[lhs_mapping[lhs_v]])
                    # )
                    # print(
                    #     "self.relative_access_level: ",
                    #     lhs.graph.vp.vertex_property[lhs_v].relative_access_level,
                    # )
                    # print(
                    #     "other.relative_access_level: ",
                    #     big_graph.vp.vertex_property[
                    #         lhs_mapping[lhs_v]
                    #     ].relative_access_level,
                    # )
                    # print(
                    #     "self.access_level: ",
                    #     lhs.graph.vp.vertex_property[lhs_v].access_level,
                    # )
                    # print(
                    #     "other.access_level: ",
                    #     big_graph.vp.vertex_property[lhs_mapping[lhs_v]].access_level,
                    # )
                    # print("----------------------")
                    rhs.graph.vp.vertex_property[rhs_v].match(
                        big_graph.vp.vertex_property[lhs_mapping[lhs_v]]
                    )
                    rhs_mapping[rhs_v] = lhs_mapping[lhs_v]
                    lhs_rhs_mapping[lhs_v] = rhs_v

        return rhs_mapping, lhs_rhs_mapping

    # --------------------------------------------------------------------------
    def _deleteMissingVertices(self, graph, lhs, rhs, lhs_mapping, lhs_rhs_mapping):
        """
        Deletes vertices from graph that appear in rule.lhs but not in
        rule.rhs. It also deletes all edges to lead to or from the
        deleted vertex.
        Inputs:
            * graph - Graph tuple for graph to which to apply the rule
            * rule - Rule to apply
            * lhs_mapping - {vid->vid} mapping between rule.lhs
                    and graph
        Outputs: None
        """
        logging.debug(">>> _deleteMissingVertices <<<")
        for lhs_v in lhs.graph.vertices():
            rhs_v = lhs_rhs_mapping[lhs_v]
            if not rhs.graph.vp.vertex_property[rhs_v]:
                graph_vertex = lhs_mapping[lhs_v]
                logging.debug(
                    "deleting vertex %s" % lhs.graph.vp.vertex_property[lhs_v]
                )
                graph.remove_vertex(graph_vertex)

        return

    # --------------------------------------------------------------------------
    def _deleteMissingEdges(self, graph, lhs, rhs, lhs_mapping, rhs_mapping):
        """
        Deletes edges from graph that appear in rule.lhs but not in
        rule.rhs. Assumes new vertices on rhs have been added to
        graph.
        Inputs:
            * graph_tuple - Graph tuple for graph to which to apply the rule
            * rule - Rule to apply
            * lhs_mapping - {vid->vid} mapping from rule.lhs
              to graph
            * rhs_mapping - {vid->vid} mapping from rule.rhs
              to graph
        Outputs: None
        """
        logging.debug(">>> _deleteMissingEdges <<<")
        for lhsEdge in lhs.graph.edges():  # [startVertex,endVertex]
            # Find the starting and ending vertices of the corresponding edge in graph.
            u, v = lhsEdge
            graphStartVID = lhs_mapping[u]
            graphEndVID = lhs_mapping[v]

            # Try to find the corresponding starting vertex in the rhs (if it
            # it exists at all). If it doesn't even exist, then the edge
            # doesn't exist either, so delete it from graph.
            rhsStart = [
                rhsID
                for rhsID, graphID in rhs_mapping.items()
                if graphID == graphStartVID
            ]
            if len(rhsStart) == 0:
                logging.debug(
                    "edge start from %s to %s does not appear in rhs"
                    % (
                        lhs.graph.vp.vertex_property[u].abbrev,
                        lhs.graph.vp.vertex_property[v].abbrev,
                    )
                )
                edge_to_be_deleted = graph.edge(graphStartVID, graphEndVID)
                graph.remove_edge(edge_to_be_deleted)
                continue

            # Try to find the corresponding ending vertex in the rhs (if it
            # it exists at all). If it doesn't even exist, then the edge
            # doesn't exist either, so delete it from graph.
            rhsEnd = [
                rhsID
                for rhsID, graphID in rhs_mapping.items()
                if graphID == graphEndVID
            ]
            if len(rhsEnd) == 0:
                logging.debug(
                    "edge start from %s to %s does not appear in rhs"
                    % (
                        lhs.graph.vp.vertex_property[u].abbrev,
                        lhs.graph.vp.vertex_property[v].abbrev,
                    )
                )
                edge_to_be_deleted = graph.edge(graphStartVID, graphEndVID)
                graph.remove_edge(edge_to_be_deleted)
                continue

            edge_exists = False
            for edge in rhs.graph.edges():
                u_r, v_r = edge
                if u_r == rhsStart[0] and v_r == rhsEnd[0]:
                    edge_exists = True
                    break

            # We found both rhs vertices, but are they connected with an
            # edge? If not, the delete the edge from graph.
            if not edge_exists:
                logging.debug(
                    "edge start from %s to %s does not appear in rhs"
                    % (
                        lhs.graph.vp.vertex_property[u].abbrev,
                        lhs.graph.vp.vertex_property[v].abbrev,
                    )
                )
                logging.debug(
                    "deleting edge from %s to %s"
                    % (
                        graph.vp.vertex_property[graphStartVID].abbrev,
                        graph.vp.vertex_property[graphEndVID].abbrev,
                    )
                )
                edge_to_be_deleted = graph.edge(graphStartVID, graphEndVID)

                graph.remove_edge(edge_to_be_deleted)

        logging.debug("graph is now %s" % graph)
        return

    # --------------------------------------------------------------------------
    def _addNewVertices(self, graph, lhs, rhs, rhs_mapping):
        """
        Adds vertices to graph that appear in rule.rhs but not in
        rule.lhs. New vertices are given a vid of the form 'vN' where
        N is the number of vertices currently in the graph. New graph vertices
        are also added to rhs_mapping.
        Inputs:
            * graph - Graph tuple for graph to which to apply the rule
            * rule - Rule to apply
            * rhs_mapping - {vid->vid} mapping from rule.rhs to graph.
              This is typically created by _mapRHSToGraph().
        Outputs: nothing
        """
        logging.debug(">>> _addNewVertices <<<")
        # TODO: If a wildcard is used, make the value of the of the graph is original graph's value not Any
        # that means whatever the lhs wildcard matches with, it should maintain that maht throughout the
        # application of the rule and it should figure out which rhs wildcard it corresponds to and do the same

        for rhs_v in rhs.graph.vertices():
            vertex_exists = False
            if rhs_v in rhs_mapping:
                # change its properties in the original graph to match that of the rhs
                v = rhs_mapping[rhs_v]
                rhs_vp = rhs.graph.vp.vertex_property[rhs_v]
                # graph.vp.abbrevs[v] = rhs.graph.vp.abbrevs[rhs_v]
                # graph.vp.labels[v] = rhs.graph.vp.labels[rhs_v]
                # print("setting vertes properties are you")
                #
                # print("original: ", graph.vp.vertex_property[v].abbrev)
                # print("rhs': ", rhs.graph.vp.vertex_property[rhs_v].abbrev)
                # print(
                #     "rhs' matched: ",
                # rhs.graph.vp.vertex_property[rhs_v].getMatched().abbrev,
                # )
                access_level_original = graph.vp.vertex_property[v].access_level
                # print(
                #     "Replacing "
                #     + str(graph.vp.vertex_property[v])
                #     + " with: "
                #     + str(rhs_vp.getMatched())
                # )
                graph.vp.vertex_property[v] = rhs_vp.getMatched()
                new_level = rhs_vp.relative_access_level + rhs_vp.access_level
                self._adjustAccessLevel(graph, v, new_level)
                access_level_after = graph.vp.vertex_property[v].access_level
                # print(
                #     "Actually ignore that, the node is actually : ",
                #     graph.vp.vertex_property[v],
                # )

                # print(
                #     "rhs relative_access_level: ",
                #     rhs.graph.vp.vertex_property[rhs_v]
                #     .getMatched()
                #     .relative_access_level,
                # )
                rhs_vp.clearMatch()
                # print("new: ", graph.vp.vertex_property[v].abbrev)
                # print(graph.vp.vertex_property[v].abbrev)
                # graph.vp.terminality[v] = rhs.graph.vp.terminality[rhs_v]
                vertex_exists = True
            else:
                logging.debug(
                    "name %s in rhs but not lhs"
                    % rhs.graph.vp.vertex_property[rhs_v].abbrev
                )
                v = graph.add_vertex()
                rhs_vp = rhs.graph.vp.vertex_property[rhs_v]
                # graph.vp.abbrevs[v] = rhs.graph.vp.abbrevs[rhs_v]
                # graph.vp.labels[v] = rhs.graph.vp.labels[rhs_v]
                # print("Adding " + str(rhs_vp.getMatched()))
                graph.vp.vertex_property[v] = rhs_vp.getMatched()
                new_level = rhs_vp.relative_access_level + rhs_vp.access_level
                self._adjustAccessLevel(graph, v, new_level)
                # graph.vp.terminality[v] = rhs.graph.vp.terminality[rhs_v]
                rhs_vp.clearMatch()
                rhs_mapping[rhs_v] = v
                logging.debug("added vertex %s" % graph.vp.vertex_property[v].label)
        logging.debug("graph is now %s" % graph)
        return

    # --------------------------------------------------------------------------
    def _addNewEdges(self, graph, lhs, rhs, rhs_mapping):
        """
        TODO:goes thorugh more edges than needed, should check against lhs.edges not graph.edges
        Adds edges to graph that appear in rule.rhs but not in
        rule.lhs. Assumes that all the new vertices in the rule.rhs
        have been added to graph, and rhs_mapping has been updated accordingly.
        Inputs:
            * graph - Graph tuple for graph to which to apply the rule
            * rule - Rule to apply
            * rhs_mapping - {vid->vid} mapping between rule.rhs
              and graph
        Outputs: None
        """
        logging.debug(">>> _addNewEdges <<<")
        for rhsEdge in rhs.graph.edges():
            u_r, v_r = rhsEdge
            edge_exists = False
            for edge in graph.edges():
                u, v = edge
                if rhs_mapping[u_r] == u and rhs_mapping[v_r] == v:
                    edge_exists = True
                    break
            if not edge_exists:
                logging.debug(
                    "adding edge (%s,%s) to graph"
                    % (
                        rhs.graph.vp.vertex_property[u_r].abbrev,
                        rhs.graph.vp.vertex_property[v_r].abbrev,
                    )
                )
                graph.add_edge(rhs_mapping[u_r], rhs_mapping[v_r])
        logging.debug("graph is now %s" % graph)

        return

    # debug, info, warning, error and critical
    # if __name__ == '__main__':
    #     if len(sys.argv) != 2:
    #         print("Usage: %s GRAMMAR_FILE" % sys.argv[0], file = sys.stderr)
    #         sys.exit(1)
    #     e = Generator()
    #     e.generateFromFile(sys.argv[1])

    # make sure new_level is  self.access_level + self.relative_access_level of the rhs
    def _adjustAccessLevel(self, graph, u, new_level):
        graph_vp = graph.vp.vertex_property[u]
        if graph_vp.access_level == new_level:
            return
        graph_vp.access_level = new_level
        if graph_vp.label != "Lever":
            children = graph.get_out_neighbors(u)
            for c in children:
                if (
                    graph.vp.vertex_property[c].label == "Lock"
                    or graph.vp.vertex_property[u].label == "Puzzle"
                ):
                    self._adjustAccessLevel(graph, c, new_level + 1)

                else:
                    self._adjustAccessLevel(graph, c, new_level)
