import pickle
class Node(object):
    def __init__(self, rule_id, rule_set):
        """
        Constructor.
        Inputs:
            * rule_id - The id of the rule chosen at this step
            * rule_set - The set of rules available at this step
        Outputs: N/A
        """
        self._rule_id = rule_id
        self._rule_set = rule_set
        self._rhs_index = None
        self._next_node = None

    def __str__(self):

        return_str = str(self.rule_id)
        return  return_str

    @property
    def rule_id(self):
        return self._rule_id
    @rule_id.setter
    def rule_id(self, value):
        self._rule_id = value

    @property
    def rule_set(self):
        return self._rule_set
    @rule_set.setter
    def rule_set(self, value):
        self._rule_set = value

    @property
    def rhs_index(self):
        return self._rhs_index
    @rhs_index.setter
    def rhs_index(self, value):
        self._rhs_index = value

    @property
    def next_node(self):
        return self._next_node
    @next_node.setter
    def next_node(self, value):
        self._next_node = value

class GenTree(object):
    def __init__(self, sample_id = None):
        """
        Constructor.
        Inputs:
        Outputs: N/A
        """
        self._sample_id = sample_id
        self._inital_node = None
        self._current_node = None


        self._rule_set_map = {}
        self._rule_map = {}
        self._prob = 1

    def __str__(self):
        return_str = ""
        node = self.inital_node
        rule_ids = []
        while node is not None:
            rule_ids.append(node.rule_id)
            node = node.next_node

        for rule_id in reversed(rule_ids):
            return_str = return_str + rule_id + ">"
        return  return_str[:-1]+"\n"

    def __repr__(self):
        return self.__str__()

    @property
    def inital_node(self):
        return self._inital_node
    @inital_node.setter
    def inital_node(self, value):
        self._inital_node = value


    @property
    def current_node(self):
        return self._inital_node
    @current_node.setter
    def current_node(self, value):
        self._inital_node = value

    @property
    def rule_map(self):
        return self._rule_map
    @rule_map.setter
    def rule_map(self, value):
        self._rule_map= value

    @property
    def rule_set_map(self):
        return self._rule_set_map
    @rule_set_map.setter
    def rule_set_map(self, value):
        self._rule_set_map= value


    @property
    def sample_id(self):
        return self._sample_id
    @sample_id.setter
    def sample_id(self, value):
        self._sample_id= value

    @property
    def prob(self):
        return self._prob
    @prob.setter
    def prob(self, value):
        self._prob= value


    def addNode(self, rule_id, rule_set):
        if not self.inital_node:
            self.inital_node = Node(rule_id,rule_set)
            self.current_node = self.inital_node
        else:
            nx = Node(rule_id, rule_set)
            nx.next_node = self.current_node
            self.current_node = nx

        if rule_set not in self.rule_set_map:
           self.rule_set_map[rule_set] = {"rules": {}, "total_count": 0}

        if rule_id in self.rule_set_map[rule_set]["rules"]:
            self.rule_set_map[rule_set]["rules"][rule_id]  =  self.rule_set_map[rule_set]["rules"][rule_id]+ 1
        else:
            self.rule_set_map[rule_set]["rules"][rule_id] = 1

        self.rule_set_map[rule_set]["total_count"] = self.rule_set_map[rule_set]["total_count"] +1
                                

    def drawTree(self, img_filename):
        return

    def saveTree(self, filename):
        outfile = open(filename, 'wb')
        pickle.dump(self,outfile)
        outfile.close()

    def addRuleApplication(self, rule, selected_rhs_index):
        if rule.name not in self.rule_map:
            self.rule_map[rule.name] = {"rhs_index": {}, "total_count": 0}

            for (i,rhs) in enumerate(rule.rhss):
                self.rule_map[rule.name]["rhs_index"][i] = 0

        self.rule_map[rule.name]["rhs_index"][selected_rhs_index]  = self.rule_map[rule.name]["rhs_index"][selected_rhs_index] + 1

        self.rule_map[rule.name]["total_count"] = self.rule_map[rule.name]["total_count"] +1

        self.current_node.rhs_index = selected_rhs_index
