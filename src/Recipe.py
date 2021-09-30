class RHS(object):
    def __init__(self, graph, prob):
        """
        Constructor.
        Inputs:
            * graph - Graph object on the LHS.
            * abbrevs - Vertex labels
            * prob - Probablity
        Outputs: N/A
        """
        self.graph = graph
        self.prob = prob

    def __str__(self):

        return_str = str(self.graph) + ': ' + str(self.prob)
        return return_str

    def __repr__(self):
        return self.__str__()


class LHS(object):
    def __init__(self, graph):
        """
        Constructor.
        Inputs:
            * graph - Graph object on the LHS.
            * abbrevs - Vertex labels
            * prob - Probablity
        Outputs: N/A
        """
        self.graph = graph

    def __str__(self):

        return_str = str(self.graph)
        return return_str

    def __repr__(self):
        return self.__str__()

class Rule(object):
    """
    Represents a rule consisting of a left-hand side and a list of right-hand
    sides and a probability.
    """

    def __init__(self, name, lhs, rhss, prob):
        """
        Constructor.
        Inputs:
            * lhs - Graph object on the LHS.
            * rhss - List of RHSs
            * prob - Probablity
        Outputs: N/A
        """
        self.name = name
        self.lhs = lhs
        self.rhss = rhss
        self.prob = prob

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

# if __name__ == '__main__':
#     print('Hi')
