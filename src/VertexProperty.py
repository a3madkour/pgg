#!/usr/bin/env python3
class VertexProperty(object):
    wildcard_labels = ["Any"]

    def __init__(
        self,
        label,
        abbrev,
        terminality,
        mark,
        access_level=0,
        relative_access_level=0,
        is_wildcard=False,
        matched_vertex=None,
    ):
        """
        Constructor.
        Inputs:
            * label - Probablity
        Outputs: N/A
        """
        self._label = label
        self._abbrev = abbrev
        self._terminality = terminality
        self._mark = mark
        self._access_level = access_level
        self._relative_access_level = relative_access_level
        self._is_wildcard = is_wildcard
        self._matched_vertex = matched_vertex
        if label in VertexProperty.wildcard_labels:
            self._is_wildcard = True

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    @property
    def access_level(self):
        return self._access_level

    @access_level.setter
    def access_level(self, value):
        self._access_level = value

    @property
    def abbrev(self):
        # return self._abbrev + str(self.access_level)
        return self._abbrev

    @abbrev.setter
    def abbrev(self, value):
        self._abbrev = value

    @property
    def relative_access_level(self):
        return self._relative_access_level

    @relative_access_level.setter
    def relative_access_level(self, value):
        self._relative_access_level = value

    @property
    def terminality(self):
        return self._terminality

    @terminality.setter
    def terminality(self, value):
        self._terminality = value

    @property
    def mark(self):
        return self._mark

    @mark.setter
    def mark(self, value):
        self._mark = value

    @property
    def is_wildcard(self):
        return self._is_wildcard

    @is_wildcard.setter
    def is_wildcard(self, value):
        self._is_wildcard = value

    def __str__(self):

        return_str = (
            str(self.label)
            + ": "
            + str(self.access_level)
            + ": "
            + str(self.relative_access_level)
        )
        # return_str = str(self.label) + ": " + str(self.access_level)
        return return_str

    def __repr__(self):
        return self.__str__()

    def __copy__(self):
        return VertexProperty(
            self._label,
            self._abbrev,
            self._terminality,
            self._mark,
            self._access_level,
            self._relative_access_level,
            self._is_wildcard,
            self._matched_vertex,
        )

    def __deepcopy__(self, memo):
        # print(memo)
        return self.__copy__()

    def _checkAccessLevel(self, other):
        if self.access_level == (other.access_level - self.relative_access_level):
            # if self.label == "Start" or other.label == "Start":
            # print("self.relative_access_level: ", self.relative_access_level)
            # print("other.relative_access_level: ", other.relative_access_level)
            # print("self.access_level: ", self.access_level)
            # print("other.access_level: ", other.access_level)
            return True
        else:
            return False

    def _checkLabel(self, other):
        if self.label == other.label:
            return True
        else:
            return False

    def __eq__(self, other):
        if self.is_wildcard or other.is_wildcard:
            return self._checkAccessLevel(other)
        else:
            return self._checkLabel(other) and self._checkAccessLevel(other)

    def match(self, other):
        if not self.is_wildcard:
            self._matched = None
            return
        self._matched = other
        return

    def clearMatch(self):
        self._matched = None

    def getMatched(self):
        if self._matched and self.is_wildcard:
            self._matched.relative_access = 0
            return self._matched
        else:
            return self
