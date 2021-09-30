from os.path import isfile, join, isdir


class RecipeEntry(object):
    def __init__(self, name, min_ap, max_ap):
        self.name = name
        self.min_ap = min_ap
        self.max_ap = max_ap

    def __str__(self):
        return self.name + "," + str(self.min_ap) + "," + str(self.max_ap) + "\n"

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def read_file(filename):
        if not isfile(filename):
            print("filename does not exist")
            return []

        input_file = open(filename, "r")
        lines = input_file.readlines()
        recipe = []
        for line in lines:
            components = line.split(",")
            if len(components) < 1:
                continue
            name = components[0].strip()
            if len(components) < 3:
                min_ap = 1
                max_ap = 1
                recipe.append(RecipeEntry(name, min_ap, max_ap))
            else:
                min_ap = int(components[1])
                max_ap = int(components[2])
                recipe.append(RecipeEntry(name, min_ap, max_ap))

        return recipe


if __name__ == "__main__":
    recipe = RecipeEntry.read_file("khalifa-grammar/graphRecipe.txt")
    print(recipe)
