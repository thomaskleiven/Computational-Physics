import math as m
import random

rand = True


def main():
	test = getData("data/test.txt")
	training = getData("data/training.txt")

	tree = dct(training, range( len(training[0])-1 ), [])
	print tree.print_Tree()
	run_tests(tree, test)
	#print "random: "
	#tree = dct(training, range( len(training[0])-1 ), [], True)
	#run_tests(tree, test)
	#print tree.print_Tree()


def dct(examples, attributes, parent_examples):
	if not examples:
		return Node(plurality_value(parent_examples))
	elif classification(examples):
		return Node(examples[0][-1])
	elif not attributes:
		return Node(plurality_value(examples))
	else:
		A = getImportance(examples, attributes)

		tree = Node(A)
		attributes.remove(A)

		for n in range(1,3):
			liste = []
			for e in examples:
				if int(e[A]) == n:
					liste.append(e)
			sub_tree = dct(liste, list(attributes), examples)
			tree.children[n] = sub_tree
	return tree


def getData(filename):
    input_file = open(filename, 'r')
    data = []
    for line in input_file.readlines():
        data.append(line.rstrip("\n").split("\t"))
    return data

def plurality_value(examples):
	ones = 0
	twos = 0

	for n in range(0, len(examples)):
		if examples[-1] == '1':
			ones += 1
		else:
			twos += 1

	if ones > twos:
		return 1
	elif twos > ones:
		return 2
	else:
		return random.randint(1,2)

def classification(examples):
	classification = examples[0][-1] #Last element
	for n in range(1, len(examples)):
		if examples[n][-1] != classification:
			return False
	return True


def getNewEntropy(q):
	if q >= 1 or q <= 0: return 0
	return -(q * np.log2(q) + ((1.0-q) * np.log2(1.0-q)))


def getImportance(data, attributes):
	if rand:
		return attributes[ random.randint(0, len(attributes)-1) ]
	attribute_entropy = {}
	for attribute in attributes:
		num = 0
		for liste in data:
			if liste[attribute] == data[0][attribute]:
				num += 1
		attribute_entropy[attribute] = getNewEntropy(num/len(data) )

	minimum = 1.1 #Sum of propabilities, has to be less than or equal to 1
	loc = None
	for n in attribute_entropy:
		if attribute_entropy[n] < minimum:
			minimum = attribute_entropy[n]
			loc = n
	return loc


def compare(root, line):
        current = root
        while current.children:
                current = current.children[int(line[current.data])]
        return current.data


def run_tests(tree, data):
        correct_tests = 0
        for line in data:
                if line[-1] == compare(tree, line):
                        correct_tests += 1
        print "Correct / number of tests = ", correct_tests, "/", len(data)


class Node():
	def __init__(self, data):
		self.data = data
		self.children = {}

	def print_Tree(self):
		if len(self.children) == 0:
			return "[" + str(self.data) + "]"
		else:
			temp = "[" + str(self.data) + " "

		for key, value in self.children.items():
			temp += self.children[key].print_Tree()

		return temp + "]"

if __name__ == "__main__":
	main()
