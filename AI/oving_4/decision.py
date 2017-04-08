#Exercise 4
#Create tree: http://ironcreek.net/phpsyntaxtree/

import math as m
import random as r
debug = True

class node():
	def __init__(self, data):
		self.data = data
		self.children = {}

	def print_Tree(self):
		if debug:
			print 'len', len(self.children)
			print 'self', self.data
			print 'children', self.children

		if len(self.children) == 0:
			return "[" + str(self.data) + "]"
		else:
			temp = "[" + str(self.data) + " "

		for key, dont_need_but_have_to_put_in_something in self.children.items():
			temp += self.children[key].print_Tree()

		if debug: print "*******************************************************"
		return temp + "]"

#Reads data sets from file
def read_from_file(name):
	liste = []
	fil = open(name, 'r')
	for line in fil.readlines():
		liste.append(line.rstrip("\n").split("\t"))
	return liste

#prints the data set. For debuging and testing during development
def print_list(l):
	for line in l:
		print line

#Finds the most commmon goal/class. If equal chooses at random(1 or 2)
def Plurality(examples):
	ones_count = 0
	twos_count = 0

	for n in xrange(0, len(examples)):
		if examples[-1] == '1':
			ones_count += 1
		else:
			twos_count += 1

	if ones_count > twos_count:
		return 1
	elif twos_count > ones_count:
		return 2
	else:
		return r.randint(1,2)

#checks if all examples has the same class
def same_classification(examples):
	classification = examples[0][-1]
	for n in xrange(1, len(examples)):
		if examples[n][-1] != classification:
			return False
	return True

#Calculate the B value from the book(page 715).
def B(q):
	#print "q", q
	if q == 0:
		return q
	else:
		return -( ( q*m.log(q, 2) ) + ( (1.0-q)*m.log((1.0-q), 2) ) )

#Chooses which attribute is the most important based on (information theory)entropy
def importance(data, attributes):
	#print "attributelistavitarinn", attributes
	attribute_entropy = {}
	#print "attribute list", attribute_entropy
	for attribute in attributes:
		count = 0
		for liste in data:
			if liste[attribute] == data[0][attribute]:
				count += 1
		#print "attribute", attribute
		attribute_entropy[attribute] = B(count/len(data) )
		if debug: print "B", attribute_entropy[attribute], count

	if debug: print "entropy", attribute_entropy


	minimum = 1.1
	index = None
	for n in attribute_entropy:
		if attribute_entropy[n] < minimum:
			minimum = attribute_entropy[n]
			index = n
	#rint "N", index
	#print '--------------------------------------------------------------------------------'
	return index

#the method that generates the tree based on the training set
def decision_tree_learinng(examples, attributes, parent_examples, random_importance):
	if debug: print"-------------------------------------------------------------------"

	if not examples:
		return node(Plurality(parent_examples))
	elif same_classification(examples):
		return node(examples[0][-1])
	elif not attributes:
		return node(Plurality(examples))
	else:
		if random_importance:
			A = attributes[ r.randint(0, len(attributes)-1) ]
		else:
			A = importance(examples, attributes)

		tree = node(A)
		attributes.remove(A)

		if debug:
			print "A", A
			print "attributes after remove", attributes

		for n in xrange(1,3):
			liste = []
			for e in examples:
				if int(e[A]) == n:
					liste.append(e)
			sub_tree = decision_tree_learinng(liste, list(attributes), examples, random_importance)
			tree.children[n] = sub_tree
	return tree


def classify(root, line):
        current = root
        while current.children:
                current = current.children[int(line[current.data])]
        return current.data


def run_tests(tree, data):
        correct_tests = 0
        for line in data:
                if line[-1] == classify(tree, line):
                        correct_tests += 1
        print "Correct / number of tests = ", correct_tests, "/", len(data)

def main():
	test = read_from_file("data/test.txt")
	training = read_from_file("data/training.txt")

	if debug:
		print_list(test)
		print "\n"+"\n"
		print_list(training)
		print "************************************************************************"

	tree = decision_tree_learinng(training, range( len(training[0])-1 ), [], False)
	print tree.print_Tree()
	run_tests(tree, test)
	print "random: "
	tree = decision_tree_learinng(training, range( len(training[0])-1 ), [], True)
	run_tests(tree, test)
	print tree.print_Tree()

main()
