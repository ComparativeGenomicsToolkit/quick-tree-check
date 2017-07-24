import newick
from argparse import ArgumentParser
from collections import deque 
from copy import deepcopy
# ^ equivalent to "import argparse" and using "argparse.ArgumentParser"


def parse_args():
	parser = ArgumentParser(description=__doc__)
	parser.add_argument('inputTree', help='newick format tree (in a file)')
	return parser.parse_args()




def post_order(tree):
	"""Return a list of nodes in tree in post-order."""
	ret = []
	for child in tree.descendants:
		ret.extend(post_order(child))
	ret.append(tree)
	return ret
def pre_order(tree):
	ret=[]
	ret.append(tree)
	for child in tree.descendants:
		ret.extend(pre_order(child))	
	return ret




def list_dict(list):
	dict = {}
	branch_val = 0
	for i in list:
		dict[i] = branch_val
		branch_val += 1
	return dict

def scroll(list):
	for i in range(len(list)-1):
		for j in range(i+1,len(list)):
			print list[i],list[j]



def scan_leaves(node):
	ret = []
	for i in range(len(node)):
		if node[i].is_leaf == True:
			ret.append(node[i])
	return ret

def all_ancestors(node):
	anc = []
	current_node = node
	while current_node.ancestor.ancestor != None:
		anc.append(current_node.ancestor)
		current_node = current_node.ancestor
	return anc

def ancestor_list(node):
	ret = []
	i=0
	leaves = scan_leaves(node)
	for i in range (i,len(leaves)):
		ret.append(all_ancestors(leaves[i]))
	return ret

def distance(anc,leaves):
	list_1 = []
	new_list = deepcopy(anc)
	for i in range(len(anc)-1):
		for j in range(i+1,len(anc)):
			ret = deque([])
			if len(anc[j]) < len(anc[i]):
				for k in range(len(anc[j])-1,-1,-1):
					if anc[i][k+len(anc[i])-len(anc[j])] == anc[j][k]:
						anc[i].pop()
						anc[j].pop()
				ret.extend(anc[i])
				ret.extend(anc[j])
			else:
				for k in range(len(anc[i])-1,-1,-1):
					if anc[j][k+len(anc[j])-len(anc[i])] == anc[i][k]:
						anc[i].pop()
						anc[j].pop()
				ret.extend(anc[i])
				ret.extend(anc[j])
			ret.appendleft(leaves[i])
			ret.append(leaves[j])
			list_1.append(ret)
			anc = deepcopy(new_list)
	return list_1



			



def main():
	opts = parse_args()
	with open(opts.inputTree) as f:
		tree = newick.load(f)

	

	for node in tree:
		po = post_order(node)
		pre = pre_order(node)
		pre_dict = list_dict(pre)
		print node.ascii_art()
		print po
		ancA = ancestor_list(po)
		"""win = naked_problem(ancA)"""
		print "ancA is: %s" % ancA
		"""print "win is: %s" % win"""
	for node in tree:
		leafs = scan_leaves(po)
		win = distance(ancA,leafs)
		print "win is: %s" % win
		

	"""for node in po:
		print node.is_leaf
	for node in po:
		print node.ancestor
	print scan(po)"""




	print 'trees is %s' % tree
	print newick.dumps(tree)
if __name__ == "__main__":
	main()

