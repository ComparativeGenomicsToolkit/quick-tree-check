import newick
from argparse import ArgumentParser
from collections import deque 
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




def scan_leaves(list):
	ret = []
	for i in range(len(list)):
		if list[i].is_leaf == True:
			ret.append(list[i])
	return ret

def all_ancestors(node):
	anc = []
	current_node = node
	while current_node.ancestor != None:
		anc.append(current_node.ancestor)
		current_node = current_node.ancestor
	return anc

def ancestor_list(list):
	ret = []
	i=1
	leaves = scan_leaves(list)
	for i in range (i,len(leaves)):
		ret.append(all_ancestors(leaves[i]))
	return ret




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
		print "ancA is: %s" % ancA
	
		

	"""for node in po:
		print node.is_leaf
	for node in po:
		print node.ancestor
	print scan(po)"""




	print 'trees is %s' % tree
	print newick.dumps(tree)
if __name__ == "__main__":
	main()

