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

def scan(list):
	ret = []
	for i in range(len(list)):
		if list[i].is_leaf == True:
			ret.append(list[i])
	return ret

def main():
	opts = parse_args()
	with open(opts.inputTree) as f:
		trees = newick.load(f)

	for tree in trees:
		print tree.ascii_art()
		po = post_order(tree)
		pre = pre_order(tree)
		pre_dict = list_dict(pre)
		print po
	for node in po:
		print node.is_leaf
	print scan(po)



	print 'trees is %s' % trees
	print newick.dumps(trees)
if __name__ == "__main__":
	main()

