import newick
from argparse import ArgumentParser
from collections import deque 
from copy import copy
from itertools import permutations
import numpy as np
from numpy.linalg import inv
# ^ equivalent to "import argparse" and using "argparse.ArgumentParser"

def read_distance_matrix(file, post_order):
        """Read a distance matrix from a file.

        The distance matrix should be in the following format, tab-separated:
        A<TAB>B<TAB>d_ab
        A<TAB>C<TAB>d_ac
        A<TAB>D<TAB>d_ad
        B<TAB>C<TAB>d_bc
        B<TAB>D<TAB>d_bd
        C<TAB>D<TAB>d_cd

        The order of the rows doesn't matter. Providing both d_ab and
        d_ba is allowed, but if they aren't the same, an exception
        will be raised.
        """
        # Get leaves
        leaf_names = [n.name for n in post_order if len(n.descendants) == 0]

        # Get leaf name -> index mapping
        name_to_index = dict((leaf, i) for i, leaf in enumerate(leaf_names))

        num_leaves = len(name_to_index)
        matrix = np.zeros((num_leaves, num_leaves))

        # Helper function to give a useful error on a missing key
        def get_leaf_index(name):
                try:
                        return name_to_index[name]
                except KeyError:
                        raise RuntimeError("Didn't find a leaf with name '%s'"
                                           "in tree" % name)

        # Parse the file
        visited_pairs = set()
        for line in file:
                if len(line.strip()) == 0:
                        # Skip blank lines
                        continue
                fields = line.strip().split('\t')
                if len(fields) != 3:
                        raise RuntimeError("Expected 3 fields in distance "
                                           "matrix")
                i_name, j_name, dist = fields[0], fields[1], float(fields[2])
                visited_pairs.add((i_name, j_name))
                visited_pairs.add((j_name, i_name))
                i = get_leaf_index(i_name)
                j = get_leaf_index(j_name)
                if matrix[i, j] != 0 and matrix[i, j] != dist:
                        raise RuntimeError("Distance between %s and %s (%s) conflicts with "
                                           "earlier entry, %s" % (i_name, j_name, dist,
                                                                  matrix[i, j]))
                matrix[i, j] = dist
                matrix[j, i] = dist

        # Verify we've got distances for all pairs
        for leaf1, leaf2 in permutations(leaf_names, 2):
                if (leaf1, leaf2) not in visited_pairs:
                        raise RuntimeError("Distance matrix doesn't have entries for"
                                           " pair (%s, %s)" % (leaf1, leaf2))
        return matrix

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
	new_list = copy([copy(sublist) for sublist in anc])
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
			anc = copy([copy(sublist) for sublist in new_list])
	return list_1

def X_matrix(paths,post,root):
	dict = {}
	post = post[:-1]
	if len(root.descendants) == 2:
		arbitrary_node = root.descendants[0]
		for path in paths:
			if arbitrary_node in path:
				path.remove(arbitrary_node)
		post.remove(arbitrary_node)
	X = np.zeros((len(paths),len(post)))
	for i in range(len(paths)):
		for j in range(len(paths[i])):
			for k in range(len(post)):
				if paths[i][j] == post[k]:
					X[i,k] = 1
					break
	return X

def v_matrix(x_matrix, d_matrix):
	X_T = np.transpose(x_matrix)
	X_T_X = X_T.dot(x_matrix)
	XXX = inv(X_T_X).dot(X_T)
	v = XXX.dot(d_matrix)
	return v

def assign_length(postorder,v):
	if len(postorder[len(postorder)-1].descendants) == 2:
		j = 0
		for i in range(len(postorder)-1):
			if postorder[i] == postorder[len(postorder)-1].descendants[0]:
				postorder[i].length = 0
			else:
				postorder[i].length = v[j]
				j+=1
	else:
		for i in range(len(postorder)-1):
			postorder[i].length = v[i]
	return postorder
		

		

def main():
	opts = parse_args()
	with open(opts.inputTree) as f:
		tree = newick.load(f)

	D_MATRIX = np.vstack(np.array([3,9,10,10,11,7]))
	print "D Matrix is: " 
	print D_MATRIX

	for node in tree:
		po = post_order(node)
		ancA = ancestor_list(po)
		leafs = scan_leaves(po)
		win = distance(ancA,leafs)
		x = X_matrix(win, po,node)
		V = v_matrix(x,D_MATRIX)
		print node.ascii_art()
		print "X Matrix is: "
		print x
		assign_length(po, V)
	



	print 'Root is %s' % tree
	print "Estimated Tree is %s" % newick.dumps(tree)
if __name__ == "__main__":
	main()

