import os
import time
import newick
from argparse import ArgumentParser
from collections import deque 
from copy import copy, deepcopy
from itertools import permutations
from tempfile import mkstemp
import numpy as np
from numpy.linalg import inv
from StringIO import StringIO
from subprocess import check_output
# ^ equivalent to "import argparse" and using "argparse.ArgumentParser"




def parse_args():
	parser = ArgumentParser(description=__doc__)
	parser.add_argument('inputTree', help='newick format tree (in a file)')
	parser.add_argument('--files', nargs='+')
	parser.add_argument('--labels', nargs='+')
        parser.add_argument('--method', choices=['mash', 'kmacs', 'spaced'],
                            default='mash')
        parser.add_argument('--tsv', action='store_true', help='Output metrics as TSV')
        parser.add_argument('--noHeader', action='store_true', help='Suppress TSV header')
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

def D_matrix(dist_matrix, post_order, paths):
	D = np.zeros((len(paths),1))
	leaf_order = [node for node in post_order if node.is_leaf]
	k=0
	for i in range(len(leaf_order)-1):
		for j in range(1+i,len(leaf_order)):
			D[k,0] = dist_matrix[i,j]
			k += 1
	return D

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
				postorder[i].length = v[j][0]
				j+=1
	else:
		for i in range(len(postorder)-1):
			postorder[i].length = v[i][0]
	return postorder



		
def performance_metric(postorder, true_tree):
	#postorder comes from assign_length so that if binary tree, the altered length is included
	metric = 0
	TruTre = post_order(true_tree)
	if len(postorder[len(postorder)-1].descendants) == 2:
		for i in range(len(postorder)-1):
			if TruTre[i] == TruTre[len(postorder)-1].descendants[0]:
				TruTre[len(TruTre)-1].descendants[1].length += TruTre[i].length
			else:
				metric += (postorder[i].length - TruTre[i].length)**2
	else:
		for i in range(len(postorder)-1):
			metric += (postorder[i].length - TruTre[i].length)**2
	return metric





def read_distance_matrix(file, post_order):
	"""
        Read a distance matrix from a file.

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
                        raise RuntimeError("Didn't find a leaf with name '%s' "
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

def run_kmacs_and_get_matrix(input_files, post_order, file_to_label, k=0):
        fasta_path = produce_concatenated_fasta(input_files, [file_to_label[f] for f in input_files])
        start_time = time.time()
        check_output(['kmacs', fasta_path, str(k)])
        runtime = time.time() - start_time
        matrix = ''
        with open('DMat') as f:
                # We make two passes. First, we just get the ordering
                # of taxa names in the distance matrix. Next, we'll go
                # through the file again, constructing the distance
                # matrix in a format we can parse.

                # 1st pass
                names = []
                # Skip first line
                f.readline()
                for line in f:
                        fields = line.split()
                        names.append(fields[0])

                # 2nd pass
                f.seek(0)
                f.readline()
                for line in f:
                        fields = line.strip().split()
                        name1 = fields[0]
                        distances = fields[1:]
                        assert len(distances) == len(names), "Incorrect number of entries in distance matrix"
                        for i, distance in enumerate(distances):
                                name2 = names[i]
                                matrix += '%s\t%s\t%s\n' % (name1, name2, distance)
        matrix_file = StringIO(matrix)
        return read_distance_matrix(matrix_file, post_order), runtime

def run_mash_and_get_matrix(input_files, post_order, file_to_label):
	ret = []
	new_output = ''

	start_time = time.time()
	check_output(['mash', 'sketch'] + input_files + ['-o', 'sketch'])
	output = check_output(['mash', 'dist', 'sketch.msh', 'sketch.msh'])
	runtime = time.time() - start_time
	for i in output.split('\n'):
		if len(i) == 0:
			# skip blank lines
			continue
		fields = i.split('\t')
		del fields[3:]
		fields[0] = file_to_label[fields[0]]
		fields[1] = file_to_label[fields[1]]
		ret.append(fields)
	for i in range(len(ret)-1):
		str1 = '\t'.join(map(str,ret[i]))
		new_output += str1
		new_output += '\n'
	mash_matrix = read_distance_matrix(StringIO(new_output),post_order)
	return mash_matrix, runtime





def run_spaced_and_get_matrix(input_files, post_order, file_to_label):
	fasta_path = produce_concatenated_fasta(input_files, [file_to_label[f] for f in input_files])
        start_time = time.time()
        check_output(['spaced', fasta_path])
        runtime = time.time() - start_time
        matrix = ''
        with open('DMat') as f:
                # We make two passes. First, we just get the ordering
                # of taxa names in the distance matrix. Next, we'll go
                # through the file again, constructing the distance
                # matrix in a format we can parse.

                # 1st pass
                names = []
                # Skip first line
                f.readline()
                for line in f:
                        fields = line.split()
                        names.append(fields[0])

                # 2nd pass
                f.seek(0)
                f.readline()
                for line in f:
                        fields = line.strip().split()
                        name1 = fields[0]
                        distances = fields[1:]
                        assert len(distances) == len(names), "Incorrect number of entries in distance matrix"
                        for i, distance in enumerate(distances):
                                name2 = names[i]
                                matrix += '%s\t%s\t%s\n' % (name1, name2, distance)
        matrix_file = StringIO(matrix)
        return read_distance_matrix(matrix_file, post_order), runtime




def produce_concatenated_fasta(input_paths, names):
        """Given a list of fasta files (representing one genome each), get a
        fasta file containing a "sequence" for each genome, which has all
        sequences in the genome concatenated together, separated by Ns."""
        # Get temporary file
        fd, path = mkstemp()
        fh = os.fdopen(fd, 'w')
        for input_path, name in zip(input_paths, names):
                seq = []
                with open(input_path) as f:
                        for line in f:
                                if line == '':
                                        # Ignore blank lines
                                        continue
                                if line[0] == '>':
                                        # fasta header. Append some Ns
                                        # to try to break up any
                                        # k-mers that may otherwise
                                        # span sequences
                                        seq.append('N' * 100)
                                else:
                                        seq.append(line.strip())
                fh.write('>' + name + '\n')
                for line in seq:
                        fh.write(line + '\n')

        return path


def main():
    opts = parse_args()
    file_to_label = {}
    for file, label in zip(opts.files, opts.labels):
        file_to_label[file] = label
    with open(opts.inputTree) as f:
        tree = newick.load(f)

    """test_D_MATRIX = np.vstack(np.array([3,9,10,10,11,7]))

    distance_mat = np.matrix([[0,3,9,10,6],
                              [3,0,10,11,7],
                              [9,10,0,7,3],
                              [10,11,7,0,4],
                              [6,7,3,4,0]])"""


    for node in tree:
        true_tree = deepcopy(node)
        po = post_order(node)
        ancA = ancestor_list(po)
        leafs = scan_leaves(po)
        win = distance(ancA,leafs)
        x = X_matrix(win, po,node)
        if opts.method == 'mash':
            matrix, runtime = run_mash_and_get_matrix(opts.files, po, file_to_label)
        elif opts.method == 'kmacs':
            matrix, runtime = run_kmacs_and_get_matrix(opts.files, po, file_to_label)
        elif opts.method == 'spaced':
            matrix, runtime = run_spaced_and_get_matrix(opts.files, po, file_to_label)
        #D_MATRIX = D_matrix(distance_mat,po,win)
        D_MATRIX = D_matrix(matrix,po,win)
        V = v_matrix(x,D_MATRIX)
        l = assign_length(po, V)
        perf=performance_metric(l, true_tree)
        if opts.tsv:
            # Print TSV-style information
            if not opts.noHeader:
                print 'TestSet\tMethod\tRuntime\tSumOfSquaredDifferences'
            print '%s\t%s\t%s\t%s' % (opts.inputTree, opts.method, runtime, perf)
        else:
            # Print other (debugging) information
            print matrix
            print "The true tree is: "
            print "the metric is "
            print perf
            print node.ascii_art()
            print "X Matrix is: "
            print x
            print "Distance Matrix is: "
            print matrix
            print "D Matrix is: "
            print D_MATRIX
            print "v matrix is: "
            print V
            print 'Root is %s' % tree
            print "Estimated Tree is %s" % newick.dumps(tree)
if __name__ == "__main__":
    main()

