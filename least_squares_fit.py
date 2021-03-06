import os
import time
from argparse import ArgumentParser
from collections import deque 
from copy import copy, deepcopy
from itertools import permutations
from tempfile import mkstemp
import numpy as np
from numpy.linalg import inv
from StringIO import StringIO
from subprocess import check_output
from cogent import LoadTree
from cogent.maths.stats.test import correlation
# ^ equivalent to "import argparse" and using "argparse.ArgumentParser"
'''
testdata
0.146424835371
1.1893932483
[[0.       0.500501 0.507581 0.386589 0.360539]
 [0.500501 0.       0.176098 0.599054 0.573004]
 [0.507581 0.176098 0.       0.606134 0.580084]
 [0.386589 0.599054 0.606134 0.       0.35211 ]
 [0.360539 0.573004 0.580084 0.35211  0.      ]]
   
[[0.       1.       1.       0.263022 0.295981]
 [1.       0.       0.159146 1.       1.      ]
 [1.       0.159146 0.       1.       1.      ]
 [0.263022 1.       1.       0.       0.295981]
 [0.295981 1.       1.       0.295981 0.      ]]

python2 least_squares_fit.py /mnt/fasta/apes/tree.nh --files /mnt/fasta/apes/hg38.fa /mnt/fasta/apes/panTro5.fa /mnt/fasta/apes/susie.fa /mnt/fasta/apes/ponAbe2.fa /mnt/fasta/apes/nomLeu3
.fa --labels hg38 panTro5 susie ponAbe2 nomLeu3
0.00678913242032
0.00031755682992
[[0.       0.01339  0.019734 0.039403 0.046204]
 [0.01339  0.       0.020024 0.039693 0.046494]
 [0.019734 0.020024 0.       0.037597 0.044398]
 [0.039403 0.039693 0.037597 0.       0.044681]
 [0.046204 0.046494 0.044398 0.044681 0.      ]]

[[0.        0.0135016 0.0168391 0.032829  0.037311 ]
 [0.0135016 0.        0.017653  0.0326168 0.0385787]
 [0.0168391 0.017653  0.        0.032829  0.0392331]
 [0.032829  0.0326168 0.032829  0.        0.0403101]
 [0.037311  0.0385787 0.0392331 0.0403101 0.       ]]

Species swap
0.418525842762
0.00185018178192
[[0.       0.01339  0.019734 0.039403 0.046204]
 [0.01339  0.       0.020024 0.039693 0.046494]
 [0.019734 0.020024 0.       0.037597 0.044398]
 [0.039403 0.039693 0.037597 0.       0.044681]
 [0.046204 0.046494 0.044398 0.044681 0.      ]]
[[0.        0.032829  0.0168391 0.0135016 0.037311 ]
 [0.032829  0.        0.032829  0.0326168 0.0403101]
 [0.0168391 0.032829  0.        0.017653  0.0392331]
 [0.0135016 0.0326168 0.017653  0.        0.0385787]
 [0.037311  0.0403101 0.0392331 0.0385787 0.       ]]

Branch length factor * 10
0.00553857372597
1.0139071292
[[0.      0.1339  0.15936 0.35605 0.42406]
 [0.1339  0.      0.16226 0.35895 0.42696]
 [0.15936 0.16226 0.      0.37597 0.44398]
 [0.35605 0.35895 0.37597 0.      0.44681]
 [0.42406 0.42696 0.44398 0.44681 0.     ]]
[[0.        0.0135016 0.0168391 0.032829  0.037311 ]
 [0.0135016 0.        0.017653  0.0326168 0.0385787]
 [0.0168391 0.017653  0.        0.032829  0.0392331]
 [0.032829  0.0326168 0.032829  0.        0.0403101]
 [0.037311  0.0385787 0.0392331 0.0403101 0.       ]]
        
python2 least_squares_fit.py /mnt/fasta/nematodes/tree.nh --files /mnt/fasta/nematodes/C_elegans.fa /mnt/fasta/nematodes/C_japonica.fa /mnt/fasta/nematodes/P_pacificus.fa --labels C_elegans C_japonica P_pacificus
0.424419128204
0.364650584073
[[0.       0.262598 0.399917]
 [0.262598 0.       0.405149]
 [0.399917 0.405149 0.      ]]
[[0.       0.263022 0.295981]
 [0.263022 0.       1.      ]
 [0.295981 1.       0.      ]]
'''

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
        leaf_names = [n.Name for n in post_order]

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
                        print line
                        fields = line.split()
                        names.append(fields[0])

                # 2nd pass
                f.seek(0)
                f.readline()
                for line in f:
                        fields = line.strip().split()
                        name1 = fields[0]
                        distances = fields[1:]
                        print distances, names
                        assert len(distances) == len(names), "Incorrect number of entries in distance matrix"
                        for i, distance in enumerate(distances):
                                name2 = names[i]
                                matrix += '%s\t%s\t%s\n' % (name1, name2, distance)
        matrix_file = StringIO(matrix)
        return read_distance_matrix(matrix_file, post_order), runtime

def run_mash_and_get_matrix(input_files, post_order, file_to_label):
	"""
	Run Mash and get out a numpy matrix.
	
	input_files: array of filenames
	post_order: ordering of leaves within the matrix
	file_to_label: correspondence between filename and leaf name
	"""
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

def distance_from_r_squared(m1, m2):
    """Estimates distance as 1-r^2: no correl = max distance"""
    return 1 - (correlation(m1.flat, m2.flat)[0])**2

def main():
    opts = parse_args()
    file_to_label = {}
    for file, label in zip(opts.files, opts.labels):
        file_to_label[file] = label
    tr = LoadTree(opts.inputTree)

    disMatrix, tip_order = tr.tipToTipDistances()
    mashMatrix, _ = run_mash_and_get_matrix(opts.files,tip_order,file_to_label)
    print distance_from_r_squared(disMatrix,mashMatrix)
  
    sqaured_difference = (disMatrix - mashMatrix) ** 2
    sum_of_squares = np.tril(sqaured_difference).sum()
    print sum_of_squares

    print disMatrix
    print mashMatrix

if __name__ == "__main__":
    main()
