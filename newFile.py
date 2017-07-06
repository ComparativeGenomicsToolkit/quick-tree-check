import newick
from argparse import ArgumentParser
# ^ equivalent to "import argparse" and using "argparse.ArgumentParser"

def parse_args():
	parser = ArgumentParser(description=__doc__)
	parser.add_argument('inputTree', help='newick format tree (in a file)')
	return parser.parse_args()

def main():
	opts = parse_args()
	with open(opts.inputTree) as f:
		trees = newick.load(f)

	for tree in trees:
		print tree.ascii_art()
	print 'trees is %s' % trees
	print newick.dumps(trees)
if __name__ == "__main__":
	main()

