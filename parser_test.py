import argparse

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('dirs',metavar='N',type=str,nargs=2)
	#parser.add_argument('num',metavar='N',type=int,nargs='+')
	parser.add_argument('--test',dest='test',action='store_true')
	parser.set_defaults(test=False)

	args = parser.parse_args()
	return args

if __name__ == "__main__":
	args = parse_args()
	print args
	print args.dirs


