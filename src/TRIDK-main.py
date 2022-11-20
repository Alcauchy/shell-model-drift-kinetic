import shellmodel_v5
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='TRIDK package v.1.0')
	parser.add_argument('-f','--filename', type=str, default = None,
                    help='parameter filename in JSON format')
	args = parser.parse_args()
	filename = args.filename
	shellmodel_v5.main(filename)

