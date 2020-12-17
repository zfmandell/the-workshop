import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('thresh',type=str,help='')
    parser.add_argument('terms',type=str,help='')
    parser.add_argument('upstream',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    if args.terms.split("_")[1].split('.')[0] == 'dA':
        with open(args.upstream,'r') as inp:
            firstline = inp.readline
        with open(args.terms,'')

    elif args.terms.split("_")[1].split('.')[0] == 'dG':

    elif args.terms.split("_")[1].split('.')[0] == 'dR':

    elif args.terms.split("_")[1].split('.')[0] == 'dAG':

    elif args.terms.split("_")[1].split('.')[0] == 'dAR':

    elif args.terms.split("_")[1].split('.')[0] == 'dGR':

    elif args.terms.split("_")[1].split('.')[0] == 'dAGR':

if __name__ == '__main__':
    main()
