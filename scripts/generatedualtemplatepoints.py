
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-A', '--observations-A', type = int, default = 50, help = 'No. A observations')
    parser.add_argument('-B', '--observations-B', type = int, default = 50, help = 'No. B observations')
    parser.add_argument('--xmin', type = float, default = -1.0, help = 'Min x')
    parser.add_argument('--ymin', type = float, default = -1.0, help = 'Min y')
    parser.add_argument('--xmax', type = float, default = 1.0, help = 'Max x')
    parser.add_argument('--ymax', type = float, default = 1.0, help = 'Max y')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    args = parser.parse_args()

    f = open(args.output, 'w')

    N = args.observations_A + args.observations_B
    
    f.write('%d\n' % N)

    for i in range(args.observations_A):

        x = numpy.random.uniform() * (args.xmax - args.xmin) + args.xmin
        y = numpy.random.uniform() * (args.ymax - args.ymin) + args.ymin

        f.write('%15.9f %15.9f 0 0.0 1.0\n' % (x, y))

    for i in range(args.observations_B):

        x = numpy.random.uniform() * (args.xmax - args.xmin) + args.xmin
        y = numpy.random.uniform() * (args.ymax - args.ymin) + args.ymin

        f.write('%15.9f %15.9f 1 0.0 1.0\n' % (x, y))

    f.close()
    

    
    
