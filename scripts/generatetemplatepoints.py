
import argparse
import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-N', '--observations', type = int, default = 100, help = 'No. observations')
    parser.add_argument('--xmin', type = float, default = -1.0, help = 'Min x')
    parser.add_argument('--ymin', type = float, default = -1.0, help = 'Min y')
    parser.add_argument('--xmax', type = float, default = 1.0, help = 'Max x')
    parser.add_argument('--ymax', type = float, default = 1.0, help = 'Max y')

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('--seed', type = int, default = 983, help = 'Random seed')

    args = parser.parse_args()

    numpy.random.seed(args.seed)
    
    f = open(args.output, 'w')

    f.write('%d\n' % args.observations)

    for i in range(args.observations):

        x = numpy.random.uniform() * (args.xmax - args.xmin) + args.xmin
        y = numpy.random.uniform() * (args.ymax - args.ymin) + args.ymin

        f.write('%15.9f %15.9f 0.0 1.0\n' % (x, y))

    f.close()
    

    
    
