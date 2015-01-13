import sys

call = []
begall = {}
endall = {}

with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.strip()
        if line[:6] == 'CONECT':
            print line
        elif line[:4] == 'ATOM':
            beg = line[:29]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            end = line[55:]
            if (x, y, z) not in call:
                call.append((x, y, z))
                print line
            else:
                while (x, y, z) in call:
                    x += 0.2
                print beg, "{0[0]:>8.3f}{0[1]:>8.3f}{0[2]:>8.3f}".format((x, y, z)), end
                call.append((x, y, z))
