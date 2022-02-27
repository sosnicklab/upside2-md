import sys
import numpy as np
import matplotlib.pyplot as plt

def compact_sigmoid(x, sharpness):
    y = x*sharpness;
    result = 0.25 * (y+2) * (y-1)**2
    result = np.where((y< 1), result, np.zeros_like(result))
    result = np.where((y>-1), result, np.ones_like (result))
    return result


def main(pro_state_file):

    fields = [ln.split() for ln in open(pro_state_file)]

    header_fields = 'residue init_state energy k width'.split()
    if [x.lower() for x in fields[0]] == header_fields:
        len_header_fields = len(header_fields)
    else:
        parser.error('First line of restraint table must be "%s"'%(" ".join(header_fields)))
    if not all(len(f)==len_header_fields for f in fields):
        parser.error('Invalid format for restraint file')

    fields = fields[1:]
    n_pro = len(fields)

    energy     = np.zeros(n_pro)
    k          = np.zeros(n_pro)
    width      = np.zeros(n_pro)
    for i,f in enumerate(fields):
        energy[i]     = float(f[2])  # the height from the trans state to the transition state
        k[i]          = float(f[3])  # k*energy is the  height from the cis state to the transition state
        width[i]      = float(f[4])  # compact_sigmoid cuts off at angle +/- width (degree)

    x1 = np.linspace(-np.pi, -0.5*np.pi, 50)
    x2 = np.linspace(-0.5*np.pi, 0., 50)
    x3 = np.linspace(0., 0.5*np.pi, 50)
    x4 = np.linspace(0.5*np.pi, np.pi, 50)

    c1 = -0.75*np.pi
    c2 = -0.25*np.pi
    c3 =  0.25*np.pi
    c4 =  0.75*np.pi

    for i,f in enumerate(fields):
        S = 1./(width[i]/180.*np.pi)

        y1 = (1.-compact_sigmoid(x1 - c1, S)) * energy[i]
        y2 = ((1-k[i]) + k[i]*compact_sigmoid(x2 - c2, S)) * energy[i]
        y3 = (1 - k[i]*compact_sigmoid(x3 - c3, S)) * energy[i]
        y4 = compact_sigmoid(x4 - c4, S) * energy[i]

        fig,ax = plt.subplots()

        ax.plot(x1,y1)
        ax.plot(x2,y2)
        ax.plot(x3,y3)
        ax.plot(x4,y4)
        ax.set_xlabel('angle')
        ax.set_ylabel('energy')
        ax.set_title(str(i))

    plt.show()

if __name__ == '__main__':
    main(sys.argv[1])
