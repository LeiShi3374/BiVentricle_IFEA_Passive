## This is used to calculate the fluid volume encompassed by a surface
## Currently a divergence theory is used
## point: the mesh/surface node position data, edge_gnode: the node index of the edge,
## surf_gnode: the node index of the surface, surf_conn: the surface connectivity

import numpy as np

def FluidVolume(point, edge_gnode, surf_gnode, surf_conn):
    center_cord = [np.sum(point[edge_gnode, 0]) / len(edge_gnode),
                      np.sum(point[edge_gnode, 1]) / len(edge_gnode),
                      np.sum(point[edge_gnode, 2]) / len(edge_gnode)]

    pa = center_cord

    NoC = len(surf_conn)

    vlm = np.zeros(NoC)
    # Cauculate the fluid volume
    for i in range(NoC):
        pb = point[surf_gnode[surf_conn[i, 0]], :]
        pc = point[surf_gnode[surf_conn[i, 1]], :]
        pd = point[surf_gnode[surf_conn[i, 2]], :]

        ab = pb - pa
        ac = pc - pa
        ad = pd - pa

        vlm[i] = 1.0 / 6.0 * np.abs(np.inner((np.cross(ab, ac)), ad))
        bcd0a = pa - 1.0 / 3.0 * (pb + pc + pd)
        bcdn = np.cross((pc - pb), (pd - pb))
        vlm[i] = np.sign(np.dot(bcd0a, bcdn)) * vlm[i]

    sum_vlm = np.sum(vlm)

    return sum_vlm