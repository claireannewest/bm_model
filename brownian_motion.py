import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.spatial import distance
import random
from itertools import chain


class Aggregation_Model:
    def __init__(self, dim, num_part, num_steps, d, dt):
        num_part = 1000
        num_steps = 202
        d = 50
        dt = .01
        scale = d**2 * dt



# fig, ax = plt.subplots(figsize=[7, 7])
# fig.dpi = 144
fig = plt.figure(figsize=[7, 7])
ax = plt.axes(projection='3d')


def dynamic_plot(parts_x, parts_y, parts_z, colors):
    ax.clear()

    # fig, ax = plt.subplots(figsize=[7, 7])
    # fig.dpi = 144

    box_limit = 1E4  # nm
    ax.set_xlim([-box_limit, box_limit])
    ax.set_ylim([-box_limit, box_limit])
    ax.set_zlim([-box_limit, box_limit])

    # ax.set_aspect(1)
    ax.set_box_aspect((num_part, num_part, num_part))

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width = bbox.width
    pixels_per_axis = width * fig.dpi
    points_per_axis = pixels_per_axis * 72 / fig.dpi

    coord_to_pixel = points_per_axis / (box_limit * 2)  # pixels per nm
    markersize = np_eff_rad * coord_to_pixel
    s = markersize**2
    ax.set_title(f"ti = {ti}")

    ax.scatter(parts_x, parts_y, parts_z,
               c=colors,
               s=s,
               cmap='plasma',
               )
    plt.pause(.05)
    return None


def random_walk(ti, ni, delta, final_colls=None):
    if delta is True:
        dx = norm.rvs(scale=scale)
        dy = norm.rvs(scale=scale)
        dz = norm.rvs(scale=scale)

    if delta is False:
        dx = (part_loc[ti + 1, final_colls[0], 0]
              - part_loc[ti, final_colls[0], 0])
        dy = (part_loc[ti + 1, final_colls[0], 1]
              - part_loc[ti, final_colls[0], 1])
        dz = (part_loc[ti + 1, final_colls[0], 2]
              - part_loc[ti, final_colls[0], 2])

    part_loc[ti + 1, ni, 0] = part_loc[ti, ni, 0] + dx
    part_loc[ti + 1, ni, 1] = part_loc[ti, ni, 1] + dy
    part_loc[ti + 1, ni, 2] = part_loc[ti, ni, 2] + dz

    return None




part_loc = np.zeros((num_steps, num_part, 3))

for ni in range(num_part):
    part_loc[0, ni, :] = np.array([random.randint(-1E4, 1E4),
                                   random.randint(-1E4, 1E4),
                                   random.randint(-1E4, 1E4),
                                   ])

np_eff_rad = 200
thresh = np_eff_rad
colors = np.arange(0, num_part)


for ti in range(num_steps-1):
    # Plot the particles at time ti

    # Find distance between each particle, and make a collision matrix.
    all_dists = distance.cdist(part_loc[ti, :, :],
                               part_loc[ti, :, :],
                               'euclidean')
    idxc = np.where((all_dists <= thresh) & (all_dists > 0))
    all_dists[idxc] = 1
    coll_mat = np.zeros((num_part, num_part))
    coll_mat[idxc] = 1

    # Iterate through all particles, find their collisions
    for ni in range(num_part):
        num_col = int(np.sum(coll_mat[:, ni]))
        if num_col == 0:
            random_walk(ti=ti, ni=ni, delta=True)
        if num_col > 0:
            all_colls = []
            idx_col = np.where(coll_mat[:, ni] == 1)[0]
            all_colls.append(list(idx_col))
            layer_i = set(list(chain.from_iterable(all_colls)))
            flag = True
            while flag is True:
                for coli in layer_i:
                    # Of the particles who have collided,
                    # who did those particles collide with?
                    all_colls.append(list(np.where(coll_mat[:, coli] == 1)[0]))
                layer_ip1 = set(list(chain.from_iterable(all_colls)))
                layer_ip1.remove(ni)
                if layer_ip1 == layer_i:
                    # you've found all the collisions
                    flag = False
                else:
                    # keep going, you've added a new collision
                    layer_i = layer_ip1
            final_colls = np.asarray(list(layer_ip1))
            final_colls.sort()
            check_small = np.where(final_colls < ni)[0]

            if len(check_small) > 0:
                # If collided with idx smaller than you, grab their dx, dy
                random_walk(ti=ti, ni=ni, delta=False, final_colls=final_colls)
                colors[ni] = colors[final_colls[0]]
            else:
                random_walk(ti=ti, ni=ni, delta=True)

    dynamic_plot(parts_x=part_loc[ti, :, 0],
                 parts_y=part_loc[ti, :, 1],
                 parts_z=part_loc[ti, :, 2],
                 colors=colors)

    if set(colors) == {0}:
        break
    # if ti == 1: break
plt.show()
