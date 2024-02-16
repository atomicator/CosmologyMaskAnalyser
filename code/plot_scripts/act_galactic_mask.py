import numpy as np
from toolkit import toolkit
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--save_path", default="act_galactic_mask_array.npy")
parser.add_argument("--raise_path", type=int, default=2)
args = parser.parse_args()


act_mask = toolkit.load_mask("act", args.raise_path)

act_galactic_mask_map = np.int8(np.ones(np.shape(act_mask.map)))

initial_point = (0, 0)

if act_mask.map[initial_point] == 1:
    raise ValueError("The initial point in unmasked")

point_queue = [initial_point]
act_galactic_mask_map[initial_point[0], initial_point[1]] = 0

while len(point_queue) != 0:
    point = point_queue.pop(0)
    points_to_check = [[point[0] - 1, point[1]], [point[0] + 1, point[1]],
                       [point[0], point[1] - 1], [point[0], point[1] + 1]]
    for point_to_check in points_to_check:
        if act_mask.map[point_to_check[0], point_to_check[1]] == 0 and act_galactic_mask_map[point_to_check[0], point_to_check[1]] == 1:
            act_galactic_mask_map[point_to_check[0], point_to_check[1]] = 0
            #print(point_to_check)
            point_queue.append([point_to_check[0], point_to_check[1]])

print(act_mask.map)
print(act_galactic_mask_map)

np.save(args.save_path, act_galactic_mask_map)
