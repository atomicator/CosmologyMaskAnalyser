import numpy as np
import toolkit.toolkit as toolkit
import os
import matplotlib.pyplot as plt

path_old = "../../data/test/"
path_new = "../../data/test2/"
files = os.listdir(path_old)
i = 0

for file_name in files:
    i += 1
    print(f"{i}/{len(files)}")
    mask = toolkit.PixellMask(path_old + file_name, suppress_warnings=True)
    mask_to_edit = mask.map.copy()
    plt.imshow(mask_to_edit)
    plt.title(f"Initial: {i}")
    plt.show()
    checked = np.zeros(mask_to_edit.shape)
    queue = [[0, 0]]

    for point in queue:
        checked[point] = 1

    while queue:
        point = queue.pop(0)
        to_check = []
        if point[0] > 0:
            to_check.append([point[0] - 1, point[1]])
        if point[0] < mask_to_edit.shape[0] - 1:
            to_check.append([point[0] + 1, point[1]])
        if point[1] > 0:
            to_check.append([point[0], point[1] - 1])
        if point[1] < mask_to_edit.shape[1] - 1:
            to_check.append([point[0], point[1] + 1])
        for point in to_check:
            #print(point)
            #print(checked[point[0], point[1]])
            if checked[point[0], point[1]] == 0 and mask_to_edit[point[0], point[1]] == 0:
                checked[point[0], point[1]] = 1
                mask_to_edit[point[0], point[1]] = 1
                queue.append(point)

    #plt.imshow(checked)
    #plt
    #plt.show()
    plt.imshow(mask_to_edit)
    plt.title(f"Final: {i}")
    plt.show()
    np.save(path_new + "point_" + file_name, mask_to_edit)
