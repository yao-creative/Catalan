from multiprocessing import Pool
import numpy as np
from simulation import *
import matplotlib.pyplot as plt
import json

N = 5 # size of object
M = 100 # number of samples
method = 2
dyck_path_objs =[]
dyck_paths = [""] * M
binary_trees = [""] * M
planar_trees =[""] * M
noncrossings = [""] * M
planar_tree_graphs = [""] * M
binary_tree_graphs = [""] * M
distribution = {}
times = {"dyck":0.,"noncrossing":0.,"plane_tree":0., "plane_tree_graph":0., "binary_tree_graph":0.}
#TODO generate M samples with each of 4 categories above in mathematica format
#TODO save to file
#TODO load and plot on mathematica
def generate_mathematica_sample(i):
    A = Catalan_Structure(N)
    A.generate(method)
    A.to_mathematica()
    dyck_path_objs.append(A.dyck_path)
    # print(f"\n\n ############################\n Generated sample {i}: \n")
    # A.show()
    # print(f"\n ############################\n\n")
    
    # A.check()
    # A.show_times()
    for time in times:
        times[time] += A.times[time]
    prefix = f"p{i} = "
    dyck_path = (prefix + A.dyck_path_mathematica)
    if A.dyck_path_mathematica in distribution:
            distribution[A.dyck_path_mathematica] += 1
    else:
        distribution[A.dyck_path_mathematica] = 1
    # binary_tree = (prefix + A.binary_tree_mathematica)
    planar_tree= (prefix + A.plane_tree_mathematica)
    noncrossing = (prefix + A.non_crossing_partition_mathematica)
    planar_tree_graph = (prefix + A.plane_tree_graph_mathematica)
    binary_tree_graph = (prefix + A.binary_tree_graph_mathematica)
    return (dyck_path, planar_tree, noncrossing,planar_tree_graph,binary_tree_graph)


    
if __name__ == "__main__":
    pool = Pool(8)
    # res = pool.map(generate_mathematica_sample, range(M))
    res = [generate_mathematica_sample(i) for i in range(M)]
        
    # pool.close()
    # pool.join()
    
    for i in range(M):
        dyck_paths[i] = res[i][0]
        # binary_trees[i] = res[i][1]
        planar_trees[i] = res[i][1]
        noncrossings[i] = res[i][2]
        planar_tree_graphs[i] = res[i][3]
        binary_tree_graphs[i] = res[i][4]
        
    dyck_paths_str = f"time: {times['dyck']}, average: {times['dyck']/M}\n"+ "\n".join(dyck_paths)
    # binary_trees_str = "\n".join(binary_trees)
    planar_trees_str = f"time: {times['plane_tree']}, average: {times['plane_tree']/M}\n"+ "\n".join(planar_trees)
    noncrossing_str = f"time: {times['noncrossing']}, average: {times['noncrossing']/M}\n"+ "\n".join(noncrossings)
    planar_tree_graphs_str = f"time: {times['plane_tree_graph']}, average: {times['plane_tree_graph']/M}\n"+ "\n".join(planar_tree_graphs)
    binary_tree_graphs_str = f"time: {times['binary_tree_graph']}, average: {times['binary_tree_graph']/M}\n"+ "\n".join(binary_tree_graphs)
    with open(f"dyck_paths_N{N}_M{M}.txt", "w") as f:
        f.write(dyck_paths_str)
    # with open(f"binary_trees_N{N}_M{M}.txt", "w") as f:
    #     f.write(binary_trees_str)
    with open(f"planar_trees_N{N}_M{M}.txt", "w") as f:
        f.write(planar_trees_str)
    with open(f"noncrossing_N{N}_M{M}.txt", "w") as f:
        f.write(noncrossing_str)
    with open(f"planar_tree_graphs_N{N}_M{M}.txt", "w") as f:
        f.write(planar_tree_graphs_str)
    with open(f"binary_tree_graphs_N{N}_M{M}.txt", "w") as f:
        f.write(binary_tree_graphs_str)
        
    print(f"Distribution of Dyck paths: {distribution}")
    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    fig.tight_layout(pad=3.0)
    ax1.set_title(f"Distribution of Dyck paths for N={N} and M={M}")
    ax1.set_ylim(0,max(distribution.values())+ M/(2*(len(distribution))))
    ax1.plot(list(distribution.values()))
    ax2.set_title(f"Dyck paths for N={N} and M={M} Method{method}")
    c = 4*np.log(M)**2
    
    for path in dyck_path_objs:
        ax2.plot(path.ys, color = "green", alpha = 1/c)
    
    fig.savefig(f"distribution_and_dyck_paths_N{N}_M{M}_Method{method}.png")
    
    fig2 = plt.figure(figsize=(20,20))
    ys = np.linspace(0,N+1,N+1)
    xs = np.linspace(0,2*N+1,2*N+1)
    zs = np.zeros((N+1,2*N+1))
    # print(f"dimensions: zs: {zs.shape} xs: {xs.shape} ys: {ys.shape}")
    for path in dyck_path_objs:
        for i in range(len(path.ys[:-1])):
            zs[path.ys[i]][i] += 1

    plt.xlabel("x")
    plt.ylabel("y")

    contour = plt.contourf(xs,ys, zs,extend="max")
    plt.colorbar(contour)
    fig2.savefig(f"contour_N{N}_M{M}_Method{method}.png")

    np.save(f"point_distribution_N{N}_M{M}_Method{method}",zs)
    
    
