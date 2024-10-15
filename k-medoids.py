import numpy as np
from itertools import combinations
from Levenshtein import distance as levenshtein_distance
import random

# 提取迹变体
trace_variants = [
    "A,B,C,D,F,E,G,H",
    "A,B,C,D,E,F,G",
    "A,B,C,H",
    "A,B,C,D,E,G,F,H",
    "A,B,C,D,H",
    "A,H",
    "A,D,E,G,H",
    "A,B,F,E,G,H",
    "A,D,F,H",
    "A,F,B,C",
    "A,C,E,F,G",
    "B,F,G"
]

# 计算Levenshtein距离矩阵
n = len(trace_variants)
distance_matrix = np.zeros((n, n))

for i, j in combinations(range(n), 2):
    dist = levenshtein_distance(trace_variants[i], trace_variants[j])
    distance_matrix[i, j] = dist
    distance_matrix[j, i] = dist


# 手动实现k-medoids算法
def k_medoids(distance_matrix, k, max_iter=100):
    m, n = distance_matrix.shape

    # Step 1: 随机选择k个数据点作为初始medoids
    current_medoids = np.random.choice(n, k, replace=False)
    clusters = None

    for iteration in range(max_iter):
        # Step 2: 分配数据点到最近的medoid
        clusters = {medoid: [] for medoid in current_medoids}

        for data_point in range(n):
            distances_to_medoids = [distance_matrix[data_point, medoid] for medoid in current_medoids]
            closest_medoid = current_medoids[np.argmin(distances_to_medoids)]
            clusters[closest_medoid].append(data_point)

        # Step 3: 更新medoids
        new_medoids = []
        for medoid in clusters.keys():
            cluster_points = clusters[medoid]
            medoid_distances = np.sum(distance_matrix[np.ix_(cluster_points, cluster_points)], axis=1)
            new_medoid = cluster_points[np.argmin(medoid_distances)]
            new_medoids.append(new_medoid)

        new_medoids = np.array(new_medoids)

        # 检查medoids是否变化
        if np.array_equal(new_medoids, current_medoids):
            print(f"收敛在第{iteration + 1}次迭代")
            break

        current_medoids = new_medoids

    return current_medoids, clusters


# 设置簇数为3
k = 3
medoids, clusters = k_medoids(distance_matrix, k)

# 输出聚类结果
for medoid, cluster in clusters.items():
    print(f"Medoid: {trace_variants[medoid]}, Cluster: {[trace_variants[idx] for idx in cluster]}")
