"""
2024.5.20

"""
from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.log import ProgressBar
from scipy.spatial import KDTree
from scipy.linalg import eigh
import warnings
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

warnings.filterwarnings('ignore')


# suppress some MDAnalysis warnings when writing PDB/GRO files


class SZ(AnalysisBase):
    def __init__(self, universe, residuesGroup: dict, k: int = None, path=None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.k = k
        self.filePath = path
        self.residues = [i for i in residuesGroup]
        self.tailSp = {}
        self.headSp = {}

        self.headAtoms = self.u.atoms[[]]
        self.tailAtoms = self.u.atoms[[]]

        for sp in residuesGroup:

            self.headSp[sp] = residuesGroup[sp][0][0]

            if residuesGroup[sp][-1][0] == 'Chain A':
                self.tailSp[sp] = '??A'
            elif residuesGroup[sp][-1][0] == 'Chain B':
                self.tailSp[sp] = '??B'
            elif residuesGroup[sp][-1][0] == 'Chain A + Chain B':
                self.tailSp[sp] = '??A ??B'

        # 将选择的头部原子和尾部原子按照残基名称进行排序
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]))

            self.tailAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.tailSp[sp]))

        # 计算每个残基尾部原子的平均位置
        self._atomTailMean = np.array([np.mean(i.positions, axis=0) for i in self.tailAtoms.split('residue')])
        # 得到每个残基头部原子和尾部平均位置的向量，用于后续进行同一区域的区分
        self._atomTailToHead = self._atomTailMean - self.headAtoms.positions
        # 计算分析的原子中所包含的残疾数目
        self._n_residues = self.headAtoms.n_residues

        self._headMask = {
            sp :self.headAtoms.resnames == sp for sp in self.residues
        }

        self._tailMask = {
            sp: self.tailAtoms.resnames == sp for sp in self.residues
        }

        self._allAtomsMask = {
            sp: self.u.atoms.resnames == sp for sp in self.residues
        }
        # 获得每种类型的残基数量
        self.numSp = {
            sp: i.atoms.n_residues for sp,i in enumerate(self.headAtoms.groupby('resnames').items())
        }

        self.results.SZ = None

    def _prepare(self):
        self.results.SZ = np.full([self._n_residues, self.n_frames],
                                  fill_value=np.NaN)

    def _single_frame(self):
        normals = self.find_k_nearest_neighbors_and_normals()

        def compute(tail_atm, residue):
            tail_positions = tail_atm.positions.reshape([self.numSp[residue], -1, 3])
            # 获得每个残基每条尾链的最后一个原子的位置
            tail_tail = tail_positions[:, -1, :]
            head_to_tail = self.headAtoms[self._headMask[sp]].positions - tail_tail
            consistent_normals = self.ensure_consistent_normals(sp_normal, head_to_tail)
            print('consistent',consistent_normals.shape)

            # /////////////////////////////////////////////////////////////////
            # 待优化，不需要每帧都计算链的数目
            chain_num = tail_positions.shape[1]

            vectors = tail_positions[:, :chain_num - 2, :] - tail_positions[:, 2:, :]
            vectors_norm = np.linalg.norm(vectors, axis=2, keepdims=True)  # 保留维度以进行广播
            arr = np.sum(vectors * consistent_normals[:, np.newaxis, :], axis=2)[:,:,np.newaxis] / vectors_norm
            # 计算theta和sz
            theta = np.mean(arr, axis=1)
            sz = (3 * theta ** 2 - 1) / 2
            return sz.reshape(-1)

        for sp in self.residues:
            # 获得每个残基他的法向量
            sp_normal = normals[self._headMask[sp]]
            # 得到需要进行分析的原子
            sp_tail = self.tailAtoms[self._tailMask[sp]]
            # 判断该种类型的残基是否选取了两条链进行分析
            if self.tailSp[sp] == '??A ??B':
                chainList = ['??A', '??B']
                spALL = np.full([self.numSp[sp], 2], fill_value=np.NaN)
                for i, chain in enumerate(chainList):
                    sp_chain = sp_tail.select_atoms('name %s' % chain)
                    spALL[:, i] = compute(sp_chain, sp)
                sz = np.mean(spALL, axis=1)
            else:
                sz = compute(sp_tail, sp)
            self.results.SZ[self._headMask[sp], self._frame_index] = sz

    @property
    def SZ(self):
        return self.results.SZ

    def fit_plane_and_get_normal(self, points):
        mean = np.mean(points, axis=0)
        centered_points = points - mean
        covariance_matrix = np.cov(centered_points.T)
        eigenvalues, eigenvectors = eigh(covariance_matrix)
        min_eigenvalue_index = np.argmin(eigenvalues)
        normal = eigenvectors[:, min_eigenvalue_index]
        normal = normal / np.linalg.norm(normal)
        return normal

    def find_k_nearest_neighbors_and_normals(self):
        kdtree = KDTree(self.headAtoms.positions)
        normals = np.zeros((self._n_residues, 3))
        # num = 0
        for i, point in enumerate(self.headAtoms.positions):
            dists, idxs = kdtree.query(point, self.k + 1)
            # vector_arr = self.vector_h_to_t[idxs[1:]]
            # if np.any(np.dot(vector_arr, self.vector_h_to_t[idxs[0]]) < 0):
            #     num += 1
            #     self.find_n_neighbors(kdtree, point, k, dists[-1], self.vector_h_to_t, self.vector_h_to_t[idxs[0]])
            nearest_neighbors = self.headAtoms.positions[idxs[1:]]
            nearest_neighbors = np.vstack((nearest_neighbors, point))
            normals[i] = self.fit_plane_and_get_normal(nearest_neighbors)
        return normals

    @staticmethod
    def find_n_neighbors(kdtree, query_point, n_neighbors, dist, data, vector, max_distance=None, max_iterations=10):
        """
        :param kdtree:
        :param query_point:
        :param n_neighbors:
        :param dist:
        :param data:  points 的头尾向量
        :param vector:point 的头尾向量
        :param max_distance:
        :param max_iterations:
        :return:
        """

        def condition_dis(vector_arr, vector_point):
            dot = np.dot(vector_arr, vector_point)
            return dot > 0

        iteration_count = 3
        for i in range(iteration_count):
            # 搜索邻居点
            distances, indices = kdtree.query(query_point, k=n_neighbors + 1, distance_upper_bound=dist)
            # 检查返回的索引是否有效
            if len(indices) == 0 or np.any(indices >= data.shape[0]):
                dist *= 1.1
                continue
            neighbors = data[indices]
            # 过滤邻居点
            filtered_neighbors = neighbors[condition_dis(data[indices], vector)]
            # 如果符合条件的邻居点数量足够，返回这些点
            if filtered_neighbors.shape[0] >= n_neighbors:
                return filtered_neighbors[:n_neighbors]
                # 如果循环退出且未找到足够的邻居点，则返回找到的所有符合条件的邻居点（如果有的话）
        return filtered_neighbors

    def ensure_consistent_normals(self, normals, vector):
        multi = normals * vector
        multi_naga = multi[:, ] < 0
        normals[multi_naga] = -normals[multi_naga]
        return normals

    def _conclude(self):
        # self.results.SZ = np.mean(self.results.SZ,axis=1)
        pass

    def run(self, start=None, stop=None, step=None, frames=None,
            verbose=None, *, progressbar_kwargs={},callBack=None):

        verbose = getattr(self, '_verbose',
                          False) if verbose is None else verbose

        self._setup_frames(self._trajectory, start=start, stop=stop,
                           step=step, frames=frames)
        self._prepare()
        for i, ts in enumerate(ProgressBar(
                self._sliced_trajectory,
                verbose=verbose,
                **progressbar_kwargs)):
            self._frame_index = i
            self._ts = ts
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            self._single_frame()
            if callBack:
                callBack((i * self.step/(self.stop - self.start) + 0.01) * 100)
        self._conclude()
        return self

    def vmd(self, array=None):
        sz_array = np.zeros([self.u.atoms.n_atoms, ])
        for sp in self._residues:
            n_sp = self._mask_head[sp].sum()
            n_sp_atoms = sz_array[self._mask_all[sp]].shape[0] / n_sp
            if array is not None:
                sz_array[self._mask_all[sp]] = np.repeat(array[self._mask_head[sp]], n_sp_atoms)
            else:
                sz_array[self._mask_all[sp]] = np.repeat(self.results.SZ[self._mask_head[sp]], n_sp_atoms)

        lines = ['proc sz {id} {',
                 'set sel [atomselect $id all]',
                 'set natom [$sel num]',
                 'set rduser [open "output.txt" r]',
                 'set userlist {}',
                 'for {set iatm 0} {$iatm<=[expr $natom-1]} {incr iatm} {',
                 'gets $rduser line',
                 'scan $line "%f" user',
                 'lappend userlist $user',
                 '}',
                 '$sel set user $userlist',
                 'close $rduser',
                 'puts "ALL DONE!"',
                 '}']
        if self.file_path:
            np.savetxt(self.file_path, sz_array, fmt='%.5f', delimiter='\n')
            with open('D:/VMD32/sz.tcl', 'w', encoding='utf-8') as file:
                for line in lines:
                    file.write(line + '\n')
        else:
            np.savetxt('D:/VMD32/output.txt', sz_array, fmt='%.5f', delimiter='\n')
            with open('D:/VMD32/sz.tcl', 'w', encoding='utf-8') as file:
                for line in lines:
                    file.write(line + '\n')

    def make_figure_3d(self, vmin=None, vmax=None, s=30, xlabel='X', ylabel='Y', zlabel='Z'):
        points_pos = self._sel_head.positions
        values = self.results.SZ
        if vmin is None:
            vmin = min(values)
            vmax = - vmin
        else:
            vmin = vmin
            vmax = vmax
        norm = Normalize(vmin=vmin, vmax=vmax)
        colors = plt.cm.bwr(norm(values))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        sc = ax.scatter(points_pos[:, 0], points_pos[:, 1], points_pos[:, 2], c=colors, s=s)
        sm = plt.cm.ScalarMappable(cmap='bwr', norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        plt.show()

    def make_figure_2d(self, method, x_axix='Resnames', y_axis='SZ'):

        if method == 'bar':
            for sp in self.residues:
                plt.bar(sp, np.mean(np.mean(self.results.SZ[self._headMask[sp]])), width=0.4, label='%s' % sp)
        elif method == 'plot':
            for sp in self.residues:
                plt.plot(np.mean(self.results.SZ[self._headMask[sp]], axis=0), label='%s' % sp)
        plt.legend()
        plt.xlabel(x_axix)
        plt.ylabel(y_axis)
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro")
    cls2 = SZ(u, {'DPPC': (['PO4'], ['Chain A + Chain B']), 'D3PC': (['PO4'], ['Chain A'])},k=11)
    cls2.run()
    cls2.make_figure_2d(method='bar')
