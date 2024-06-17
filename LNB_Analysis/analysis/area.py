"""
2024.6.17
未修改完成
1.无法修复点的位置在边界的情况
2.没有修改搜索k近邻问题

"""

import pandas as pd
from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.log import ProgressBar
from scipy.spatial import KDTree, Voronoi
from scipy.linalg import eigh
import warnings
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cProfile
warnings.filterwarnings('ignore')
# suppress some MDAnalysis warnings when writing PDB/GRO files


class Area(AnalysisBase):
    def __init__(self, universe, residueGroup: dict, k: int, file_path=None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.k = k
        self.filePath = file_path
        self.residues = [i for i in residueGroup]
        self.headSp = {
            sp: residueGroup[sp][0] for sp in self.residues
        }

        self.headAtoms = self.u.atoms[[]]
        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]))

        self._n_residues = self.headAtoms.n_residues

        # self._mask_head = {
        #     sp: self._sel_head.resnames == sp for sp in self._residues
        # }
        # self._mask_all = {
        #     sp: self.u.atoms.resnames == sp for sp in self._residues
        # }

        self.results.Area = None

    @property
    def Area(self):
        return self.results.Area

    def _prepare(self):
        self.results.Area = np.full([self._n_residues, self.n_frames],
                                    fill_value=np.NaN)

    def _single_frame(self):
        kdtree = KDTree(self.headAtoms.positions)
        area_arr = np.zeros([self._n_residues])
        for i, point in enumerate(self.headAtoms.positions):
            _, idxs = kdtree.query(point, self.k + 1)
            nearest_neighbors = np.vstack((self.headAtoms.positions[idxs[1:]], point))
            area_arr[i] = self.get_voronoi_area(nearest_neighbors)
        self.results.Area[:, self._frame_index] = area_arr

    # def _single_frame(self):
    #     kdtree = KDTree(self._sel_head.positions)
    #     area_arr = np.zeros([len(self._sel_head), ])
    #     with ThreadPoolExecutor(max_workers=5) as executor:
    #         # 为每个点创建一个 future
    #         futures = [executor.submit(self.calculate_voronoi_area, point, kdtree, self.k)
    #                    for point in self._sel_head.positions]
    #         # 收集结果
    #         for i, future in enumerate(concurrent.futures.as_completed(futures)):
    #             area_arr[i] = future.result()
    #     self.results.Area[:, self._frame_index] = area_arr
    #
    # def calculate_voronoi_area(self, point, kdtree, k):
    #     dists, idxs = kdtree.query(point, k + 1)
    #     nearest_neighbors = self._sel_head.positions[idxs[1:]]
    #     nearest_neighbors = np.vstack((nearest_neighbors, point))
    #     area = self.get_voronoi_area(nearest_neighbors)
    #     return area

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

    def _conclude(self):
        # self.writeExcel()
        pass
    def writeExcel(self):
        columnFrame = [i*self.step for i in range(self.n_frames)]

        columnHead = ['resid'] + ['resname'] + columnFrame
        resids = self._sel_head.resids
        resnames = self._sel_head.resnames
        resPos = np.column_stack((resids, resnames, self.results.Area))
        df = pd.DataFrame(resPos)

        df.to_excel('%s/area2.xlsx' % self.file_path, header=columnHead, index=False)

    @property
    def Area(self):
        return self.results.Area

    def get_voronoi_area(self, points):
        mean = np.mean(points, axis=0)
        # 将每个点减去中心点得到中心化点
        centered_points = points - mean
        # 计算中心化点的协方差矩阵
        # 计算协方差矩阵的特征值和特征向量
        eigenvalues, eigenvectors = eigh(np.cov(centered_points.T))
        # 找到最小特征值对应的特征向量作为平面的法向量
        normal = eigenvectors[:, np.argmin(eigenvalues)]
        # 计算点到平面的距离
        # 将点投影到平面上
        projected_p = points - np.dot(centered_points, normal)[:, np.newaxis] * normal

        # 选择最后一个点作为参考点
        reference_point = projected_p[-1]
        # 计算x轴方向，这里我们选择最后一个点与中心点的向量
        x_axis = (reference_point - mean) / np.linalg.norm(reference_point - mean)
        # 计算y轴方向，叉乘法向量和x轴得到y轴
        y_axis = np.cross(normal, x_axis)
        y_axis = y_axis / np.linalg.norm(y_axis)  # 标准化y轴

        # 将投影点转换到局部坐标系
        local_coordinates = projected_p - mean
        local_coordinates = np.column_stack((np.dot(local_coordinates, x_axis), np.dot(local_coordinates, y_axis)))

        # 创建Voronoi图
        vor = Voronoi(local_coordinates)

        # 找到最后一个点的Voronoi区域
        region_index = vor.point_region[len(points) - 1]

        # 获取Voronoi区域的顶点
        region_vertices = vor.regions[region_index]
        if -1 in region_vertices:
            return 0

        # 过滤掉无限区域的顶点
        finite_vertices_indices = [index for index in region_vertices if index != -1]
        vertices = vor.vertices[finite_vertices_indices]

        # 计算多边形面积
        return self.polygon_area(vertices)

    @staticmethod
    def polygon_area(vertices):
        n = vertices.shape[0]
        area = 0.0
        for i in range(n):
            x1, y1 = vertices[i]
            x2, y2 = vertices[(i+1) % n]
            area += x1 * y2 - y1 * x2
        return 0.5 * abs(area)

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


if __name__ == "__main__":
    # u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc",all_coordinates=False)
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro")
    cls2 = Area(u, {'DPPC': ['PO4'],'D3PC':['PO4'],'CHOL':['ROH']}, 10,file_path='E:/excel')

    # t1 = time.time()
    # cProfile.run('cls2.run()')
    # t2 = time.time()
    # print('运行时间：' ,t2- t1)
    # cls2.make_figure_2d(method='plot')
    # print(cls2.Area.shape)

    # u = mda.Universe("E:/awork/TEST/ach.gro")
    # cls1 = Area(u, {'DPPC':['PO4'],'DIPC':['PO4']}, k=12, file_path='E:/excel')
    # cls1.run()
    # cls1.make_figure_2d('bar')
    # cls2.vmd(arr1)

    # atm1 = u.select_atoms('resid 1 and name PO4').positions
    # atm2 = u.select_atoms('resname DPPC DAPC and name PO4').positions
    #
    # cls3 = Mocartoo(atm1[0],atm2,arr1)
    # cls3.make_figure()
