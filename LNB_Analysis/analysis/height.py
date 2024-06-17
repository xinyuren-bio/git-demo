"""
2024.5.20
2024.6.17
修改：
1.统一了格式问题
2.使用了新的mask方式（不同于lip）
问题：
同area
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
import pandas as pd

warnings.filterwarnings('ignore')


class Height(AnalysisBase):

    def __init__(self, universe, residuesGroup: dict, k=None, file_path=None):
        super().__init__(universe.trajectory)

        self.u = universe
        self.residues = [i for i in residuesGroup]
        self.k = k
        self.file_path = file_path

        self.headSp = {sp: ' '.join(residuesGroup[sp][0]) for sp in residuesGroup}
        self.tailSp = {sp: ' '.join(residuesGroup[sp][-1]) for sp in residuesGroup}

        self.headAtoms = self.u.atoms[[]]
        self.tailAtoms = self.u.atoms[[]]

        self.numSp = {}

        self._headBool = {}
        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]))

            self.tailAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.tailSp[sp]))

        self._n_residues = self.headAtoms.n_residues

        self._tailMask = {
            sp: self.tailAtoms.resnames == sp for sp in self.residues
        }

        self.numSp = {
            sp: i.n_residues for sp, i in self.headAtoms.groupby('resnames').items()
        }
        # ///////////////////////////////////////////////////////////////////////////
        # 待优化
        self._headMask = {}
        for x in range(len(self.residues)):
            if x == 0:
                num = 0
            sp = self.residues[x]
            arrBool = np.full(self._n_residues, False)
            arrBool[num:num+self.numSp[sp]] = True
            self._headMask[sp] = arrBool
            num += self.numSp[sp]

        self.results.Height = None

    def _prepare(self):
        self.results.Height = np.full([self._n_residues, self.n_frames],
                                      fill_value=np.NaN)

    def _single_frame(self):
        centerHeadSp = self.headAtoms.center_of_geometry(compound='residues')
        normals = self.find_k_nearest_neighbors_and_normals(centerHeadSp)
        centerTailSp = self.tailAtoms.center_of_geometry(compound='residues')
        for sp in self.residues:
            sp_normal = normals[self._headMask[sp]]
            head_to_tail = centerHeadSp[self._headMask[sp]] - centerTailSp[self._headMask[sp]]
            distance = (np.abs(np.einsum('ij,ij->i', sp_normal, head_to_tail)))
            self.results.Height[self._headMask[sp], self._frame_index] = distance

    @property
    def Height(self):
        return self.results.Height

    def fit_plane_and_get_normal(self, points):
        centered_points = points - np.mean(points, axis=0)
        covariance_matrix = np.cov(centered_points.T)
        eigenvalues, eigenvectors = eigh(covariance_matrix)
        min_eigenvalue_index = np.argmin(eigenvalues)
        normal = eigenvectors[:, min_eigenvalue_index]
        return normal / np.linalg.norm(normal)

    def find_k_nearest_neighbors_and_normals(self, particles):
        kdtree = KDTree(particles)
        normals = np.zeros((particles.shape[0], 3))
        for i, point in enumerate(particles):
            dists, idxs = kdtree.query(point, self.k + 1)
            nearest_neighbors = particles[idxs[1:]]
            if nearest_neighbors.shape[0] >= self.k:
                nearest_neighbors = np.vstack((nearest_neighbors, point))
                normal = self.fit_plane_and_get_normal(nearest_neighbors)
                normals[i] = normal
            else:
                normals[i] = np.zeros(3)
        return normals

    def ensure_consistent_normals(self, normals, vector):
        multi = normals * vector
        multi_naga = multi[:, ] < 0
        normals[multi_naga] = -normals[multi_naga]
        return normals

    def _conclude(self):
            # self.writeExcel()
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


    def writeExcel(self):
        columnFrame = [i*self.step for i in range(self.n_frames)]
        columnHead = ['resid'] + columnFrame
        resids = self.headAtoms.resids
        resPos = np.column_stack((resids, self.results.Height))
        df = pd.DataFrame(resPos)

        df.to_excel('%s/area.xlsx' % self.file_path, header=columnHead, index=False)


    def vmd(self, array=None):
        sz_array = np.zeros([self.u.atoms.n_atoms, ])
        for sp in self._residues:
            n_sp = self._mask_head[sp].sum()
            n_sp_atoms = sz_array[self._mask_all[sp]].shape[0] / n_sp
            if array is True:

                array /= 10
                array -= 0.7
                sz_array[self._mask_all[sp]] = np.repeat(array[self._mask_head[sp]], n_sp_atoms)
            else:
                array_max_value = np.max(self.results.Height)
                array = self.results.Height / 10
                array -= 0.5
                sz_array[self._mask_all[sp]] = np.repeat(array[self._mask_head[sp]], n_sp_atoms)

        lines = ['proc height {id} {',
                 'set sel [atomselect $id all]',
                 'set natom [$sel num]',
                 'set rduser [open "output_height.txt" r]',
                 'set userlist {}',
                 'for {set iatm 0} {$iatm<=[expr $natom-1]} {incr iatm} {',
                 'gets $rduser line',
                 'scan $line "%f" user',
                 'lappend userlist $user',
                 '}',
                 '$sel set user3 $userlist',
                 'close $rduser',
                 'puts "ALL DONE!"',
                 '}']
        if self.file_path:
            np.savetxt(self.file_path, sz_array, fmt='%.5f', delimiter='\n')
            with open('D:/VMD32/height.tcl', 'w', encoding='utf-8') as file:
                for line in lines:
                    file.write(line + '\n')
        else:
            np.savetxt('D:/VMD32/output_height.txt', sz_array, fmt='%.5f', delimiter='\n')
            with open('D:/VMD32/height.tcl', 'w', encoding='utf-8') as file:
                for line in lines:
                    file.write(line + '\n')

    def make_figure_3d(self, vmin=None, vmax=None, s=30, xlabel='X', ylabel='Y', zlabel='Z'):
        points_pos = self._sel_head.positions
        values = self.results.Height
        vmax_value = np.max(values)
        values /= vmax_value
        if vmin is None:
            vmin = min(values)
            vmax = max(values)
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

    def make_figure_2d(self, method):
        if method == 'bar':
            for sp in self.residues:
                plt.bar(sp, np.mean(np.mean(self.results.Height[self._headMask[sp]],axis=-1)), width=0.4, label='%s' % sp)
        elif method == 'plot':
            for sp in self.residues:
                plt.plot(np.mean(self.results.Height[self._headMask[sp]], axis=0), label='%s' % sp)
        plt.legend()
        plt.xlabel('Frames')
        plt.ylabel('Lipid Height(A)')
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc", all_coordinates=False)
    dict_residue = {'DPPC': (['PO4'], ['C4A']),'D3PC': (['PO4'], ['C4A']),'CHOL':(['ROH'],['R5'])}
    cls2 = Height(u,dict_residue,k=12)
    cls2.run(500,1000,1,verbose=True)
    cls2.make_figure_2d('bar')