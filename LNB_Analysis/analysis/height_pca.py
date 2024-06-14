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
from mocarto import Mocartoo
import pandas as pd

warnings.filterwarnings('ignore')


class Height(AnalysisBase):

    def __init__(self, universe, residue_group: dict, file_path, k=10, method=None, callBack=None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues_dict = residue_group
        self.k = k
        self.file_path = file_path
        self.figureMethod = method
        self.call = callBack
        self._residues = list(self.residues_dict.keys())
        self._head = [' '.join(i[0]) for i in self.residues_dict.values()]
        self._tail = [' '.join(i[1]) for i in self.residues_dict.values()]

        self._sel_head = self.u.atoms[[]]
        self._sel_tail = self.u.atoms[[]]
        self._mask = self.u.atoms[[]]

        # for i in range(len(self._residues)):

        for i in range(len(self._residues)):
            self._sel_head += self.u.select_atoms('resname %s and name %s'
                                                  % (self._residues[i], self._head[i]))

        for i in range(len(self._residues)):
            self._sel_tail += self.u.select_atoms('resname %s and name %s'
                                                  % (self._residues[i], self._tail[i]))


        self._n_residues = self._sel_head.n_residues

        self._mask_head = {
            sp: self._sel_head.resnames == sp for sp in self._residues
        }

        self._mask_tail = {
            sp: self._sel_tail.resnames == sp for sp in self._residues
        }

        self._mask_all = {
            sp: self.u.atoms.resnames == sp for sp in self._residues
        }

        self.results.SZ = None

    def _prepare(self):
        self.results.Height = np.full([self._sel_head.n_residues, self.n_frames],
                                      fill_value=np.NaN)
        self.meanHeadPositions = np.full(shape=[self._sel_head.n_atoms, 3],
                                         fill_value=np.NaN)

    def _single_frame(self):
        for sp in self._residues:
            n_sp1 = self._mask_head[sp].sum()
            sp_head = self._sel_head[self._mask_head[sp]]
            spHeadPos = sp_head.positions.reshape([n_sp1, -1, 3])
            self.meanHeadPositions[self._mask_head[sp]] = np.mean(spHeadPos, axis=1)
        normals = self.find_k_nearest_neighbors_and_normals(self.meanHeadPositions, self.k)
        for sp in self._residues:
            sp_normal = normals[self._mask_head[sp]]
            n_sp = self._mask_head[sp].sum()
            sp_tail = self._sel_tail[self._mask_tail[sp]]

            tail_positions = sp_tail.positions.reshape([n_sp, -1, 3])
            tail_center = np.mean(tail_positions, axis=1)
            head_to_tail = self.meanHeadPositions[self._mask_head[sp]] - tail_center
            distance = (np.abs(np.einsum('ij,ij->i', sp_normal, head_to_tail)))
            self.results.Height[self._mask_head[sp], self._frame_index] = distance

    @property
    def Height(self):
        return self.results.Height

    def fit_plane_and_get_normal(self, points):
        mean = np.mean(points, axis=0)
        centered_points = points - mean
        covariance_matrix = np.cov(centered_points.T)
        eigenvalues, eigenvectors = eigh(covariance_matrix)
        min_eigenvalue_index = np.argmin(eigenvalues)
        normal = eigenvectors[:, min_eigenvalue_index]
        normal = normal / np.linalg.norm(normal)
        return normal

    def find_k_nearest_neighbors_and_normals(self, particles, k):
        kdtree = KDTree(particles)
        normals = np.zeros((particles.shape[0], 3))
        for i, point in enumerate(particles):
            dists, idxs = kdtree.query(point, k + 1)
            nearest_neighbors = particles[idxs[1:]]
            if nearest_neighbors.shape[0] >= k:
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
            self.writeExcel()

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
            callBack((i * self.step/(self.stop - self.start) + 0.01) * 100)
        self._conclude()
        return self


    def writeExcel(self):
        columnFrame = [i*self.step for i in range(self.n_frames)]
        columnHead = ['resid'] + columnFrame
        resids = self._sel_head.resids
        resPos = np.column_stack((resids, self.results.Area))
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
        data = self.results.Height
        if method == 'bar':
            for sp in self._residues:
                plt.bar(sp, np.mean(np.mean(data[self._mask_head[sp]],axis=-1)), width=0.4, label='%s' % sp)
        elif method == 'plot':
            for sp in self._residues:
                plt.plot(np.mean(data[self._mask_head[sp]], axis=0), label='%s' % sp)
        plt.legend()
        plt.xlabel('Frames')
        plt.ylabel('Lipid Height(A)')
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc",all_coordinates=False)
    dict_residue = {'DPPC': (['PO4'], ['C4A', 'C4B']), 'CHOL': (['ROH'], ['R5']), 'D3PC': (['PO4'], ['C5A'])}
    cls2 = Height(u,dict_residue,k=12,method='bar')
    cls2.run(1,100,1,verbose=True)
