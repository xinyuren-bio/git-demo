import pandas as pd
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
from MDAnalysis.lib.log import ProgressBar
from scipy.linalg import eig
import warnings
warnings.filterwarnings('ignore')


class PCA(AnalysisBase):
    def __init__(self, u, residues, path):
        super().__init__(u.trajectory)
        self.u = u
        self._residue = residues
        self.filePath = path
        self.results.pca = None
        self.atmSel = self.u.select_atoms('resname %s' % ' '.join(self._residue))

    def _prepare(self):
        self.results.pca = np.full([self.n_frames], fill_value=np.NaN)

    def _single_frame(self):
        atmCenterPos = self.atmSel.center_of_mass(compound='residues')
        atmsCenter = self.atmSel.center_of_mass()
        atmToCenter = atmCenterPos - atmsCenter
        covMatix = np.cov(atmToCenter, rowvar=False)
        eigenValues, _ = eig(covMatix)
        sorted_indices = np.argsort(eigenValues)[::-1]
        eigenvalues = eigenValues[sorted_indices]
        varianceRatio = eigenvalues[-1] / eigenvalues[0]
        self.results.pca[self._frame_index] = varianceRatio

    def _conclude(self):
        self._writeExcel()

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

    def _writeExcel(self):
        columnFrame = [self.start + (i*self.step) for i in range(self.n_frames)]
        res = np.column_stack((columnFrame, self.results.pca))

        df = pd.DataFrame(res)
        df.to_excel('%s/pca.xlsx' % self.filePath, header=['Frame', '形状比'],index=False)


if __name__ == '__main__':
    import MDAnalysis as mda
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc",all_coordinates=False)
    cls2 = PCA(u, ['DPPC','D3PC','CHOL'], path='E:/excel/pca.xlsx')
    cls2.run(500,1500,1,verbose=True)