import pandas as pd
import matplotlib.pyplot as plt


class Figure:
    def __init__(self, filePath, methodFigure):
        self.filepath = filePath
        self.methodFigure = methodFigure

        self.resultBar = {}
        self.resultPlot = {}

        self.df = pd.read_excel(self.filepath)

    def readFile(self):
        groups = self.df.groupby('resname')
        for type_name, group_df in groups:
            data_columns = group_df.iloc[1:, 2:]
            if self.methodFigure == 'bar':
                # data_columns
                self.resultBar[type_name] = data_columns.mean(axis=1).mean()
            if self.methodFigure == 'plot':
                self.resultPlot[type_name] = data_columns.mean(axis=0)

    def figure(self):
        if self.methodFigure == 'bar':
            for sp in self.resultBar.keys():
                plt.bar(sp, self.resultBar[sp], width=0.4, label=sp)
        elif self.methodFigure == 'plot':
            for sp in self.resultPlot.keys():
                print(self.resultPlot[sp].shape)
                plt.plot( self.resultPlot[sp], label=sp)
        plt.legend()
        # plt.xlabel(x_axix)
        # plt.ylabel(y_axis)
        plt.grid(True)
        plt.show()

    def run(self):
        self.readFile()
        self.figure()


if __name__ == '__main__':
    cls = Figure('E:/excel/area.xlsx', 'plot')
    cls.run()






