# ///////////////////////////////////////////////////////////////
#
# BY: RenXinYu
# PROJECT MADE WITH: Qt Designer and PySide6
# V: 0.0.0
#
#
#
# ///////////////////////////////////////////////////////////////
# MAIN FILE
# ///////////////////////////////////////////////////////////////

from threading import Thread
from main import *
# WITH ACCESS TO MAIN WINDOW WIDGETS
# ///////////////////////////////////////////////////////////////


class AppFunctions(MainWindow):
    def setThemeHack(self):
        Settings.BTN_LEFT_BOX_COLOR = "background-color: #495474;"
        Settings.BTN_RIGHT_BOX_COLOR = "background-color: #495474;"
        Settings.MENU_SELECTED_STYLESHEET = MENU_SELECTED_STYLESHEET = """
        border-left: 22px solid qlineargradient(spread:pad, x1:0.034, y1:0, x2:0.216, y2:0, stop:0.499 rgba(255, 121, 198, 255), stop:0.5 rgba(85, 170, 255, 0));
        background-color: #566388;
        """

        # SET MANUAL STYLES
        self.ui.lineEdit.setStyleSheet("background-color: #6272a4;")
        self.ui.pushButton.setStyleSheet("background-color: #6272a4;")
        self.ui.plainTextEdit.setStyleSheet("background-color: #6272a4;")
        self.ui.tableWidget.setStyleSheet(
            "QScrollBar:vertical { background: #6272a4; } QScrollBar:horizontal { background: #6272a4; }")
        self.ui.scrollArea.setStyleSheet(
            "QScrollBar:vertical { background: #6272a4; } QScrollBar:horizontal { background: #6272a4; }")
        self.ui.comboBox.setStyleSheet("background-color: #6272a4;")
        self.ui.horizontalScrollBar.setStyleSheet("background-color: #6272a4;")
        self.ui.verticalScrollBar.setStyleSheet("background-color: #6272a4;")
        self.ui.commandLinkButton.setStyleSheet("color: #ff79c6;")

    # Height Edit Display
    # ///////////////////////////////////////////////////////////////
    def btnUploadClick(self, btn_sender):
        btnName = btn_sender.objectName()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        if btnName == 'btn_gro_height':
            filter = "GRO Files (*.gro);;TPR Files (*.tpr);;PDB Files (*.pdb)"
            filename, _ = QFileDialog.getOpenFileName(widgets, "选择文件", "",
                                                      filter)
            if filename:
                self.gro_path = filename
                self.edit_gro_height.setText(self.gro_path)
            else:
                self.gro_path = None
                self.edit_gro_height.setText('导入文件失败！')

        if btnName == 'btn_xtc_height':
            filter = "xtc Files (*.xtc)"
            filename, _ = QFileDialog.getOpenFileName(widgets, "选择文件", "",
                                                      filter)
            if filename:
                # 设置文件路径
                self.xtc_path = filename
                self.edit_xtc_height.setText(self.xtc_path)
            else:
                self.xtc_path = None
                self.edit_xtc_height.setText('导入文件失败！')

        if btnName == 'btnSaveHeight':
            options = QFileDialog.Options()
            options |= QFileDialog.ShowDirsOnly
            folder_selected = QFileDialog.getExistingDirectory(widgets, "选择文件夹", QDir.currentPath(),
                                                               options=options)

            if folder_selected:
                self.save_path = folder_selected
                self.edit_save_height.setText(self.save_path)
            else:
                self.edit_save_height.setText('选取路径失败！')
    # Next Click
    # ///////////////////////////////////////////////////////////////


class NextBtnClick:
    def __init__(self, ui):
        self.ui = ui

        self.firstFrame = int(self.ui.edit_first_height.text())
        self.lastFrame = int(self.ui.lineEdit_11.text())
        self.stepFrame = int(self.ui.lineEdit_10.text())
        self.k = int(self.ui.lineEdit_12.text())
        self.runMethod = self.ui.comboBox_2.currentText()
        self.path = self.ui.edit_save_height.text()

        self.u = mda.Universe(self.ui.gro_path, self.ui.xtc_path, all_coordinates=False, update=False)
        self.residues = np.unique(self.u.atoms.resnames)

        self.checkBoxs = []
        self.residueClick = []
        self.residueAtoms = {}
        self.checkAtoms = []
        self.heightHeadAtoms = {}
        self.heightTailAtoms = {}
        self.radioAtoms = {}
        # 函数执行
        self.layOut()

    def compute(self):
        if self.runMethod == 'Height':
            self.Height()
        if self.runMethod == 'SZ':
            self.SZ()
        if self.runMethod == 'Area':
            self.Area()
        if self.runMethod == 'PCA':
            self.PCA()
        if self.runMethod == 'MeanCurvature':
            self.MeanCurvature()
        if self.runMethod == 'GussainCurvature':
            self.GussainCurvature()


    def layOut(self):
        """
        进行布局
        进行需要分析的残基的选择
        """
        self.ui.VLayoutRightAll = QVBoxLayout(self.ui.extraRightBox)
        self.ui.frameResidue = QFrame()
        self.ui.VLayoutRightResidue = QVBoxLayout(self.ui.frameResidue)
        for sp in self.residues:
            checkbox = QCheckBox(sp)
            checkbox.setAutoFillBackground(False)
            checkbox.setStyleSheet(u"")
            self.ui.VLayoutRightResidue.addWidget(checkbox)
            self.checkBoxs.append(checkbox)
        if self.runMethod == 'PCA':
            self.ui.btnNextResidue = QPushButton('RUN!')
            self.ui.btnNextResidue.clicked.connect(lambda: self.compute())
        else:
            self.ui.btnNextResidue = QPushButton('Next')
            self.ui.btnNextResidue.clicked.connect(lambda: self.getClickResidue())
        self.ui.btnNextResidue.setCursor(QCursor(Qt.PointingHandCursor))
        self.ui.btnNextResidue.setStyleSheet("background-color: #6272a4;")
        self.ui.VLayoutRightResidue.addWidget(self.ui.btnNextResidue)
        self.ui.VLayoutRightAll.addWidget(self.ui.frameResidue)
        self.ui.extraRightBox.setLayout(self.ui.VLayoutRightAll)

    def getClickResidue(self):
        """
        获得勾选的残基名称
        """
        self.residueClick.clear()
        for box in self.checkBoxs:
            if box.isChecked():
                self.residueClick.append(box.text())
        self.residueNext()

    def residueNext(self):
        """
        获得勾选的残基的原子
        :return:
        """
        self.residueAtoms.clear()

        def getAtom():
            for sp in self.residueClick:
                self.residueAtoms[sp] = np.unique(self.u.select_atoms('resname %s' % sp).names)

        if self.runMethod == 'Height':
            getAtom()
            self.DisplayAtom('选择头部原子(Head Atom)', 'Select Tail', QCheckBox, method='Height')
        if self.runMethod == 'SZ':
            getAtom()
            self.DisplayAtom('请选择头部原子(Head Atom)', 'Next', QRadioButton, method='SZ')
        if self.runMethod == 'Area':
            getAtom()
            self.DisplayAtom('请选择代表性原子', 'RUN!', QRadioButton, method='Area')

    @staticmethod
    def clearLayout(layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)

    def DisplayAtom(self, label: str, btnName, btnType, method):
        self.clearLayout(self.ui.VLayoutRightAll)
        self.checkAtoms = []
        widgetAtoms = QWidget()
        widgetAtoms.setStyleSheet(u"background-color: black;")
        self.ui.VLayoutRightAtoms = QVBoxLayout(widgetAtoms)
        self.ui.labelHead = QLabel(label)
        self.ui.VLayoutRightAtoms.addWidget(self.ui.labelHead)

        def makeGroup(Check_Ratio):
            for sp in self.residueClick:
                groupbox = QGroupBox(sp)
                groupbox.setStyleSheet("""QGroupBox {background-color: black;}
                                                   QGroupBox::title {color: white;}
                                                   QGroupBox QCheckBox {color: white;}""")
                grouplayout = QVBoxLayout(groupbox)
                if self.runMethod == 'SZ' and btnName == 'RUN!':
                    for chain in ['Chain A', 'Chain B', 'Chain A + Chain B']:
                        checkbox = Check_Ratio(chain)
                        grouplayout.addWidget(checkbox)
                else:
                    for atm in self.residueAtoms[sp]:
                        checkbox = Check_Ratio(atm)
                        grouplayout.addWidget(checkbox)
                self.checkAtoms.append(groupbox)
                self.ui.VLayoutRightAtoms.addWidget(groupbox)

        self.ui.btnNextCircle = QPushButton(btnName)
        self.ui.btnNextCircle.setStyleSheet(u"color: rgb(52, 59, 72);")
        makeGroup(btnType)
        scrollArea = QScrollArea()
        scrollArea.setWidget(widgetAtoms)


        scrollArea.setWidgetResizable(True)
        self.ui.VLayoutRightAll.addWidget(scrollArea)
        self.ui.VLayoutRightAll.addWidget(self.ui.btnNextCircle)
        self.ui.extraRightBox.setLayout(self.ui.VLayoutRightAll)

        if method == 'Height':
            if btnName == 'Select Tail':
                self.ui.btnNextCircle.clicked.connect(
                    lambda: self.checkGroup(self.checkAtoms, self.heightHeadAtoms, QCheckBox))
                self.ui.btnNextCircle.clicked.connect(lambda: self.DisplayAtom(
                    '请选择尾部原子(Tail Atom)', 'RUN!', QCheckBox, 'Height'))
            if btnName == 'RUN!':
                self.ui.btnNextCircle.clicked.connect(
                    lambda: self.checkGroup(self.checkAtoms, self.heightTailAtoms, QCheckBox))
                self.ui.btnNextCircle.clicked.connect(lambda: self.compute())

        if method == 'Area':
            self.ui.btnNextCircle.clicked.connect(lambda :self.checkGroup(self.checkAtoms, self.radioAtoms, QRadioButton))
            self.ui.btnNextCircle.clicked.connect(lambda :self.compute())

        if method == 'SZ':
            if btnName == 'Next':
                self.ui.btnNextCircle.clicked.connect(
                    lambda: self.checkGroup(self.checkAtoms, self.heightHeadAtoms, QRadioButton)
                )
                self.ui.btnNextCircle.clicked.connect(lambda: self.DisplayAtom(
                    '请选择需要分析的链', 'RUN!', QRadioButton, 'SZ'))
            if btnName == 'RUN!':
                self.ui.btnNextCircle.clicked.connect(lambda :self.checkGroup(self.checkAtoms, self.radioAtoms, QRadioButton))
                self.ui.btnNextCircle.clicked.connect(lambda :self.compute())

    def addProgressBar(self):
        self.ui.progressBar = QProgressBar()
        self.ui.progressBar.setStyleSheet("""
                            QProgressBar {
                                border: 2px solid grey; /* 设置进度条边框 */
                                border-radius: 5px;    /* 设置进度条圆角 */
                            }
                            QProgressBar::chunk {
                                background-color: #6272a4 ; /* 设置进度条填充颜色 */
                                width: 20px;           /* 设置进度条块的宽度，从而加粗进度条 */
                            }
                        """)
        self.ui.VLayoutRightAll.addWidget(self.ui.progressBar)

    def updateProgressBar(self, value):
        self.ui.progressBar.setValue(value)

    @staticmethod
    def checkGroup(groupbox, save, childrenType):
        save.clear()
        for group in groupbox:
            listAtoms = []
            for box in group.findChildren(childrenType):
                if box.isChecked():
                    listAtoms.append(box.text())
                save[group.title()] = listAtoms

    @staticmethod
    def keyValue(dict1, dict2):
        dict3 = {}
        for key in dict1.keys() & dict2.keys():
            atoms_1 = dict1[key]
            atoms_2 = dict2[key]
            atomsAdd = (atoms_1, atoms_2)
            dict3[key] = atomsAdd
        return dict3

    def Height(self):
        self.addProgressBar()
        dictAtomsHeadTail = self.keyValue(self.heightHeadAtoms, self.heightTailAtoms)
        self.cls = Height(self.u, dictAtomsHeadTail, self.path, k=self.k)

        worker = Worker(self.cls, self.firstFrame, self.lastFrame, self.stepFrame)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)

        self.thread.start()


    def SZ(self):
        self.addProgressBar()
        print(self.keyValue(self.heightHeadAtoms, self.radioAtoms))
        self.cls = SZ(self.u, self.keyValue(self.heightHeadAtoms, self.radioAtoms), self.path, k=self.k)

        worker = Worker(self.cls, self.firstFrame, self.lastFrame, self.stepFrame)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)

        self.thread.start()



    def Area(self):
        self.addProgressBar()
        self.cls = Area(universe=self.u, residue_group=self.radioAtoms, file_path=self.path, k=self.k)
        worker = Worker(self.cls, self.firstFrame, self.lastFrame, self.stepFrame)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)
        self.thread.start()

    def PCA(self):
        self.addProgressBar()
        self.cls = PCA(u=self.u, residues=self.residues,path=self.path)
        worker = Worker(self.cls, self.firstFrame, self.lastFrame, self.stepFrame)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)
        self.thread.start()

    def MeanCurvature(self):
        pass

    def GussainCurvature(self):
        pass

class Worker(QObject):
    progressValueChanged = Signal(int)

    def __init__(self, cls_instance, *args, **kwargs):
        super().__init__()
        self.cls_instance = cls_instance
        self.args = args
        self.kwargs = kwargs

    def run(self):
        self.cls_instance.run(*self.args, **self.kwargs,callBack=self.update_progress)

    def update_progress(self, value):
        # 发出进度更新信号
        self.progressValueChanged.emit(value)
