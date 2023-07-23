#!/usr/bin/env python
from PyQt5 import QtWidgets, QtCore, QtGui
import os
import shutil
import PyDP4
import queue
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
from scipy.stats import norm
from PyQt5.QtSvg import QSvgWidget
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import rdMolDraw2D
import pickle
from matplotlib import cm

from pathlib import Path

class Window(QtWidgets.QMainWindow):
    # create signal to start background job

    signal_start_PyDP4 = QtCore.pyqtSignal()

    # signal_finished_PyDP4 = QtCore.pyqtSignal()

    def __init__(self):
        super(Window, self).__init__()

        self.table_widget = TabWidget(self)
        self.setCentralWidget(self.table_widget)

        self.resize(self.table_widget.Tab1.Overall_widget.sizeHint())
        #self.resize(self.table_widget.sizeHint())
        self.show()

        # check if DP4outfile is there


class TabWidget(QtWidgets.QWidget):

    def __init__(self, parent):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.layouttabs = QtWidgets.QVBoxLayout()

        # Initialize tab screen
        self.tabs = QtWidgets.QTabWidget()

        self.tab1 = QtWidgets.QWidget()

        self.Tab1 = CalculationTab()

        self.tabs.addTab(self.Tab1, "Calculation")
        self.tab1.layout = QtWidgets.QVBoxLayout()

        # Add tabs to widget
        self.layouttabs.addWidget(self.tabs)

        self.setLayout(self.layouttabs)

    def addplottabs(self):

        self.tab2 = QtWidgets.QWidget()
        self.tab3 = QtWidgets.QWidget()
        self.tab4 = QtWidgets.QWidget()
        self.tab5 = QtWidgets.QWidget()
        self.tab6 = QtWidgets.QWidget()

        self.tabs.addTab(self.tab5, "Conformers")

        #only add proton/carbon plot tabs if the raw data is there

        NMR_list = [str(self.Tab1.NMR_list.item(i).text()) for i in range(self.Tab1.NMR_list.count())]

        if 'Proton' in NMR_list or 'Proton.dx' in NMR_list or 'proton' in NMR_list or 'proton.dx' in NMR_list:

            self.tabs.addTab(self.tab2, "Proton Plot")
            self.Tab2 = ProtonPlotTab()
            self.tab2.layout = QtWidgets.QVBoxLayout(self)
            self.tab2.layout.addWidget(self.Tab2)
            self.tab2.setLayout(self.tab2.layout)

        if 'Carbon' in NMR_list or 'Carbon.dx' in NMR_list or 'carbon' in NMR_list or 'carbon.dx' in NMR_list:

            self.tabs.addTab(self.tab3, "Carbon Plot")
            self.Tab3 = CarbonPlotTab()
            self.tab3.layout = QtWidgets.QVBoxLayout(self)
            self.tab3.layout.addWidget(self.Tab3)
            self.tab3.setLayout(self.tab3.layout)

        self.Tab5 = ConformerTab()
        self.tab5.layout = QtWidgets.QVBoxLayout(self)
        self.tab5.layout.addWidget(self.Tab5)
        self.tab5.setLayout(self.tab5.layout)

        self.Tab6 = DP5Tab()

        if 's' in ui.table_widget.Tab1.settings.Workflow:
            self.tabs.addTab(self.tab6, "DP5")
            self.tab6.layout = QtWidgets.QVBoxLayout(self)
            self.tab6.layout.addWidget(self.Tab6)
            self.tab6.setLayout(self.tab6.layout)

        self.Tab4 = StatsTab()
        if 's' in ui.table_widget.Tab1.settings.Workflow and len(ui.table_widget.Tab1.settings.InputFiles) > 1:
            self.tabs.addTab(self.tab4, "DP4")
            self.tab4.layout = QtWidgets.QVBoxLayout(self)
            self.tab4.layout.addWidget(self.Tab4)
            self.tab4.setLayout(self.tab4.layout)




        self.layouttabs.addWidget(self.tabs)
        self.setLayout(self.layouttabs)

        ui.update()


class DP5Tab(QtWidgets.QWidget):

    def __init__(self):
        super(DP5Tab, self).__init__()

        # main horizontal widget

        self.main_layout = QtWidgets.QHBoxLayout()
        self.setLayout(self.main_layout)
        self.Isomers = ui.table_widget.Tab1.worker.Isomers
        self.DP5data = ui.table_widget.Tab1.worker.DP5Data
        self.DP4data = ui.table_widget.Tab1.worker.DP4Data

        self.Isomer_DP5_table = QtWidgets.QTableWidget(self)

        self.Isomer_DP5_table.setColumnCount(3)

        self.Isomer_DP5_table.setHorizontalHeaderLabels(["Isomer", "DP5", "DP4"])

        self.Isomer_DP5_table.setRowCount(len(self.DP5data.DP5scaledprobs))

        c = 0

        for scaled_DP5, DP4 in zip(self.DP5data.DP5scaledprobs, self.DP4data.DP4probs):

            self.Isomer_DP5_table.setItem(c, 0, QtWidgets.QTableWidgetItem(str(c + 1)))

            self.Isomer_DP5_table.setItem(c, 1, QtWidgets.QTableWidgetItem(str(round(scaled_DP5, 2))))

            self.Isomer_DP5_table.setItem(c, 2, QtWidgets.QTableWidgetItem(str(round(DP4 , 3))))

            c += 1

        self.Isomer_DP5_table.selectRow(0)

        # add error table

        self.DP5_table = QtWidgets.QTableWidget(self)

        self.DP5_table.setColumnCount(4)

        self.DP5_table.setHorizontalHeaderLabels(
            ["Atom Label", "Atomic Probability",  "Scaled error (ppm)", "Unscaled error (ppm)"])

        self.DP5_table.setColumnWidth(0, 70)
        self.DP5_table.setColumnWidth(1, 70)
        self.DP5_table.setColumnWidth(2, 70)
        self.DP5_table.setColumnWidth(3, 70)

        #########

        # add isomer canvas

        self.table_widget = QtWidgets.QWidget(self)

        self.image = QSvgWidget()

        self.image.setMinimumWidth(self.image.height())

        # organise layouts

        self.table_layout = QtWidgets.QVBoxLayout()

        self.table_layout.addWidget(self.Isomer_DP5_table)

        self.table_layout.addWidget(self.DP5_table)

        #self.image_layout.addWidget(self.image)

        self.table_widget.setLayout(self.table_layout)

        self.main_layout.addWidget(self.table_widget)

        self.main_layout.addWidget(self.image)

        # connections

        self.Isomer_DP5_table.currentCellChanged.connect(self.Mol_drawer)

        self.Isomer_DP5_table.currentCellChanged.connect(self.populatetable)

        self.DP5_table.currentCellChanged.connect(self.Mol_drawer)

    def populatetable(self):

        isomerindex = self.Isomer_DP5_table.currentRow()

        # which isomer is selected

        self.DP5_table.setRowCount(0)

        self.labels, self.shifts, self.exps, self.scaleds = self.DP5data.Clabels[isomerindex], \
                                                            self.DP5data.Cshifts[isomerindex], \
                                                            self.DP5data.Cexp[isomerindex], \
                                                            self.DP5data.Cscaled[isomerindex], \

        self.scaled_probs = self.DP5data.BScaledAtomProbs[isomerindex]

        self.DP5_table.setRowCount(len(self.labels))

        # set info in rows and columns

        c = 0

        for label, shift, exp, scaled_shift,  scaled_prob in zip(self.labels, self.shifts, self.exps,
                                                                               self.scaleds,
                                                                               self.scaled_probs):
            self.DP5_table.setItem(c, 0, QtWidgets.QTableWidgetItem(label))
            self.DP5_table.setItem(c, 1, QtWidgets.QTableWidgetItem(str(round(scaled_prob, 2))))

            self.DP5_table.setItem(c, 2, QtWidgets.QTableWidgetItem(str(round(abs(scaled_shift - exp), 2))))
            self.DP5_table.setItem(c, 3, QtWidgets.QTableWidgetItem(str(round(abs(shift - exp), 2))))

            c += 1

        self.DP5_table.selectRow(0)

        # self.plot()

    def Isomer_number(self):
        Isomer_list = []

        for c, i in enumerate(self.Isomers):
            Isomer_list.append("Isomer " + str(c + 1))

        return Isomer_list

    def Mol_drawer(self):

        isomerindex = self.Isomer_DP5_table.currentRow()

        highlight = {}

        inds = []

        for label in self.DP5data.Clabels[isomerindex]:
            inds.append(int(label.split("C")[1]) - 1)

        for i, s in zip(inds, self.DP5data.BScaledAtomProbs[isomerindex]):

            highlight[i] = cm.RdYlGn(s)

        # pick selected atom

        if self.DP5_table.item(self.DP5_table.currentRow(), 0):

            error_ind = int(self.DP5_table.currentRow())

            atom_ind = inds[error_ind]

            highlight[atom_ind] = cm.viridis(  self.DP5data.BScaledAtomProbs[isomerindex][error_ind])

        #m = Chem.MolFromMolFile(str(ui.table_widget.Tab1.worker.settings.InputFiles[0]).split('.sdf')[0] + '.sdf', removeHs=False)

        m = Chem.MolFromMolFile(str(ui.table_widget.Tab1.worker.settings.InputFilesPaths[isomerindex]).split('.sdf')[0] + '.sdf',removeHs=False)

        Chem.Compute2DCoords(m)

        drawer = rdMolDraw2D.MolDraw2DSVG(self.image.width(), self.image.height())

        drawer.DrawMolecule(m, highlightAtoms=inds, highlightAtomColors=highlight, highlightBonds=False)

        # ,addStereoAnnotation = True,addAtomIndices = True

        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '')

        svg_bytes = bytearray(svg, encoding='utf-8')

        self.image.renderer().load(svg_bytes)

        ui.update()


class StatsTab(QtWidgets.QWidget):

    def __init__(self):
        super(StatsTab, self).__init__()

        #main horizontal widget

        self.main_layout = QtWidgets.QHBoxLayout()
        self.setLayout(self.main_layout)
        self.Isomers = ui.table_widget.Tab1.worker.Isomers
        self.dp4data = ui.table_widget.Tab1.worker.DP4Data

        # add isomer selection
        self.IsomerSelect = QtWidgets.QComboBox(self)
        self.IsomerSelect.addItems(self.Isomer_number())

        # add isomer canvas

        self.image = QSvgWidget()

        #make vertical isomer layout

        self.isomer_layout = QtWidgets.QVBoxLayout()
        self.isomer_layout.addWidget(self.IsomerSelect)
        self.isomer_layout.addWidget(self.image)

        self.isomer_widget = QtWidgets.QWidget(self)
        self.isomer_widget.setLayout(self.isomer_layout)

        # make error table

        # populate table based on the Isomer selected

        self.Hplot = plotstats('H')

        #check if there is a DP4 result for this atomt type:

        if not np.isnan(ui.table_widget.Tab1.worker.DP4Data.HDP4probs).any() :

            self.IsomerSelect.currentIndexChanged.connect(self.Hplot.populatetable)

            self.IsomerSelect.currentIndexChanged.connect(self.RenderImage)

            self.Hplot.errortable.itemSelectionChanged.connect(self.RenderImage)

        self.Cplot = plotstats('C')

        if not np.isnan(ui.table_widget.Tab1.worker.DP4Data.CDP4probs).any() :

            self.IsomerSelect.currentIndexChanged.connect(self.Cplot.populatetable)

            self.IsomerSelect.currentIndexChanged.connect(self.RenderImage)

            self.Cplot.errortable.itemSelectionChanged.connect(self.RenderImage)

        #make error table layout

        self.error_layout = QtWidgets.QVBoxLayout()

        self.error_layout.addWidget(self.Cplot.errortable)
        self.error_layout.addWidget(self.Hplot.errortable)

        self.error_widget = QtWidgets.QWidget(self)
        self.error_widget.setLayout(self.error_layout)

        #make plot layout

        self.plot_layout = QtWidgets.QVBoxLayout()

        self.plot_layout.addWidget(self.Cplot.statscanvas)
        self.plot_layout.addWidget(self.Hplot.statscanvas)

        self.plot_widget = QtWidgets.QWidget(self)
        self.plot_widget.setLayout(self.plot_layout)

        # add these widgets to main horizontal layout

        self.isomer_widget.setMinimumWidth(400)

        self.error_widget.setFixedWidth(400)

        self.main_layout.addWidget(self.isomer_widget)
        self.main_layout.addWidget(self.error_widget)
        self.main_layout.addWidget(self.plot_widget)

    def Isomer_number(self):
        Isomer_list = []

        for c, i in enumerate(self.Isomers):
            Isomer_list.append("Isomer " + str(c + 1))

        return Isomer_list

    def RenderImage(self):

        '''colors = [(0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17), (0.84, 0.15, 0.16), (0.58, 0.4, 0.74),
                  (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.5, 0.5, 0.5), (0.74, 0.74, 0.13), (0.09, 0.75, 0.81)]'''



        isomerindex = self.IsomerSelect.currentIndex()

        atom = []

        highlight = {}

        if self.Cplot.errortable.item(self.Cplot.errortable.currentRow(), 0):

            C_atom = int(self.Cplot.errortable.item(self.Cplot.errortable.currentRow(), 0).text().split("C")[1])

            atom.append(C_atom - 1 )

            highlight[C_atom - 1] = (0.12, 0.47, 0.71)

        if self.Hplot.errortable.item(self.Hplot.errortable.currentRow(), 0):

            H_atom = int(self.Hplot.errortable.item(self.Hplot.errortable.currentRow(), 0).text().split("H")[1])

            atom.append(H_atom - 1)

            highlight[H_atom - 1] = (0.84, 0.15, 0.16)

        m = Chem.MolFromMolFile(str(ui.table_widget.Tab1.worker.settings.InputFilesPaths[isomerindex]).split('.sdf')[0] + '.sdf',removeHs=False)
        Chem.Compute2DCoords(m)

        drawer = rdMolDraw2D.MolDraw2DSVG(self.image.width(), self.image.height())

        drawer.DrawMolecule(m, highlightAtoms=atom, highlightAtomColors=highlight,highlightBonds = False)

        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '')

        svg_bytes = bytearray(svg, encoding='utf-8')

        self.image.renderer().load(svg_bytes)

        #self.image.setGeometry(QtCore.QRect(0, 50, 500, 500))

        ui.update()


class plotstats(QtWidgets.QWidget):

    def __init__(self, Atom):

        super(plotstats, self).__init__()

        # when a cell is selected signal to plot stats

        self.atom = Atom

        self.findmeans()

        self.dp4data = ui.table_widget.Tab1.worker.DP4Data

        self.Isomers = ui.table_widget.Tab1.worker.Isomers

        self.errortable = QtWidgets.QTableWidget(self)

        self.errortable.setColumnCount(5)

        self.errortable.setHorizontalHeaderLabels(["Atom Label", "Calc Shift", "Scaled", "Exp", "Error"])

        self.errortable.setColumnWidth(0, 70)
        self.errortable.setColumnWidth(1, 70)
        self.errortable.setColumnWidth(2, 70)
        self.errortable.setColumnWidth(3, 70)
        self.errortable.setColumnWidth(4, 70)

        self.statsfigure = Figure()

        self.statscanvas = FigureCanvas(self.statsfigure)

        # self.errortable.cellClicked.connect(self.plot)

        self.errortable.itemSelectionChanged.connect(self.plot)

        self.statscanvas.mpl_connect('button_press_event', self.selectpoint)

    def plot(self):

        r = int(self.errortable.currentRow())

        if r >= 0:

            error = float(self.errortable.item(r, 4).text())

            # decide which stats params have been used

            self.statsfigure.clear()

            self.statsfig = self.statsfigure.add_subplot(111)

            if self.atom == 'H':

                m = abs(max([item for sublist in self.dp4data.Herrors for item in sublist], key=abs))

                p_color = 3

            elif self.atom == 'C':

                m = abs(max([item for sublist in self.dp4data.Cerrors for item in sublist], key=abs))

                p_color = 0

            # plot all errors at low transparency

            for e in self.errors:

                self.statsfig.plot([e, e], [0, self.multipdf([float(e)])], color='C'+ str(p_color), alpha=0.25)

                self.statsfig.plot(e, self.multipdf([float(e)]), 'o', color='C'+ str(p_color), alpha=0.25)

            x = np.linspace(- 2 * m, 2 * m, 1000)

            self.statsfig.plot(x, self.multipdf(x),color = "C1")


            self.statsfig.plot([error, error], [0, self.multipdf([float(error)])], color='C'+ str(p_color), alpha=1)

            self.statsfig.plot(error, self.multipdf([float(error)]), 'o', color='C'+ str(p_color), alpha= 1)


            self.statsfig.set_xlim([x[0], x[-1]])

            self.statsfig.set_xlabel("error (ppm)")

            self.statsfig.set_ylabel("probability density")

            self.statscanvas.draw()

    def multipdf(self, x):

        y = np.zeros(len(x))

        if self.atom == 'H':

            for mean, std in zip(self.Hmeans, self.Hstdevs):
                y += norm.pdf(x, loc=mean, scale=std)

        elif self.atom == 'C':

            for mean, std in zip(self.Cmeans, self.Cstdevs):
                y += norm.pdf(x, loc=mean, scale=std)

        return y

    def findmeans(self):

        if ui.table_widget.Tab1.worker.settings.StatsParamFile == 'none':

            self.Hmeans = [0]

            self.Hstdevs = [0.18731058105269952]

            self.Cmeans = [0]

            self.Cstdevs = [2.269372270818724]

        else:

            self.Cmeans, self.Cstdevs, self.Hmeans, self.Hstdevs = ReadParamFile(
                ui.table_widget.Tab1.worker.settings.StatsParamFile, 'g')

    def populatetable(self):

        # which isomer is selected

        self.isomerindex = int(str(ui.table_widget.Tab4.IsomerSelect.currentText())[-1]) - 1

        self.errortable.setRowCount(0)

        if self.atom == 'H':

            self.labels, self.shifts, self.exps, self.scaleds, self.errors = self.dp4data.Hlabels[self.isomerindex], \
                                                                             self.dp4data.Hshifts[self.isomerindex], \
                                                                             self.dp4data.Hexp[self.isomerindex], \
                                                                             self.dp4data.Hscaled[self.isomerindex], \
                                                                             self.dp4data.Herrors[self.isomerindex]

        elif self.atom == 'C':

            self.labels, self.shifts, self.exps, self.scaleds, self.errors = self.dp4data.Clabels[self.isomerindex], \
                                                                             self.dp4data.Cshifts[self.isomerindex], \
                                                                             self.dp4data.Cexp[self.isomerindex], \
                                                                             self.dp4data.Cscaled[self.isomerindex], \
                                                                             self.dp4data.Cerrors[self.isomerindex]

        self.errortable.setRowCount(len(self.labels))

        # set info in rows and columns

        c = 0

        for label, shift, exp, scaled, error in zip(self.labels, self.shifts, self.exps, self.scaleds, self.errors):
            self.errortable.setItem(c, 0, QtWidgets.QTableWidgetItem(label))
            self.errortable.setItem(c, 1, QtWidgets.QTableWidgetItem(str(round(shift, 2))))
            self.errortable.setItem(c, 2, QtWidgets.QTableWidgetItem(str(round(scaled, 2))))
            self.errortable.setItem(c, 3, QtWidgets.QTableWidgetItem(str(round(exp, 2))))
            self.errortable.setItem(c, 4, QtWidgets.QTableWidgetItem(str(round(error, 2))))

            c += 1

        self.errortable.selectRow(0)

        self.plot()

    def selectpoint(self, event):

        self.xpos = event.xdata

        self.ypos = event.ydata

        # find the closest point to the click

        coords = self.statsfig.transData.transform((self.xpos, self.ypos))

        coordinates = np.array(self.statsfig.transData.transform(list(zip(self.errors, self.multipdf(self.errors)))))

        mindis = np.argmin((coordinates[:, 0] - coords[0]) ** 2 + (coordinates[:, 1] - coords[1]) ** 2)

        self.errortable.selectRow(mindis)


class CalculationTab(QtWidgets.QWidget):

    signal_start_PyDP4 = QtCore.pyqtSignal()

    def __init__(self):
        self.log_file = ''
        super(CalculationTab, self).__init__()
        self.cwd = Path(os.getcwd())

        self.Overall_verticallayout =QtWidgets.QVBoxLayout(self)
        self.Overall_widget = QtWidgets.QWidget(self)

        self.title = QtWidgets.QLabel(self)
        self.title.setText("DP5 GUI")
        self.title.setObjectName("Title")
        self.Overall_verticallayout.addWidget(self.title)

        ########### input

        self.input_widget = QtWidgets.QWidget(self)
        self.input_layout = QtWidgets.QHBoxLayout(self)

        self.structure_widget = QtWidgets.QWidget(self)
        self.structure_layout = QtWidgets.QVBoxLayout(self)

        self.Add_structure = QtWidgets.QPushButton(self)
        self.Add_structure.setObjectName("Add_structure")
        self.Add_structure.setText("Add structure")

        self.remove_structure = QtWidgets.QPushButton(self)
        self.remove_structure.setObjectName("remove_structure")
        self.remove_structure.setText("Remove selected")

        self.structure_layout.addWidget(self.Add_structure)
        self.structure_layout.addWidget(self.remove_structure)

        self.structure_widget.setLayout(self.structure_layout)

        self.input_layout.addWidget(self.structure_widget)

        self.Structure_list = QtWidgets.QListWidget(self)
        self.Structure_list.setMaximumHeight(self.Add_structure.sizeHint().height() * 2)

        self.Structure_list.setObjectName("Structure_list")
        self.input_layout.addWidget(self.Structure_list)

        self.NMR_widget = QtWidgets.QWidget(self)
        self.NMR_layout = QtWidgets.QVBoxLayout(self)

        self.Add_NMR = QtWidgets.QPushButton(self)
        self.Add_NMR.setObjectName("Add_NMR")
        self.Add_NMR.setText("Add raw NMR")

        self.Add_NMR_desc = QtWidgets.QPushButton(self)
        self.Add_NMR_desc.setObjectName("Add_NMR")
        self.Add_NMR_desc.setText("Add NMR desc")

        self.remove_NMR = QtWidgets.QPushButton(self)
        self.remove_NMR.setObjectName("remove_NMR")
        self.remove_NMR.setText("Remove selected")

        self.NMR_layout.addWidget(self.Add_NMR)
        self.NMR_layout.addWidget(self.Add_NMR_desc)
        self.NMR_layout.addWidget(self.remove_NMR)

        self.NMR_widget.setLayout(self.NMR_layout)
        self.input_layout.addWidget(self.NMR_widget)

        self.NMR_list = QtWidgets.QListWidget(self)
        self.NMR_list.setObjectName("NMR_list")
        self.NMR_list.setMaximumHeight(self.Add_structure.sizeHint().height() * 3)
        self.input_layout.addWidget(self.NMR_list)

        self.Output_add = QtWidgets.QPushButton(self)
        self.Output_add.setObjectName("Out add")
        self.Output_add.setText("Output Folder")
        self.input_layout.addWidget(self.Output_add)

        self.Output_list = QtWidgets.QListWidget(self)
        self.Output_list.setObjectName("Out list")
        self.Output_list.setMaximumHeight(self.Add_structure.sizeHint().height() * 2)
        self.input_layout.addWidget(self.Output_list)

        self.input_widget.setLayout(self.input_layout)
        self.Overall_verticallayout.addWidget(self.input_widget)

        ############ workflow

        self.workflow_layout = QtWidgets.QGridLayout(self)
        self.workflow_widget = QtWidgets.QWidget(self)

        self.workflow_title = QtWidgets.QLabel(self)
        self.workflow_title.setObjectName("workflow_title")
        self.workflow_title.setText("Workflow")
        self.workflow_layout.addWidget(self.workflow_title, 0, 0)

        self.CleanUp_yn = QtWidgets.QCheckBox(self)
        self.CleanUp_yn.setObjectName("CleanUp_yn")
        self.CleanUp_yn.setText("Clean\n"
                                          "Structures")

        self.workflow_layout.addWidget(self.CleanUp_yn,1,0)

        self.Gen_diastereomers_yn = QtWidgets.QCheckBox(self)
        self.Gen_diastereomers_yn.setObjectName("Gen_diastereomers_yn")
        self.Gen_diastereomers_yn.setText("Generate\n"
                                          "Diastereomers")

        self.workflow_layout.addWidget(self.Gen_diastereomers_yn,1,1)

        self.Solvent_yn = QtWidgets.QCheckBox(self)
        self.Solvent_yn.setObjectName("Solvent_yn")
        self.Solvent_yn.setText("Solvent")
        self.workflow_layout.addWidget(self.Solvent_yn,1,2)

        self.solvent_drop = QtWidgets.QComboBox(self)
        self.solvent_drop.setObjectName("solvent_drop")
        solvents = ['chloroform', 'dimethylsulfoxide', 'benzene', 'methanol', 'pyridine', 'acetone']
        self.solvent_drop.addItems(solvents)

        self.workflow_layout.addWidget(self.solvent_drop,2,2)

        self.MM_yn = QtWidgets.QCheckBox(self)
        self.MM_yn.setObjectName("MM_yn")
        self.MM_yn.setText("Molecular\nMechanics")
        self.workflow_layout.addWidget(self.MM_yn,1,3)

        self.MM_drop = QtWidgets.QComboBox(self)
        self.MM_drop.setObjectName("MM_drop")
        MM = ["MacroModel","Tinker"]
        self.MM_drop.addItems(MM)
        self.workflow_layout.addWidget(self.MM_drop,2,3)

        self.MM_advanced = QtWidgets.QPushButton(self)
        self.MM_advanced.setText("Advanced Settings")
        self.workflow_layout.addWidget(self.MM_advanced,3,3)

        self.DFT_yn = QtWidgets.QCheckBox(self)
        self.DFT_yn.setObjectName("DFT_yn")
        self.DFT_yn.setText("DFT\nCalculations")
        self.workflow_layout.addWidget(self.DFT_yn, 1, 4)

        self.DFT_drop = QtWidgets.QComboBox(self)
        self.DFT_drop.setObjectName("DFT_drop")
        DFT = ["Gaussian", "NWChem"]
        self.DFT_drop.addItems(DFT)
        self.workflow_layout.addWidget(self.DFT_drop, 2, 4)

        self.DFT_advanced = QtWidgets.QPushButton(self)
        self.DFT_advanced.setText("Advanced Settings")
        self.workflow_layout.addWidget(self.DFT_advanced,3,4)

        self.Assignment_yn = QtWidgets.QCheckBox(self)
        self.Assignment_yn.setObjectName("Assignment_yn")
        self.Assignment_yn.setText("NMR\nAssignment")
        self.workflow_layout.addWidget(self.Assignment_yn, 1, 5)

        self.DP5_stat_yn = QtWidgets.QCheckBox(self)
        self.DP5_stat_yn.setObjectName("DP5_stat_yn")
        self.DP5_stat_yn.setText("DP5 Statistics")
        self.workflow_layout.addWidget(self.DP5_stat_yn, 1, 6)

        self.DP4_stat_yn = QtWidgets.QCheckBox(self)
        self.DP4_stat_yn.setObjectName("DP4_stat_yn")
        self.DP4_stat_yn.setText("DP4 Statistics")
        self.workflow_layout.addWidget(self.DP4_stat_yn, 1, 7)

        self.Add_stats_model = QtWidgets.QPushButton(self)
        self.Add_stats_model.setObjectName("Add_stats_model")
        self.Add_stats_model.setText("Add stats model")
        self.workflow_layout.addWidget(self.Add_stats_model, 2, 7)

        self.Stats_list = QtWidgets.QListWidget(self)
        self.Stats_list.setMaximumHeight(self.Add_structure.sizeHint().height())
        self.Stats_list.setObjectName("Stats_list")

        self.workflow_layout.addWidget(self.Stats_list, 3, 7)
        self.workflow_widget.setLayout(self.workflow_layout)
        self.Overall_verticallayout.addWidget(self.workflow_widget)

        ############ calc

        self.Gobutton = QtWidgets.QPushButton(self)
        self.Gobutton.setObjectName("Gobutton")
        self.Gobutton.setText("Calculate")
        self.Overall_verticallayout.addWidget(self.Gobutton)

        self.Output_box = QtWidgets.QTextEdit(self)
        self.Output_box.setObjectName("Output_box")
        self.Overall_verticallayout.addWidget(self.Output_box)

        self.Overall_widget.setLayout(self.Overall_verticallayout)

        #################################################################################################buttons methods

        # adding structures to the gui

        self.NMR_paths = []

        self.Structure_paths = []

        self.Output_folder = Path.cwd()

        self.Add_structure.clicked.connect(self.addstructure)

        self.remove_structure.clicked.connect(self.removestructure)

        self.Add_NMR.clicked.connect(self.addNMR)

        self.Add_NMR_desc.clicked.connect(self.remove_all_NMR)

        self.Add_NMR_desc.clicked.connect(self.addNMRdesc)

        self.remove_NMR.clicked.connect(self.removeNMR)

        self.Output_add.clicked.connect(self.addoutputfolder)

        # selecting solvent check box

        self.solvent_drop.setEnabled(False)

        self.Solvent_yn.stateChanged.connect(self.solventtoggle)

        self.MM_drop.setEnabled(False)
        self.MM_advanced.setEnabled(False)

        self.MM_yn.stateChanged.connect(self.MMtoggle)

        self.DFT_drop.setEnabled(False)
        self.DFT_advanced.setEnabled(False)

        self.DFT_yn.stateChanged.connect(self.DFTtoggle)

        self.Assignment_yn.stateChanged.connect(self.Assignment_toggle)

        # adding a stats model

        self.Add_stats_model.setEnabled(False)

        self.DP4_stat_yn.stateChanged.connect(self.DP4toggle)

        self.DP5_stat_yn.stateChanged.connect(self.DP5toggle)

        self.Add_stats_model.clicked.connect(self.addstats)

        #####################################################################################running PyDP4 in new thread

        # self.Gobutton.clicked.connect(self.Gotoggle)

        self.Gobutton.clicked.connect(self.get_current_values)

        # make worker object

        self.worker = PyDP4WorkerObject()

        # make thread

        self.thread = QtCore.QThread()

        # move the worker object to the ne thread

        self.worker.moveToThread(self.thread)

        # connect workers finished signal to quit thread

        self.worker.finished.connect(self.thread.quit)

        # connect workers finished signal to enable go button

        self.worker.finished.connect(self.Gotoggle)

        self.worker.finished.connect(self.enabletabs)

        # connect start background signal to background job slot

        self.Gobutton.clicked.connect(self.start_PyDP4)

        # start pydp4 give start signal

        # go button press disables go button

        self.Gobutton.clicked.connect(self.Gotoggle)

        # start signal starts pydp4

        self.signal_start_PyDP4.connect(self.worker.runPyDP4)

        self.DFT_settings = DFT_advanced_settings()

        self.DFT_advanced.clicked.connect(self.DFT_pop)

        self.MM_settings = MM_advanced_settings()

        self.MM_advanced.clicked.connect(self.MM_pop)

    def DFT_pop(self):

        self.DFT_settings.show()


    def MM_pop(self):

        self.MM_settings.show()

    def enabletabs(self):

        ui.table_widget.addplottabs()

    def addoutputfolder(self):

        # filename = QtWidgets.QFileDialog.getOpenFileName()
        filename = QtWidgets.QFileDialog.getExistingDirectory()

        # self.NMR_list.addItem(str(filename[0].split("/")[-1]))
        self.Output_list.clear()
        self.Output_list.addItem(filename)
        self.Output_folder = Path(filename)

    def get_current_values(self):

        import PyDP4
        import datetime
        import getpass

        self.settings = PyDP4.Settings()

        self.settings.GUIRunning = True

        now = datetime.datetime.now()
        self.settings.StartTime = now.strftime('%d%b%H%M')

        self.settings.user = getpass.getuser()
        self.settings.DarwinScrDir.replace('/u/', self.settings.user)

        # Read config file and fill in settings in from that
        self.settings = PyDP4.ReadConfig(self.settings)

        # add output folder

        self.settings.OutputFolder = self.Output_folder

        self.log_file = open(self.Output_folder / "DP4_log.log", "w+")

        self.settings.InputFilesPaths = self.Structure_paths

        # copy structures to output folder

        for f in self.Structure_paths:

            if not Path(self.Output_folder / f).exists():
                shutil.copyfile(f, self.settings.OutputFolder / f)

        # add structures

        Smiles = []
        Smarts = []
        InChIs = []

        for index in range(self.Structure_list.count()):

            list_text = self.Structure_list.item(index).text()

            if list_text != '':

                # check if file is text or sdf

                if list_text.endswith('.sdf'):

                    self.settings.InputFiles.append(list_text[:-4])

                # else its a text file - work out if it contains smiles smarts or inchis

                elif (list_text.endswith('Smiles')) or (list_text.endswith('smiles')):

                    Smiles.append(list_text)

                elif (list_text.endswith('Smarts')) or (list_text.endswith('smarts')):

                    Smarts.append(list_text)

                elif (list_text.endswith('InChI')) or (list_text.endswith('InChi')) or (list_text.endswith('inchi')):

                    InChIs.append(list_text)

                else:

                    print(
                        "file type not recognised to use smiles, smarts or inchi input please use .smiles, .smarts or .inchi file extension respectively")

            if len(Smiles) == 1:
                self.settings.Smiles = Smiles[0]

            elif len(Smiles) > 1:

                # if the user has added more than one Smiles string make a new file and concatenate them

                SmilesFile = open(self.settings.OutputFolder / "Smiles_Input.smiles", "w+")

                AllSmiles = []

                for f in Smiles:

                    for line in open(self.settings.OutputFolder / f).readlines():
                        AllSmiles.append(line.strip())

                for s in AllSmiles[:-1]:
                    SmilesFile.write(s + "\n")

                SmilesFile.write(AllSmiles[-1])

                SmilesFile.close()

                self.settings.Smiles = "Smiles_Input.smiles"

            if len(Smarts) == 1:
                self.settings.Smarts = Smarts[0]

            elif len(Smarts) > 1:

                # if the user has added more than one Smiles string make a new file and concatenate them

                SmartsFile = open(self.settings.OutputFolder / "Smarts_Input.smarts", "w+")

                AllSmarts = []

                for f in Smarts:

                    for line in open(self.settings.OutputFolder / f).readlines():
                        AllSmarts.append(line.strip())

                for s in AllSmarts[:-1]:
                    SmartsFile.write(s + "\n")

                SmartsFile.write(AllSmarts[-1])

                SmartsFile.close()

                self.settings.Smarts = "Smarts_Input.smarts"

            if len(InChIs) == 1:
                self.settings.InChIs = InChIs[0]

            elif len(InChIs) > 1:

                # if the user has added more than one Smiles string make a new file and concatenate them

                InchIsFile = open(self.settings.OutputFolder / "InchIs_Input.inchi", "w+")

                AllInchIs = []

                for f in InChIs:

                    for line in open(self.settings.OutputFolder / f).readlines():
                        AllInchIs.append(line.strip())

                for s in AllInchIs[:-1]:
                    InchIsFile.write(s + "\n")

                InchIsFile.write(AllInchIs[-1])

                InchIsFile.close()

                self.settings.InChIs = "InchIs_Input.smarts"


        for f in self.Structure_paths:

            if not Path(self.Output_folder / f).exists():
                shutil.copyfile(f, self.settings.OutputFolder / f)

        # add NMR

        self.settings.NMRsource = self.NMR_paths

        # get workflow information

        # solvent

        if self.Solvent_yn.isChecked() == 1:
            self.settings.Solvent = self.solvent_drop.currentText()

        # generate diastereomers

        self.settings.Workflow = ''

        if self.CleanUp_yn.isChecked() == 1:
            self.settings.Workflow += 'c'

        if self.Gen_diastereomers_yn.isChecked() == 1:
            self.settings.Workflow += 'g'

        # molecular mechanics

        if self.MM_yn.isChecked() == 1:

            self.settings.Workflow += 'm'

            if self.MM_drop.currentText() == "MacroModel":

                self.settings.MM = 'm'

            else:

                self.settings.MM = 't'

        if self.DFT_yn.isChecked() == 1:

            if self.DFT_drop.currentText() == "Gaussian":
                self.settings.DFT = 'g'

            else:
                self.settings.DFT = 'n'

            # Split single point

            if self.DFT_settings.Energy_yn.isChecked() == 1:
                self.settings.Workflow += 'e'
                self.settings.eBasisSet = self.DFT_settings.Energy_basis_drop.currentText()
                self.settings.eFunctional = self.DFT_settings.Energy_functional_drop.currentText()

            if self.DFT_settings.DFTGeom_yn.isChecked() == 1:
                self.settings.Workflow += 'o'
                self.settings.oBasisSet = self.DFT_settings.DFT_geom_basis_drop.currentText()
                self.settings.oFunctional = self.DFT_settings.DFT_geom_functional_drop.currentText()

            if self.DFT_settings.NMR_calc_yn.isChecked() == 1:
                self.settings.Workflow += 'n'
                self.settings.nBasisSet = self.DFT_settings.NMR_basis_drop.currentText()
                self.settings.nFunctional = self.DFT_settings.NMR_functional_drop.currentText()

        if self.DP4_stat_yn.isChecked():

            self.settings.Workflow += 's'

            if self.Stats_list.item(0) != None:
                self.settings.StatsParamFile = self.Stats_list.item(0).text()
                self.settings.StatsModel = 'm'

        if self.DP5_stat_yn.isChecked():

            self.settings.Workflow += 'w'

        elif self.Assignment_yn.isChecked():

            self.settings.Workflow += 'a'

        self.settings.ScriptDir = os.path.dirname(os.path.realpath(sys.argv[0]))

    def start_PyDP4(self):

        # No harm in calling thread.start() after the thread is already started.

        # start the thread

        self.thread.start()

        self.thread.finished.connect(self.Gotoggle)

        # send background job start signal

        self.signal_start_PyDP4.emit()

        self.thread.disconnect()

    def addstructure(self):

        f = QtWidgets.QFileDialog.getOpenFileName()[0]

        if f:
            filename = Path(f)

            self.Structure_list.addItem(filename.name)
            self.Structure_paths.append(filename)

    def removestructure(self):

        item = self.Structure_list.selectedItems()
        if not item:
            return
        for i in item:
            self.Structure_list.takeItem(self.Structure_list.row(i))
            self.Structure_paths.pop(self.Structure_list.row(i))

    def addNMR(self):

        # filename = QtWidgets.QFileDialog.getOpenFileName()

        i = QtWidgets.QFileDialog.getExistingDirectory()

        if i:

            filename = Path(i)

            p_switch = 0

            c_switch = 0

            for f in filename.iterdir():

                if f.name == "Carbon" or f.name == "carbon" or f.name == "Carbon.dx" or f.name == "carbon.dx":
                    self.NMR_list.addItem(f.name)
                    self.NMR_paths.append(f)
                    c_switch = 1

                elif f.name == "Proton" or f.name == "proton" or f.name == "Proton.dx" or f.name == "proton.dx":
                    self.NMR_list.addItem(f.name)
                    self.NMR_paths.append(f)
                    p_switch = 1

                if p_switch == 1 and c_switch == 1:
                    break

            # self.NMR_list.addItem(str(filename[0].split("/")[-1]))

            if p_switch == 0 and c_switch == 0:
                self.NMR_list.addItem(filename.name)
                self.NMR_paths.append(filename.name)

            self.Add_NMR_desc.setEnabled(False)

    def addNMRdesc(self):

        f = QtWidgets.QFileDialog.getOpenFileName()[0]

        if f:
            filename = Path(f)

            self.NMR_list.addItem(filename.name)
            self.NMR_paths.append(filename)

            self.Add_NMR.setEnabled(False)
            self.Add_NMR_desc.setEnabled(False)


    def remove_all_NMR(self):

        self.NMR_list.clear()
        self.NMR_paths = []
        self.Add_NMR.setEnabled(True)
        self.Add_NMR_desc.setEnabled(True)

    def removeNMR(self):

        item = self.NMR_list.selectedItems()
        if not item:
            return
        for i in item:
            self.NMR_list.takeItem(self.NMR_list.row(i))
            self.NMR_paths.pop(self.NMR_list.row(i))

        if len(self.NMR_paths) == 0:
            self.Add_NMR.setEnabled(True)
            self.Add_NMR_desc.setEnabled(True)


    def addstats(self):
        filename = QtWidgets.QFileDialog.getOpenFileName()
        self.Stats_list.addItem(str(filename[0]))

    def solventtoggle(self, state):

        if state > 0:
            self.solvent_drop.setEnabled(True)
            #self.Gen_diastereomers_yn.setChecked(True)

        else:
            self.solvent_drop.setEnabled(False)
            self.MM_yn.setChecked(False)

    def MMtoggle(self, state):

        if state > 0:

            #self.Gen_diastereomers_yn.setChecked(True)

            self.Solvent_yn.setChecked(True)
            self.solvent_drop.setEnabled(True)
            self.MM_drop.setEnabled(True)
            self.MM_advanced.setEnabled(True)


        else:
            self.solvent_drop.setEnabled(False)
            self.Solvent_yn.setChecked(False)
            self.MM_drop.setEnabled(False)
            self.MM_advanced.setEnabled(False)
            self.DFT_yn.setChecked(False)

    def DFTtoggle(self,state):

        if state > 0:

            # self.Gen_diastereomers_yn.setChecked(True)

            self.Solvent_yn.setChecked(True)
            self.solvent_drop.setEnabled(True)
            self.MM_yn.setChecked(True)
            self.MM_drop.setEnabled(True)
            self.DFT_drop.setEnabled(True)
            self.DFT_advanced.setEnabled(True)
            self.DFT_settings.NMR_calc_yn.setChecked(True)

        else:
            self.DFT_advanced.setEnabled(False)
            self.DFT_drop.setEnabled(False)
            self.solvent_drop.setEnabled(False)
            self.Solvent_yn.setChecked(False)
            self.MM_yn.setChecked(False)
            self.MM_drop.setEnabled(False)
            self.Assignment_yn.setChecked(False)
            self.DP4_stat_yn.setChecked(False)

            self.DFT_settings.DFTGeom_yn.setChecked(False)
            self.DFT_settings.Energy_yn.setChecked(False)
            self.DFT_settings.NMR_calc_yn.setChecked(False)

    def Assignment_toggle(self, state):

        if state > 0:

            self.DFT_yn.setChecked(True)
            self.Solvent_yn.setChecked(True)
            self.solvent_drop.setEnabled(True)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(True)

        else:
            self.DP4_stat_yn.setChecked(False)
            self.Add_stats_model.setEnabled(False)

            self.DFT_yn.setChecked(False)
            self.Solvent_yn.setChecked(False)
            self.solvent_drop.setEnabled(False)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(False)

    def DP4toggle(self, state):

        if state > 0:
            self.Add_stats_model.setEnabled(True)

            self.DFT_yn.setChecked(True)
            self.Solvent_yn.setChecked(True)
            self.solvent_drop.setEnabled(True)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(True)
            self.Assignment_yn.setChecked(True)

        else:
            self.Add_stats_model.setEnabled(False)
            self.DFT_yn.setChecked(False)
            self.Solvent_yn.setChecked(False)
            self.solvent_drop.setEnabled(False)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(False)
            self.Assignment_yn.setChecked(False)

    def DP5toggle(self, state):

        if state > 0:
            self.DFT_yn.setChecked(True)
            self.Solvent_yn.setChecked(True)
            self.solvent_drop.setEnabled(True)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(True)
            self.Assignment_yn.setChecked(True)

        else:
            self.DFT_yn.setChecked(False)
            self.Solvent_yn.setChecked(False)
            self.solvent_drop.setEnabled(False)
            #self.Gen_diastereomers_yn.setChecked(True)
            self.MM_yn.setChecked(False)
            self.Assignment_yn.setChecked(False)

    def append_text(self, text):
        self.Output_box.moveCursor(QtGui.QTextCursor.End)
        self.Output_box.insertPlainText(text)
        # self.log_file.write(text)

        self.log_file.write(text)

    def Gotoggle(self):

        if self.Gobutton.isEnabled() == True:

            self.Gobutton.setEnabled(False)

        else:

            self.Gobutton.setEnabled(True)


class DFT_advanced_settings(QtWidgets.QWidget):

    def __init__(self):
        super(DFT_advanced_settings, self).__init__()
        self.layout = QtWidgets.QVBoxLayout(self)

        self.title = QtWidgets.QLabel(self)
        self.title.setText("DFT Advanced Options")
        self.layout.addWidget(self.title)

        ##################################### DFT Geom Optimisation

        self.DFTGeom_yn = QtWidgets.QCheckBox(self)
        self.DFTGeom_yn.setText("DFT Geometry Optimisation")
        self.layout.addWidget(self.DFTGeom_yn)

        self.DFT_geom_functional_drop_label = QtWidgets.QLabel(self)
        self.DFT_geom_functional_drop_label.setText("Functional")
        self.layout.addWidget(self.DFT_geom_functional_drop_label)

        self.DFT_geom_functional_drop = QtWidgets.QComboBox(self)
        DFTopt_functional = ['B3LYP', 'm062x', 'mPW1PW91']
        self.DFT_geom_functional_drop.addItems(DFTopt_functional)
        self.layout.addWidget(self.DFT_geom_functional_drop)

        self.DFT_geom_basis_drop_label = QtWidgets.QLabel(self)
        self.DFT_geom_basis_drop_label.setText("Basis Set")
        self.layout.addWidget(self.DFT_geom_basis_drop_label)

        self.DFT_geom_basis_drop = QtWidgets.QComboBox(self)
        DFTopt_basis = ['6-31g(d)', '6-311g(d)', 'def2svp', 'def2tzvp']
        self.DFT_geom_basis_drop.addItems(DFTopt_basis)
        self.layout.addWidget(self.DFT_geom_basis_drop)

        ##################################### DFT Energy

        self.Energy_yn = QtWidgets.QCheckBox(self)
        self.Energy_yn.setText("Split Single Point DFT")
        self.layout.addWidget(self.Energy_yn)

        self.Energy_functional_drop_label = QtWidgets.QLabel(self)
        self.Energy_functional_drop_label.setText("Functional")
        self.layout.addWidget(self.Energy_functional_drop_label)

        self.Energy_functional_drop = QtWidgets.QComboBox(self)
        Energy_functional = ['m062x', 'mPW1PW91', 'B3LYP']
        self.Energy_functional_drop.addItems(Energy_functional)
        self.layout.addWidget(self.Energy_functional_drop)

        self.Energy_basis_drop_label = QtWidgets.QLabel(self)
        self.Energy_basis_drop_label.setText("Basis Set")
        self.layout.addWidget(self.Energy_basis_drop_label)

        self.Energy_basis_drop = QtWidgets.QComboBox(self)
        Energy_basis = ['6-31g(d)', '6-311g(d)', 'def2svp', 'def2tzvp']
        self.Energy_basis_drop.addItems(Energy_basis)
        self.layout.addWidget(self.Energy_basis_drop)

        ##################################### DFT NMR

        self.NMR_calc_yn = QtWidgets.QCheckBox(self)
        self.NMR_calc_yn.setText("NMR Calculations")
        self.layout.addWidget(self.NMR_calc_yn)

        self.NMR_functional_drop_label = QtWidgets.QLabel(self)
        self.NMR_functional_drop_label.setText("Functional")
        self.layout.addWidget(self.NMR_functional_drop_label)

        self.NMR_functional_drop = QtWidgets.QComboBox(self)
        NMR_functional = ['mPW1PW91', 'B3LYP', 'M062X']
        self.NMR_functional_drop.addItems(NMR_functional)
        self.layout.addWidget(self.NMR_functional_drop)

        self.NMR_basis_drop_label = QtWidgets.QLabel(self)
        self.NMR_basis_drop_label.setText("Basis Set")
        self.layout.addWidget(self.NMR_basis_drop_label)

        self.NMR_basis_drop = QtWidgets.QComboBox(self)
        NMR_basis = ['6-311g(d)', '6-31g(d)', 'def2svp', 'def2tzvp']
        self.NMR_basis_drop.addItems(NMR_basis)
        self.layout.addWidget(self.NMR_basis_drop)

        ###################################### methods

        # selecting DFT geometry opt box

        self.DFT_geom_functional_drop.setEnabled(False)

        self.DFT_geom_basis_drop.setEnabled(False)

        self.DFTGeom_yn.stateChanged.connect(self.DFTopttoggle)

        # selecting Split Single point box

        self.Energy_functional_drop.setEnabled(False)

        self.Energy_basis_drop.setEnabled(False)

        self.Energy_yn.stateChanged.connect(self.Energytoggle)

        # selecting NMR box

        self.NMR_functional_drop.setEnabled(False)

        self.NMR_basis_drop.setEnabled(False)

        self.NMR_calc_yn.stateChanged.connect(self.NMRtoggle)


    def Energytoggle(self, state):

        if state > 0:
            self.Energy_functional_drop.setEnabled(True)
            self.Energy_basis_drop.setEnabled(True)

        else:
            self.Energy_functional_drop.setEnabled(False)
            self.Energy_basis_drop.setEnabled(False)

    def DFTopttoggle(self, state):

        if state > 0:
            self.DFT_geom_functional_drop.setEnabled(True)
            self.DFT_geom_basis_drop.setEnabled(True)

        else:
            self.DFT_geom_functional_drop.setEnabled(False)
            self.DFT_geom_basis_drop.setEnabled(False)

    def NMRtoggle(self, state):

        if state > 0:
            self.NMR_functional_drop.setEnabled(True)
            self.NMR_basis_drop.setEnabled(True)

        else:
            self.NMR_functional_drop.setEnabled(False)
            self.NMR_basis_drop.setEnabled(False)


class MM_advanced_settings(QtWidgets.QWidget):

    def __init__(self):
        super(MM_advanced_settings, self).__init__()

        self.layout = QtWidgets.QVBoxLayout(self)

        self.title = QtWidgets.QLabel(self)
        self.title.setText("MM Advanced Options")
        self.layout.addWidget(self.title)


class ProtonPlotTab(QtWidgets.QWidget):

    def __init__(self):

        super(ProtonPlotTab, self).__init__()

        self.cwd = Path(os.getcwd())

        # create a figure instance

        self.figure = Figure()

        # create the canvas widget thats displays figure

        self.canvas = FigureCanvas(self.figure)

        # make the navigation widget - this takes the canvas widget as the parent

        self.toolbar = NavigationToolbar(self.canvas, self)

        self.IsomerSelect = QtWidgets.QComboBox(self)

        self.IsomerSelect.addItems(self.Isomer_number())

        self.IsomerSelect.currentIndexChanged.connect(self.PlotProton)
        self.image = QSvgWidget()

        self.widget1 = QtWidgets.QWidget(self)

        self.widget2 = QtWidgets.QWidget(self)

        self.hl1 = QtWidgets.QVBoxLayout(self.widget1)

        self.hl2 = QtWidgets.QVBoxLayout(self.widget2)

        self.hl1.addWidget(self.toolbar)

        self.hl1.addWidget(self.canvas)

        self.hl2.addWidget(self.IsomerSelect)

        self.hl2.addWidget(self.image)

        self.widget1.setGeometry(QtCore.QRect(300, 0, 1200, 875))

        # self.IsomerSelect.setGeometry(QtCore.QRect(0,0,100,50))

        # self.image.setGeometry(QtCore.QRect(0,50,300,300))

        self.widget2.setGeometry(QtCore.QRect(0, 0, 300, 400))

        self.canvas.mpl_connect('button_press_event', self.selectpoint)

        ################################################################################################################

        if ui.table_widget.Tab1.settings.OutputFolder == '':

            pdir = self.cwd / "Pickles"

        else:
            pdir = ui.table_widget.Tab1.settings.OutputFolder / "Pickles"

        if os.path.isfile(pdir / ui.table_widget.Tab1.settings.InputFiles[0] / "protondata"):
            self.protondata = pickle.load(
                Path(pdir / ui.table_widget.Tab1.settings.InputFiles[0] / "protondata").open(mode="rb"))

        self.xdata = self.protondata["xdata"]

        self.ydata = self.protondata["ydata"]

        self.centres = self.protondata["centres"]

        self.exp_peaks = self.protondata["exppeaks"]

        self.peak_regions = self.protondata["peakregions"]

        self.cummulative_vectors = self.protondata["cummulativevectors"]

        self.integral_sum = self.protondata["integralsum"]

        self.integrals = self.protondata["integrals"]

        self.sim_regions = self.protondata["sim_regions"]

        self.PlotProton()

        ################################################################################################################

    def Isomer_number(self):

        Isomer_list = []

        for c, i in enumerate(ui.table_widget.Tab1.worker.Isomers):
            Isomer_list.append("Isomer " + str(c + 1))

        return Isomer_list

    def PlotProton(self):

        self.RenderImage([], None)

        self.isomerindex = int(str(self.IsomerSelect.currentText())[-1]) -1

        self.isomer = ui.table_widget.Tab1.worker.Isomers[self.isomerindex]

        self.assigned_shifts = self.isomer.Hshifts

        self.assigned_peaks = []

        for peak in self.isomer.Hexp:

            if peak != '':
                self.assigned_peaks.append(peak)

        self.assigned_labels = self.isomer.Hlabels

        # check if pickle and DP4output files are in place

        if ui.table_widget.Tab1.worker.settings.OutputFolder == '':

            pdir = Path.cwd() / "Pickles"

        else:

            pdir = ui.table_widget.Tab1.worker.settings.OutputFolder / "Pickles"

        if Path(pdir / ui.table_widget.Tab1.worker.settings.InputFiles[0] / "protondata").exists():

            self.figure.clear()

            fig = self.figure.add_subplot(111)

            self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.07, top=0.95, wspace=0.05, hspace=0.05)

            #################################### will probs need to fix sorting here

            fig.set_xlim([10, 0])

            fig.set_xlabel("ppm")

            fig.plot(self.xdata, self.ydata, label='data', color='grey')

            set_exp = sorted(list(set(self.exp_peaks)))[::-1]

            # simulate_spectrum(xdata, assigned_shifts, assigned_peaks, set_exp)

            ##############################

            for ind, shift in enumerate(self.assigned_shifts):
                exp_p = self.assigned_peaks[ind]

                ind2 = set_exp.index(exp_p)
                y = lorentzian(self.xdata, 0.001, shift, 0.2)

                fig.plot(self.xdata, y + 1.05, color='C' + str(ind2 % 10))

            ##############################

            fig.axhline(1.05, color='grey')

            # plt integral information

            prev = 15

            count = 0

            for index in range(0, len(self.peak_regions)):

                if abs(prev - self.xdata[self.centres[index]]) < 0.45:
                    count += 1
                else:
                    count = 0
                    prev = self.xdata[self.centres[index]]

                fig.annotate(str(self.integrals[index]) + ' Hs',
                             xy=(self.xdata[self.centres[index]], -(0.1) - 0.1 * count), color='C' + str(index % 10))

                fig.annotate(str(round(self.xdata[self.centres[index]], 3)) + ' ppm',
                             xy=(self.xdata[self.centres[index]], -(0.15) - 0.1 * count), color='C' + str(index % 10))

                fig.plot(self.xdata[self.peak_regions[index]],
                         self.cummulative_vectors[index] + self.integral_sum[index],
                         color='C' + str(index % 10),
                         linewidth=2)

            for index in range(0, len(self.peak_regions) - 1):
                fig.plot([self.xdata[self.peak_regions[index][-1]], self.xdata[self.peak_regions[index + 1][0]]],
                         [self.integral_sum[index + 1], self.integral_sum[index + 1]], color='grey')

            for index, region in enumerate(self.peak_regions):
                fig.plot(self.xdata[region], self.sim_regions[index], color='C' + str(index % 10))

            # plt.legend()

            ### plotting assignment

            fig.set_yticks([], [])

            # fig.set_title(str(ui.table_widget.Tab1.settings.InputFiles[0]) +
            #             "\nProton NMR of Isomer " + str(isomerindex + 1) + "\n Number of Peaks Found = " + str(
            #  len(exp_peaks)))

            # plot assignments

            for ind1, peak in enumerate(self.assigned_peaks):
                fig.plot([peak, self.assigned_shifts[ind1]],
                         [1, 1.05], linewidth=0.5, color='cyan')

            # annotate peak locations

            for x, txt in enumerate(self.exp_peaks):

                if self.exp_peaks[x] in self.assigned_peaks:

                    color = 'C1'

                else:

                    color = 'grey'

                fig.plot(txt, -0.02, 'o', color=color)

            # annotate shift positions

            prev = 0

            count = 0

            s = np.argsort(np.array(self.assigned_shifts))

            s_assigned_shifts = np.array(self.assigned_shifts)[s]

            s_assigned_labels = np.array(self.assigned_labels)[s]

            s_assigned_peaks = np.array(self.assigned_peaks)[s]

            for x, txt in enumerate(s_assigned_labels[::-1]):

                w = np.where(set_exp == s_assigned_peaks[::-1][x])[0][0]

                color = w % 10

                if abs(prev - s_assigned_shifts[::-1][x]) < 0.2:
                    count += 1

                else:
                    count = 0
                    prev = s_assigned_shifts[::-1][x]

                fig.annotate(txt, (s_assigned_shifts[::-1][x], + 1.25 + 0.05 * count), color='C' + str(color))

            # ax1.plot(picked_peaks_ppm,ydata[picked_peaks],
            #        'co', label='Picked Peaks')

            fig.set_ylim([-0.5, 2.0])

            self.canvas.draw()

        else:
            pass

    def selectpoint(self, event):

        if event.dblclick:

            self.xpos = event.xdata

            self.ypos = event.ydata

            # find the cloest point to the click

            if self.xpos is None and self.xpos is None:

                self.PlotProton()

            elif self.ypos < 1:

                mindis = np.argmin(abs(self.xdata[self.centres] - self.xpos))

                self.PlotProtonSelected(mindis)

                # find which protons this peak has been assigned to

                p = np.where([round(i, 4) for i in self.assigned_peaks] == round(self.xdata[self.centres[mindis]], 4))[
                    0]

                la = [int(self.assigned_labels[j][1:]) - 1 for j in p]

                self.RenderImage(la, mindis)

            else:

                # find closest assigned shift

                m = np.argmin(abs(np.array(self.assigned_shifts) - self.xpos))

                # find which peak this is assigned to

                p = self.assigned_peaks[m]

                mindis = np.argmin(abs(p - self.xdata[self.centres]))

                p = np.where(self.assigned_peaks == round(self.xdata[self.centres[mindis]], 4))[0]

                la = [int(self.assigned_labels[j][1:]) - 1 for j in p]

                self.PlotProtonSelected(mindis)

                self.RenderImage(la, mindis)

    def PlotProtonSelected(self, mindis):

        # check if pickle and DP4output files are in place

        if ui.table_widget.Tab1.worker.settings.OutputFolder == '':

            pdir = Path.cwd() / "Pickles"

        else:

            pdir = ui.table_widget.Tab1.worker.settings.OutputFolder / "Pickles"

        if Path(pdir / ui.table_widget.Tab1.worker.settings.InputFiles[0] / "protondata").exists():

            self.figure.clear()

            fig = self.figure.add_subplot(111)

            self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.07, top=0.95, wspace=0.05, hspace=0.05)

            assigned_shifts = self.isomer.Hshifts

            assigned_peaks = []

            for peak in self.isomer.Hexp:

                if peak != '':
                    assigned_peaks.append(peak)

            assigned_labels = self.isomer.Hlabels

            #################################### will probs need to fix sorting here

            fig.set_xlim([10, 0])

            fig.set_xlabel("ppm")

            fig.plot(self.xdata, self.ydata, label='data', color='grey', alpha=0.5)

            set_exp = sorted(list(set(self.exp_peaks)))[::-1]

            # simulate_spectrum(xdata, assigned_shifts, assigned_peaks, set_exp)

            ##############################

            for ind, shift in enumerate(assigned_shifts):

                exp_p = assigned_peaks[ind]

                ind2 = set_exp.index(exp_p)

                y = lorentzian(self.xdata, 0.001, shift, 0.2)

                if ind2 == mindis:

                    fig.plot(self.xdata, y + 1.05, color='C' + str(ind2 % 10))

                else:

                    fig.plot(self.xdata, y + 1.05, color='grey', alpha=0.5)

            ##############################

            fig.axhline(1.05, color='grey')

            # plt integral information

            prev = 15

            count = 0

            for index, region in enumerate(self.peak_regions):
                fig.plot(self.xdata[region], self.sim_regions[index], color='grey')

            fig.plot(self.xdata[self.peak_regions[mindis]], self.sim_regions[mindis], color='C' + str(mindis))

            fig.annotate(str(self.integrals[mindis]) + ' Hs',
                         xy=(self.xdata[self.centres[mindis]], -(0.1) - 0.1 * count), color='C' + str(mindis))

            fig.annotate(str(round(self.xdata[self.centres[mindis]], 3)) + ' ppm',
                         xy=(self.xdata[self.centres[mindis]], -(0.15) - 0.1 * count), color='C' + str(mindis))

            fig.plot(self.xdata[self.peak_regions[mindis]],
                     self.cummulative_vectors[mindis] + self.integral_sum[mindis],
                     color='C' + str(mindis % 10),
                     linewidth=2)
            # plt.legend()

            ### plotting assignment

            fig.set_yticks([], [])

            # fig.set_title(str(ui.table_widget.Tab1.settings.InputFiles[0]) +
            #             "\nProton NMR of Isomer " + str(isomerindex + 1) + "\n Number of Peaks Found = " + str(
            #  len(exp_peaks)))

            # plot assignments

            # annotate shift positions

            prev = 0

            count = 0

            s = np.argsort(np.array(assigned_shifts))

            s_assigned_shifts = np.array(assigned_shifts)[s]

            s_assigned_labels = np.array(assigned_labels)[s]

            s_assigned_peaks = np.array(assigned_peaks)[s]

            for x, txt in enumerate(s_assigned_labels[::-1]):

                w = np.where(set_exp == s_assigned_peaks[::-1][x])[0][0]

                if w == mindis:

                    color = w % 10

                    if abs(prev - s_assigned_shifts[::-1][x]) < 0.2:
                        count += 1

                    else:
                        count = 0
                        prev = s_assigned_shifts[::-1][x]

                    fig.annotate(txt, (s_assigned_shifts[::-1][x], + 1.25 + 0.05 * count), color='C' + str(color))

            # ax1.plot(picked_peaks_ppm,ydata[picked_peaks],
            #        'co', label='Picked Peaks')

            fig.set_ylim([-0.5, 2.0])

            for x, txt in enumerate(self.exp_peaks):

                if self.exp_peaks[x] in self.assigned_peaks:

                    color = 'C1'

                else:

                    color = 'grey'

                fig.plot(txt, -0.02, 'o', color=color)

            self.canvas.draw()

        else:
            pass

    def RenderImage(self, atom, color):

        colors = [(0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17), (0.84, 0.15, 0.16), (0.58, 0.4, 0.74),
                  (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.5, 0.5, 0.5), (0.74, 0.74, 0.13), (0.09, 0.75, 0.81)]

        highlight = {}

        for i in atom:
            highlight[i] = colors[color % 10]

        isomerindex = int(str(self.IsomerSelect.currentText())[-1]) -1

        m = Chem.MolFromMolFile(str(ui.table_widget.Tab1.worker.settings.InputFilesPaths[isomerindex]).split('.sdf')[0] + '.sdf',removeHs=False)

        Chem.Compute2DCoords(m)

        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)

        drawer.DrawMolecule(m, highlightAtoms=atom, highlightAtomColors=highlight)

        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '')

        svg_bytes = bytearray(svg, encoding='utf-8')

        self.image.renderer().load(svg_bytes)

        self.image.setGeometry(QtCore.QRect(0, 50, 300, 300))

        ui.update()

        # f = open("f.svg", "w+")

        # f.write(str(svg))


class CarbonPlotTab(QtWidgets.QWidget):

    def __init__(self):

        super(CarbonPlotTab, self).__init__()

        # create a figure instance

        self.figure = Figure()

        # create the canvas widget thats displays figure

        self.canvas = FigureCanvas(self.figure)

        # make the navigation widget - this takes the canvas widget as the parent

        self.toolbar = NavigationToolbar(self.canvas, self)

        self.IsomerSelect = QtWidgets.QComboBox(self)

        self.IsomerSelect.addItems(self.Isomer_number())

        self.IsomerSelect.currentIndexChanged.connect(self.PlotCarbon)

        self.image = QSvgWidget()

        self.widget1 = QtWidgets.QWidget(self)

        self.widget2 = QtWidgets.QWidget(self)

        self.hl1 = QtWidgets.QVBoxLayout(self.widget1)

        self.hl2 = QtWidgets.QVBoxLayout(self.widget2)

        self.hl1.addWidget(self.toolbar)

        self.hl1.addWidget(self.canvas)

        self.hl2.addWidget(self.IsomerSelect)

        self.hl2.addWidget(self.image)

        self.widget1.setGeometry(QtCore.QRect(300, 0, 1200, 875))

        self.widget2.setGeometry(QtCore.QRect(0, 0, 300, 400))

        self.canvas.mpl_connect('button_press_event', self.selectpoint)

        #############################

        if ui.table_widget.Tab1.settings.OutputFolder == '':

            pdir = Path.cwd() / "Pickles"

        else:
            pdir = ui.table_widget.Tab1.settings.OutputFolder / "Pickles"

        if Path(pdir / ui.table_widget.Tab1.settings.InputFiles[0] / "carbondata").exists():
            self.carbondata = pickle.load(
                Path(pdir / ui.table_widget.Tab1.settings.InputFiles[0] / "carbondata").open(mode="rb"))

        self.xdata = self.carbondata["xdata"]

        self.ydata = self.carbondata["ydata"]

        self.exppeaks = self.carbondata["exppeaks"]

        self.simulated_ydata = self.carbondata["simulated_ydata"]

        self.removed = self.carbondata["removed"]

        self.PlotCarbon()

        #############################

    def Isomer_number(self):

        Isomer_list = []

        for c, i in enumerate(ui.table_widget.Tab1.worker.Isomers):
            Isomer_list.append("Isomer " + str(c + 1))

        return Isomer_list

    def PlotCarbon(self):

        self.RenderImage([])

        if ui.table_widget.Tab1.worker.settings.OutputFolder == '':

            pdir = Path.cwd() / "Pickles"

        else:

            pdir = ui.table_widget.Tab1.worker.settings.OutputFolder / "Pickles"

        if Path(pdir / ui.table_widget.Tab1.worker.settings.InputFiles[0] / "carbondata").exists():

            self.figure.clear()

            fig = self.figure.add_subplot(111)

            self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.07, top=0.95, wspace=0.05, hspace=0.05)

            self.isomerindex = int(str(self.IsomerSelect.currentText())[-1]) -1

            self.isomer = ui.table_widget.Tab1.worker.Isomers[self.isomerindex]
            self.assigned_shifts = self.isomer.Cshifts

            self.assigned_peaks = []

            for peak in self.isomer.Cexp:

                if peak != '':
                    self.assigned_peaks.append(peak)

            self.assigned_labels = self.isomer.Clabels

            ###################will probs need to fix sorting here

            exppeaks_ppm = self.xdata[self.exppeaks].tolist()

            shiftlist = self.assigned_shifts

            totallist = exppeaks_ppm + shiftlist

            fig.set_xlim([max(totallist) + 10, min(totallist) - 10])

            fig.plot(self.xdata, self.ydata, color='grey', linewidth=0.75, label='experimental spectrum')
            fig.plot(self.xdata, self.simulated_ydata, label='simulated spectrum')

            fig.set_xlabel('PPM')  # axis labels

            # plt.yticks([], [])
            # fig.set_title(str(ui.table_widget.Tab1.worker.settings.InputFiles[0]) +
            #             "\nCarbon NMR of Isomer " + str(self.isomerindex + 1) + "\n Number of Peaks Found = " + str(
            #  len(self.exppeaks)))

            # plot assignments

            for ind1, peak in enumerate(self.assigned_peaks):
                wh = np.argmin(abs(self.xdata - peak))

                fig.plot([peak, self.assigned_shifts[ind1]],
                         [self.ydata[wh], 1.1], linewidth=0.5, color='cyan')

            prev = round(exppeaks_ppm[0], 2)

            count = 0

            # annotate peak locations

            for x, txt in enumerate([round(i, 2) for i in exppeaks_ppm]):

                if abs(prev - txt) < 5:

                    count += 1
                else:
                    count = 0
                    prev = txt

                if exppeaks_ppm[x] in self.assigned_peaks:
                    color = 'C1'
                else:
                    color = 'grey'

                fig.annotate(txt, (exppeaks_ppm[x], -0.06 - 0.025 * count), color=color)

                fig.plot(exppeaks_ppm[x], self.ydata[self.exppeaks[x]], 'o', color=color)

            if len(self.removed) > 0:
                fig.plot(self.xdata[self.removed],
                         self.simulated_ydata[self.removed], "ro")

            # annotate shift positions

            count = 0

            ####some quick sorting

            argss = np.argsort(self.assigned_shifts)
            sortshifts = np.sort(self.assigned_shifts)[::-1]
            slabels = np.array(self.assigned_labels)[argss][::-1]

            prev = sortshifts[0]

            for x, txt in enumerate(slabels):

                if abs(prev - sortshifts[x]) < 4:
                    count += 1
                else:
                    count = 0
                    prev = sortshifts[x]

                fig.annotate(txt, (sortshifts[x], + 2.05 + 0.05 * count))

            ##########

            simulated_calc_ydata = np.zeros(len(self.xdata))

            for peak in self.assigned_shifts:
                y = np.exp(-0.5 * ((self.xdata - peak) / 0.002) ** 2)
                simulated_calc_ydata += y

            scaling_factor = np.amax(self.simulated_ydata) / np.amax(simulated_calc_ydata)

            simulated_calc_ydata = simulated_calc_ydata * scaling_factor

            #########

            fig.plot(self.xdata, simulated_calc_ydata + 1.1, label='calculated spectrum')

            self.canvas.draw()

    def selectpoint(self, event):

        if event.dblclick:

            self.xpos = event.xdata

            self.ypos = event.ydata

            # find the cloest point to the click

            if self.xpos is None and self.xpos is None:

                self.PlotCarbon()

            elif self.ypos < 1:

                mindis = np.argmin(abs(self.xdata[self.exppeaks] - self.xpos))

                self.PlotCarbonSelected(mindis)

                # find which protons this peak has been assigned to

                p = np.where(self.assigned_peaks == self.xdata[self.exppeaks[mindis]])[0]

                la = [int(self.assigned_labels[j][1:]) - 1 for j in p]

                self.RenderImage(la)

            else:

                # find closest assigned shift

                m = np.argmin(abs(np.array(self.assigned_shifts) - self.xpos))

                # find which peak this is assigned to

                p = self.assigned_peaks[m]

                mindis = np.argmin(abs(p - self.xdata[self.exppeaks]))

                p = np.where(self.assigned_peaks == self.xdata[self.exppeaks[mindis]])[0]

                la = [int(self.assigned_labels[j][1:]) - 1 for j in p]

                self.PlotCarbonSelected(mindis)

                self.RenderImage(la)

    def PlotCarbonSelected(self, mindis):

        if ui.table_widget.Tab1.worker.settings.OutputFolder == '':

            pdir = Path.cwd() / "Pickles"

        else:

            pdir = ui.table_widget.Tab1.worker.settings.OutputFolder / "Pickles"

        if Path(pdir / ui.table_widget.Tab1.worker.settings.InputFiles[0] / "carbondata").exists():

            self.figure.clear()

            fig = self.figure.add_subplot(111)

            self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.07, top=0.95, wspace=0.05, hspace=0.05)

            exppeaks_ppm = self.xdata[self.exppeaks].tolist()

            shiftlist = self.assigned_shifts

            totallist = exppeaks_ppm + shiftlist

            fig.set_xlim([max(totallist) + 10, min(totallist) - 10])

            fig.plot(self.xdata, self.ydata, color='grey', linewidth=0.75, label='experimental spectrum')

            fig.plot(self.xdata, self.simulated_ydata, label='simulated spectrum')

            fig.set_xlabel('PPM')  # axis labels

            assigned_peak = exppeaks_ppm[mindis]

            s = np.where(np.round(self.assigned_peaks, 8) == np.round(assigned_peak, 8))[0]

            assigned_shift = np.array(self.assigned_shifts)[s]

            assigned_label = np.array(self.assigned_labels)[s]

            for i in assigned_shift:
                fig.plot([assigned_peak, i], [self.ydata[self.exppeaks[mindis]], 1.1], color='cyan')

            # plot assignments

            # for ind1, peak in enumerate(self.assigned_peaks):
            #   wh = np.argmin(abs(self.xdata - peak))

            #  fig.plot([peak, self.assigned_shifts[ind1]],
            #          [self.ydata[wh], 1.1], linewidth=0.5, color='cyan')

            prev = round(exppeaks_ppm[0], 2)

            count = 0

            # annotate peak locations

            for x, txt in enumerate([round(i, 2) for i in exppeaks_ppm]):

                if exppeaks_ppm[x] == assigned_peak:

                    color = 'cyan'

                elif exppeaks_ppm[x] in self.assigned_peaks:

                    color = 'C1'

                else:
                    color = 'grey'

                fig.plot(exppeaks_ppm[x], self.ydata[self.exppeaks[x]], 'o', color=color)

            fig.annotate(str(round(assigned_peak, 2)), (assigned_peak, -0.06), color='cyan')

            if len(self.removed) > 0:
                fig.plot(self.xdata[self.removed],
                         self.simulated_ydata[self.removed], "ro")

            for l, s in zip(assigned_label, assigned_shift):
                fig.annotate(l, (s, 2.05), color="cyan")

            ##########

            simulated_calc_ydata = np.zeros(len(self.xdata))

            for peak in self.assigned_shifts:
                y = np.exp(-0.5 * ((self.xdata - peak) / 0.002) ** 2)

                simulated_calc_ydata += y

            scaling_factor = np.amax(self.simulated_ydata) / np.amax(simulated_calc_ydata)

            for peak in self.assigned_shifts:
                y = np.exp(-0.5 * ((self.xdata - peak) / 0.002) ** 2)

                fig.plot(self.xdata, y * scaling_factor + 1.1, label='calculated spectrum', color='grey', alpha=0.5)

            for peak in assigned_shift:
                y = np.exp(-0.5 * ((self.xdata - peak) / 0.002) ** 2)

                fig.plot(self.xdata, y * scaling_factor + 1.1, label='calculated spectrum', color='C1')

            #########

            self.canvas.draw()

    def RenderImage(self, atom):

        highlight = {}

        for i in atom:
            highlight[i] = (0, 1, 1)

        isomerindex = int(str(self.IsomerSelect.currentText())[-1]) -1

        m = Chem.MolFromMolFile(str(ui.table_widget.Tab1.worker.settings.InputFilesPaths[isomerindex]).split('.sdf')[0] + '.sdf',removeHs=False)

        Chem.Compute2DCoords(m)

        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)

        drawer.DrawMolecule(m, highlightAtoms=atom, highlightAtomColors=highlight)

        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace('svg:', '')

        svg_bytes = bytearray(svg, encoding='utf-8')

        self.image.renderer().load(svg_bytes)

        self.image.setGeometry(QtCore.QRect(0, 50, 300, 300))

        ui.update()

        # f = open("f.svg", "w+")

        # f.write(str(svg))


class ConformerTab(QtWidgets.QWidget):
    def __init__(self):

        super(ConformerTab, self).__init__()

        self.layout = QtWidgets.QGridLayout()

        self.setLayout(self.layout)

        self.IsomerSelect = QtWidgets.QComboBox(self)

        self.Isomers = ui.table_widget.Tab1.worker.Isomers

        self.IsomerSelect.addItems(self.Isomer_number())

        self.conformertable = QtWidgets.QTableWidget(self)

        self.layout.addWidget(self.conformertable, 0, 0)

        self.conformertable.setColumnCount(3)

        self.conformertable.setHorizontalHeaderLabels(["Conformer", "Energy", "Population"])

        self.IsomerSelect.currentIndexChanged.connect(self.populate_table)

        self.layout.setContentsMargins(0, 50, 0, 0)

        self.conffigure = Figure()

        self.confcanvas = FigureCanvas(self.conffigure)

        self.conformertable.itemSelectionChanged.connect(self.plot_conformers)

        self.layout.addWidget(self.confcanvas, 0, 1)

        self.confcanvas.mpl_connect('button_press_event', self.selectpoint)

    def selectpoint(self, event):

        self.xpos = event.xdata

        self.ypos = event.ydata

        # find the cloest point to the click

        coords = self.conffig.transData.transform((self.xpos, self.ypos))

        coordinates = np.array(self.conffig.transData.transform(list(zip(self.energies, self.populations))))

        mindis = np.argmin((coordinates[:, 0] - coords[0]) ** 2 + (coordinates[:, 1] - coords[1]) ** 2)

        self.conformertable.selectRow(mindis)

    def Isomer_number(self):

        Isomer_list = []

        for c, i in enumerate(self.Isomers):
            Isomer_list.append("Isomer " + str(c + 1))

        return Isomer_list

    def populate_table(self):

        c = 0

        self.isomerindex = int(str(self.IsomerSelect.currentText())[-1]) - 1

        self.conformertable.setRowCount(0)

        self.conformertable.setRowCount(len(self.Isomers[self.isomerindex].Energies))

        for energy, population in zip(self.Isomers[self.isomerindex].Energies,
                                      self.Isomers[self.isomerindex].Populations):
            self.conformertable.setItem(c, 0, QtWidgets.QTableWidgetItem(str(c)))

            self.conformertable.setItem(c, 1, QtWidgets.QTableWidgetItem(str(energy)))

            self.conformertable.setItem(c, 2, QtWidgets.QTableWidgetItem(str(population)))

            c += 1

        self.conformertable.selectRow(0)

        self.plot_conformers()

    def plot_conformers(self):

        self.conffigure.clear()

        self.conffig = self.conffigure.add_subplot(111)

        self.isomerindex = int(str(self.IsomerSelect.currentText())[-1]) - 1

        self.energies = np.array(self.Isomers[self.isomerindex].Energies)

        self.populations = np.array(self.Isomers[self.isomerindex].Populations)

        s = np.argsort(self.energies)

        self.conffig.plot(self.energies[s], self.populations[s],color = "C1")

        self.conffig.plot(self.energies, self.populations, 'o', color='C0', alpha=0.5)

        self.conffig.set_xlabel("Energy (Kcal)")

        self.conffig.set_ylabel("Population")

        selected_conformer = int(self.conformertable.currentRow())

        E = self.energies[selected_conformer]

        P = self.populations[selected_conformer]

        self.conffig.plot(E, P, 'o', color='red', alpha=0.75)

        self.conffig.plot([0, E], [P, P], color='red', alpha=0.75)

        self.conffig.plot([E, E], [0, P], color='red', alpha=0.75)

        self.confcanvas.draw()


class PyDP4WorkerObject(QtCore.QObject):

    finished = QtCore.pyqtSignal()

    def runPyDP4(self):

        launchdir = Path.cwd()

        print(ui.table_widget.Tab1.settings.OutputFolder)
        os.chdir(ui.table_widget.Tab1.settings.OutputFolder)

        self.NMRData, self.Isomers, self.settings, self.DP4Data,self.DP5Data = PyDP4.main(ui.table_widget.Tab1.settings)
        os.chdir(launchdir)
        self.finished.emit()


class WriteStream(object):
    def __init__(self, queue):
        self.queue = queue

    def write(self, text):
        self.queue.put(text)


class MyReceiver(QtCore.QObject):
    mysignal = QtCore.pyqtSignal(str)

    def __init__(self, queue, *args, **kwargs):
        QtCore.QObject.__init__(self, *args, **kwargs)
        self.queue = queue

    @QtCore.pyqtSlot()
    def run(self):
        while True:
            text = self.queue.get()
            self.mysignal.emit(text)


def lorentzian(p, w, p0, A):
    x = (p0 - p) / (w / 2)
    L = A / (1 + x ** 2)

    return L


def ReadParamFile(f, t):
    infile = open(f, 'r')
    inp = infile.readlines()
    infile.close()

    if t not in inp[0]:
        print("Wrong parameter file type, exiting...")
        quit()

    if t == 'm':
        Cmeans = [float(x) for x in inp[1].split(',')]
        Cstdevs = [float(x) for x in inp[2].split(',')]
        Hmeans = [float(x) for x in inp[3].split(',')]
        Hstdevs = [float(x) for x in inp[4].split(',')]

        return Cmeans, Cstdevs, Hmeans, Hstdevs


q = queue.Queue()

sys.stdout = WriteStream(q)

app = QtWidgets.QApplication(sys.argv)

ui = Window()

ui.show()

thread = QtCore.QThread()
my_receiver = MyReceiver(q)
my_receiver.mysignal.connect(ui.table_widget.Tab1.append_text)
my_receiver.moveToThread(thread)
thread.started.connect(my_receiver.run)
thread.start()

sys.exit(app.exec_())


