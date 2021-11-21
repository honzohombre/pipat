# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtGui, QtWidgets
from tkinter import *
from tkinter import messagebox

import re, sys, os
import glob, shutil #for removing folders
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import datetime #to see how long the program takes to run

import pandas as pd #to read MyResult file from PDBbind
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB.PDBParser import PDBParser

"""
UPDATE NOTES
1-3
- created log files to show progress
- updated output file format to csv, showing details of each match
- searchdb will remove segmened files once blast search is complete
- bug fixes with changing directories on filterhits
- working on linux compatibility
"""

# DECLARING GLOBAL VARIABLES
dic_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
default_directory = os.getcwd()

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(591, 801)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("icon.png"), QtGui.QIcon.Selected, QtGui.QIcon.On)
        MainWindow.setWindowIcon(icon)
        MainWindow.setSizePolicy(sizePolicy)
        #MainWindow.setMinimumSize(QtCore.QSize(591, 801))
        #MainWindow.setMaximumSize(QtCore.QSize(591, 811))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tab_select = QtWidgets.QTabWidget(self.centralwidget)
        self.tab_select.setGeometry(QtCore.QRect(10, 10, 571, 581))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.tab_select.setFont(font)
        self.tab_select.setAutoFillBackground(False)
        self.tab_select.setObjectName("tab_select")
        self.tab0 = QtWidgets.QWidget()
        self.tab0.setObjectName("tab0")
        self.lbl0_overview = QtWidgets.QLabel(self.tab0)
        self.lbl0_overview.setGeometry(QtCore.QRect(10, 10, 551, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl0_overview.setFont(font)
        self.lbl0_overview.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.lbl0_overview.setWordWrap(True)
        self.lbl0_overview.setObjectName("lbl0_overview")
        self.lbl0_selectdb = QtWidgets.QLabel(self.tab0)
        self.lbl0_selectdb.setGeometry(QtCore.QRect(20, 80, 221, 31))
        self.lbl0_selectdb.setObjectName("lbl0_selectdb")
        self.cmb0_selectdb = QtWidgets.QComboBox(self.tab0)
        self.cmb0_selectdb.setGeometry(QtCore.QRect(250, 80, 301, 31))
        self.cmb0_selectdb.setObjectName("cmb0_selectdb")
        self.cmb0_selectdb.addItem("")
        self.grb0_parameters = QtWidgets.QGroupBox(self.tab0)
        self.grb0_parameters.setGeometry(QtCore.QRect(20, 120, 521, 221))
        self.grb0_parameters.setObjectName("grb0_parameters")
        self.lbl0_dist = QtWidgets.QLabel(self.grb0_parameters)
        self.lbl0_dist.setGeometry(QtCore.QRect(30, 40, 151, 61))
        self.lbl0_dist.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl0_dist.setWordWrap(True)
        self.lbl0_dist.setObjectName("lbl0_dist")
        self.spn0_dist = QtWidgets.QSpinBox(self.grb0_parameters)
        self.spn0_dist.setGeometry(QtCore.QRect(80, 111, 61, 31))
        self.spn0_dist.setMinimum(1)
        self.spn0_dist.setMaximum(20)
        self.spn0_dist.setProperty("value", 5)
        self.spn0_dist.setObjectName("spn0_dist")
        self.lbl0_lenfrag = QtWidgets.QLabel(self.grb0_parameters)
        self.lbl0_lenfrag.setGeometry(QtCore.QRect(210, 40, 141, 61))
        self.lbl0_lenfrag.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl0_lenfrag.setWordWrap(True)
        self.lbl0_lenfrag.setObjectName("lbl0_lenfrag")
        self.lbl0_overlap = QtWidgets.QLabel(self.grb0_parameters)
        self.lbl0_overlap.setGeometry(QtCore.QRect(370, 40, 131, 61))
        self.lbl0_overlap.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl0_overlap.setWordWrap(True)
        self.lbl0_overlap.setObjectName("lbl0_overlap")
        self.spn0_overlap = QtWidgets.QSpinBox(self.grb0_parameters)
        self.spn0_overlap.setGeometry(QtCore.QRect(410, 111, 61, 31))
        self.spn0_overlap.setObjectName("spn0_overlap")
        self.spn0_overlap.setProperty("value", 10)
        self.spn0_lenfrag = QtWidgets.QSpinBox(self.grb0_parameters)
        self.spn0_lenfrag.setGeometry(QtCore.QRect(250, 111, 61, 31))
        self.spn0_lenfrag.setMinimum(19)
        self.spn0_lenfrag.setValue(20)
        self.spn0_lenfrag.setMaximum(100)
        self.spn0_lenfrag.setObjectName("spn0_lenfrag")
        self.lbl0_calcdist = QtWidgets.QLabel(self.grb0_parameters)
        self.lbl0_calcdist.setGeometry(QtCore.QRect(20, 170, 221, 31))
        self.lbl0_calcdist.setObjectName("lbl0_calcdist")
        self.cmb0_calcdist = QtWidgets.QComboBox(self.grb0_parameters)
        self.cmb0_calcdist.setGeometry(QtCore.QRect(250, 170, 251, 31))
        self.cmb0_calcdist.setObjectName("cmb0_calcdist")
        self.cmb0_calcdist.addItem("")
        self.cmb0_calcdist.addItem("")
        self.grb0_advanced = QtWidgets.QGroupBox(self.tab0)
        self.grb0_advanced.setGeometry(QtCore.QRect(20, 350, 521, 171))
        self.grb0_advanced.setObjectName("grb0_advanced")
        self.lbl0_omit = QtWidgets.QLabel(self.grb0_advanced)
        self.lbl0_omit.setGeometry(QtCore.QRect(390, 30, 91, 41))
        self.lbl0_omit.setObjectName("lbl0_omit")
        self.spn0_hits = QtWidgets.QSpinBox(self.grb0_advanced)
        self.spn0_hits.setEnabled(False)
        self.spn0_hits.setGeometry(QtCore.QRect(310, 35, 61, 31))
        self.spn0_hits.setMinimum(1)
        self.spn0_hits.setMaximum(100)
        self.spn0_hits.setProperty("value", 1)
        self.spn0_hits.setObjectName("spn0_hits")
        self.chk_omit = QtWidgets.QCheckBox(self.grb0_advanced)
        self.chk_omit.setGeometry(QtCore.QRect(20, 30, 291, 41))
        self.chk_omit.setObjectName("chk_omit")
        self.chk_checkbinding = QtWidgets.QCheckBox(self.grb0_advanced)
        self.chk_checkbinding.setGeometry(QtCore.QRect(20, 60, 381, 41))
        self.chk_checkbinding.setObjectName("chk_checkbinding")
        self.chk_writedetails = QtWidgets.QCheckBox(self.grb0_advanced)
        self.chk_writedetails.setGeometry(QtCore.QRect(20, 90, 411, 41))
        self.chk_writedetails.setObjectName("chk_writedetails")
        self.chk_createdb = QtWidgets.QCheckBox(self.grb0_advanced)
        self.chk_createdb.setGeometry(QtCore.QRect(20, 120, 471, 41))
        self.chk_createdb.setObjectName("chk_createdb")
        self.tab_select.addTab(self.tab0, "")
        self.tab3 = QtWidgets.QWidget()
        self.tab3.setObjectName("tab3")
        self.lbl3_ = QtWidgets.QLabel(self.tab3)
        self.lbl3_.setGeometry(QtCore.QRect(10, 10, 551, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl3_.setFont(font)
        self.lbl3_.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.lbl3_.setWordWrap(True)
        self.lbl3_.setObjectName("lbl3_")
        self.lbl3_1 = QtWidgets.QLabel(self.tab3)
        self.lbl3_1.setGeometry(QtCore.QRect(20, 70, 221, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl3_1.setFont(font)
        self.lbl3_1.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl3_1.setWordWrap(True)
        self.lbl3_1.setObjectName("lbl3_1")
        self.cmb3_db = QtWidgets.QComboBox(self.tab3)
        self.cmb3_db.setGeometry(QtCore.QRect(250, 70, 291, 31))
        self.cmb3_db.setObjectName("cmb3_db")
        self.lbl3_2 = QtWidgets.QLabel(self.tab3)
        self.lbl3_2.setGeometry(QtCore.QRect(20, 120, 251, 31))
        self.lbl3_2.setObjectName("lbl3_2")
        self.spn3_hits = QtWidgets.QSpinBox(self.tab3)
        self.spn3_hits.setGeometry(QtCore.QRect(280, 120, 61, 31))
        self.spn3_hits.setMinimum(2)
        self.spn3_hits.setObjectName("spn3_hits")
        self.lbl3_3 = QtWidgets.QLabel(self.tab3)
        self.lbl3_3.setGeometry(QtCore.QRect(360, 120, 91, 31))
        self.lbl3_3.setObjectName("lbl3_3")
        self.lbl3_4 = QtWidgets.QLabel(self.tab3)
        self.lbl3_4.setGeometry(QtCore.QRect(20, 170, 201, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl3_4.setFont(font)
        self.lbl3_4.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl3_4.setWordWrap(True)
        self.lbl3_4.setObjectName("lbl3_4")
        self.chk3_createdb = QtWidgets.QCheckBox(self.tab3)
        self.chk3_createdb.setGeometry(QtCore.QRect(20, 220, 471, 41))
        self.chk3_createdb.setObjectName("chk3_createdb")
        self.txt3_newdb = QtWidgets.QLineEdit(self.tab3)
        self.txt3_newdb.setEnabled(False)
        self.txt3_newdb.setGeometry(QtCore.QRect(200, 170, 341, 31))
        self.txt3_newdb.setObjectName("txt3_newdb")
        self.tab_select.addTab(self.tab3, "")
        self.tab1 = QtWidgets.QWidget()
        self.tab1.setObjectName("tab1")
        self.lbl1_ = QtWidgets.QLabel(self.tab1)
        self.lbl1_.setGeometry(QtCore.QRect(10, 10, 541, 51))
        self.lbl1_.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.lbl1_.setObjectName("lbl1_")
        self.grb1_selectfile = QtWidgets.QGroupBox(self.tab1)
        self.grb1_selectfile.setGeometry(QtCore.QRect(20, 70, 521, 80))
        self.grb1_selectfile.setAcceptDrops(True)
        self.grb1_selectfile.setObjectName("grb1_selectfile")
        self.lbl1_1 = QtWidgets.QLabel(self.grb1_selectfile)
        self.lbl1_1.setGeometry(QtCore.QRect(30, 30, 131, 31))
        self.lbl1_1.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl1_1.setObjectName("lbl1_1")
        self.txt1_query = QtWidgets.QLineEdit(self.grb1_selectfile)
        self.txt1_query.setGeometry(QtCore.QRect(170, 31, 211, 31))
        self.txt1_query.setObjectName("txt1_query")
        self.btn1_browse = QtWidgets.QPushButton(self.grb1_selectfile)
        self.btn1_browse.setGeometry(QtCore.QRect(392, 30, 111, 31))
        self.btn1_browse.setObjectName("btn1_browse")
        self.grb1_parameters = QtWidgets.QGroupBox(self.tab1)
        self.grb1_parameters.setGeometry(QtCore.QRect(20, 159, 521, 191))
        self.grb1_parameters.setObjectName("grb1_parameters")
        self.lbl1_3 = QtWidgets.QLabel(self.grb1_parameters)
        self.lbl1_3.setGeometry(QtCore.QRect(270, 39, 131, 61))
        self.lbl1_3.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl1_3.setWordWrap(True)
        self.lbl1_3.setObjectName("lbl1_3")
        self.spn1_overlap = QtWidgets.QSpinBox(self.grb1_parameters)
        self.spn1_overlap.setGeometry(QtCore.QRect(310, 110, 61, 31))
        self.spn1_overlap.setObjectName("spn1_overlap")
        self.lbl1_2 = QtWidgets.QLabel(self.grb1_parameters)
        self.lbl1_2.setGeometry(QtCore.QRect(110, 39, 141, 61))
        self.lbl1_2.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl1_2.setWordWrap(True)
        self.lbl1_2.setObjectName("lbl1_2")
        self.spn1_lenfrag = QtWidgets.QSpinBox(self.grb1_parameters)
        self.spn1_lenfrag.setGeometry(QtCore.QRect(150, 110, 61, 31))
        self.spn1_lenfrag.setMinimum(20)
        self.spn1_lenfrag.setMaximum(100)
        self.spn1_lenfrag.setObjectName("spn1_lenfrag")
        self.grb1_outputfile = QtWidgets.QGroupBox(self.tab1)
        self.grb1_outputfile.setGeometry(QtCore.QRect(20, 360, 521, 80))
        self.grb1_outputfile.setObjectName("grb1_outputfile")
        self.lbl1_4 = QtWidgets.QLabel(self.grb1_outputfile)
        self.lbl1_4.setGeometry(QtCore.QRect(30, 30, 131, 31))
        self.lbl1_4.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl1_4.setObjectName("lbl1_4")
        self.txt1_outputfile = QtWidgets.QLineEdit(self.grb1_outputfile)
        self.txt1_outputfile.setGeometry(QtCore.QRect(170, 31, 331, 31))
        self.txt1_outputfile.setObjectName("txt1_outputfile")
        self.txt1_outputfile.setEnabled(False)
        self.tab_select.addTab(self.tab1, "")
        self.tab2 = QtWidgets.QWidget()
        self.tab2.setObjectName("tab2")
        self.lbl2_ = QtWidgets.QLabel(self.tab2)
        self.lbl2_.setGeometry(QtCore.QRect(10, 10, 541, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_.setFont(font)
        self.lbl2_.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.lbl2_.setWordWrap(True)
        self.lbl2_.setObjectName("lbl2_")
        self.lbl2_1 = QtWidgets.QLabel(self.tab2)
        self.lbl2_1.setGeometry(QtCore.QRect(20, 50, 221, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_1.setFont(font)
        self.lbl2_1.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl2_1.setWordWrap(True)
        self.lbl2_1.setObjectName("lbl2_1")
        self.grb2_segment = QtWidgets.QGroupBox(self.tab2)
        self.grb2_segment.setGeometry(QtCore.QRect(20, 90, 521, 281))
        self.grb2_segment.setAcceptDrops(True)
        self.grb2_segment.setObjectName("grb2_segment")
        self.rdb2_matchdb = QtWidgets.QRadioButton(self.grb2_segment)
        self.rdb2_matchdb.setGeometry(QtCore.QRect(40, 60, 401, 41))
        self.rdb2_matchdb.setChecked(True)
        self.rdb2_matchdb.setObjectName("rdb2_matchdb")
        self.rdb2_halfdb = QtWidgets.QRadioButton(self.grb2_segment)
        self.rdb2_halfdb.setGeometry(QtCore.QRect(40, 90, 391, 41))
        self.rdb2_halfdb.setObjectName("rdb2_halfdb")
        self.rdb2_manualsegment = QtWidgets.QRadioButton(self.grb2_segment)
        self.rdb2_manualsegment.setGeometry(QtCore.QRect(40, 120, 161, 41))
        self.rdb2_manualsegment.setObjectName("rdb2_manualsegment")
        self.spn2_querylength = QtWidgets.QSpinBox(self.grb2_segment)
        self.spn2_querylength.setEnabled(False)
        self.spn2_querylength.setGeometry(QtCore.QRect(200, 125, 42, 31))
        self.spn2_querylength.setProperty("value", 20)
        self.spn2_querylength.setObjectName("spn2_querylength")
        self.lbl2_2 = QtWidgets.QLabel(self.grb2_segment)
        self.lbl2_2.setGeometry(QtCore.QRect(270, 120, 81, 41))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_2.setFont(font)
        self.lbl2_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl2_2.setWordWrap(True)
        self.lbl2_2.setObjectName("lbl2_2")
        self.spn2_queryoverlap = QtWidgets.QSpinBox(self.grb2_segment)
        self.spn2_queryoverlap.setEnabled(False)
        self.spn2_queryoverlap.setGeometry(QtCore.QRect(350, 125, 42, 31))
        self.spn2_queryoverlap.setProperty("value", 10)
        self.spn2_queryoverlap.setObjectName("spn2_queryoverlap")
        self.chk2_segment = QtWidgets.QCheckBox(self.grb2_segment)
        self.chk2_segment.setGeometry(QtCore.QRect(20, 30, 341, 31))
        self.chk2_segment.setObjectName("chk2_segment")
        self.chk2_segment.setChecked(True)
        self.lbl2_5 = QtWidgets.QLabel(self.grb2_segment)
        self.lbl2_5.setGeometry(QtCore.QRect(270, 250, 231, 21))
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.lbl2_5.setFont(font)
        self.lbl2_5.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lbl2_5.setWordWrap(True)
        self.lbl2_5.setObjectName("lbl2_5")
        self.txt2_query1 = QtWidgets.QLineEdit(self.grb2_segment)
        self.txt2_query1.setGeometry(QtCore.QRect(150, 170, 231, 31))
        self.txt2_query1.setObjectName("txt2_query1")
        self.lbl2_4 = QtWidgets.QLabel(self.grb2_segment)
        self.lbl2_4.setGeometry(QtCore.QRect(20, 210, 121, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_4.setFont(font)
        self.lbl2_4.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl2_4.setWordWrap(True)
        self.lbl2_4.setObjectName("lbl2_4")
        self.lbl2_3 = QtWidgets.QLabel(self.grb2_segment)
        self.lbl2_3.setGeometry(QtCore.QRect(20, 170, 121, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_3.setFont(font)
        self.lbl2_3.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl2_3.setWordWrap(True)
        self.lbl2_3.setObjectName("lbl2_3")
        self.txt2_query2 = QtWidgets.QLineEdit(self.grb2_segment)
        self.txt2_query2.setGeometry(QtCore.QRect(150, 210, 231, 31))
        self.txt2_query2.setObjectName("txt2_query2")
        self.btn2_browse1 = QtWidgets.QPushButton(self.grb2_segment)
        self.btn2_browse1.setGeometry(QtCore.QRect(390, 170, 111, 31))
        self.btn2_browse1.setObjectName("btn2_browse1")
        self.btn2_browse2 = QtWidgets.QPushButton(self.grb2_segment)
        self.btn2_browse2.setGeometry(QtCore.QRect(390, 210, 111, 31))
        self.btn2_browse2.setObjectName("btn2_browse2")
        self.cmb2_db = QtWidgets.QComboBox(self.tab2)
        self.cmb2_db.setGeometry(QtCore.QRect(250, 50, 291, 31))
        self.cmb2_db.setObjectName("cmb2_db")
        self.grb2_advanced = QtWidgets.QGroupBox(self.tab2)
        self.grb2_advanced.setGeometry(QtCore.QRect(20, 380, 521, 141))
        self.grb2_advanced.setTitle("")
        self.grb2_advanced.setObjectName("grb2_advanced")
        self.lbl2_6 = QtWidgets.QLabel(self.grb2_advanced)
        self.lbl2_6.setGeometry(QtCore.QRect(20, 10, 221, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.lbl2_6.setFont(font)
        self.lbl2_6.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbl2_6.setWordWrap(True)
        self.lbl2_6.setObjectName("lbl2_6")
        self.txt2_evalue = QtWidgets.QLineEdit(self.grb2_advanced)
        self.txt2_evalue.setGeometry(QtCore.QRect(250, 10, 251, 31))
        self.txt2_evalue.setObjectName("txt2_evalue")
        self.chk2_makeblastdb = QtWidgets.QCheckBox(self.grb2_advanced)
        self.chk2_makeblastdb.setGeometry(QtCore.QRect(20, 50, 501, 31))
        #self.chk2_makeblastdb.setChecked(True)
        self.chk2_makeblastdb.setObjectName("chk2_makeblastdb")
        self.chk2_overwritelog = QtWidgets.QCheckBox(self.grb2_advanced)
        self.chk2_overwritelog.setGeometry(QtCore.QRect(20, 80, 341, 31))
        #self.chk2_overwritelog.setChecked(True)
        self.chk2_overwritelog.setObjectName("chk2_overwritelog")
        self.tab_select.addTab(self.tab2, "")
        self.btn_run = QtWidgets.QPushButton(self.centralwidget)
        self.btn_run.setGeometry(QtCore.QRect(290, 740, 141, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.btn_run.setFont(font)
        self.btn_run.setObjectName("btn_run")
        self.btn_cancel = QtWidgets.QPushButton(self.centralwidget)
        self.btn_cancel.setGeometry(QtCore.QRect(440, 740, 141, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.btn_cancel.setFont(font)
        self.btn_cancel.setObjectName("btn_cancel")
        self.txt_progress = QtWidgets.QTextEdit(self.centralwidget)
        self.txt_progress.setEnabled(False)
        self.txt_progress.setGeometry(QtCore.QRect(10, 600, 571, 81))
        self.txt_progress.setObjectName("txt_progress")
        self.prg_progress = QtWidgets.QProgressBar(self.centralwidget)
        self.prg_progress.setGeometry(QtCore.QRect(10, 690, 571, 23))
        self.prg_progress.setProperty("value", 0)
        self.prg_progress.setObjectName("prg_progress")
        self.btn_help = QtWidgets.QPushButton(self.centralwidget)
        self.btn_help.setGeometry(QtCore.QRect(10, 740, 141, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.btn_help.setFont(font)
        self.btn_help.setObjectName("btn_help")
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionNew = QtWidgets.QAction(MainWindow)
        self.actionNew.setObjectName("actionNew")

        self.setConnections(MainWindow)
        self.retranslateUi(MainWindow)
        #self.tab_select.setCurrentIndex(3)

        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "PIPAT"))
        self.lbl0_overview.setText(_translate("MainWindow", "Use a selection of PDB files containing protein interactions to generate a local BLAST database."))
        self.lbl0_selectdb.setText(_translate("MainWindow", "Select PDB collection folder: "))
        self.cmb0_selectdb.setItemText(0, _translate("MainWindow", "PP (default)"))
        self.grb0_parameters.setTitle(_translate("MainWindow", "Choose Parameters"))
        self.lbl0_dist.setText(_translate("MainWindow", "Maximum distance between objects (angstroms)"))
        self.lbl0_lenfrag.setText(_translate("MainWindow", "Segment Length (# of residues)"))
        self.lbl0_overlap.setText(_translate("MainWindow", "Overlap \n(# of residues)"))
        self.lbl0_calcdist.setText(_translate("MainWindow", "Calculate distance between:"))
        self.cmb0_calcdist.setItemText(0, _translate("MainWindow", "geometric center of side chains"))
        self.cmb0_calcdist.setItemText(1, _translate("MainWindow", "alpha carbons"))
        self.grb0_advanced.setTitle(_translate("MainWindow", "Advanced Configuration"))
        self.lbl0_omit.setText(_translate("MainWindow", "match(es)"))
        self.chk_omit.setText(_translate("MainWindow", "Omit pairings containing less than"))
        self.chk_checkbinding.setText(_translate("MainWindow", "Check MyResult.xls for binding affinity"))
        self.chk_writedetails.setText(_translate("MainWindow", "Print details of each hit to runtime log"))
        self.chk_createdb.setText(_translate("MainWindow", "Continue to create BLAST database when run finishes"))
        self.tab_select.setTabText(self.tab_select.indexOf(self.tab0), _translate("MainWindow", "Create Database"))
        self.lbl3_.setText(_translate("MainWindow", "Having already created a local database, you can generate a smaller replicate by filtering entries with low interaction counts."))
        self.lbl3_1.setText(_translate("MainWindow", "Select local database folder:"))
        self.lbl3_2.setText(_translate("MainWindow", "Filter out pairings with less than"))
        self.lbl3_3.setText(_translate("MainWindow", "matches."))
        self.lbl3_4.setText(_translate("MainWindow", "New database file:"))
        self.chk3_createdb.setText(_translate("MainWindow", "Continue to create BLAST database after filtering"))
        self.tab_select.setTabText(self.tab_select.indexOf(self.tab3), _translate("MainWindow", "Filter Hits"))
        self.lbl1_.setText(_translate("MainWindow", "Segment an amino acid sequence into smaller peptide sequences,\nwhere the segmented peptides are saved in a single FASTA file."))
        self.grb1_selectfile.setTitle(_translate("MainWindow", "Select File"))
        self.lbl1_1.setText(_translate("MainWindow", "Query filename: "))
        self.btn1_browse.setText(_translate("MainWindow", "Browse..."))
        self.grb1_parameters.setTitle(_translate("MainWindow", "Choose Parameters"))
        self.lbl1_3.setText(_translate("MainWindow", "Overlap \n(# of residues)"))
        self.lbl1_2.setText(_translate("MainWindow", "Segment Length (# of residues)"))
        self.grb1_outputfile.setTitle(_translate("MainWindow", "Output File"))
        self.lbl1_4.setText(_translate("MainWindow", "Output filename: "))
        self.tab_select.setTabText(self.tab_select.indexOf(self.tab1), _translate("MainWindow", "Segment Query"))
        self.lbl2_.setText(_translate("MainWindow", "Search against a local database to predict interactions between queries."))
        self.lbl2_1.setText(_translate("MainWindow", "Select local database folder:"))
        self.grb2_segment.setTitle(_translate("MainWindow", "Choose Parameters"))
        self.rdb2_matchdb.setText(_translate("MainWindow", "match segment length / overlap with database"))
        self.rdb2_halfdb.setText(_translate("MainWindow", "use half the length / overlap of the database"))
        self.rdb2_manualsegment.setText(_translate("MainWindow", "segment length:"))
        self.lbl2_2.setText(_translate("MainWindow", "/ overlap: "))
        self.chk2_segment.setText(_translate("MainWindow", "Segment query sequences."))
        self.lbl2_5.setText(_translate("MainWindow", "*You can drag & drop a file here."))
        self.lbl2_4.setText(_translate("MainWindow", "Select query 2:"))
        self.lbl2_3.setText(_translate("MainWindow", "Select query 1:"))
        self.btn2_browse2.setText(_translate("MainWindow", "Browse..."))
        self.btn2_browse1.setText(_translate("MainWindow", "Browse..."))
        self.lbl2_6.setText(_translate("MainWindow", "BLAST search e-value cutoff:"))
        self.txt2_evalue.setText(_translate("MainWindow", "3200"))
        self.chk2_makeblastdb.setText(_translate("MainWindow", "Create (or recreate) BLAST database before search."))
        self.chk2_overwritelog.setText(_translate("MainWindow", "Overwrite existing run summaries."))
        self.tab_select.setTabText(self.tab_select.indexOf(self.tab2), _translate("MainWindow", "Search Database"))
        self.btn_run.setText(_translate("MainWindow", "&Run"))
        self.btn_cancel.setText(_translate("MainWindow", "Close"))
        self.btn_help.setText(_translate("MainWindow", "Help"))

    #CONNECTIONS are set here - move items from this function back to SetupUi when finished
    def setConnections(self, MainWindow):
        #CONNECTING BUTTONS TO FUNCTIONS
        self.btn_cancel.clicked.connect(self.quitProgram)
        self.btn_run.clicked.connect(self.runProgram)

        # we use lambda to pass an argument into a function
        self.btn1_browse.clicked.connect(lambda: self.fetchFile("btn1_browse"))
        self.btn2_browse1.clicked.connect(lambda: self.fetchFile("btn2_browse1"))
        self.btn2_browse2.clicked.connect(lambda: self.fetchFile("btn2_browse2"))

        #showing database names when tab indexes changes
        self.tab_select.currentChanged.connect(self.loadBoxes)

        self.chk_omit.stateChanged.connect(self.enableSpin)
        self.chk2_segment.stateChanged.connect(self.changeSegmentOption)
        self.rdb2_manualsegment.toggled.connect(self.changeSegmentOption)
        #self.cmb2_db.currentIndexChanged.connect(self.setMaxLength)
        #self.cmb3_db.currentIndexChanged.connect(self.setMaxLength)

        self.cmb3_db.currentIndexChanged.connect(self.showOutputfiles)
        self.spn3_hits.valueChanged.connect(self.showOutputfiles)

#PERIPHERAL FUNCTIONS - REACT TO USER INPUT AS FORM IS BEING FILLED
    #When a combobox is selected, show available databases in the folder
    def loadBoxes(self):
        self.cmb2_db.clear()
        self.cmb3_db.clear()
        lst_directories = next(os.walk('.'))[1]
        for d in lst_directories:
            if d[0] == "d" and d[1].isdigit():
                self.cmb2_db.addItem(d)
                lastparam = d.split('-')[-1][0]
                if lastparam != "h": #because we don't want to edit already-edited databases
                    self.cmb3_db.addItem(d)

    def changeSegmentOption(self):
        if self.chk2_segment.isChecked() == True:
            self.rdb2_halfdb.setEnabled(True)
            self.rdb2_manualsegment.setEnabled(True)
            self.rdb2_matchdb.setEnabled(True)
            if self.rdb2_manualsegment.isChecked() == True:
                self.spn2_queryoverlap.setEnabled(True)
                self.spn2_querylength.setEnabled(True)
            else:
                self.spn2_queryoverlap.setEnabled(False)
                self.spn2_querylength.setEnabled(False)
        else:
            self.rdb2_halfdb.setEnabled(False)
            self.rdb2_manualsegment.setEnabled(False)
            self.rdb2_matchdb.setEnabled(False)
            self.spn2_querylength.setEnabled(False)
            self.spn2_queryoverlap.setEnabled(False)

    #When user selects an item from the combobox, show them what the output looks like
    def showOutputfiles(self):
        str_inputfile = self.cmb3_db.currentText()
        int_minhits = self.spn3_hits.value()
        str_outputfile = str_inputfile.split('.f')[0]+"-h"+str(int_minhits)
        self.txt3_newdb.setText(str_outputfile)

    #When the browse button is pressed, take the file and put it in the textbox next to the button
    def fetchFile(self, btnName):
        cwd = os.getcwd()
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", cwd, "FASTA Files (*.fa *.fasta)")
        if fileName:
            if btnName == "btn1_browse":
                self.txt1_query.setText(fileName)
                
            elif btnName == "btn2_browse1":
                self.txt2_query1.setText(fileName)
            elif btnName == "btn2_browse2":
                self.txt2_query2.setText(fileName)
    
    #Does not allow the user to increase the segment length beyond the database (maybe there's no need for this?)
    def setMaxLength(self):
        try:
            str_dbname = self.cmb2_db.currentText()
            res = re.search(r'd(\d.*)-l(\d+)-o(\d+)(-h\d+)?', str_dbname)
            
            int_maxlength = int(res.group(2))
            int_maxoverlap = int(res.group(3))-1
            self.spn2_querylength.setMaximum(int_maxlength)
            self.spn2_queryoverlap.setMaximum(int_maxlength-1)
        except:
            print("Cannot set max length because no db folders are loaded")
            print("db name: "+str_dbname)


    def enableSpin(self):
        if self.chk_omit.isChecked: self.spn0_hits.setEnabled(True)
        else: self.spn0_hits.setEnabled(False)

    def quitProgram(self):
        quit()

# MAIN FUNCTION
    def runProgram(self):
        global default_directory
        os.chdir(default_directory)
        int_index = self.tab_select.currentIndex()
        self.prg_progress.setValue(0)

        time_start = datetime.datetime.now()

        if int_index == 0: #Create Database
            self.txt_progress.setText("Creating local database... ")
            str_logfile = "log_createdb.txt"
            write_logfile = open(str_logfile,"wt")
            write_logfile.write("Gathering information to create database...\nStart time: "+str(time_start)+"\n")

            #Loading default parameters
            input_files = str(Path("PP") / "*.ent.pdb")
            num_files = len(glob.glob(input_files))

            flt_dist = self.spn0_dist.value()
            int_overlap = self.spn0_overlap.value()
            int_lenfrag = self.spn0_lenfrag.value()
            if self.cmb0_calcdist.currentIndex() == 0: tf_sidechain = True
            elif self.cmb0_calcdist.currentIndex() == 1: tf_sidechain = False

            str_filename = "d"+str(flt_dist)+"-l"+str(int_lenfrag)+"-o"+str(int_overlap)
            int_minhits = self.spn0_hits.value()
            if int_minhits > 1: str_filename += "-h"+str(int_minhits)
            write_logfile.write("Database name: "+str_filename+"\n")
            
            if self.chk_checkbinding.isChecked(): 
                tf_checkbinding = True
                write_logfile.write("- checking MyResult.xls for binding affinities\n")
            else: tf_checkbinding = False
            if self.chk_createdb.isChecked(): 
                tf_makedb = True
                write_logfile.write("- proceed to make BLAST database upon finishing\n")
            else: tf_makedb = False
            if self.chk_writedetails.isChecked():
                tf_details = True
                write_logfile.write("- print details of each hit to the logbook\n")
            else: tf_details = False

            outputfile = Path(str_filename) / (str_filename+".fasta")

            dic_parameters = {"d":flt_dist, "l":int_lenfrag, "o":int_overlap, "sc":tf_sidechain, "h":int_minhits, \
                "binding":tf_checkbinding, "output":outputfile, "log":write_logfile, "log_details":tf_details}
            #seqread: number of segments (grand total), filesread: completed files, files: total files in the archive
            dic_info = {"num_seqsread":0, "num_filesread":0, "num_files":num_files}

            x = reset_files(str_filename)
            #todo - msgbox does not work
            x = True
            if not x: 
                #do not proceed to running program - show status in textbox
                self.txt_progress.setText("File: "+str_filename+" already exists.")
                write_logfile.write(("File: "+str_filename+" already exists.\n"))
            else:
                counter = 0
                for pdb_filename in glob.glob(input_files):
                    dic_info = writelines(pdb_filename, dic_parameters, dic_info, self)
                    counter += 1
                str_showinfo = "Finished reading "+str(dic_info['num_filesread'])+" of "+str(dic_info['num_files'])+" files."
                if tf_makedb:
                    self.txt_progress.setText(str_showinfo+"\nMaking BLAST database...")
                    makedb(str_filename)
                    str_showinfo += "\nCreated BLAST database: "+str_filename

            time_end = datetime.datetime.now()
            time_difference = time_end - time_start

            self.txt_progress.setText("Complete.\n"+str_showinfo)

        elif int_index == 1:
            self.txt_progress.setText("Filtering low hit counts... ")

            str_logfile = "log_filterhits.txt"
            write_logfile = open(str_logfile, "wt")
            write_logfile.write("Filtering low hit counts...\nStart time: "+str(time_start)+"\n")

            str_inputfile = self.cmb3_db.currentText()
            str_outputfile = self.txt3_newdb.text()
            str_hitcount = self.spn3_hits.value()

            if Path(str_outputfile).exists(): shutil.rmtree(str_outputfile)
            os.mkdir(str_outputfile)
            pth_outputfile = Path(str_outputfile) / (str_outputfile+".fasta")
            pth_inputfile = Path(str_inputfile) / (str_inputfile+".fasta")
            int_seqsread, int_seqstotal = filterhits(pth_inputfile, pth_outputfile, str_hitcount, self)

            str_output =  "Extracted "+str(int_seqsread)+" of "+str(int_seqstotal)+" sequences.\n" \
                "Destination: "+str_outputfile+"\n"
            self.txt_progress.setText("Filtering complete.\n"+str_output)
            write_logfile.write(str_output)

            if self.chk3_createdb.isChecked():
                self.txt_progress.setText("Filtering complete.\n"+str_output+"\nMaking BlastDB...")
                write_logfile.write("Continuing to make BLAST database: "+str_outputfile)
                makedb(str_outputfile) #making db from filterhits

            self.txt_progress.setText("Complete.")
            
        elif int_index == 2: #Segments the fasta and saves result in same directory
            self.txt_progress.setText("Segmenting query sequences...")

            str_logfile = "log_segmenter.txt"
            write_logfile = open(str_logfile, "wt")
            write_logfile.write("Running Segmenter...\nStart time: "+str(time_start)+"\n")

            #retrieve information that we need from the form
            str_query = self.txt1_query.text()
            int_lenfrag = self.spn1_lenfrag.value()
            int_overlap = self.spn1_overlap.value()

            #make the changes necessary
            str_output = str_query.split('.')[0]+"_l"+str(int_lenfrag)+"-o"+str(int_overlap)+".fasta"
            self.txt1_outputfile.setText(str_output)
            str_status = "Segmenting query sequences...\n" \
                "File name: "+str_query+"\n" \
                "Segment length: "+str(int_lenfrag)+"\n" \
                "Overlap distance: "+str(int_overlap)
            self.txt_progress.setText(str_status)
            write_logfile.write(str_status)

            #Call the function "segmenter"
            int_numpeps = segmenter(str_query, str_output, int_lenfrag, int_overlap)
            self.txt_progress.setText("Complete.\n Contents saved in: "+str(str_output))
            self.prg_progress.setValue(100)
            write_logfile.write("Created "+str_output+" containing "+str(int_numpeps)+" segmented sequences.")

        elif int_index == 3:
            self.txt_progress.setText("Performing database search... ")
            str_logfile = "log_searchdb.txt"
            str_dbname = self.cmb2_db.currentText() #e.g. "d5-l20-o10"
            str_query1 = self.txt2_query1.text() #file address of the query
            str_query2 = self.txt2_query2.text()

            write_logfile = open(str_logfile,"wt")
            write_logfile.write("Performing database search...\nStart time: "+str(time_start)+"\n")

            if self.chk2_makeblastdb.isChecked(): #creating new BLAST database (if checked)
                print("Creating BLAST database... ")
                self.txt_progress.setText("Creating BLAST database...")
                makedb(str_dbname)
                write_logfile.write("Created new BLAST database: "+str_dbname)
                self.txt_progress.setText("Finished creating BLAST database: "+str_dbname)

            flt_evalue = float(self.txt2_evalue.text())
            tf_overwrite = self.chk2_overwritelog.isChecked()

            if self.chk2_segment.isChecked():
                #Segmenting query sequences by first retrieving database information
                res = re.search(r'd(\d.*)-l(\d+)-o(\d+)(-h\d+)?', str_dbname)
                int_lenfrag = int(res.group(2))
                int_overlap = int(res.group(3))
                if self.rdb2_matchdb.isChecked():
                    print("matchdb")

                elif self.rdb2_halfdb.isChecked():
                    int_lenfrag = int_lenfrag // 2
                    int_overlap = int_overlap // 2

                elif self.rdb2_manualsegment.isChecked():
                    int_lenfrag = self.spn2_querylength.value()
                    int_overlap = self.spn2_queryoverlap.value()
                    
                write_logfile.write("- Segmented peptide length (residues): "+str(int_lenfrag)+"\t")
                write_logfile.write("- Overlapping distance (residues): "+str(int_overlap)+"\n")
                str_segmentedquery1 = str_query1.split(".f")[0]+"_l"+str(int_lenfrag)+"-o"+str(int_overlap)+".fasta"
                str_segmentedquery2 = str_query2.split(".f")[0]+"_l"+str(int_lenfrag)+"-o"+str(int_overlap)+".fasta"

                #Segmenting the queries, saving to file so we can perform BLAST search
                int_newseqs = segmenter(str_query1, str_segmentedquery1, int_lenfrag, int_overlap)
                write_logfile.write(str(int_newseqs)+" peptides created from "+str_query1+", \nsaved to"+str_segmentedquery1+"\n")
                int_newseqs = segmenter(str_query2, str_segmentedquery2, int_lenfrag, int_overlap)
                write_logfile.write(str(int_newseqs)+" peptides created from "+str_query2+", \nsaved to"+str_segmentedquery2+"\n")

                int_matches = blastsearch(str_segmentedquery1, str_segmentedquery2, str_dbname, flt_evalue, write_logfile, self)
                #opt - remove segmented query once search is complete
                os.remove(str_segmentedquery1)
                os.remove(str_segmentedquery2)

            else: #no segmenting selected
                write_logfile.write("Query 1: "+str_query1+"\n")
                write_logfile.write("Query 2: "+str_query2+"\n")
                int_matches = blastsearch(str_query1, str_query2, str_dbname, flt_evalue, write_logfile, self)

            self.txt_progress.setText("Search Complete.\nTotal matches found: "+str(int_matches))
            lst_outputline = []
            write_logfile.write("Total matches found: "+str(int_matches))

        time_end = datetime.datetime.now()
        time_difference = time_end - time_start
        write_logfile.write("Program complete. Total runtime: "+str(time_difference)+"\n")
        print("Program complete. Total runtime: "+str(time_difference)+"\n")
        write_logfile.close()



### 0 - FUNCTIONS FOR CREATING DB
def reset_files(dbname):
    #maybe just delete and remake the whole folder?
    print(dbname) #testing
    if os.path.exists(dbname):
        #Create a message box telling the user the file already exists
        root = Tk() 
        root.withdraw() #to hide another dialogue that shows up with the message box
        #root.geometry("300x200") 
        w = Label(root, text ='File already exists')          
        res = messagebox.askyesno("askyesno", "Overwrite?") 

        if not res: return False
        else: shutil.rmtree(dbname)

    os.mkdir(dbname)
    fastaname = Path(dbname) / (dbname+".fasta")
    open(fastaname,"w").close()
    return True

def get_affinity(df_interactions, pdb_code):
    affinity_line = str(df_interactions[df_interactions["PDB code"]==pdb_code]["Affinity Data"])
    kvalue = affinity_line.split()[1]
    return kvalue

def parse_pdb_ca(pdb_filename):
    #PARSING PDB FILE, returns a dictionary containing information for the whole peptide
    #whereby the keys are the name of the chain, the values are a list of all amino acids
    all_seqs = {} #Dictionary 
    chain_names = []
    parser = PDBParser()
    structure = parser.get_structure('myfile',pdb_filename)

    for model in structure.get_list():
        for chain in model.get_list():
            current_chain = str(chain.get_id()) #getting the name of the chain i.e. 'A'
            for residue in chain.get_list():
                #calculating by alpha-carbon
                res = residue.get_resname()
                if res == "HOH" or len(res) < 3: continue
                if residue.has_id("CA"):
                    ca = residue["CA"] #don't know why we have to do this
                    #lst_info: ['D696','D',[-6.401 14.699 33.043]]
                    #res: extracting 3-letter amino acid, from the header
                    try: #gives error sometimes, for some files - decided to move on if RE does not match
                        #res = re.search(r'<* ([A-Z]+) *',str(residue)).group(1)
                        lst_info = [current_chain+str(residue.get_id()[1]),dic_aa[res],ca.get_coord()]           
                        if current_chain not in all_seqs:
                            all_seqs[current_chain] = []
                            chain_names.append(current_chain)
                        all_seqs[current_chain].append(lst_info)
                    except:
                        print("cannot read: "+pdb_filename+": "+str(residue))                
    return all_seqs, chain_names

def parse_pdb_sc(pdb_filename):
    #PARSING PDB FILE, returns a dictionary containing information for the whole peptide
    #whereby the keys are the name of the chain, the values are a list of all amino acids
    all_seqs = {} #Dictionary 
    chain_names = []
    parser = PDBParser()
    structure = parser.get_structure('myfile',pdb_filename)

    for model in structure.get_list():
        for chain in model.get_list():
            current_chain = str(chain.get_id()) #getting the name of the chain i.e. 'A'
            for residue in chain.get_list():
                x_coord = 0.0; y_coord = 0.0; z_coord = 0.0; count_sc = 0
                res = residue.get_resname()
                if res == "HOH" or len(res) < 3: continue #go to next iteration in for-loop
                for atom in residue.get_list()[3:]:
                    str_id = atom.get_id()
                    if str_id == "N" or str_id == "O" or str_id == "C": continue 
                    else:
                        x_coord += atom.get_coord()[0]
                        y_coord += atom.get_coord()[1]
                        z_coord += atom.get_coord()[2]
                        count_sc += 1
                    tpl_center = (round(x_coord/count_sc,3),round(y_coord/count_sc,3),round(z_coord/count_sc,3))

                try: #gives error sometimes, for some files - decided to move on if RE does not match
                    lst_info = [current_chain+str(residue.get_id()[1]),dic_aa[res],tpl_center]           
                    if current_chain not in all_seqs:
                        all_seqs[current_chain] = []
                        chain_names.append(current_chain)
                    all_seqs[current_chain].append(lst_info)
                except:
                    print("cannot read: "+pdb_filename+": "+str(residue))                
    return all_seqs, chain_names

def cut_peptide(peptide, pep_length, overlap):
    #Cuts a peptide, returns a list of fragments, taking into account some overlaps
    #I think that this function can work on both lists and strings
    startpos = 0
    endpos = pep_length
    frags = []
    while True:
        #print("HERE: ",peptide[startpos:endpos])
        frags.append(peptide[startpos:endpos])
        startpos = endpos-overlap
        endpos = startpos+pep_length
        if startpos > len(peptide): break
    return frags

def count_hits(p_1,p_2, distance_cutoff, params):
    hit_count = 0
    #we are assuming 
    #print(len(p_1),len(p_2)) #p_1, p_2 are the peptide fragments
    for a in range(len(p_1)):
        dist = 0 #if two residues are very far apart, skip the next few calculations
        for b in range(len(p_2)):
            if dist > 20.0:
                dist -= 6.0
            else:
                pos_1 = p_1[a][2]
                pos_2 = p_2[b][2]
                #Calculating the pythagorean molecular distance between two carbon alphas
                dist = ((pos_2[0]-pos_1[0])**2+(pos_2[1]-pos_1[1])**2+(pos_2[2]-pos_1[2])**2)**0.5
                if round(dist,3) <= distance_cutoff: 
                    hit_count += 1
                    if params['log_details']:
                        params['log'].write('hit-'+str(hit_count)+': ')
                        params['log'].write(str(p_1[a][2])+" "+str(p_1[a][0])+'--')
                        params['log'].write(str(p_2[b][0])+" "+str(p_2[b][2])+' --- ')
                        params['log'].write('dist:'+str(round(dist,3))+'\t\n')
    return hit_count

def writelines(pdb_filename, params, dic_info, self): #THIS IS THE MASTER FUNCTION
    dic_seqs = {}
    chain_names = []
    
    #extracting parameters from dictionary: params
    writeinfo_fasta = open(params['output'],"at")
    writeinfo_logfile = params['log']
    distance_cutoff = params['d']
    pep_length = params['l']
    overlap = params['o']
    min_hits = params['h']

    dic_info['num_filesread'] += 1

    #dic_seqs holds information for each amino acid in the form ['D696','D',[-6.401 14.699 33.043]]
    #populate dictionary by parsing the PDBs
    if params["sc"]: dic_seqs, chain_names = parse_pdb_sc(pdb_filename)
    else: dic_seqs, chain_names = parse_pdb_ca(pdb_filename)
    
    writeinfo_logfile.write(str(pdb_filename)+" chains: "+str(chain_names)+"\n")
    pdb_code = re.search(r"^(PP(\\|\/))?([a-z0-9]+)(\.ent)?\.pdb",pdb_filename).group(3)
    
    if params["binding"]:   
        xls_filename = "MyResult.xls"
        readxls = pd.ExcelFile(xls_filename)
        df_interactions = readxls.parse("MyResult")
        affinity = get_affinity(df_interactions, pdb_code)
    else:
        affinity = "na"

    for c1 in range(len(chain_names)-1): #0,1,2
        for c2 in range(c1+1,len(chain_names)): #(1,2,3),(2,3),(3)
            #at this point, we have the sequences from chain1 and chain2 ready for comparison
            peps_1 = cut_peptide(dic_seqs[chain_names[c1]],pep_length,overlap)
            peps_2 = cut_peptide(dic_seqs[chain_names[c2]],pep_length,overlap)
            #now, we have two lists (peps_1, peps_2), each containing fragmented peptide sequences
            writeinfo_logfile.write("Comparing chains "+chain_names[c1]+" with "+chain_names[c2]+"...\n")
            self.txt_progress.setText("Reading file "+str(dic_info['num_filesread'])+" of "+str(dic_info['num_files'])+"\n" \
                "Comparing chains "+chain_names[c1]+" with "+chain_names[c2]+"...\n")

            int_numpeps = 0 #how many peps segs we have read so far (for progress bar purposes)
            int_totalpeps = len(peps_1)*len(peps_2) #total number of peptides segs to read (for progress bar)
            for p1 in peps_1: #each peptide in a list of peptides
                for p2 in peps_2:
                    int_numpeps+=1 
                    if len(p1)>19 and len(p2)>19: #because blastdb does not take sequences less than 19
                        num_hits = count_hits(p1, p2, distance_cutoff, params)
                        if num_hits > 0:
                            str_aa1 = ''
                            for aa in p1:
                                str_aa1 += aa[1]
                            str_aa2 = ''
                            for aa in p2:
                                str_aa2 += aa[1]
                            #str_line is the what we will print into the output file
                            #instead of binary(1/0) for contact/no contact, we are counting number of hits 
                            if num_hits >= min_hits:
                                dic_info['num_seqsread'] += 1
                                writeinfo_fasta.write(">A|"+str(dic_info['num_seqsread'])+"|"+str(num_hits)+"|"+affinity+"|"+pdb_code+"\n")
                                writeinfo_fasta.write(str_aa1+"\n")
                                writeinfo_fasta.write(">B|"+str(dic_info['num_seqsread'])+"|"+str(num_hits)+"|"+affinity+"|"+pdb_code+"\n")
                                writeinfo_fasta.write(str_aa2+"\n")
                    self.prg_progress.setValue(int(int_numpeps/int_totalpeps*100))

    writeinfo_fasta.close()
    return dic_info

### 1 - HIT FILTERING FUNCTIONS
def filterhits(str_inputfile, str_outputfile, int_minhits, self): #input and output files from base directory
    int_current = 0 #number of sequences already read
    int_newseqs = 0 #number of sequences after filtering

    #Check how many lines are in a file
    with open(str_inputfile) as my_file:
        int_oldseqs = sum(1 for _ in my_file)
    int_oldseqs = int_oldseqs / 4 #because def1, seq1, def2, seq2

    readfile = open(str_inputfile, "rt")
    writefile = open(str_outputfile, "wt")
    while True:	
        tagline_1 = readfile.readline()
        if not tagline_1:
            break
        seqline_1 = readfile.readline()
        tagline_2 = readfile.readline()
        seqline_2 = readfile.readline()
        int_current += 1
        self.prg_progress.setValue(int_current/int_oldseqs*100)

        int_hits = int(tagline_1.split("|")[2])
        if int_hits >= int_minhits:
            int_newseqs += 1
            writefile.write(tagline_1)
            writefile.write(seqline_1)
            writefile.write(tagline_2)
            writefile.write(seqline_2)

    readfile.close()
    writefile.close()
    #os.chdir("..")
    return int_newseqs, int_oldseqs

# 2 - FUNCTIONS FOR SEGMENTING PEPTIDES
def segmenter(str_inputfile, str_outputfile, int_lenfrag, int_overlap):
    seqFile=open(str_inputfile,"r")
    seqIOObject=SeqIO.parse(seqFile,"fasta")
    writefile = open(str_outputfile,"wt")

    for record in seqIOObject:
        temp_name = record.description
        print("...reading record \'"+temp_name+"\'")
        seq = str(record.seq)
        seq_cut = []
        startpos = 0
        endpos = int_lenfrag
        num_pep = 1
        while True:
            writefile.write(">"+temp_name+"|frag-"+str(num_pep)+"\n")
            writefile.write(seq[startpos:endpos]+"\n")
            num_pep += 1
            startpos = endpos-int_overlap
            endpos = startpos+int_lenfrag
            #don't want to enter sequences that are too short (19 is the limit for creating blastdb)
            if startpos + int_lenfrag > len(seq):
                 break
    writefile.close()
    seqFile.close()
    return num_pep

# 3 - FUNCTIONS FOR BLAST SEARCHING
#Takes the basename of the db as a string, creates a BLAST database
def makedb(str_filename):
    os.chdir(str_filename)
    db_name = str_filename+".fasta"
    #user needs to have makeblastdb added to system path
    cmd_makedb = "makeblastdb -in "+db_name+" -blastdb_version 5 -title "+db_name[:-6]+" -dbtype prot"
    os.system(cmd_makedb)
    os.chdir("..")

def blastsearch(query_1, query_2, db_name, ev, write_logfile, self):

    dic_searchresults = {}
    db_folder = Path(db_name)
    db_fasta = db_folder / str(db_name+".fasta")
    blastres_1 = str(db_folder / str("blastres_"+db_name+"_1.xml"))
    blastres_2 = db_folder / str("blastres_"+db_name+"_2.xml")

    int_totalmatches = 0
    #Searching through the first query
    self.txt_progress.setText("...performing BLAST search for "+str(query_1))
    write_logfile.write("Performing BLAST search for "+str(query_1)+"\n")
    # name of the resulting blast search (format 5: .xml)
    print(blastres_1, query_1, db_fasta, ev)
    blastx_cline_1 = NcbiblastxCommandline(cmd='blastp', out=blastres_1, outfmt=5, query=query_1, db=db_fasta, evalue=ev)
    os.system(str(blastx_cline_1))

    print("First BLAST search complete.")
    self.prg_progress.setValue(20)
    #Parsing results from the first query
    self.txt_progress.setText("1st BLAST search complete. Now parsing results...")    
    int_hits = parserecords(blastres_1, 1, dic_searchresults, write_logfile)
    self.prg_progress.setValue(40)
    if int_hits < 1: 
        self.txt_progress.setText("Parsing complete.\nNo hits found for: "+query_1)
        write_logfile.write("Parsing complete.\nNo hits found for: "+query_1+"\n")
        self.prg_progress.setValue(100)
        return -1
    write_logfile.write("Total hits found: "+str(int_hits)+"\n")

    #Searching through the second query
    self.txt_progress.setText("Parsing complete. Hits found: "+str(int_hits) \
        +"...performing BLAST search for "+str(query_2))
    blastx_cline_2 = NcbiblastxCommandline(cmd='blastp', out=blastres_2, outfmt=5, query=query_2, db=db_fasta, evalue=ev)
    os.system(str(blastx_cline_2))

    print("Second BLAST search complete.")
    self.prg_progress.setValue(60)
    self.txt_progress.setText("2nd BLAST search complete. Now parsing results...")
    int_hits = parserecords(blastres_2, 2, dic_searchresults, write_logfile)
    self.prg_progress.setValue(80)
    if int_hits < 1: 
        self.txt_progress.setText("Parsing complete.\nNo hits found for: "+query_2)
        write_logfile.write("Parsing complete.\nNo hits found for: "+query_2+"\n")
        self.prg_progress.setValue(100)
        return -1
    write_logfile.write("Total hits found: "+str(int_hits)+"\n")
    
    self.txt_progress.setText("Parsing complete. Hits found: "+str(int_hits)+ \
        "\n...checking for matches between two queries")
    #passing "self" through function because it contains the checkbox for whether or not to overwrite output file
    #may want to find an easier way to do this
    int_totalmatches = checkmatch(dic_searchresults, write_logfile, self)
    self.prg_progress.setValue(100)
    return int_totalmatches

def parserecords(xmloutput, searchtype, dic_results, write_logfile):
    #blast results are already sorted from best hit
    result_handle = open(xmloutput)
    blast_records = NCBIXML.parse(result_handle)
    lst_hits = []
    int_hits = 0
    #creating a holding variable to skip duplicates
    ex_query = ""; ex_sbjct = ""

    try:
        for blast_record in blast_records:
            #since we used blast with an e-value already set, we don't need this if-statement
            #E_VALUE_THRESH = 0.04
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #if hsp.expect<E_VALUE_THRESH:
                    my_title = str(alignment.title).split()[1]
                    matchtype = my_title.split("|")[0] #should give A or B
                    
                    #depending on the search type, you only want to put 1A and 2B into the dictionary
                    if (searchtype == 1 and matchtype == "A") or (searchtype == 2 and matchtype == "B"):
                        tup_line = my_title, alignment.length, hsp.expect, hsp.query, hsp.match, hsp.sbjct
                        ex_query = hsp.query
                        ex_sbjct = hsp.sbjct
                        lst_hits.append(tup_line)
                        int_hits += 1
                        write_logfile.write("****Alignment-"+matchtype+"****\n")
                        write_logfile.write("sequence:"+str(alignment.title)+", length:"+str(alignment.length)+", evalue:"+str(hsp.expect)+"\n")
                        write_logfile.write("query: "+hsp.query[0:75]+"\n")
                        write_logfile.write("match: "+hsp.match[0:75]+"\n")
                        write_logfile.write("sbjct: "+hsp.sbjct[0:75]+"\n")

    except Exception as errormessage:
        write_logfile.write("Error reading file: "+str(xmloutput)+"\n")
        print(errormessage)
        int_hits = -1
    
    dic_results[searchtype] = lst_hits
    result_handle.close()
    return int_hits
    
def checkmatch(dic_results, write_logfile, self):
    int_matchcount = 0
    indices_1 = {}

    ex_query1 = ""
    ex_query2 = ""

    try:
        if self.chk2_overwritelog.isChecked():
            write_outfile = open("searchdb_output.csv","wt")
            write_outfile.write("time,database,order,query1,seq1,ev1,query2,seq2,ev2,origin,avg\n")
        else:
            write_outfile = open("searchdb_output.csv","at")
    except Exception as errormessage:
        print("Cannot open file: searchdb_output.csv\n"+errormessage)
        return -1


    for alignment_1 in dic_results[1]:
        current_index = (alignment_1[0].split("|")[1])
        indices_1[current_index] = alignment_1

    for alignment_2 in dic_results[2]:
        current_index = (alignment_2[0].split("|")[1])
        if current_index in indices_1:
            #print(indices_a[current_index],case)
            #title, evalue, query
            #tup_line = my_title, alignment.length, hsp.expect, hsp.query, hsp.match, hsp.sbjct
            if (ex_query1 == indices_1[current_index][3]) and (ex_query2 == alignment_2[3]):
                print("repeated match: "+ex_query1+" and "+ex_query2)
                #int_matchcount += 1
                continue
            ex_query1 = indices_1[current_index][3]
            ex_query2 = alignment_2[3]
            int_matchcount += 1

            #old way of writing lines
            line_1 = indices_1[current_index][0]+","+str(indices_1[current_index][2])+",query:"+indices_1[current_index][3] \
                +",match:"+indices_1[current_index][4]+ ",sbjct:"+indices_1[current_index][5]+"\n"
            line_2 = alignment_2[0]+","+str(alignment_2[2])+",query:"+alignment_2[3] \
                +",match:"+alignment_2[4]+",sbjct:"+alignment_2[5]+"\n"

            write_logfile.write("Match found ("+str(int_matchcount)+"):\n")
            write_logfile.write(line_1.replace(",",", ")+"\n")
            write_logfile.write(line_2.replace(",",", ")+"\n")

            #New way of writing lines
            time_current = datetime.datetime.now()
            str_dbname = self.cmb2_db.currentText()
            str_queryname1 = os.path.basename(self.txt2_query1.text())
            str_queryname2 = os.path.basename(self.txt2_query2.text())
            flt_avgevalue = round((indices_1[current_index][2]+alignment_2[2]) / 2, 3)
            lst_line = [time_current, str_dbname, int_matchcount, \
                str_queryname1, indices_1[current_index][3], indices_1[current_index][2],\
                str_queryname2, alignment_2[3], alignment_2[2], alignment_2[0].split("|",1)[1], flt_avgevalue]
            str_line = ",".join(str(item).replace(",","-") for item in lst_line)
            write_outfile.write(str_line+"\n")

    write_outfile.close()

    return int_matchcount

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
