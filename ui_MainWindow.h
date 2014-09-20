/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.1.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QGridLayout *s_mainWindowGridLayout;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QGroupBox *s_transformGB;
    QGridLayout *gridLayout;
    QLabel *label;
    QComboBox *m_selected;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_2;
    QDoubleSpinBox *m_timeStep;
    QLabel *label_7;
    QLabel *label_6;
    QSpinBox *m_timerUpdate;
    QLabel *label_2;
    QSpinBox *m_simIter;
    QPushButton *m_stepForward;
    QPushButton *m_simBtn;
    QPushButton *m_recordBtn;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_3;
    QLabel *label_3;
    QComboBox *m_minimizationMethod;
    QLabel *label_4;
    QDoubleSpinBox *m_minTolerance;
    QLabel *label_5;
    QSpinBox *m_minMaxIter;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout_6;
    QDoubleSpinBox *m_stiction;
    QLabel *label_15;
    QLabel *label_14;
    QCheckBox *m_collisions;
    QCheckBox *m_selfInteractions;
    QDoubleSpinBox *m_repulsion;
    QGroupBox *groupBox_4;
    QGridLayout *gridLayout_4;
    QDoubleSpinBox *m_drag;
    QLabel *label_11;
    QDoubleSpinBox *m_bendStiffness;
    QDoubleSpinBox *m_twistStiffness;
    QDoubleSpinBox *m_maxForce;
    QLabel *label_9;
    QLabel *label_12;
    QLabel *label_10;
    QLabel *label_13;
    QSpinBox *m_pbdIter;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_5;
    QPushButton *m_selectDirBtn;
    QPushButton *m_expBtn;
    QLineEdit *m_filePrefix;
    QLabel *label_8;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(1464, 1028);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        s_mainWindowGridLayout = new QGridLayout(centralwidget);
        s_mainWindowGridLayout->setObjectName(QStringLiteral("s_mainWindowGridLayout"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        s_mainWindowGridLayout->addItem(horizontalSpacer, 0, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        s_mainWindowGridLayout->addItem(verticalSpacer, 8, 1, 1, 1);

        s_transformGB = new QGroupBox(centralwidget);
        s_transformGB->setObjectName(QStringLiteral("s_transformGB"));
        gridLayout = new QGridLayout(s_transformGB);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label = new QLabel(s_transformGB);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        m_selected = new QComboBox(s_transformGB);
        m_selected->setObjectName(QStringLiteral("m_selected"));
        m_selected->setEditable(false);

        gridLayout->addWidget(m_selected, 0, 1, 1, 1);


        s_mainWindowGridLayout->addWidget(s_transformGB, 0, 1, 1, 1);

        groupBox_2 = new QGroupBox(centralwidget);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        gridLayout_2 = new QGridLayout(groupBox_2);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        m_timeStep = new QDoubleSpinBox(groupBox_2);
        m_timeStep->setObjectName(QStringLiteral("m_timeStep"));
        m_timeStep->setDecimals(6);
        m_timeStep->setMinimum(1e-06);
        m_timeStep->setMaximum(0.3);
        m_timeStep->setSingleStep(0.001);
        m_timeStep->setValue(0.01);

        gridLayout_2->addWidget(m_timeStep, 0, 1, 1, 1);

        label_7 = new QLabel(groupBox_2);
        label_7->setObjectName(QStringLiteral("label_7"));

        gridLayout_2->addWidget(label_7, 2, 0, 1, 1);

        label_6 = new QLabel(groupBox_2);
        label_6->setObjectName(QStringLiteral("label_6"));

        gridLayout_2->addWidget(label_6, 0, 0, 1, 1);

        m_timerUpdate = new QSpinBox(groupBox_2);
        m_timerUpdate->setObjectName(QStringLiteral("m_timerUpdate"));
        m_timerUpdate->setMaximum(1000);
        m_timerUpdate->setValue(33);

        gridLayout_2->addWidget(m_timerUpdate, 2, 1, 1, 1);

        label_2 = new QLabel(groupBox_2);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout_2->addWidget(label_2, 1, 0, 1, 1);

        m_simIter = new QSpinBox(groupBox_2);
        m_simIter->setObjectName(QStringLiteral("m_simIter"));
        m_simIter->setMinimum(1);
        m_simIter->setMaximum(100);

        gridLayout_2->addWidget(m_simIter, 1, 1, 1, 1);

        m_stepForward = new QPushButton(groupBox_2);
        m_stepForward->setObjectName(QStringLiteral("m_stepForward"));

        gridLayout_2->addWidget(m_stepForward, 3, 1, 1, 1);

        m_simBtn = new QPushButton(groupBox_2);
        m_simBtn->setObjectName(QStringLiteral("m_simBtn"));
        m_simBtn->setCheckable(true);
        m_simBtn->setChecked(false);

        gridLayout_2->addWidget(m_simBtn, 3, 0, 1, 1);

        m_recordBtn = new QPushButton(groupBox_2);
        m_recordBtn->setObjectName(QStringLiteral("m_recordBtn"));
        m_recordBtn->setCheckable(true);

        gridLayout_2->addWidget(m_recordBtn, 4, 0, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox_2, 5, 1, 1, 1);

        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_3 = new QGridLayout(groupBox);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout_3->addWidget(label_3, 0, 0, 1, 1);

        m_minimizationMethod = new QComboBox(groupBox);
        m_minimizationMethod->setObjectName(QStringLiteral("m_minimizationMethod"));

        gridLayout_3->addWidget(m_minimizationMethod, 0, 1, 1, 1);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QStringLiteral("label_4"));

        gridLayout_3->addWidget(label_4, 1, 0, 1, 1);

        m_minTolerance = new QDoubleSpinBox(groupBox);
        m_minTolerance->setObjectName(QStringLiteral("m_minTolerance"));
        m_minTolerance->setDecimals(6);
        m_minTolerance->setMaximum(10);
        m_minTolerance->setSingleStep(0.01);

        gridLayout_3->addWidget(m_minTolerance, 1, 1, 1, 1);

        label_5 = new QLabel(groupBox);
        label_5->setObjectName(QStringLiteral("label_5"));

        gridLayout_3->addWidget(label_5, 2, 0, 1, 1);

        m_minMaxIter = new QSpinBox(groupBox);
        m_minMaxIter->setObjectName(QStringLiteral("m_minMaxIter"));
        m_minMaxIter->setMinimum(0);
        m_minMaxIter->setMaximum(100);

        gridLayout_3->addWidget(m_minMaxIter, 2, 1, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox, 3, 1, 1, 1);

        groupBox_5 = new QGroupBox(centralwidget);
        groupBox_5->setObjectName(QStringLiteral("groupBox_5"));
        gridLayout_6 = new QGridLayout(groupBox_5);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        m_stiction = new QDoubleSpinBox(groupBox_5);
        m_stiction->setObjectName(QStringLiteral("m_stiction"));
        m_stiction->setDecimals(10);
        m_stiction->setMaximum(1);
        m_stiction->setSingleStep(1e-06);

        gridLayout_6->addWidget(m_stiction, 3, 0, 1, 1);

        label_15 = new QLabel(groupBox_5);
        label_15->setObjectName(QStringLiteral("label_15"));

        gridLayout_6->addWidget(label_15, 4, 0, 1, 1);

        label_14 = new QLabel(groupBox_5);
        label_14->setObjectName(QStringLiteral("label_14"));

        gridLayout_6->addWidget(label_14, 2, 0, 1, 1);

        m_collisions = new QCheckBox(groupBox_5);
        m_collisions->setObjectName(QStringLiteral("m_collisions"));

        gridLayout_6->addWidget(m_collisions, 0, 0, 1, 1);

        m_selfInteractions = new QCheckBox(groupBox_5);
        m_selfInteractions->setObjectName(QStringLiteral("m_selfInteractions"));

        gridLayout_6->addWidget(m_selfInteractions, 1, 0, 1, 1);

        m_repulsion = new QDoubleSpinBox(groupBox_5);
        m_repulsion->setObjectName(QStringLiteral("m_repulsion"));
        m_repulsion->setDecimals(10);
        m_repulsion->setMaximum(1);
        m_repulsion->setSingleStep(1e-06);

        gridLayout_6->addWidget(m_repulsion, 5, 0, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox_5, 4, 1, 1, 1);

        groupBox_4 = new QGroupBox(centralwidget);
        groupBox_4->setObjectName(QStringLiteral("groupBox_4"));
        gridLayout_4 = new QGridLayout(groupBox_4);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        m_drag = new QDoubleSpinBox(groupBox_4);
        m_drag->setObjectName(QStringLiteral("m_drag"));
        m_drag->setDecimals(10);
        m_drag->setMaximum(1);
        m_drag->setSingleStep(1e-06);

        gridLayout_4->addWidget(m_drag, 3, 1, 1, 1);

        label_11 = new QLabel(groupBox_4);
        label_11->setObjectName(QStringLiteral("label_11"));

        gridLayout_4->addWidget(label_11, 2, 0, 1, 1);

        m_bendStiffness = new QDoubleSpinBox(groupBox_4);
        m_bendStiffness->setObjectName(QStringLiteral("m_bendStiffness"));
        m_bendStiffness->setDecimals(6);
        m_bendStiffness->setMaximum(2);
        m_bendStiffness->setSingleStep(1e-06);

        gridLayout_4->addWidget(m_bendStiffness, 0, 1, 1, 1);

        m_twistStiffness = new QDoubleSpinBox(groupBox_4);
        m_twistStiffness->setObjectName(QStringLiteral("m_twistStiffness"));
        m_twistStiffness->setDecimals(6);
        m_twistStiffness->setMaximum(2);
        m_twistStiffness->setSingleStep(1e-06);

        gridLayout_4->addWidget(m_twistStiffness, 1, 1, 1, 1);

        m_maxForce = new QDoubleSpinBox(groupBox_4);
        m_maxForce->setObjectName(QStringLiteral("m_maxForce"));
        m_maxForce->setDecimals(6);
        m_maxForce->setMinimum(1);
        m_maxForce->setMaximum(1000);

        gridLayout_4->addWidget(m_maxForce, 2, 1, 1, 1);

        label_9 = new QLabel(groupBox_4);
        label_9->setObjectName(QStringLiteral("label_9"));

        gridLayout_4->addWidget(label_9, 0, 0, 1, 1);

        label_12 = new QLabel(groupBox_4);
        label_12->setObjectName(QStringLiteral("label_12"));

        gridLayout_4->addWidget(label_12, 3, 0, 1, 1);

        label_10 = new QLabel(groupBox_4);
        label_10->setObjectName(QStringLiteral("label_10"));

        gridLayout_4->addWidget(label_10, 1, 0, 1, 1);

        label_13 = new QLabel(groupBox_4);
        label_13->setObjectName(QStringLiteral("label_13"));

        gridLayout_4->addWidget(label_13, 4, 0, 1, 1);

        m_pbdIter = new QSpinBox(groupBox_4);
        m_pbdIter->setObjectName(QStringLiteral("m_pbdIter"));
        m_pbdIter->setMinimum(0);
        m_pbdIter->setMaximum(1000);

        gridLayout_4->addWidget(m_pbdIter, 4, 1, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox_4, 1, 1, 1, 1);

        groupBox_3 = new QGroupBox(centralwidget);
        groupBox_3->setObjectName(QStringLiteral("groupBox_3"));
        gridLayout_5 = new QGridLayout(groupBox_3);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        m_selectDirBtn = new QPushButton(groupBox_3);
        m_selectDirBtn->setObjectName(QStringLiteral("m_selectDirBtn"));

        gridLayout_5->addWidget(m_selectDirBtn, 5, 0, 1, 1);

        m_expBtn = new QPushButton(groupBox_3);
        m_expBtn->setObjectName(QStringLiteral("m_expBtn"));
        m_expBtn->setCheckable(false);

        gridLayout_5->addWidget(m_expBtn, 5, 1, 1, 1);

        m_filePrefix = new QLineEdit(groupBox_3);
        m_filePrefix->setObjectName(QStringLiteral("m_filePrefix"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(m_filePrefix->sizePolicy().hasHeightForWidth());
        m_filePrefix->setSizePolicy(sizePolicy);

        gridLayout_5->addWidget(m_filePrefix, 3, 1, 1, 1);

        label_8 = new QLabel(groupBox_3);
        label_8->setObjectName(QStringLiteral("label_8"));

        gridLayout_5->addWidget(label_8, 3, 0, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox_3, 7, 1, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1464, 20));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        MainWindow->setStatusBar(statusbar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "ElasticRods", 0));
        s_transformGB->setTitle(QApplication::translate("MainWindow", "Select Object", 0));
        label->setText(QApplication::translate("MainWindow", "SelectedObject", 0));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Animation", 0));
        label_7->setText(QApplication::translate("MainWindow", "timer (ms)", 0));
        label_6->setText(QApplication::translate("MainWindow", "time step (s)", 0));
        label_2->setText(QApplication::translate("MainWindow", "iter per step", 0));
        m_stepForward->setText(QApplication::translate("MainWindow", ">>", 0));
        m_simBtn->setText(QApplication::translate("MainWindow", "sim", 0));
        m_recordBtn->setText(QApplication::translate("MainWindow", "record", 0));
        groupBox->setTitle(QApplication::translate("MainWindow", "Minimization", 0));
        label_3->setText(QApplication::translate("MainWindow", "method", 0));
        label_4->setText(QApplication::translate("MainWindow", "tolerance", 0));
        label_5->setText(QApplication::translate("MainWindow", "max iterations", 0));
        groupBox_5->setTitle(QApplication::translate("MainWindow", "Collisions and Self-interactions", 0));
        label_15->setText(QApplication::translate("MainWindow", "repulsion", 0));
        label_14->setText(QApplication::translate("MainWindow", "stiction", 0));
        m_collisions->setText(QApplication::translate("MainWindow", "resolve collisions", 0));
        m_selfInteractions->setText(QApplication::translate("MainWindow", "resolve self-interactions", 0));
        groupBox_4->setTitle(QApplication::translate("MainWindow", "Hair Props", 0));
        label_11->setText(QApplication::translate("MainWindow", "maxElasticForce", 0));
        label_9->setText(QApplication::translate("MainWindow", "bendStiffness", 0));
        label_12->setText(QApplication::translate("MainWindow", "drag", 0));
        label_10->setText(QApplication::translate("MainWindow", "twistStiffness", 0));
        label_13->setText(QApplication::translate("MainWindow", "PBD iterations", 0));
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Export", 0));
        m_selectDirBtn->setText(QApplication::translate("MainWindow", "export dir", 0));
        m_expBtn->setText(QApplication::translate("MainWindow", "export", 0));
        m_filePrefix->setText(QApplication::translate("MainWindow", "default_", 0));
        label_8->setText(QApplication::translate("MainWindow", "file prefix", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
