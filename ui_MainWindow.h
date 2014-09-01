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
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
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
    QGroupBox *groupBox;
    QGridLayout *gridLayout_3;
    QLabel *label_3;
    QDoubleSpinBox *m_k;
    QLabel *label_4;
    QDoubleSpinBox *m_b;
    QLabel *label_5;
    QDoubleSpinBox *m_move;
    QGroupBox *s_transformGB;
    QGridLayout *gridLayout;
    QLabel *label;
    QComboBox *m_selected;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_2;
    QPushButton *simBtn;
    QDoubleSpinBox *m_dt;
    QSpinBox *m_timerValue;
    QLabel *label_7;
    QLabel *label_6;
    QPushButton *expBtn;
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

        s_mainWindowGridLayout->addItem(verticalSpacer, 3, 1, 1, 1);

        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_3 = new QGridLayout(groupBox);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout_3->addWidget(label_3, 0, 0, 1, 1);

        m_k = new QDoubleSpinBox(groupBox);
        m_k->setObjectName(QStringLiteral("m_k"));
        m_k->setMinimum(-200);
        m_k->setMaximum(200);
        m_k->setSingleStep(0.01);
        m_k->setValue(0);

        gridLayout_3->addWidget(m_k, 0, 1, 1, 1);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QStringLiteral("label_4"));

        gridLayout_3->addWidget(label_4, 1, 0, 1, 1);

        m_b = new QDoubleSpinBox(groupBox);
        m_b->setObjectName(QStringLiteral("m_b"));
        m_b->setMinimum(-200);
        m_b->setMaximum(200);
        m_b->setSingleStep(0.01);
        m_b->setValue(0);

        gridLayout_3->addWidget(m_b, 1, 1, 1, 1);

        label_5 = new QLabel(groupBox);
        label_5->setObjectName(QStringLiteral("label_5"));

        gridLayout_3->addWidget(label_5, 2, 0, 1, 1);

        m_move = new QDoubleSpinBox(groupBox);
        m_move->setObjectName(QStringLiteral("m_move"));
        m_move->setMinimum(-10);
        m_move->setMaximum(10);
        m_move->setSingleStep(0.5);
        m_move->setValue(0);

        gridLayout_3->addWidget(m_move, 2, 1, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox, 1, 1, 1, 1);

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
        simBtn = new QPushButton(groupBox_2);
        simBtn->setObjectName(QStringLiteral("simBtn"));
        simBtn->setCheckable(true);
        simBtn->setChecked(false);

        gridLayout_2->addWidget(simBtn, 2, 1, 1, 1);

        m_dt = new QDoubleSpinBox(groupBox_2);
        m_dt->setObjectName(QStringLiteral("m_dt"));
        m_dt->setDecimals(3);
        m_dt->setMinimum(0.001);
        m_dt->setMaximum(4);
        m_dt->setSingleStep(0.001);
        m_dt->setValue(0.033);

        gridLayout_2->addWidget(m_dt, 0, 1, 1, 1);

        m_timerValue = new QSpinBox(groupBox_2);
        m_timerValue->setObjectName(QStringLiteral("m_timerValue"));
        m_timerValue->setMaximum(1000);
        m_timerValue->setValue(33);

        gridLayout_2->addWidget(m_timerValue, 1, 1, 1, 1);

        label_7 = new QLabel(groupBox_2);
        label_7->setObjectName(QStringLiteral("label_7"));

        gridLayout_2->addWidget(label_7, 1, 0, 1, 1);

        label_6 = new QLabel(groupBox_2);
        label_6->setObjectName(QStringLiteral("label_6"));

        gridLayout_2->addWidget(label_6, 0, 0, 1, 1);

        expBtn = new QPushButton(groupBox_2);
        expBtn->setObjectName(QStringLiteral("expBtn"));

        gridLayout_2->addWidget(expBtn, 3, 1, 1, 1);


        s_mainWindowGridLayout->addWidget(groupBox_2, 2, 1, 1, 1);

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
        groupBox->setTitle(QApplication::translate("MainWindow", "Hair Properties", 0));
        label_3->setText(QString());
        label_4->setText(QString());
        label_5->setText(QApplication::translate("MainWindow", "move", 0));
        s_transformGB->setTitle(QApplication::translate("MainWindow", "Select Object", 0));
        label->setText(QApplication::translate("MainWindow", "Selected Object", 0));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Simulation", 0));
        simBtn->setText(QApplication::translate("MainWindow", "sim", 0));
        label_7->setText(QApplication::translate("MainWindow", "timer (ms)", 0));
        label_6->setText(QApplication::translate("MainWindow", "dt step", 0));
        expBtn->setText(QApplication::translate("MainWindow", "export", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
