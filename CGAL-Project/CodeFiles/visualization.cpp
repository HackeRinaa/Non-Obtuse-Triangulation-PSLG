#include "../HeaderFiles/visualization.h"
#include <QVBoxLayout>

VisualizationWindow::VisualizationWindow(QWidget *parent) : QMainWindow(parent)
{
    setWindowTitle("CGAL Visualization");
    textEdit = new QTextEdit(this);
    textEdit->setReadOnly(true);

    QVBoxLayout *layout = new QVBoxLayout();
    layout->addWidget(textEdit);

    QWidget *centralWidget = new QWidget(this);
    centralWidget->setLayout(layout);
    setCentralWidget(centralWidget);
}

void VisualizationWindow::displayTriangulation(const CDT &cdt)
{
    QString output = "Initial triangulation:\n";
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
    {
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        output += QString("Edge (%1, %2) - (%3, %4)\n")
                      .arg(v1->point().x())
                      .arg(v1->point().y())
                      .arg(v2->point().x())
                      .arg(v2->point().y());
    }
    textEdit->setPlainText(output);
}
